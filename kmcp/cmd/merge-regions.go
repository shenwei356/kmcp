// Copyright © 2020-2022 Wei Shen <shenwei356@gmail.com>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

package cmd

import (
	"fmt"
	"os"
	"path/filepath"
	"regexp"
	"runtime"
	"strings"
	"sync"
	"time"

	"github.com/shenwei356/breader"
	"github.com/spf13/cobra"
	"github.com/zeebo/wyhash"
)

var regionsCmd = &cobra.Command{
	Use:   "merge-regions",
	Short: "Merge species/assembly-specific regions",
	Long: `Merge species/assembly-specific regions

Steps:
  # 1. Simulating reads and searching on one or more databases.
  seqkit sliding --step 10 --window 100 ref.fna.gz \
      | kmcp search -d db1.kmcp -o ref.fna.gz.kmcp@db1.tsv.gz
  seqkit sliding --step 10 --window 100 ref.fna.gz \
      | kmcp search -d db2.kmcp -o ref.fna.gz.kmcp@db2.tsv.gz
  
  # 2. Merging and filtering searching results
  kmcp merge ref.fna.gz.kmcp@*.tsv.gz \
      | kmcp utils filter -X taxdump -T taxid.map \
            -o ref.fna.gz.kmcp.uniq.tsv.gz
  
  # 3. Merging regions.
  # Here the value of --min-overlap should be k-1.
  kmcp utils merge-regions --min-overlap 20 ref.fna.gz.kmcp.uniq.tsv.gz \
      -o ref.fna.gz.kmcp.uniq.tsv.gz.bed

Output (BED6 format):
  1. chrom      - chromosome name
  2. chromStart - starting position (0-based)
  3. chromEnd   - ending position (0-based)
  4. name       - "species-specific" or "assembly-specific"
  5. score      - 0-1000, 1000 for "assembly-specific", others for ""species-specific"
  6. strand     - "."

Performance notes:
  1. Searching results are parsed in parallel, and the number of
     lines proceeded by a thread can be set by the flag --line-chunk-size.
  2. However using a lot of threads does not always accelerate
     processing, 4 threads with a chunk size of 500-5000 is fast enough.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)

		var fhLog *os.File
		if opt.Log2File {
			fhLog = addLog(opt.LogFile, opt.Verbose)
		}
		timeStart := time.Now()
		defer func() {
			if opt.Verbose || opt.Log2File {
				log.Info()
				log.Infof("elapsed time: %s", time.Since(timeStart))
				log.Info()
			}
			if opt.Log2File {
				fhLog.Close()
			}
		}()

		var err error

		outFile := getFlagString(cmd, "out-file")

		maxFPR := getFlagPositiveFloat64(cmd, "max-fpr")
		minQcov := getFlagNonNegativeFloat64(cmd, "min-query-cov")

		chunkSize := getFlagPositiveInt(cmd, "line-chunk-size")
		if opt.NumCPUs > 4 {
			if opt.Verbose || opt.Log2File {
				log.Infof("using a lot of threads does not always accelerate processing, 4-threads is fast enough")
			}
			opt.NumCPUs = 4
			runtime.GOMAXPROCS(opt.NumCPUs)
		}

		maxGap := getFlagNonNegativeInt(cmd, "max-gap")
		limitGap := maxGap > 0
		minOverlap := getFlagPositiveInt(cmd, "min-overlap")
		if minOverlap == 1 {
			log.Warningf("you may set -l/--min-overlap as k-1, where k is the k-mer size")
		}

		reQueryStr := getFlagString(cmd, "regexp")
		reQuery, err := regexp.Compile(reQueryStr)
		checkError(err)

		nameMultiple := getFlagString(cmd, "name-species")
		nameSingle := getFlagString(cmd, "name-assembly")

		ignoreType := getFlagBool(cmd, "ignore-type")

		// ---------------------------------------------------------------

		if opt.Verbose || opt.Log2File {
			log.Infof("kmcp v%s", VERSION)
			log.Info("  https://github.com/shenwei356/kmcp")
			log.Info()

			log.Info("checking input files ...")
		}
		files := getFileListFromArgsAndFile(cmd, args, true, "infile-list", true)
		if opt.Verbose || opt.Log2File {
			if len(files) == 1 && isStdin(files[0]) {
				log.Info("no files given, reading from stdin")
			} else {
				log.Infof("  %d input file(s) given", len(files))
			}
		}

		outFileClean := filepath.Clean(outFile)
		for _, file := range files {
			if isStdin(file) {
				// checkError(fmt.Errorf("stdin not supported"))
			} else if filepath.Clean(file) == outFileClean {
				checkError(fmt.Errorf("out file should not be one of the input file"))
			}
		}

		// ---------------------------------------------------------------

		// output
		outfh, gw, w, err := outStream(outFile, strings.HasSuffix(strings.ToLower(outFile), ".gz"), opt.CompressionLevel)
		checkError(err)
		defer func() {
			outfh.Flush()
			if gw != nil {
				gw.Close()
			}
			w.Close()
		}()

		if opt.Verbose || opt.Log2File {
			log.Info()
			log.Infof("-------------------- [main parameters] --------------------")

			log.Infof("match filtration: ")
			log.Infof("  maximum false positive rate: %f", maxFPR)
			log.Infof("  minimum query coverage: %4f", minQcov)
			log.Info()

			log.Infof("region merging: ")
			log.Infof("  minimum overlap of two adjacent regions: %d", minOverlap)
			log.Infof("  maximum distance of starting positions of two adjacent regions: %d", maxGap)

			log.Infof("-------------------- [main parameters] --------------------")
			log.Info()

		}
		// ---------------------------------------------------------------

		numFields := 13

		// ---------------------------------------------------------------

		pool := &sync.Pool{New: func() interface{} {
			tmp := make([]string, numFields)
			return &tmp
		}}

		fn := func(line string) (interface{}, bool, error) {
			if line == "" || line[0] == '#' { // ignoring blank line and comment line
				return "", false, nil
			}

			items := pool.Get().(*[]string)

			match, ok := parseMatchResult3(line, numFields, items, maxFPR, minQcov, reQuery)
			if !ok {
				pool.Put(items)
				return nil, false, nil
			}

			pool.Put(items)
			return match, true, nil
		}

		var nReads int
		var nRegions int

		// ---------------------------------------------------------------
		if opt.Verbose || opt.Log2File {
			log.Info("analyzing ...")
		}

		for _, file := range files {
			if opt.Verbose || opt.Log2File {
				log.Infof("  parsing file: %s", file)
			}

			var matches map[uint64]*[]*MatchResult3 // target -> match result
			var ms *[]*MatchResult3
			var ok bool
			var hTarget uint64
			var prevQuery string
			var match *MatchResult3

			var ref string
			var begin, end int
			var name string
			var score float64
			var first bool

			var ref0 string
			var begin0, end0 int
			var name0 string
			var score0 float64
			var extend bool

			var begin1, end1 int

			reader, err := breader.NewBufferedReader(file, opt.NumCPUs, chunkSize, fn)
			checkError(err)
			var data interface{}

			matches = make(map[uint64]*[]*MatchResult3)
			for chunk := range reader.Ch {
				checkError(chunk.Err)

				for _, data = range chunk.Data {
					match = data.(*MatchResult3)

					if prevQuery != match.Query { // new query
						if len(matches) > 0 { // not the first query
							nReads++

							first = true
							for _, ms = range matches {
								if first {
									ref = (*ms)[0].Ref
									begin, end = (*ms)[0].Begin, (*ms)[0].End
									score = (*ms)[0].QCov
									first = false
								} else {
									score += (*ms)[0].QCov
								}

								poolMatchResults.Put(ms)
							}
							if len(matches) == 1 {
								name = nameSingle
							} else {
								name = nameMultiple
								score /= float64(len(matches))
							}

							if begin0 > 0 {
								extend = false

								// 1. the same chromesome; 2. has overlap; 3 the same type
								if ref == ref0 && begin+minOverlap-1 <= end1 && (ignoreType || name == name0) {
									if limitGap {
										if begin-begin1 <= maxGap {
											extend = true
										}
									} else {
										extend = true
									}
								}

								if extend {
									end0 = end
									if name0 != name {
										name0 = nameMultiple
									}
									if name0 == nameMultiple {
										score0 = (score0 + score) / 2
									}
								} else { // print previous region
									nRegions++
									fmt.Fprintf(outfh, "%s\t%d\t%d\t%s\t%.0f\t.\n", ref0, begin0-1, end0, name0, score0*1000)

									ref0, begin0, end0, name0, score0 = ref, begin, end, name, score
								}
							} else { // first record
								ref0, begin0, end0, name0, score0 = ref, begin, end, name, score
							}

							begin1, end1 = begin, end
						}

						matches = make(map[uint64]*[]*MatchResult3)
					}

					hTarget = wyhash.HashString(match.Target, 1)
					if ms, ok = matches[hTarget]; !ok {
						// tmp := []*MatchResult{match}
						tmp := poolMatchResults3.Get().(*[]*MatchResult3)
						*tmp = (*tmp)[:0]
						*tmp = append(*tmp, match)
						matches[hTarget] = tmp
					} else {
						*ms = append(*ms, match)
					}

					prevQuery = match.Query
				}
			}

			if len(matches) > 0 {
				nReads++

				first = true
				for _, ms = range matches {
					if first {
						ref = (*ms)[0].Ref
						begin, end = (*ms)[0].Begin, (*ms)[0].End
						score = (*ms)[0].QCov
						first = false
					} else {
						score += (*ms)[0].QCov
					}

					poolMatchResults.Put(ms)
				}
				if len(matches) == 1 {
					name = nameSingle
				} else {
					name = nameMultiple
					score /= float64(len(matches))
				}

				extend = false

				// 1. the same chromesome; 2. has overlap; 3 the same type
				if ref == ref0 && begin+minOverlap-1 <= end1 && (ignoreType || name == name0) {
					if limitGap {
						if begin-begin1 <= maxGap {
							extend = true
						}
					} else {
						extend = true
					}
				}

				if extend {
					end0 = end
					if name0 != name {
						name0 = nameMultiple
					}
					if name0 == nameMultiple {
						score0 = (score0 + score) / 2
					}

					nRegions++
					fmt.Fprintf(outfh, "%s\t%d\t%d\t%s\t%.0f\t.\n", ref0, begin0-1, end0, name0, score0*1000)
				} else {
					nRegions++
					fmt.Fprintf(outfh, "%s\t%d\t%d\t%s\t%.0f\t.\n", ref0, begin0-1, end0, name0, score0*1000)

					nRegions++
					fmt.Fprintf(outfh, "%s\t%d\t%d\t%s\t%.0f\t.\n", ref, begin-1, end, name, score*1000)
				}
			}
		}

		if opt.Verbose || opt.Log2File {
			log.Info()
			log.Infof("%d regions merged to %d", nReads, nRegions)
		}
	},
}

func init() {
	utilsCmd.AddCommand(regionsCmd)

	regionsCmd.Flags().IntP("line-chunk-size", "", 5000, formatFlagUsage(`Number of lines to process for each thread, and 4 threads is fast enough. Type "kmcp utils merge-regions -h" for details.`))
	regionsCmd.Flags().StringP("out-file", "o", "-", formatFlagUsage(`Out file, supports and recommends a ".gz" suffix ("-" for stdout).`))

	// for single read
	regionsCmd.Flags().Float64P("max-fpr", "f", 0.05, formatFlagUsage(`Maximum false positive rate of a read in search result.`))
	regionsCmd.Flags().Float64P("min-query-cov", "t", 0.55, formatFlagUsage(`Minimum query coverage of a read in search result.`))

	// for merge
	regionsCmd.Flags().IntP("max-gap", "g", 0, formatFlagUsage(`Maximum distance of starting positions of two adjacent regions, 0 for no limitation, 1 for no merging.`))
	regionsCmd.Flags().IntP("min-overlap", "l", 1, formatFlagUsage(`Minimum overlap of two adjacent regions, recommend K-1.`))
	regionsCmd.Flags().BoolP("ignore-type", "I", false, formatFlagUsage("Merge species and assembly-specific regions."))

	regionsCmd.Flags().StringP("regexp", "r", `^(.+)_sliding:(\d+)\-(\d+)$`, formatFlagUsage(`Regular expression for extract reference name and query locations.`))

	regionsCmd.Flags().StringP("name-species", "s", "species-specific", formatFlagUsage(`Name of species-specific regions.`))
	regionsCmd.Flags().StringP("name-assembly", "a", "assembly-specific", formatFlagUsage(`Name of assembly-specific regions.`))
}
