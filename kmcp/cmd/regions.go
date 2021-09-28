// Copyright © 2020-2021 Wei Shen <shenwei356@gmail.com>
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
	Use:   "regions",
	Short: "Extract species/assembly-specific regions",
	Long: `Extract species/assembly-specific regions

Input:
  # 1. Simulating reads and searching on one or more databases.
  seqkit sliding -s 10 -W 100 ref.fna.gz \
      | kmcp search -d db1.kmcp -o ref.fna.gz.kmcp@db1.tsv.gz
  seqkit sliding -s 10 -W 100 ref.fna.gz \
      | kmcp search -d db2.kmcp -o ref.fna.gz.kmcp@db2.tsv.gz
  
  # 2. Merging and filtering searching results
  kmcp merge ref.fna.gz.kmcp@*.tsv.gz \
      | kmcp filter -X taxdump -T taxid.map \
	    -o ref.fna.gz.kmcp.uniq.tsv.gz

Performance notes:
  1. Searching results are parsed in parallel, and the number of
     lines proceeded by a thread can be set by the flag --chunk-size.
  2. However using a lot of threads does not always accelerate
     processing, 4 threads with chunk size of 500-5000 is fast enough.

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

		outFile := getFlagString(cmd, "out-prefix")

		maxFPR := getFlagPositiveFloat64(cmd, "max-fpr")
		minQcov := getFlagNonNegativeFloat64(cmd, "min-query-cov")

		chunkSize := getFlagPositiveInt(cmd, "chunk-size")
		if opt.NumCPUs > 4 {
			if opt.Verbose || opt.Log2File {
				log.Infof("using a lot of threads does not always accelerate processing, 4-threads is fast enough")
			}
			opt.NumCPUs = 4
			runtime.GOMAXPROCS(opt.NumCPUs)
		}

		reQueryStr := getFlagString(cmd, "regexp")
		reQuery, err := regexp.Compile(reQueryStr)
		checkError(err)

		nameMultiple := getFlagString(cmd, "name-species")
		nameSingle := getFlagString(cmd, "name-assembly")

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
			log.Infof("  maximal false positive rate: %f", maxFPR)
			log.Infof("  minimal query coverage: %4f", minQcov)
			log.Info()

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

			reader, err := breader.NewBufferedReader(file, opt.NumCPUs, chunkSize, fn)
			checkError(err)
			var data interface{}

			matches = make(map[uint64]*[]*MatchResult3)
			for chunk := range reader.Ch {
				checkError(chunk.Err)

				for _, data = range chunk.Data {
					match = data.(*MatchResult3)

					if prevQuery != match.Query { // new query
						nReads++

						if len(matches) > 0 { // not the first query
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
								if begin <= end0 && name == name0 { // has overlap, update region
									end0 = end
									if name0 == nameMultiple {
										score0 = (score0 + score) / 2
									}
								} else { // print previous region
									nRegions++
									fmt.Fprintf(outfh, "%s\t%d\t%d\t%s\t%.0f\t.\n", ref0, begin0, end0, name0, score0*100)

									ref0, begin0, end0, name0, score0 = ref, begin, end, name, score
								}
							} else { // first record
								ref0, begin0, end0, name0, score0 = ref, begin, end, name, score
							}
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

			nReads++
			if begin <= end0 && name == name0 { // has overlap, update region
				end0 = end
				if name0 == nameMultiple {
					score0 = (score0 + score) / 2
				}
			}
			nRegions++
			fmt.Fprintf(outfh, "%s\t%d\t%d\t%s\t%.0f\t.\n", ref0, begin0, end0, name0, score0*100)
		}

		if opt.Verbose || opt.Log2File {
			log.Info()
			log.Infof("%d regions merged to %d", nReads, nRegions)
		}
	},
}

func init() {
	utilsCmd.AddCommand(regionsCmd)

	regionsCmd.Flags().IntP("chunk-size", "", 5000, `number of lines to process for each thread, and 4 threads is fast enough. Type "kmcp profile -h" for details`)
	regionsCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)

	// for single read
	regionsCmd.Flags().Float64P("max-fpr", "f", 0.05, `maximal false positive rate of a read in search result`)
	regionsCmd.Flags().Float64P("min-query-cov", "t", 0.55, `minimal query coverage of a read in search result`)

	regionsCmd.Flags().StringP("regexp", "r", `^(.+)_sliding:(\d+)\-(\d+)$`, `regular expression for extract reference name and query locations`)
	regionsCmd.Flags().StringP("name-species", "s", "species", `name of species-specific regions`)
	regionsCmd.Flags().StringP("name-assembly", "a", "assembly", `name of assembly-specific regions`)

}