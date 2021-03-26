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
	"bufio"
	"bytes"
	"fmt"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/pkg/errors"
	"github.com/shenwei356/breader"
	"github.com/shenwei356/util/cliutil"
	"github.com/spf13/cobra"
	"github.com/twotwotwo/sorts"
	"github.com/zeebo/xxh3"
)

var profileCmd = &cobra.Command{
	Use:   "profile",
	Short: "Generate taxonomic profile from search result",
	Long: `Generate taxonomic profile from search result

Performance Note:
  1. This command parses searching result in parallel, and the number of
     lines proceeded by a thread can be set by the flag --chunk-size).
  2. However using a lot of threads does not always accelerate processing,
     4 threads with chunk size 0f 500-5000 is fast enough.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)

		timeStart := time.Now()
		defer func() {
			if opt.Verbose {
				log.Info()
				log.Infof("elapsed time: %s", time.Since(timeStart))
			}
		}()

		var err error

		outFile := getFlagString(cmd, "out-prefix")

		maxFPR := getFlagPositiveFloat64(cmd, "max-fpr")
		minQcov := getFlagNonNegativeFloat64(cmd, "min-qcov")

		minReads := float64(getFlagPositiveInt(cmd, "min-reads"))
		minUReads := float64(getFlagPositiveInt(cmd, "min-uniq-reads"))
		minFragsProp := getFlagPositiveFloat64(cmd, "min-frags-prop")

		minDReadsProp := getFlagPositiveFloat64(cmd, "min-dreads-prop")
		maxMismatchErr := getFlagPositiveFloat64(cmd, "max-mismatch-err")

		nameMappingFiles := getFlagStringSlice(cmd, "name-map")

		chunkSize := getFlagPositiveInt(cmd, "chunk-size")
		if opt.NumCPUs > 4 {
			log.Infof("using a lot of threads does not always accelerate processing, 4-threads is fast enough")
			opt.NumCPUs = 4
			runtime.GOMAXPROCS(opt.NumCPUs)
		}

		// ---------------------------------------------------------------

		if opt.Verbose {
			log.Info("checking input files ...")
		}
		files := getFileListFromArgsAndFile(cmd, args, true, "infile-list", true)
		if opt.Verbose {
			if len(files) == 1 && isStdin(files[0]) {
				// log.Info("no files given, reading from stdin")
				checkError(fmt.Errorf("stdin not supported"))
			} else {
				log.Infof("%d input file(s) given", len(files))
			}
		}

		outFileClean := filepath.Clean(outFile)
		for _, file := range files {
			if isStdin(file) {
				checkError(fmt.Errorf("stdin not supported"))
			} else if filepath.Clean(file) == outFileClean {
				checkError(fmt.Errorf("out file should not be one of the input file"))
			}
		}

		// ---------------------------------------------------------------
		// name mapping files

		var namesMap map[string]string
		mappingNames := len(nameMappingFiles) != 0
		if mappingNames {
			if opt.Verbose {
				log.Infof("loading name mapping file ...")
			}
			var nameMappingFile string
			nameMappingFile = nameMappingFiles[0]
			namesMap, err = cliutil.ReadKVs(nameMappingFile, false)
			if err != nil {
				checkError(errors.Wrap(err, nameMappingFile))
			}

			if len(nameMappingFiles) > 1 {
				for _, _nameMappingFile := range nameMappingFiles[1:] {
					_namesMap, err := cliutil.ReadKVs(_nameMappingFile, false)
					if err != nil {
						checkError(errors.Wrap(err, nameMappingFile))
					}
					for _k, _v := range _namesMap {
						namesMap[_k] = _v
					}
				}
			}

			if opt.Verbose {
				log.Infof("%d pairs of name mapping values from %d file(s) loaded", len(namesMap), len(nameMappingFiles))
			}

			mappingNames = len(namesMap) > 0
		}

		// ---------------------------------------------------------------

		numFields := 12

		profile := make(map[uint64]*Target, 128)

		floatOne := float64(1)

		// ---------------------------------------------------------------

		pool := &sync.Pool{New: func() interface{} {
			return make([]string, numFields)
		}}

		fn := func(line string) (interface{}, bool, error) {
			if line == "" || line[0] == '#' { // ignoring blank line and comment line
				return "", false, nil
			}

			items := pool.Get().([]string)

			match, ok := parseMatchResult(line, numFields, &items, maxFPR, minQcov)
			if !ok {
				pool.Put(items)
				return nil, false, nil
			}

			pool.Put(items)
			return match, true, nil
		}

		// ---------------------------------------------------------------
		// stage 1/3
		if opt.Verbose {
			log.Infof("stage 1/3: counting matches and unique matches for filtering out low-confidence references")
		}
		timeStart1 := time.Now()

		for _, file := range files {
			infh, r, _, err := inStream(file)
			checkError(err)

			scanner := bufio.NewScanner(infh)

			matches := make(map[uint64]*[]MatchResult) // target -> match result
			var m MatchResult
			var ms *[]MatchResult
			var t *Target
			var ok bool
			var hTarget, h uint64
			var prevQuery string
			var floatMsSize float64
			var match MatchResult
			var first bool

			reader, err := breader.NewBufferedReader(file, opt.NumCPUs, chunkSize, fn)
			checkError(err)
			var data interface{}

			for chunk := range reader.Ch {
				checkError(chunk.Err)

				for _, data = range chunk.Data {
					match = data.(MatchResult)

					if prevQuery != match.Query { // new query
						for h, ms = range matches {
							floatMsSize = float64(len(*ms))
							first = true
							for _, m = range *ms {
								if t, ok = profile[h]; !ok {
									t0 := Target{
										Name:       m.Target,
										GenomeSize: m.GSize,
										Match:      make([]float64, m.IdxNum),
										UniqMatch:  make([]float64, m.IdxNum),
										QLen:       make([]float64, m.IdxNum),
									}
									profile[h] = &t0
									t = &t0
								}

								if first { // count once
									if len(matches) == 1 { // record unique match
										t.UniqMatch[m.FragIdx]++
									}
									t.QLen[m.FragIdx] += float64(m.QLen)
									first = false
								}

								// for a read matching multiple regions of a reference, distribute count to multiple regions,
								// the sum is still the one.
								t.Match[m.FragIdx] += floatOne / floatMsSize
							}
						}

						matches = make(map[uint64]*[]MatchResult)
					}

					hTarget = xxh3.HashString(match.Target)
					if ms, ok = matches[hTarget]; !ok {
						tmp := []MatchResult{match}
						matches[hTarget] = &tmp
					} else {
						*ms = append(*ms, match)
					}

					prevQuery = match.Query
				}
			}

			for h, ms = range matches {
				floatMsSize = float64(len(*ms))
				first = true
				for _, m = range *ms {
					if t, ok = profile[h]; !ok {
						t0 := Target{
							Name:       m.Target,
							GenomeSize: m.GSize,
							Match:      make([]float64, m.IdxNum),
							UniqMatch:  make([]float64, m.IdxNum),
							QLen:       make([]float64, m.IdxNum),
						}
						profile[h] = &t0
						t = &t0
					}

					if first { // count once
						if len(matches) == 1 { // record unique match
							t.UniqMatch[m.FragIdx]++
						}
						t.QLen[m.FragIdx] += float64(m.QLen)
						first = false
					}

					// for a read matching multiple regions of a reference, distribute count to multiple regions,
					// the sum is still the one.
					t.Match[m.FragIdx] += floatOne / floatMsSize
				}
			}

			checkError(scanner.Err())
			r.Close()
		}

		// --------------------
		// sum up

		log.Infof("  number of references in search result: %d", len(profile))

		var c float64
		var c1 float64
		var c2 float64
		var hs []uint64
		hs = make([]uint64, 0, 8) // list to delete
		for h, t := range profile {
			for _, c = range t.Match {
				if c > minReads {
					t.FragsProp++
				}
				t.SumMatch += c
			}
			t.FragsProp = t.FragsProp / float64(len(t.Match))
			if t.FragsProp < minFragsProp { // low coverage
				hs = append(hs, h)
				continue
			}

			for _, c1 = range t.UniqMatch {
				t.SumUniqMatch += c1
			}
			if t.SumUniqMatch < minUReads { // no enough unique match
				hs = append(hs, h)
				continue
			}

			for _, c2 = range t.QLen {
				t.Qlens += c2
			}

			t.Coverage = float64(t.Qlens) / float64(t.GenomeSize)
		}
		for _, h := range hs {
			delete(profile, h)
		}

		log.Infof("  number of estimated references: %d", len(profile))
		log.Infof("  elapsed time: %s", time.Since(timeStart1))

		// ---------------------------------------------------------------
		// stage 2/3, counting ambiguous reads/matches
		if opt.Verbose {
			log.Infof("stage 2/3: counting ambiguous matches for correcting matches")
		}
		timeStart1 = time.Now()

		// hashA -> hashB -> count
		ambMatch := make(map[uint64]map[uint64]float64, 128)

		for _, file := range files {
			infh, r, _, err := inStream(file)
			checkError(err)

			scanner := bufio.NewScanner(infh)

			matches := make(map[uint64]*[]MatchResult) // target -> match result
			var ok bool
			var ms *[]MatchResult
			var hTarget, h, h1, h2 uint64
			var prevQuery string
			hs := make([]uint64, 0, 256)
			var match MatchResult
			var amb map[uint64]float64
			var i, j int
			var n, np1 int

			reader, err := breader.NewBufferedReader(file, opt.NumCPUs, chunkSize, fn)
			checkError(err)
			var data interface{}

			for chunk := range reader.Ch {
				checkError(chunk.Err)

				for _, data = range chunk.Data {
					match = data.(MatchResult)

					if prevQuery != match.Query { // new query
						if len(matches) > 1 { // skip uniq match
							hs = hs[:0]
							for h = range matches {
								if _, ok = profile[h]; !ok { // skip low-confidence match
									continue
								}
								hs = append(hs, h)
							}

							if len(hs) >= 2 {
								sorts.Quicksort(Uint64Slice(hs))

								n = len(hs)
								np1 = len(hs) - 1
								for i = 0; i < np1; i++ {
									for j = i + 1; j < n; j++ {
										h1, h2 = hs[i], hs[j]
										if amb, ok = ambMatch[h1]; !ok {
											ambMatch[h1] = make(map[uint64]float64, 128)
											ambMatch[h1][h2]++
										} else {
											amb[h2]++
										}
									}
								}
							}
						}
						matches = make(map[uint64]*[]MatchResult)
					}

					hTarget = xxh3.HashString(match.Target)
					if ms, ok = matches[hTarget]; !ok {
						tmp := []MatchResult{match}
						matches[hTarget] = &tmp
					} else {
						*ms = append(*ms, match)
					}
					prevQuery = match.Query
				}
			}

			if len(matches) > 1 {
				hs = hs[:0]
				for h = range matches {
					if _, ok = profile[h]; !ok { // skip low-confidence match
						continue
					}
					hs = append(hs, h)
				}

				if len(hs) >= 2 {
					sorts.Quicksort(Uint64Slice(hs))

					n = len(hs)
					np1 = len(hs) - 1
					for i = 0; i < np1; i++ {
						for j = i + 1; j < n; j++ {
							h1, h2 = hs[i], hs[j]
							if amb, ok = ambMatch[h1]; !ok {
								ambMatch[h1] = make(map[uint64]float64, 128)
								ambMatch[h1][h2]++
							} else {
								amb[h2]++
							}
						}
					}
				}
			}

			matches = make(map[uint64]*[]MatchResult)

			checkError(scanner.Err())
			r.Close()
		}

		log.Infof("  elapsed time: %s", time.Since(timeStart1))

		// ---------------------------------------------------------------
		// stage 3/3
		if opt.Verbose {
			log.Infof("stage 3/3: correcting matches and computing profile")
		}
		timeStart1 = time.Now()

		profile2 := make(map[uint64]*Target, 128)
		var nReads float64
		var nUnassignedReads float64

		for _, file := range files {
			infh, r, _, err := inStream(file)
			checkError(err)

			scanner := bufio.NewScanner(infh)

			matches := make(map[uint64]*[]MatchResult) // target -> match result
			var m MatchResult
			var ms *[]MatchResult
			var t *Target
			var ok bool
			var hTarget, h, h0, h1, h2 uint64
			var prevQuery string
			var floatMsSize float64
			var uniqMatch bool
			var first bool
			var sumUReads, prop float64
			hs := make([]uint64, 0, 256)
			var match MatchResult

			reader, err := breader.NewBufferedReader(file, opt.NumCPUs, chunkSize, fn)
			checkError(err)
			var data interface{}

			for chunk := range reader.Ch {
				checkError(chunk.Err)

				for _, data = range chunk.Data {
					match = data.(MatchResult)

					if prevQuery != match.Query {
						nReads++
						uniqMatch = false
						if len(matches) > 1 { // multiple match
							// skip low-confidence match
							hs = hs[:0] // list to delete
							for h = range matches {
								if _, ok = profile[h]; ok { // skip low-confidence match
									continue
								}
								hs = append(hs, h)
							}
							for _, h = range hs {
								delete(matches, h)
							}

							if len(matches) > 1 {
								hs = hs[:0] // list to delete
								first = true
								for h = range matches {
									if first {
										h0 = h // previous one
										first = false
										continue
									}
									h1, h2 = sortTowUint64s(h0, h)

									// decide which to keep
									if profile[h0].SumMatch*(1-minDReadsProp) >= ambMatch[h1][h2] &&
										profile[h].SumUniqMatch < profile[h0].SumUniqMatch*maxMismatchErr {
										hs = append(hs, h)
									} else {
										h0 = h
									}
								}
								for _, h = range hs {
									delete(matches, h)
								}

								if len(matches) > 1 { // redistribute matches
									sumUReads = 0
									for h = range matches {
										sumUReads += profile[h].SumUniqMatch
									}
									for h, ms = range matches {
										floatMsSize = float64(len(*ms))
										first = true
										for _, m = range *ms {
											if t, ok = profile2[h]; !ok {
												t0 := Target{
													Name:       m.Target,
													GenomeSize: m.GSize,
													Match:      make([]float64, m.IdxNum),
													UniqMatch:  make([]float64, m.IdxNum),
													QLen:       make([]float64, m.IdxNum),
												}
												profile2[h] = &t0
												t = &t0
											}

											prop = profile[h].SumUniqMatch / sumUReads

											if first { // count once
												t.QLen[m.FragIdx] += float64(m.QLen) * prop
												first = false
											}

											t.Match[m.FragIdx] += prop / floatMsSize
										}
									}
								} else { // len(matches) == 1
									uniqMatch = true
								}
							} else if len(matches) == 1 {
								uniqMatch = true
							} else { // 0
								nUnassignedReads++
							}
						} else if len(matches) == 1 {
							uniqMatch = true
						}

						if uniqMatch {
							for h, ms = range matches {
								floatMsSize = float64(len(*ms))
								first = true
								for _, m = range *ms {
									if t, ok = profile2[h]; !ok {
										t0 := Target{
											Name:       m.Target,
											GenomeSize: m.GSize,
											Match:      make([]float64, m.IdxNum),
											UniqMatch:  make([]float64, m.IdxNum),
											QLen:       make([]float64, m.IdxNum),
										}
										profile2[h] = &t0
										t = &t0
									}

									if first { // count once
										t.UniqMatch[m.FragIdx] += floatOne
										t.QLen[m.FragIdx] += float64(m.QLen)
										first = false
									}

									t.Match[m.FragIdx] += floatOne / floatMsSize
								}
							}
						}

						matches = make(map[uint64]*[]MatchResult)
					}

					hTarget = xxh3.HashString(match.Target)
					if ms, ok = matches[hTarget]; !ok {
						tmp := []MatchResult{match}
						matches[hTarget] = &tmp
					} else {
						*ms = append(*ms, match)
					}

					prevQuery = match.Query
				}
			}

			nReads++
			uniqMatch = false
			if len(matches) > 1 { // multiple match
				// skip low-confidence match
				hs = hs[:0] // list to delete
				for h = range matches {
					if _, ok = profile[h]; ok { // skip low-confidence match
						continue
					}
					hs = append(hs, h)
				}
				for _, h = range hs {
					delete(matches, h)
				}

				if len(matches) > 1 {
					hs = hs[:0] // list to delete
					first = true
					for h = range matches {
						if first {
							h0 = h // previous one
							first = false
							continue
						}
						h1, h2 = sortTowUint64s(h0, h)

						// decide which to keep
						if profile[h0].SumMatch*(1-minDReadsProp) >= ambMatch[h1][h2] &&
							profile[h].SumUniqMatch < profile[h0].SumUniqMatch*maxMismatchErr {
							hs = append(hs, h)
						} else {
							h0 = h
						}
					}
					for _, h = range hs {
						delete(matches, h)
					}

					if len(matches) > 1 { // redistribute matches
						sumUReads = 0
						for h = range matches {
							sumUReads += profile[h].SumUniqMatch
						}
						for h, ms = range matches {
							floatMsSize = float64(len(*ms))
							first = true
							for _, m = range *ms {
								if t, ok = profile2[h]; !ok {
									t0 := Target{
										Name:       m.Target,
										GenomeSize: m.GSize,
										Match:      make([]float64, m.IdxNum),
										UniqMatch:  make([]float64, m.IdxNum),
										QLen:       make([]float64, m.IdxNum),
									}
									profile2[h] = &t0
									t = &t0
								}

								prop = profile[h].SumUniqMatch / sumUReads

								if first { // count once
									t.QLen[m.FragIdx] += float64(m.QLen) * prop
									first = false
								}

								t.Match[m.FragIdx] += prop / floatMsSize
							}
						}
					} else { // len(matches) == 1
						uniqMatch = true
					}
				} else if len(matches) == 1 {
					uniqMatch = true
				} else { // 0
					nUnassignedReads++
				}
			} else if len(matches) == 1 {
				uniqMatch = true
			}

			if uniqMatch {
				for h, ms = range matches {
					floatMsSize = float64(len(*ms))
					first = true
					for _, m = range *ms {
						if t, ok = profile2[h]; !ok {
							t0 := Target{
								Name:       m.Target,
								GenomeSize: m.GSize,
								Match:      make([]float64, m.IdxNum),
								UniqMatch:  make([]float64, m.IdxNum),
								QLen:       make([]float64, m.IdxNum),
							}
							profile2[h] = &t0
							t = &t0
						}

						if first { // count once
							t.UniqMatch[m.FragIdx] += floatOne
							t.QLen[m.FragIdx] += float64(m.QLen)
							first = false
						}

						t.Match[m.FragIdx] += floatOne / floatMsSize
					}
				}
			}

			checkError(scanner.Err())
			r.Close()
		}

		// --------------------
		// sum up

		targets := make([]*Target, 0, 128)

		for _, t := range profile2 {
			for _, c = range t.Match {
				if c > minReads {
					t.FragsProp++
				}
				t.SumMatch += c
			}
			t.FragsProp = t.FragsProp / float64(len(t.Match))
			if t.FragsProp < minFragsProp {
				continue
			}

			for _, c1 = range t.UniqMatch {
				t.SumUniqMatch += c1
			}
			if t.SumUniqMatch < minUReads {
				continue
			}

			for _, c2 = range t.QLen {
				t.Qlens += c2
			}

			t.Coverage = float64(t.Qlens) / float64(t.GenomeSize)

			targets = append(targets, t)
		}

		log.Infof("  number of estimated references: %d", len(targets))
		log.Infof("  elapsed time: %s", time.Since(timeStart1))

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

		if mappingNames {
			outfh.WriteString(fmt.Sprint("target\tabundance\tfragsProp\treads\tureads\tannotation\n"))
		} else {
			outfh.WriteString(fmt.Sprint("target\tabundance\tfragsProp\treads\tureads\n"))
		}

		sorts.Quicksort(Targets(targets))

		var name2 string
		var totalCoverage float64

		for _, t := range targets {
			totalCoverage += t.Coverage
		}
		for _, t := range targets {
			if mappingNames {
				name2 = namesMap[t.Name]
				outfh.WriteString(fmt.Sprintf("%s\t%.6f\t%.2f\t%.0f\t%.0f\t%s\n",
					t.Name, t.Coverage/totalCoverage*100, t.FragsProp, t.SumMatch, t.SumUniqMatch, name2))
			} else {
				outfh.WriteString(fmt.Sprintf("%s\t%0.6f\t%.2f\t%.0f\t%.0f\n",
					t.Name, t.Coverage/totalCoverage*100, t.FragsProp, t.SumMatch, t.SumUniqMatch))
			}
		}

		log.Infof("#input reads: %.0f, #reads belonging to references in profile: %0.f, proportion: %.6f%%", nReads, (nReads - nUnassignedReads), (nReads-nUnassignedReads)/nReads*100)
	},
}

func init() {
	RootCmd.AddCommand(profileCmd)

	profileCmd.Flags().IntP("chunk-size", "", 5000, `number of lines to process for each thread, and 4 threads is fast enough. Type "kmcp profile -h" for details`)

	profileCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)

	// for single read
	profileCmd.Flags().Float64P("max-fpr", "f", 0.01, `maximal false positive rate of a read`)
	profileCmd.Flags().Float64P("min-qcov", "t", 0.7, `minimal query coverage of a read`)

	// for ref fragments
	profileCmd.Flags().IntP("min-reads", "r", 50, `minimal number of reads for a reference fragment`)
	profileCmd.Flags().IntP("min-uniq-reads", "u", 5, `minimal number of unique matched reads for a reference fragment`)
	profileCmd.Flags().Float64P("min-frags-prop", "p", 0.3, `minimal proportion of matched fragments`)

	// for the two-stage taxonomy assignment algorithm in MagaPath
	profileCmd.Flags().Float64P("min-dreads-prop", "D", 0.05, `minimal proportion of distinct reads, for determing the right reference for ambigous reads`)
	profileCmd.Flags().Float64P("max-mismatch-err", "R", 0.05, `maximal error rate of a read being matched to a wrong reference, for determing the right reference for ambigous reads`)

	// name mapping
	profileCmd.Flags().StringSliceP("name-map", "M", []string{}, `tabular two-column file(s) mapping names to user-defined values`)
}

type MatchResult struct {
	Query   string
	QLen    int
	QKmers  int
	FPR     float64
	Hits    int
	Target  string
	FragIdx int
	IdxNum  int
	GSize   uint64
	MKmers  int
	QCov    float64
}

func parseMatchResult(line string, numFields int, items *[]string, maxPFR float64, minQcov float64) (MatchResult, bool) {
	stringSplitN(line, "\t", numFields, items)
	if len(*items) < numFields {
		checkError(fmt.Errorf("invalid kmcp search result format"))
	}

	var m MatchResult

	var err error

	// too slow
	m.FPR, err = strconv.ParseFloat((*items)[3], 64)
	if err != nil {
		checkError(fmt.Errorf("failed to parse FPR: %s", (*items)[3]))
	}
	if m.FPR > maxPFR {
		return m, false
	}

	// too slow
	m.QCov, err = strconv.ParseFloat((*items)[10], 64)
	if err != nil {
		checkError(fmt.Errorf("failed to parse qCov: %s", (*items)[10]))
	}
	if m.QCov < minQcov {
		return m, false
	}

	// -----------

	m.Query = (*items)[0]

	m.QLen, err = strconv.Atoi((*items)[1])
	if err != nil {
		checkError(fmt.Errorf("failed to parse qLen: %s", (*items)[1]))
	}

	m.QKmers, err = strconv.Atoi((*items)[2])
	if err != nil {
		checkError(fmt.Errorf("failed to parse qKmers: %s", (*items)[2]))
	}

	m.Hits, err = strconv.Atoi((*items)[4])
	if err != nil {
		checkError(fmt.Errorf("failed to parse hits: %s", (*items)[4]))
	}

	m.Target = (*items)[5]

	m.FragIdx, err = strconv.Atoi((*items)[6])
	if err != nil {
		checkError(fmt.Errorf("failed to parse fragIdx: %s", (*items)[6]))
	}

	m.IdxNum, err = strconv.Atoi((*items)[7])
	if err != nil {
		checkError(fmt.Errorf("failed to parse IdxNum: %s", (*items)[7]))
	}

	m.GSize, err = strconv.ParseUint((*items)[8], 10, 64)
	if err != nil {
		checkError(fmt.Errorf("failed to parse genomeSize: %s", (*items)[8]))
	}

	m.MKmers, err = strconv.Atoi((*items)[9])
	if err != nil {
		checkError(fmt.Errorf("failed to parse mKmers: %s", (*items)[9]))
	}

	return m, true
}

type Target struct {
	Name string

	GenomeSize uint64

	// Counting matches in all frags
	// some reads match multiple sites in the same genome,
	// the count should be divided by number of sites.
	Match []float64

	// sum of read (query) length
	QLen []float64

	// unique match
	UniqMatch []float64

	SumMatch     float64
	SumUniqMatch float64
	FragsProp    float64
	Coverage     float64
	Qlens        float64
}

func (t Target) String() string {
	var buf bytes.Buffer
	buf.WriteString(t.Name)
	for i := range t.Match {
		buf.WriteString(fmt.Sprintf(", %d: %.0f(%.0f)", i, t.Match[i], t.UniqMatch[i]))
	}
	buf.WriteString("\n")
	return buf.String()
}

type Targets []*Target

func (t Targets) Len() int { return len(t) }
func (t Targets) Less(i, j int) bool {
	return t[i].Coverage > t[j].Coverage || t[i].FragsProp > t[j].FragsProp
}
func (t Targets) Swap(i, j int) {
	t[i], t[j] = t[j], t[i]
}
