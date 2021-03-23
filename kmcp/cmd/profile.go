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
	"strconv"
	"strings"
	"time"

	"github.com/pkg/errors"
	"github.com/shenwei356/util/cliutil"
	"github.com/spf13/cobra"
	"github.com/twotwotwo/sorts"
	"github.com/zeebo/xxh3"
)

var profileCmd = &cobra.Command{
	Use:   "profile",
	Short: "Generate taxonomic profile from search result",
	Long: `Generate taxonomic profile from search result

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

		minReads := getFlagPositiveInt(cmd, "min-reads")
		minUReads := float64(getFlagPositiveInt(cmd, "min-uniq-reads"))
		minFragsProp := getFlagPositiveFloat64(cmd, "min-frags-prop")

		minDReadsProp := getFlagPositiveFloat64(cmd, "min-dreads-prop")
		maxMismatchErr := getFlagPositiveFloat64(cmd, "max-mismatch-err")

		nameMappingFiles := getFlagStringSlice(cmd, "name-map")

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
		for _, file := range files {
			if isStdin(file) {
				checkError(fmt.Errorf("stdin not supported"))
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
		items := make([]string, numFields)

		profile := make(map[uint64]*Target, 128)

		floatOne := float64(1)

		// ---------------------------------------------------------------
		// stage 1/3
		if opt.Verbose {
			log.Infof("stage 1/3: basic counting for filtering out low-confidence references")
		}

		for _, file := range files {
			infh, r, _, err := inStream(file)
			checkError(err)

			scanner := bufio.NewScanner(infh)

			matches := make(map[uint64][]MatchResult) // target -> match result
			var m MatchResult
			var ms []MatchResult
			var t *Target
			var ok bool
			var hTarget, h uint64
			var prevQuery string
			firtLine := true
			var floatMsSize float64
			var i int
			for scanner.Scan() {
				if firtLine {
					firtLine = false
					continue
				}

				var match MatchResult
				match, ok = parseMatchResult(scanner.Text(), numFields, &items, maxFPR, minQcov)
				if !ok {
					continue
				}

				if prevQuery != match.Query { // new query
					for h, ms = range matches {
						floatMsSize = float64(len(ms))
						for i, m = range ms {
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

							if i == 0 { // count once
								if len(matches) == 1 { // record unique match
									t.UniqMatch[m.FragIdx]++
								}
								t.QLen[m.FragIdx] += float64(m.QLen)
							}

							// for a read matching multiple regions of a reference, distribute count to multiple regions,
							// the sum is still the one.
							t.Match[m.FragIdx] += floatOne / floatMsSize
						}
					}

					matches = make(map[uint64][]MatchResult)
				}

				hTarget = xxh3.HashString(match.Target)
				if _, ok = matches[hTarget]; !ok {
					matches[hTarget] = make([]MatchResult, 0, 1)
				}
				matches[hTarget] = append(matches[hTarget], match)

				prevQuery = match.Query
			}

			for h, ms = range matches {
				floatMsSize = float64(len(ms))
				for i, m = range ms {
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

					if i == 0 { // count once
						if len(matches) == 1 { // record unique match
							t.UniqMatch[m.FragIdx]++
						}
						t.QLen[m.FragIdx] += float64(m.QLen)
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

		var c float64
		var c1 float64
		var c2 float64
		var hs []uint64
		hs = make([]uint64, 0, 8)
		for h, t := range profile {
			for _, c = range t.Match {
				if c > float64(minReads) {
					t.FragsProp++
				}
				t.SumMatch += c
			}
			t.FragsProp = t.FragsProp / float64(len(t.Match))
			if t.FragsProp < minFragsProp {
				hs = append(hs, h)
				continue
			}

			for _, c1 = range t.UniqMatch {
				t.SumUniqMatch += c1
			}
			if t.SumUniqMatch < minUReads {
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

		// ---------------------------------------------------------------
		// stage 2/3, counting ambiguous reads/matches
		if opt.Verbose {
			log.Infof("stage 2/3: counting ambiguous matches for correcting matches")
		}

		// hashA -> hashB -> count
		ambMatch := make(map[uint64]map[uint64]float64, 128)
		// ambMatchDecision := make(map[uint64]map[uint64]uint64, 128)

		for _, file := range files {
			infh, r, _, err := inStream(file)
			checkError(err)

			scanner := bufio.NewScanner(infh)

			matches := make(map[uint64][]MatchResult) // target -> match result
			var ok bool
			var hTarget, h, h1, h2 uint64
			var prevQuery string
			hs := make([]uint64, 0, 256)
			var comb []uint64
			firtLine := true
			for scanner.Scan() {
				if firtLine {
					firtLine = false
					continue
				}

				var match MatchResult
				match, ok = parseMatchResult(scanner.Text(), numFields, &items, maxFPR, minQcov)
				if !ok {
					continue
				}

				if prevQuery != match.Query { // new query
					if len(matches) > 1 { // skip uniq match
						hs = hs[:0]
						for h = range matches {
							_, ok = profile[h]
							if !ok { // skip low-confidence match
								continue
							}
							hs = append(hs, h)
						}

						if len(hs) >= 2 {
							for _, comb = range Combinations(hs, 2) {
								h1, h2 = sortTowUint64s(comb[0], comb[1])
								if _, ok = ambMatch[h1]; !ok {
									ambMatch[h1] = make(map[uint64]float64, 128)
								}
								ambMatch[h1][h2]++
							}
						}
					}
					matches = make(map[uint64][]MatchResult)
				}

				hTarget = xxh3.HashString(match.Target)
				if _, ok = matches[hTarget]; !ok {
					matches[hTarget] = make([]MatchResult, 0, 1)
				}
				matches[hTarget] = append(matches[hTarget], match)
				prevQuery = match.Query
			}

			if len(matches) > 1 {
				hs = hs[:0]
				for h = range matches {
					_, ok = profile[h]
					if !ok { // skip low-confidence match
						continue
					}
					hs = append(hs, h)
				}

				if len(hs) >= 2 {
					for _, comb = range Combinations(hs, 2) {
						h1, h2 = sortTowUint64s(comb[0], comb[1])
						if _, ok = ambMatch[h1]; !ok {
							ambMatch[h1] = make(map[uint64]float64, 128)
						}
						ambMatch[h1][h2]++
					}
				}
			}

			matches = make(map[uint64][]MatchResult)

			checkError(scanner.Err())
			r.Close()
		}

		// ---------------------------------------------------------------
		// stage 3/3
		if opt.Verbose {
			log.Infof("stage 3/3: computing profile")
		}

		profile2 := make(map[uint64]*Target, 128)
		var nReads float64
		var nUnassignedReads float64

		for _, file := range files {
			infh, r, _, err := inStream(file)
			checkError(err)

			scanner := bufio.NewScanner(infh)

			matches := make(map[uint64][]MatchResult) // target -> match result
			var m MatchResult
			var ms []MatchResult
			var t *Target
			var ok bool
			var hTarget, h, h0, h1, h2 uint64
			var prevQuery string
			var floatMsSize float64
			var i int
			var uniqMatch bool
			var first bool
			var sumUReads, prop float64
			hs := make([]uint64, 0, 256)
			firtLine := true
			for scanner.Scan() {
				if firtLine {
					firtLine = false
					continue
				}

				var match MatchResult
				match, ok = parseMatchResult(scanner.Text(), numFields, &items, maxFPR, minQcov)
				if !ok {
					continue
				}

				if prevQuery != match.Query {
					nReads++
					uniqMatch = false
					if len(matches) > 1 { // multiple match
						// skip low-confidence match
						hs = hs[:0] // list to delete
						for h = range matches {
							_, ok = profile[h]
							if ok { // skip low-confidence match
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

							if len(matches) > 1 {
								sumUReads = 0
								for h = range matches {
									sumUReads += profile[h].SumUniqMatch
								}
								for h, ms = range matches {
									floatMsSize = float64(len(ms))
									for i, m = range ms {
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

										if i == 0 { // count once
											t.QLen[m.FragIdx] += float64(m.QLen) * prop
										}

										t.Match[m.FragIdx] += prop / floatMsSize
									}
								}
							} else { // len(matches) == 1
								uniqMatch = true
							}
						} else if len(matches) == 1 {
							uniqMatch = true
						} else {
							nUnassignedReads++
						}
					} else if len(matches) == 1 {
						uniqMatch = true
					}

					if uniqMatch {
						for h, ms = range matches {
							floatMsSize = float64(len(ms))

							for i, m = range ms {
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

								if i == 0 { // count once
									t.UniqMatch[m.FragIdx] += floatOne
									t.QLen[m.FragIdx] += float64(m.QLen)
								}

								t.Match[m.FragIdx] += floatOne / floatMsSize
							}
						}
					}

					matches = make(map[uint64][]MatchResult)
				}

				hTarget = xxh3.HashString(match.Target)
				if _, ok = matches[hTarget]; !ok {
					matches[hTarget] = make([]MatchResult, 0, 1)
				}

				matches[hTarget] = append(matches[hTarget], match)
				prevQuery = match.Query
			}

			uniqMatch = false
			if len(matches) > 1 { // multiple match
				// skip low-confidence match
				hs = hs[:0] // list to delete
				for h = range matches {
					_, ok = profile[h]
					if ok { // skip low-confidence match
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

					if len(matches) > 1 {
						sumUReads = 0
						for h = range matches {
							sumUReads += profile[h].SumUniqMatch
						}

						for h, ms = range matches {
							floatMsSize = float64(len(ms))
							for i, m = range ms {
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
								if i == 0 { // count once
									t.QLen[m.FragIdx] += float64(m.QLen) * prop
								}

								t.Match[m.FragIdx] += prop / floatMsSize
							}
						}
					} else if len(matches) == 1 {
						uniqMatch = true
					} else {
						nUnassignedReads++
					}
				} else if len(matches) == 1 {
					uniqMatch = true
				}
			} else if len(matches) == 1 {
				uniqMatch = true
			}
			if uniqMatch {
				for h, ms = range matches {
					floatMsSize = float64(len(ms))

					for i, m = range ms {
						if t, ok = profile2[h]; !ok {
							t0 := Target{
								Name:       m.Target,
								GenomeSize: m.GSize,
								Match:      make([]float64, m.IdxNum),
								UniqMatch:  make([]float64, m.IdxNum),
								QLen:       make([]float64, m.IdxNum),
							}
							profile2[h2] = &t0
							t = &t0
						}

						if i == 0 { // count once
							t.UniqMatch[m.FragIdx] += floatOne
							t.QLen[m.FragIdx] += float64(m.QLen)
						}

						t.Match[m.FragIdx] += floatOne / floatMsSize
					}
				}
			}
			matches = make(map[uint64][]MatchResult)

			checkError(scanner.Err())
			r.Close()
		}

		// --------------------
		// sum up

		targets := make([]*Target, 0, 128)

		for _, t := range profile2 {
			for _, c = range t.Match {
				if c > float64(minReads) {
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

		log.Infof("#reads: %.0f, #reads of profile: %0.f, proportion: %.6f%%", nReads, (nReads - nUnassignedReads), (nReads-nUnassignedReads)/nReads*100)
	},
}

func init() {
	RootCmd.AddCommand(profileCmd)

	profileCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)

	// for single read
	profileCmd.Flags().Float64P("max-fpr", "f", 0.01, `maximal false positive rate of a read`)
	profileCmd.Flags().Float64P("min-qcov", "t", 0.7, `minimal query coverage of a read`)

	// for ref fragments
	profileCmd.Flags().IntP("min-reads", "r", 50, `minimal number of reads for a reference fragment`)
	profileCmd.Flags().IntP("min-uniq-reads", "u", 10, `minimal number of unique matched reads for a reference fragment`)
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

	m.FPR, err = strconv.ParseFloat((*items)[3], 64)
	if err != nil {
		checkError(fmt.Errorf("failed to parse FPR: %s", (*items)[3]))
	}
	if m.FPR > maxPFR {
		return m, false
	}

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
