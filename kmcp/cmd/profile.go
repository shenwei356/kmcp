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
		minQcov := getFlagNonNegativeFloat64(cmd, "max-qcov")

		minReads := getFlagPositiveInt(cmd, "min-reads")
		minUReads := getFlagPositiveInt(cmd, "min-uniq-reads")
		minFragsProp := getFlagPositiveFloat64(cmd, "min-frags-prop")

		nameMappingFiles := getFlagStringSlice(cmd, "name-map")

		if opt.Verbose {
			log.Info("checking input files ...")
		}
		files := getFileListFromArgsAndFile(cmd, args, true, "infile-list", true)
		if opt.Verbose {
			if len(files) == 1 && isStdin(files[0]) {
				log.Info("no files given, reading from stdin")
			} else {
				log.Infof("%d input file(s) given", len(files))
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
			outfh.WriteString(fmt.Sprint("target\tfragsProp\tabundance\treads\tureads\tannotation\n"))
		} else {
			outfh.WriteString(fmt.Sprint("target\tfragsProp\tabundance\treads\tureads\n"))
		}

		numFields := 12
		items := make([]string, numFields)

		profile := make(map[uint64]*Target, 128)

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
			var floatOne, floatMsSize float64
			floatOne = float64(1)
			for scanner.Scan() {
				if firtLine {
					firtLine = false
					continue
				}

				match, ok := parseMatchResult(scanner.Text(), numFields, &items, maxFPR, minQcov)
				if !ok {
					continue
				}

				if prevQuery != match.Query { // new query
					for h, ms = range matches {
						floatMsSize = float64(len(ms))
						for _, m = range ms {
							if t, ok = profile[h]; !ok {
								t0 := Target{
									Name:      m.Target,
									Match:     make([]float64, m.IdxNum),
									UniqMatch: make([]int, m.IdxNum),
									QLen:      make([]uint64, m.IdxNum),
								}
								profile[h] = &t0
								t = &t0
							}

							t.Name = m.Target
							t.GenomeSize = m.GSize

							// for a read matching multiple regions of a reference, distribute count to multiple regions,
							// the sum is still the one.
							t.Match[m.FragIdx] += floatOne / floatMsSize

							// record unique match
							if len(matches) == 1 {
								t.UniqMatch[m.FragIdx] += 1
							}
							t.QLen[m.FragIdx] += uint64(m.QLen)
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
				for _, m = range ms {
					if t, ok = profile[h]; !ok {
						t0 := Target{
							Name:      m.Target,
							Match:     make([]float64, m.IdxNum),
							UniqMatch: make([]int, m.IdxNum),
							QLen:      make([]uint64, m.IdxNum),
						}
						profile[h] = &t0
						t = &t0
					}

					t.Name = m.Target
					t.GenomeSize = m.GSize

					// for a read matching multiple regions of a reference, distribute count to multiple regions,
					// the sum is still the one.
					t.Match[m.FragIdx] += floatOne / floatMsSize

					// record unique match
					if len(matches) == 1 {
						t.UniqMatch[m.FragIdx] += 1
					}
					t.QLen[m.FragIdx] += uint64(m.QLen)
				}
			}

			matches = make(map[uint64][]MatchResult)

			checkError(scanner.Err())
			r.Close()
		}

		targets := make([]*Target, 0, 128)
		var c float64
		var c1 int
		var c2 uint64
		for h, t := range profile {
			for _, c = range t.Match {
				if c > float64(minReads) {
					t.FragsProp++
				}
				t.SumMatch += c
			}
			t.FragsProp = t.FragsProp / float64(len(t.Match))
			if t.FragsProp < minFragsProp {
				delete(profile, h)
				continue
			}

			for _, c1 = range t.UniqMatch {
				t.SumUniqMatch += c1
			}
			if t.SumUniqMatch < minUReads {
				delete(profile, h)
				continue
			}

			for _, c2 = range t.QLen {
				t.Qlens += c2
			}

			t.Coverage = float64(t.Qlens) / float64(t.GenomeSize)

			targets = append(targets, t)

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
				outfh.WriteString(fmt.Sprintf("%s\t%.6f\t%.2f\t%.0f\t%d\t%s\n",
					t.Name, t.Coverage/totalCoverage*100, t.FragsProp, t.SumMatch, t.SumUniqMatch, name2))
			} else {
				outfh.WriteString(fmt.Sprintf("%s\t%0.6f\t%.2f\t%.0f\t%d\n",
					t.Name, t.Coverage/totalCoverage*100, t.FragsProp, t.SumMatch, t.SumUniqMatch))
			}
		}
	},
}

func init() {
	RootCmd.AddCommand(profileCmd)

	profileCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)

	// for single read
	profileCmd.Flags().Float64P("max-fpr", "f", 0.01, `maximum false positive rate of a read`)
	profileCmd.Flags().Float64P("max-qcov", "t", 0.7, `maximum query coverage of a read`)

	// for ref fragments
	profileCmd.Flags().IntP("min-reads", "r", 50, `minimum number of reads for a reference fragment`)
	profileCmd.Flags().IntP("min-uniq-reads", "u", 10, `minimum number of unique matched reads for a reference fragment`)
	profileCmd.Flags().Float64P("min-frags-prop", "p", 0.3, `minimum proportion of matched fragments`)

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
	QLen []uint64

	// unique match
	UniqMatch []int

	SumMatch     float64
	SumUniqMatch int
	FragsProp    float64
	Coverage     float64
	Qlens        uint64
}

func (t Target) String() string {
	var buf bytes.Buffer
	buf.WriteString(t.Name)
	for i := range t.Match {
		// buf.WriteString(fmt.Sprintf(", %d: %.1f(%d)-%d", i, t.Match[i], t.UniqMatch[i], t.QLen[i]))
		buf.WriteString(fmt.Sprintf(", %d: %.0f(%d)", i, t.Match[i], t.UniqMatch[i]))
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
