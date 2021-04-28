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
	"fmt"
	"path/filepath"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/pkg/errors"
	"github.com/shenwei356/breader"
	"github.com/shenwei356/unikmer"
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
     lines proceeded by a thread can be set by the flag --chunk-size.
  2. However using a lot of threads does not always accelerate processing,
     4 threads with chunk size 0f 500-5000 is fast enough.

Profiling output format
  1. kmcp
  2. CAMI         (-M/--metaphlan-report)
  3. MetaPhlAn v2 (-C/--cami-report)

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
		topNScore := getFlagNonNegativeInt(cmd, "keep-top-scores")
		keepFullMatch := getFlagBool(cmd, "keep-full-match")

		minReads := float64(getFlagPositiveInt(cmd, "min-reads"))
		minUReads := float64(getFlagPositiveInt(cmd, "min-uniq-reads"))
		minFragsProp := getFlagPositiveFloat64(cmd, "min-frags-prop")

		minDReadsProp := getFlagPositiveFloat64(cmd, "min-dreads-prop")
		maxMismatchErr := getFlagPositiveFloat64(cmd, "max-mismatch-err")

		nameMappingFiles := getFlagStringSlice(cmd, "name-map")

		taxidMappingFiles := getFlagStringSlice(cmd, "taxid-map")
		taxonomyDataDir := getFlagString(cmd, "taxonomy-dir")
		if len(taxidMappingFiles) > 0 && taxonomyDataDir == "" {
			checkError(fmt.Errorf("flag -X/--taxonomy-dir is needed when -N/--taxid-map given"))
		}
		if len(taxidMappingFiles) == 0 && taxonomyDataDir != "" {
			checkError(fmt.Errorf("flag -N/--taxid-map is needed when -X/--taxonomy-dir given"))
		}
		separator := getFlagString(cmd, "separator")
		if separator == "" {
			log.Warningf("value of -s/--separator better not be empty")
		}

		chunkSize := getFlagPositiveInt(cmd, "chunk-size")
		if opt.NumCPUs > 4 {
			log.Infof("using a lot of threads does not always accelerate processing, 4-threads is fast enough")
			opt.NumCPUs = 4
			runtime.GOMAXPROCS(opt.NumCPUs)
		}

		camiReportFile := getFlagString(cmd, "cami-report")
		outputCamiReport := camiReportFile != ""
		if outputCamiReport && !strings.HasSuffix(camiReportFile, ".profile") {
			camiReportFile = camiReportFile + ".profile"
		}

		metaphlanReportFile := getFlagString(cmd, "metaphlan-report")
		outputMetaphlanReport := metaphlanReportFile != ""
		if outputMetaphlanReport && !strings.HasSuffix(metaphlanReportFile, ".profile") {
			metaphlanReportFile = metaphlanReportFile + ".profile"
		}

		showRanks := getFlagStringSlice(cmd, "show-rank")
		rankPrefixes := getFlagStringSlice(cmd, "rank-prefix")
		if outputMetaphlanReport && len(showRanks) != len(rankPrefixes) {
			checkError(fmt.Errorf("number of ranks to show and ther prefixes should match"))
		}
		rankOrder := make(map[string]int, len(showRanks))
		for _i, _r := range showRanks {
			rankOrder[_r] = _i
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
				log.Infof("  %d input file(s) given", len(files))
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
				log.Infof("  %d pairs of name mapping values from %d file(s) loaded", len(namesMap), len(nameMappingFiles))
			}

			mappingNames = len(namesMap) > 0
		}

		// ---------------------------------------------------------------
		// taxid mapping files

		mappingTaxids := len(taxidMappingFiles) != 0
		var taxdb *unikmer.Taxonomy
		var taxidMap map[string]uint32

		if mappingTaxids {
			if opt.Verbose {
				log.Infof("loading TaxId mapping file ...")
			}
			var taxidMappingFile string
			taxidMappingFile = taxidMappingFiles[0]
			taxidMapStr, err := cliutil.ReadKVs(taxidMappingFile, false)
			if err != nil {
				checkError(errors.Wrap(err, taxidMappingFile))
			}
			taxidMap = make(map[string]uint32, len(taxidMapStr))
			var taxid uint64
			for k, s := range taxidMapStr {
				taxid, err = strconv.ParseUint(s, 10, 32)
				if err != nil {
					checkError(fmt.Errorf("invalid TaxId: %s", s))
				}
				taxidMap[k] = uint32(taxid)
			}

			if len(taxidMappingFiles) > 1 {
				for _, taxidMappingFile := range taxidMappingFiles[1:] {
					_taxidMapStr, err := cliutil.ReadKVs(taxidMappingFile, false)
					if err != nil {
						checkError(errors.Wrap(err, taxidMappingFile))
					}
					for _k, _v := range _taxidMapStr {
						taxid, err = strconv.ParseUint(_v, 10, 32)
						if err != nil {
							checkError(fmt.Errorf("invalid TaxId: %s", _v))
						}
						taxidMap[_k] = uint32(taxid)
					}
				}
			}

			if opt.Verbose {
				log.Infof("  %d pairs of TaxId mapping values from %d file(s) loaded", len(taxidMap), len(taxidMappingFiles))
			}

			mappingTaxids = len(taxidMap) > 0

			if mappingTaxids {
				taxdb = loadTaxonomy(opt, taxonomyDataDir)
			} else {
				checkError(fmt.Errorf("no valid TaxIds found in TaxId mapping file: %s", strings.Join(taxidMappingFiles, ", ")))
			}
		}

		// ---------------------------------------------------------------

		numFields := 13

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

			var matches map[uint64]*[]MatchResult // target -> match result
			var m MatchResult
			var ms *[]MatchResult
			var t *Target
			var ok bool
			var hTarget, h uint64
			var prevQuery string
			var floatMsSize float64
			var match MatchResult
			var first bool

			onlyTopNScore := topNScore > 0
			var nScore int
			var pScore float64
			var processThisMatch bool
			pScore = 1024
			nScore = 0
			processThisMatch = true

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
						pScore = 1024
						nScore = 0
						processThisMatch = true
					}

					if onlyTopNScore {
						if match.QCov < pScore { // match with a smaller score
							nScore++
							if nScore > topNScore {
								processThisMatch = false
							}
							pScore = match.QCov
						}
						if !processThisMatch {
							prevQuery = match.Query
							continue
						}
						if keepFullMatch && match.QCov == 1 {
							processThisMatch = false
							prevQuery = match.Query
							continue
						}
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

			var matches map[uint64]*[]MatchResult // target -> match result
			var ok bool
			var ms *[]MatchResult
			var hTarget, h, h1, h2 uint64
			var prevQuery string
			hs := make([]uint64, 0, 256)
			var match MatchResult
			var amb map[uint64]float64
			var i, j int
			var n, np1 int

			onlyTopNScore := topNScore > 0
			var nScore int
			var pScore float64
			var processThisMatch bool
			pScore = 1024
			nScore = 0
			processThisMatch = true

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
						pScore = 1024
						nScore = 0
						processThisMatch = true
					}

					if onlyTopNScore {
						if match.QCov < pScore { // match with a smaller score
							nScore++
							if nScore > topNScore {
								processThisMatch = false
							}
							pScore = match.QCov
						}
						if !processThisMatch {
							prevQuery = match.Query
							continue
						}
						if keepFullMatch && match.QCov == 1 {
							processThisMatch = false
							prevQuery = match.Query
							continue
						}
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

			var matches map[uint64]*[]MatchResult // target -> match result
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

			onlyTopNScore := topNScore > 0
			var nScore int
			var pScore float64
			var processThisMatch bool
			pScore = 1024
			nScore = 0
			processThisMatch = true

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
						pScore = 1024
						nScore = 0
						processThisMatch = true
					}

					if onlyTopNScore {
						if match.QCov < pScore { // match with a smaller score
							nScore++
							if nScore > topNScore {
								processThisMatch = false
							}
							pScore = match.QCov
						}
						if !processThisMatch {
							prevQuery = match.Query
							continue
						}
						if keepFullMatch && match.QCov == 1 {
							processThisMatch = false
							prevQuery = match.Query
							continue
						}
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
		log.Infof("#input reads: %.0f, #reads belonging to references in profile: %0.f, proportion: %.6f%%",
			nReads, (nReads - nUnassignedReads), (nReads-nUnassignedReads)/nReads*100)

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

		sorts.Quicksort(Targets(targets))

		var totalCoverage float64
		for _, t := range targets {
			totalCoverage += t.Coverage
		}

		var taxid uint32
		var ok bool

		// for limit ranks to show
		showRanksMap := make(map[string]interface{}, 128)
		for _, _rank := range showRanks {
			showRanksMap[_rank] = struct{}{}
		}

		rankPrefixesMap := make(map[string]string, len(rankPrefixes))
		for _i, _r := range showRanks {
			rankPrefixesMap[_r] = rankPrefixes[_i]
		}

		outfh.WriteString(fmt.Sprint("ref\tpercentage\tfragsProp\treads\tureads\trefsize\trefname\ttaxid\trank\ttaxname\ttaxpath\ttaxpathsn\n"))

		for _, t := range targets {
			if mappingNames {
				t.RefName = namesMap[t.Name]
			}

			if mappingTaxids {
				if taxid, ok = taxidMap[t.Name]; !ok {
					log.Warningf("%s is not mapped to any TaxId", t.Name)
				} else {
					t.AddTaxonomy(taxdb, showRanksMap, taxid)
				}
			}

			t.Percentage = t.Coverage / totalCoverage * 100

			outfh.WriteString(fmt.Sprintf("%s\t%.6f\t%.2f\t%.0f\t%.0f\t%d\t%s\t%d\t%s\t%s\t%s\t%s\n",
				t.Name, t.Percentage, t.FragsProp, t.SumMatch, t.SumUniqMatch, t.GenomeSize,
				t.RefName,
				taxid, t.Rank, t.TaxonName,
				strings.Join(t.LineageNames, separator),
				strings.Join(t.LineageTaxids, separator)))
		}

		// ---------------------------------------------------------------
		// more output

		var profile3 map[uint32]*ProfileNode
		var nodes []*ProfileNode

		if outputCamiReport || outputMetaphlanReport {
			profile3 = generateProfile(taxdb, targets)

			nodes = make([]*ProfileNode, 0, len(profile3))
			for _, node := range profile3 {
				nodes = append(nodes, node)
			}

			sort.Slice(nodes, func(i, j int) bool {
				if rankOrder[nodes[i].Rank] < rankOrder[nodes[j].Rank] {
					return true
				}

				if rankOrder[nodes[i].Rank] == rankOrder[nodes[j].Rank] {
					return nodes[i].Percentage > nodes[j].Percentage
				}

				return false
			})
		}

		// metaphlan format

		if outputMetaphlanReport {
			outfh2, gw2, w2, err := outStream(metaphlanReportFile, strings.HasSuffix(strings.ToLower(outFile), ".gz"), opt.CompressionLevel)
			checkError(err)
			defer func() {
				outfh2.Flush()
				if gw2 != nil {
					gw2.Close()
				}
				w2.Close()
			}()

			outfh2.WriteString("#SampleID\tKMCP\n")

			var lineageNames string
			filterByRank := len(showRanksMap) > 0
			names := make([]string, 0, 8)
			for _, node := range nodes {
				if filterByRank {
					if _, ok = showRanksMap[taxdb.Rank(node.Taxid)]; !ok {
						continue
					}

					names = names[:0]
					for i, taxid := range node.LineageTaxids {
						if _, ok = showRanksMap[taxdb.Rank(taxid)]; ok {
							names = append(names, rankPrefixesMap[taxdb.Rank(taxid)]+node.LineageNames[i])
						}
					}
					lineageNames = strings.Join(names, "|")
				} else {
					lineageNames = strings.Join(node.LineageNames, "|")
				}

				outfh2.WriteString(fmt.Sprintf("%s\t%.6f\n", lineageNames, node.Percentage))
			}
		}

		// cami format
		// https://github.com/bioboxes/rfc/blob/master/data-format/profiling.mkd

		if outputCamiReport {
			outfh3, gw3, w3, err := outStream(camiReportFile, strings.HasSuffix(strings.ToLower(outFile), ".gz"), opt.CompressionLevel)
			checkError(err)
			defer func() {
				outfh3.Flush()
				if gw3 != nil {
					gw3.Close()
				}
				w3.Close()
			}()

			outfh3.WriteString("@SampleID:\n")
			outfh3.WriteString("@Version:0.10.0\n")
			outfh3.WriteString("@Ranks:superkingdom|phylum|class|order|family|genus|species|strain\n")
			outfh3.WriteString("@TaxonomyID:ncbi-taxonomy\n")
			outfh3.WriteString("@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n")

			var lineageTaxids, lineageNames string
			filterByRank := len(showRanksMap) > 0
			names := make([]string, 0, 8)
			taxids := make([]string, 0, 8)
			for _, node := range nodes {
				if filterByRank {
					if _, ok = showRanksMap[taxdb.Rank(node.Taxid)]; !ok {
						continue
					}

					names = names[:0]
					taxids = taxids[:0]
					for i, taxid := range node.LineageTaxids {
						if _, ok = showRanksMap[taxdb.Rank(taxid)]; ok {
							taxids = append(taxids, strconv.Itoa(int(taxid)))
							names = append(names, node.LineageNames[i])
						}
					}
					lineageTaxids = strings.Join(taxids, "|")
					lineageNames = strings.Join(names, "|")
				} else {
					taxids = taxids[:0]
					for _, taxid := range node.LineageTaxids {
						taxids = append(taxids, strconv.Itoa(int(taxid)))
					}
					lineageTaxids = strings.Join(taxids, "|")
					lineageNames = strings.Join(node.LineageNames, "|")
				}

				outfh3.WriteString(fmt.Sprintf("%d\t%s\t%s\t%s\t%.6f\n",
					node.Taxid, node.Rank, lineageTaxids, lineageNames, node.Percentage))
			}
		}

	},
}

func init() {
	RootCmd.AddCommand(profileCmd)

	profileCmd.Flags().IntP("chunk-size", "", 5000, `number of lines to process for each thread, and 4 threads is fast enough. Type "kmcp profile -h" for details`)

	profileCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)

	// for single read
	profileCmd.Flags().Float64P("max-fpr", "f", 0.01, `maximal false positive rate of a read`)
	profileCmd.Flags().Float64P("min-qcov", "t", 0.7, `minimal query coverage of a read`)
	profileCmd.Flags().IntP("keep-top-scores", "n", 0, `keep matches with the top N score, 0 for all`)
	profileCmd.Flags().BoolP("keep-full-match", "F", false, `only keep the full matches if have`)

	// for ref fragments
	profileCmd.Flags().IntP("min-reads", "r", 50, `minimal number of reads for a reference fragment`)
	profileCmd.Flags().IntP("min-uniq-reads", "u", 5, `minimal number of unique matched reads for a reference fragment`)
	profileCmd.Flags().Float64P("min-frags-prop", "p", 0.5, `minimal proportion of matched fragments`)

	// for the two-stage taxonomy assignment algorithm in MagaPath
	profileCmd.Flags().Float64P("min-dreads-prop", "D", 0.05, `minimal proportion of distinct reads, for determing the right reference for ambigous reads`)
	profileCmd.Flags().Float64P("max-mismatch-err", "R", 0.05, `maximal error rate of a read being matched to a wrong reference, for determing the right reference for ambigous reads`)

	// name mapping
	profileCmd.Flags().StringSliceP("name-map", "N", []string{}, `tabular two-column file(s) mapping reference IDs to reference names`)

	// taxonomy
	profileCmd.Flags().StringSliceP("taxid-map", "T", []string{}, `tabular two-column file(s) mapping reference IDs to TaxIds`)
	profileCmd.Flags().StringP("taxonomy-dir", "X", "", `directory of NCBI taxonomy dump files: names.dmp, nodes.dmp, optional with merged.dmp and delnodes.dmp`)
	profileCmd.Flags().StringP("separator", "s", ";", `separator of TaxIds and taxonomy names`)
	profileCmd.Flags().StringSliceP("show-rank", "", []string{"superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain"}, "only show TaxIds and names of these ranks")
	profileCmd.Flags().StringSliceP("rank-prefix", "", []string{"k__", "p__", "c__", "o__", "f__", "g__", "s__", "t__"}, "prefixes of taxon name in certain ranks, used with --metaphlan-report ")

	// other output formats
	profileCmd.Flags().StringP("metaphlan-report", "M", "", `save extra metaphlan-like report`)
	profileCmd.Flags().StringP("cami-report", "C", "", `save extra CAMI-like report`)
}
