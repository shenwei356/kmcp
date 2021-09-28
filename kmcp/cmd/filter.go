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
	"runtime"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/pkg/errors"
	"github.com/shenwei356/breader"
	"github.com/shenwei356/unikmer"
	"github.com/shenwei356/util/cliutil"
	"github.com/spf13/cobra"
	"github.com/zeebo/wyhash"
)

var filterCmd = &cobra.Command{
	Use:   "filter",
	Short: "Filter search results and find species/assembly-specific queries",
	Long: `Filter search results and find species/assembly-specific queries

Taxonomy data:
  1. Mapping references IDs to TaxIds: -T/--taxid-map
  2. NCBI taxonomy dump files: -X/--taxdump

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

		noHeaderRow := getFlagBool(cmd, "no-header-row")

		outFile := getFlagString(cmd, "out-prefix")

		maxFPR := getFlagPositiveFloat64(cmd, "max-fpr")
		minQcov := getFlagNonNegativeFloat64(cmd, "min-query-cov")

		level := strings.ToLower(getFlagString(cmd, "level"))
		var levelSpecies bool
		switch level {
		case "species":
			levelSpecies = true
		case "strain", "assembly":
			levelSpecies = false
		default:
			checkError(fmt.Errorf("invalid value for --level, available values: species, strain/assembly"))
		}

		// -----

		nameMappingFiles := getFlagStringSlice(cmd, "name-map")

		taxidMappingFiles := getFlagStringSlice(cmd, "taxid-map")
		taxonomyDataDir := getFlagString(cmd, "taxdump")

		if len(taxidMappingFiles) > 0 && taxonomyDataDir == "" {
			checkError(fmt.Errorf("flag -X/--taxonomy-dir is needed when -T/--taxid-map given"))
		}
		if len(taxidMappingFiles) == 0 && taxonomyDataDir != "" {
			checkError(fmt.Errorf("flag -T/--taxid-map is needed when -X/--taxonomy-dir given"))
		}
		if (len(taxidMappingFiles) == 0 || taxonomyDataDir == "") && levelSpecies {
			checkError(fmt.Errorf("TaxID mapping files (-T/--taxid-map) and taxonomy dump files are needed for --level species"))
		}

		mappingTaxids := len(taxidMappingFiles) != 0
		if levelSpecies && !mappingTaxids {
			checkError(fmt.Errorf("-T/--taxid-map needed for --level species"))
		}

		chunkSize := getFlagPositiveInt(cmd, "chunk-size")
		if opt.NumCPUs > 4 {
			if opt.Verbose || opt.Log2File {
				log.Infof("using a lot of threads does not always accelerate processing, 4-threads is fast enough")
			}
			opt.NumCPUs = 4
			runtime.GOMAXPROCS(opt.NumCPUs)
		}

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
		// name mapping files

		var namesMap map[string]string
		mappingNames := len(nameMappingFiles) != 0
		if mappingNames {
			if opt.Verbose || opt.Log2File {
				log.Infof("loading name mapping file ...")
			}
			nameMappingFile := nameMappingFiles[0]
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

			if opt.Verbose || opt.Log2File {
				log.Infof("  %d pairs of name mapping values from %d file(s) loaded", len(namesMap), len(nameMappingFiles))
			}

			mappingNames = len(namesMap) > 0
		}

		// ---------------------------------------------------------------
		// taxid mapping files

		var taxdb *unikmer.Taxonomy
		var taxidMap map[string]uint32

		if mappingTaxids {
			if opt.Verbose || opt.Log2File {
				log.Infof("loading TaxId mapping file ...")
			}
			taxidMappingFile := taxidMappingFiles[0]
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

			if opt.Verbose || opt.Log2File {
				log.Infof("  %d pairs of TaxId mapping values from %d file(s) loaded", len(taxidMap), len(taxidMappingFiles))
			}

			mappingTaxids = len(taxidMap) > 0

			if mappingTaxids {
				taxdb = loadTaxonomy(opt, taxonomyDataDir)
				taxdb.CacheLCA()
			} else {
				checkError(fmt.Errorf("no valid TaxIds found in TaxId mapping file: %s", strings.Join(taxidMappingFiles, ", ")))
			}
		}

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

		if !noHeaderRow {
			outfh.WriteString("#query\tqLen\tqKmers\tFPR\thits\ttarget\tfragIdx\tfrags\ttLen\tkSize\tmKmers\tqCov\ttCov\tjacc\tqueryIdx\n")
		}

		if opt.Verbose || opt.Log2File {
			log.Info()
			log.Infof("-------------------- [main parameters] --------------------")

			log.Infof("match filtration: ")
			log.Infof("  maximal false positive rate: %f", maxFPR)
			log.Infof("  minimal query coverage: %4f", minQcov)
			log.Info()

			if mappingTaxids {
				log.Infof("taxonomy data:")
				log.Infof("  taxdump directory: %s", taxonomyDataDir)
				log.Infof("  mapping reference IDs to TaxIds: %s", taxidMappingFiles)
				log.Info()
			}

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

			match, ok := parseMatchResult2(line, numFields, items, maxFPR, minQcov)
			if !ok {
				pool.Put(items)
				return nil, false, nil
			}

			pool.Put(items)
			return match, true, nil
		}

		var nReads int
		var nPassed int

		// ---------------------------------------------------------------
		if opt.Verbose || opt.Log2File {
			log.Info("filtering ...")
		}

		for _, file := range files {
			if opt.Verbose || opt.Log2File {
				log.Infof("  parsing file: %s", file)
			}

			var matches map[uint64]*[]*MatchResult2 // target -> match result
			var m *MatchResult2
			var ms *[]*MatchResult2
			var ok bool
			var hTarget uint64
			var prevQuery string
			var match *MatchResult2

			taxids := make([]uint32, 0, 128)
			var taxid1, taxid2 uint32
			var theSameSpecies bool

			reader, err := breader.NewBufferedReader(file, opt.NumCPUs, chunkSize, fn)
			checkError(err)
			var data interface{}

			matches = make(map[uint64]*[]*MatchResult2)
			for chunk := range reader.Ch {
				checkError(chunk.Err)

				for _, data = range chunk.Data {
					match = data.(*MatchResult2)

					if prevQuery != match.Query { // new query
						nReads++

						if len(matches) > 0 { // not the first query
							if levelSpecies {
								taxids = taxids[:0]
								for _, ms = range matches {
									taxids = append(taxids, taxidMap[(*ms)[0].Target])
								}

								theSameSpecies = false
								taxid1 = taxids[0]
								for _, taxid2 = range taxids[1:] {
									taxid1 = taxdb.LCA(taxid1, taxid2)
								}
								if taxdb.AtOrBelowRank(taxid1, "species") {
									theSameSpecies = true
								}
							}

							if len(matches) == 1 || theSameSpecies {
								nPassed++
								for _, ms = range matches {
									for _, m = range *ms { // multiple matches in different fragments
										outfh.WriteString(*m.Line)
									}
									poolMatchResults.Put(ms)
								}
							}
						}

						matches = make(map[uint64]*[]*MatchResult2)
					}

					hTarget = wyhash.HashString(match.Target, 1)
					if ms, ok = matches[hTarget]; !ok {
						// tmp := []*MatchResult{match}
						tmp := poolMatchResults2.Get().(*[]*MatchResult2)
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

			if len(matches) == 1 || theSameSpecies {
				nPassed++
				for _, ms = range matches {
					for _, m = range *ms { // multiple matches in different fragments
						outfh.WriteString(*m.Line)
					}
					poolMatchResults.Put(ms)
				}
			}
		}

		if opt.Verbose || opt.Log2File {
			log.Info()
			log.Infof("%.4f%% (%d/%d) queries filtered", float64(nPassed)/float64(nReads)*100, nPassed, nReads)
		}
	},
}

func init() {
	RootCmd.AddCommand(filterCmd)

	filterCmd.Flags().IntP("chunk-size", "", 5000, `number of lines to process for each thread, and 4 threads is fast enough. Type "kmcp profile -h" for details`)
	filterCmd.Flags().BoolP("no-header-row", "H", false, `do not print header row`)
	filterCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)

	// for single read
	filterCmd.Flags().Float64P("max-fpr", "f", 0.05, `maximal false positive rate of a read in search result`)
	filterCmd.Flags().Float64P("min-query-cov", "t", 0.55, `minimal query coverage of a read in search result`)

	// taxonomy
	filterCmd.Flags().StringSliceP("taxid-map", "T", []string{}, `tabular two-column file(s) mapping reference IDs to TaxIds`)
	filterCmd.Flags().StringP("taxdump", "X", "", `directory of NCBI taxonomy dump files: names.dmp, nodes.dmp, optional with merged.dmp and delnodes.dmp`)

	// name mapping
	filterCmd.Flags().StringSliceP("name-map", "N", []string{}, `tabular two-column file(s) mapping reference IDs to reference names`)

	filterCmd.Flags().StringP("level", "", "species", `level to estimate abundance at. available values: species, strain/assembly`)
}
