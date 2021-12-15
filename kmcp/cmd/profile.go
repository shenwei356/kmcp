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
	"io"
	"math"
	"os"
	"path/filepath"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/pkg/errors"
	"github.com/shenwei356/bio/taxdump"
	"github.com/shenwei356/breader"
	"github.com/shenwei356/util/cliutil"
	"github.com/shenwei356/util/stats"
	"github.com/spf13/cobra"
	"github.com/twotwotwo/sorts"
	"github.com/zeebo/wyhash"
)

var profileCmd = &cobra.Command{
	Use:   "profile",
	Short: "Generate taxonomic profile from search results",
	Long: `Generate taxonomic profile from search results

Methods:
  1. Reference genomes can be splitted into fragments when computing
     k-mers (sketches), which could help to increase the specificity
     via a threshold, i.e., the minimal proportion of matched fragments
     (-p/--min-frags-fraction). (***highly recommended***)
     Another flag -d/--max-frags-depth-stdev further reduces false positives.
  2. We require part of the uniquely matched reads of a reference
     having high similarity, i.e., with high confidence for decreasing
     the false positive rate.
  3. We also use the two-stage taxonomy assignment algorithm in MegaPath
     to reduce the false positive of ambiguous matches.
  4. Multi-aligned queries are proportionally assigned to references
     with a similar strategy in Metalign.
  5. Input files are parsed fours times, therefore STDIN is not supported.

Reference:
  1. MegaPath: https://doi.org/10.1186/s12864-020-06875-6
  2. Metalign: https://doi.org/10.1186/s13059-020-02159-0

Accuracy notes:
  *. Smaller -t/--min-qcov increase sensitivity in cost of higher false
     positive rate (-f/--max-fpr) of a query.
  *. We require part of the uniquely matched reads of a reference
     having high similarity, i.e., with high confidence for decreasing
     the false positive rate.
     E.g., -H >= 0.8 and -P >= 0.1 equals to 90th percentile >= 0.8
     *. -U/--min-hic-ureads,      minimal number, >= 1
     *. -H/--min-hic-ureads-qcov, minimal query coverage, >= -t/--min-qcov
     *. -P/--min-hic-ureads-prop, minimal proportion, higher values
        increase precision in cost of sensitivity.
  *. -R/--max-mismatch-err and -D/--min-dreads-prop is for determing
     the right reference for ambigous reads.
  *. --keep-perfect-match is not recommended, which decreases sensitivity. 
  *. --keep-main-match is not recommended, which affects accuracy of
     abundance estimation.
  *. -n/--keep-top-qcovs  is not recommended, which affects accuracy of
     abundance estimation.

Profiling modes:
  We preset five profiling modes, availabe with the flag -m/--mode:
    - 0 (for pathogen detection)
    - 1 (highest recall)
    - 2 (high recall)
    - 3 (default)
    - 4 (higher precision)
    - 5 (highest precision)
  Using this flag will overide the relevant options.

    options                      m=0    m=1   m=2   m=3    m=4   m=5
    --------------------------   ----   ---   ---   ----   ---   ----
    -r/--min-frags-reads         1      20    30    50     100   100
    -p/--min-frags-fraction      0.01   0.5   0.7   0.8    1     1
    -d/--max-frags-depth-stdev   100    10    3     2      2     1.5
    -u/--min-uniq-reads          1      20    20    20     50    50
    -U/--min-hic-ureads          1      5     5     5      10    10
    -H/--min-hic-ureads-qcov     0.55   0.7   0.7   0.75   0.8   0.8
    -P/--min-hic-ureads-prop     0.01   0.1   0.2   0.1    0.1   0.15


Taxonomy data:
  1. Mapping references IDs to TaxIds: -T/--taxid-map
  2. NCBI taxonomy dump files: -X/--taxdump

Performance notes:
  1. Searching results are parsed in parallel, and the number of
     lines proceeded by a thread can be set by the flag --chunk-size.
  2. However using a lot of threads does not always accelerate
     processing, 4 threads with chunk size of 500-5000 is fast enough.

Profiling output formats:
  1. KMCP      (-o/--out-prefix)
  2. CAMI      (-M/--metaphlan-report)
  3. MetaPhlAn (-C/--cami-report)

Taxonomic binning formats:
  1. CAMI      (-B/--binning-result)
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

		// ---------------- debug ---------

		// debugFile := getFlagString(cmd, "debug")
		debugFile := ""

		debug := debugFile != ""
		var outfhD *bufio.Writer
		var gwD io.WriteCloser
		var wD *os.File
		if debug {
			outfhD, gwD, wD, err = outStream(debugFile, strings.HasSuffix(strings.ToLower(debugFile), ".gz"), opt.CompressionLevel)
			checkError(err)
			defer func() {
				outfhD.Flush()
				if gwD != nil {
					gwD.Close()
				}
				wD.Close()
			}()
		}
		// ---------------- debug ---------

		outFile := getFlagString(cmd, "out-prefix")

		maxFPR := getFlagPositiveFloat64(cmd, "max-fpr")
		minQcov := getFlagNonNegativeFloat64(cmd, "min-query-cov")
		topNScore := getFlagNonNegativeInt(cmd, "keep-top-qcovs")
		keepFullMatch := getFlagBool(cmd, "keep-perfect-match")
		keepMainMatch := getFlagBool(cmd, "keep-main-match")
		maxScoreGap := getFlagFloat64(cmd, "max-qcov-gap")

		var minReads float64
		var minFragsProp float64
		var maxFragsDepthStdev float64
		var minUReads float64
		var minHicUreads float64
		var hicUreadsMinQcov float64
		var HicUreadsMinProp float64

		minReads = float64(getFlagPositiveInt(cmd, "min-frags-reads"))
		minUReads = float64(getFlagPositiveInt(cmd, "min-uniq-reads"))
		minFragsProp = getFlagNonNegativeFloat64(cmd, "min-frags-fraction")
		maxFragsDepthStdev = getFlagPositiveFloat64(cmd, "max-frags-depth-stdev")

		minHicUreads = float64(getFlagPositiveInt(cmd, "min-hic-ureads"))
		if minHicUreads > minUReads {
			minUReads = minHicUreads
		}
		hicUreadsMinQcov = getFlagPositiveFloat64(cmd, "min-hic-ureads-qcov")
		if hicUreadsMinQcov < minQcov {
			hicUreadsMinQcov = minQcov
		}
		HicUreadsMinProp = getFlagPositiveFloat64(cmd, "min-hic-ureads-prop")

		minDReadsProp := getFlagPositiveFloat64(cmd, "min-dreads-prop")
		if minDReadsProp > 1 {
			checkError(fmt.Errorf("the value of -D/--min-dreads-prop (%f) should be in range of (0, 1]", minDReadsProp))
		}
		maxMismatchErr := getFlagPositiveFloat64(cmd, "max-mismatch-err")
		if maxMismatchErr >= 1 {
			checkError(fmt.Errorf("the value of -R/--max-mismatch-err (%f) should be in range of (0, 1)", maxMismatchErr))
		}

		mode := getFlagNonNegativeInt(cmd, "mode")
		switch mode {
		case 3:
		case 0:
			minReads = 1
			minFragsProp = 0.01
			maxFragsDepthStdev = 100
			minUReads = 1
			minHicUreads = 1
			hicUreadsMinQcov = 0.55
			HicUreadsMinProp = 0.01
		case 1:
			minReads = 20
			minFragsProp = 0.5
			maxFragsDepthStdev = 10
			minUReads = 20
			minHicUreads = 5
			hicUreadsMinQcov = 0.7
			HicUreadsMinProp = 0.1
		case 2:
			minReads = 30
			minFragsProp = 0.7
			maxFragsDepthStdev = 3
			minUReads = 20
			minHicUreads = 5
			hicUreadsMinQcov = 0.7
			HicUreadsMinProp = 0.2
		case 4:
			minReads = 100
			minFragsProp = 1
			maxFragsDepthStdev = 2
			minUReads = 50
			minHicUreads = 10
			hicUreadsMinQcov = 0.8
			HicUreadsMinProp = 0.1
		case 5:
			minReads = 100
			minFragsProp = 1
			maxFragsDepthStdev = 1.5
			minUReads = 50
			minHicUreads = 10
			hicUreadsMinQcov = 0.8
			HicUreadsMinProp = 0.15
		default:
			checkError(fmt.Errorf("invalid mode: %d", mode))
		}

		lowAbcPct := getFlagNonNegativeFloat64(cmd, "filter-low-pct")
		if lowAbcPct >= 100 {
			checkError(fmt.Errorf("the value of -F/--filter-low-pct (%f) should be in range of [0, 100)", lowAbcPct))
		} else if lowAbcPct > 10 {
			log.Warningf("the value of -F/--filter-low-pct (%v) may be too big", lowAbcPct)
		}
		fileterLowAbc := lowAbcPct > 0

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

		// UregionPropFiles := getFlagStringSlice(cmd, "uregion-prop-map")
		UregionPropFiles := []string{}
		considerUregionProp := len(UregionPropFiles) > 0

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
		if len(taxidMappingFiles) == 0 || taxonomyDataDir == "" {
			log.Warningf("TaxID mapping files (-T/--taxid-map) and taxonomy dump files are recommended to add taxonomy information")
		}

		mappingTaxids := len(taxidMappingFiles) != 0
		if levelSpecies && !mappingTaxids {
			checkError(fmt.Errorf("-T/--taxid-map needed for --level species"))
		}

		separator := getFlagString(cmd, "separator")
		if separator == "" {
			log.Warningf("value of -s/--separator better not be empty")
		}

		chunkSize := getFlagPositiveInt(cmd, "chunk-size")
		if opt.NumCPUs > 4 {
			if opt.Verbose || opt.Log2File {
				log.Infof("using a lot of threads does not always accelerate processing, 4-threads is fast enough")
			}
			opt.NumCPUs = 4
			runtime.GOMAXPROCS(opt.NumCPUs)
		}

		sampleID := getFlagString(cmd, "sample-id")
		taxonomyID := getFlagString(cmd, "taxonomy-id")
		binningFile := getFlagString(cmd, "binning-result")
		outputBinningResult := binningFile != ""
		if outputBinningResult && !(strings.HasSuffix(binningFile, ".binning") || strings.HasSuffix(binningFile, ".binning.gz")) {
			binningFile = binningFile + ".binning.gz"
		}
		if outputBinningResult && (len(taxidMappingFiles) == 0 || taxonomyDataDir == "") {
			checkError(fmt.Errorf("flag -T/--taxid-map and -T/--taxid-map needed when -B/--binning-result given"))
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

		if (outputBinningResult || outputCamiReport || outputMetaphlanReport) && !mappingTaxids {
			log.Warningf("TaxID mapping files (-T/--taxid-map) and taxonomy dump files are needed to output CAMI/MetaPhlan/binning report")
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

		normAbund := getFlagString(cmd, "norm-abund")
		switch normAbund {
		case "mean", "min", "max":
		default:
			checkError(fmt.Errorf("invalid value of --norm-abund: %s. available: mean, min, max", normAbund))
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

		var taxdb *taxdump.Taxonomy
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

		// ---------------------------------------------------------------

		var uregionPropMap map[string]float64

		if considerUregionProp {
			if opt.Verbose || opt.Log2File {
				log.Infof("loading reference unique regions proportion mapping file ...")
			}
			uregionPropFile := UregionPropFiles[0]
			uregionPropStr, err := cliutil.ReadKVs(uregionPropFile, false)
			if err != nil {
				checkError(errors.Wrap(err, uregionPropFile))
			}
			uregionPropMap = make(map[string]float64, len(uregionPropStr))
			var p float64
			for k, s := range uregionPropStr {
				p, err = strconv.ParseFloat(s, 64)
				if err != nil {
					checkError(fmt.Errorf("invalid proportion: %s", s))
				}
				uregionPropMap[k] = p
			}

			if len(UregionPropFiles) > 1 {
				for _, uregionPropFile := range taxidMappingFiles[1:] {
					uregionPropStr, err := cliutil.ReadKVs(uregionPropFile, false)
					if err != nil {
						checkError(errors.Wrap(err, uregionPropFile))
					}
					uregionPropMap = make(map[string]float64, len(uregionPropStr))
					var p float64
					for k, s := range uregionPropStr {
						p, err = strconv.ParseFloat(s, 64)
						if err != nil {
							checkError(fmt.Errorf("invalid proportion: %s", s))
						}
						uregionPropMap[k] = p
					}
				}
			}

			if opt.Verbose || opt.Log2File {
				log.Infof("  %d pairs of reference unique regions proportion mapping values from %d file(s) loaded", len(uregionPropMap), len(UregionPropFiles))
			}

			considerUregionProp = len(uregionPropMap) > 0

			if !considerUregionProp {
				checkError(fmt.Errorf("no valid reference unique regions proportions found in mapping file: %s", strings.Join(UregionPropFiles, ", ")))
			}
		}

		if opt.Verbose || opt.Log2File {
			log.Info()
			log.Infof("-------------------- [main parameters] --------------------")

			log.Infof("match filtration: ")
			log.Infof("  maximal false positive rate: %f", maxFPR)
			log.Infof("  minimal query coverage: %4f", minQcov)
			log.Infof("  keep matches with the top N scores: N=%d", topNScore)
			log.Infof("  only keep the full matches: %v", keepFullMatch)
			log.Infof("  only keep main matches: %v, maximal score gap: %f", keepMainMatch, maxScoreGap)
			log.Info()

			log.Infof("deciding the existence of a reference:")
			log.Infof("  minimal number of reads per reference fragment: %.0f", minReads)
			log.Infof("  minimal number of uniquely matched reads: %.0f", minUReads)
			log.Infof("  minimal proportion of matched reference fragments: %f", minFragsProp)
			log.Infof("  maximal standard deviation of relative depths of all fragments: %f", maxFragsDepthStdev)
			log.Info()

			log.Infof("  minimal number of high-confidence uniquely matched reads: %.0f", minHicUreads)
			log.Infof("  minimal query coverage of high-confidence uniquely matched reads: %f", hicUreadsMinQcov)
			log.Infof("  minimal proportion of high-confidence uniquely matched reads: %f", HicUreadsMinProp)
			log.Info()

			if mappingTaxids {
				log.Infof("taxonomy data:")
				log.Infof("  taxdump directory: %s", taxonomyDataDir)
				log.Infof("  mapping reference IDs to TaxIds: %s", taxidMappingFiles)
				log.Info()
			}

			log.Infof("reporting:")
			if mappingNames {
				log.Infof("  mapping reference IDs to names: %s", nameMappingFiles)
			}
			if fileterLowAbc {
				log.Infof("  filter out predictions with the smallest relative abundances summing up %d%%", lowAbcPct)
			}
			log.Infof("  default format  : %s", outFile)
			if outputCamiReport {
				log.Infof("  CAMI format     : %s", camiReportFile)
				log.Infof("    Sample ID     : %s", sampleID)
			}
			if outputMetaphlanReport {
				log.Infof("  MetaPhlan format: %s", metaphlanReportFile)
				log.Infof("    Sample ID     : %s", sampleID)
				log.Infof("    Taxonomy ID   : %s", taxonomyID)
			}
			if outputBinningResult {
				log.Infof("  Binning result  : %s", binningFile)
			}

			log.Infof("-------------------- [main parameters] --------------------")
			log.Info()

		}
		// ---------------------------------------------------------------

		numFields := 13

		profile := make(map[uint64]*Target, 128)

		floatOne := float64(1)

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

			match, ok := parseMatchResult(line, numFields, items, maxFPR, minQcov)
			if !ok {
				pool.Put(items)
				return nil, false, nil
			}

			pool.Put(items)
			return match, true, nil
		}

		var nReads float64

		// ---------------------------------------------------------------
		// stage 1/4
		if opt.Verbose || opt.Log2File {
			log.Infof("stage 1/4: counting matches and unique matches for filtering out low-confidence references")
		}
		timeStart1 := time.Now()

		for _, file := range files {
			if opt.Verbose || opt.Log2File {
				log.Infof("  parsing file: %s", file)
			}

			var matches map[uint64]*[]*MatchResult // target -> match result
			var m *MatchResult
			var ms *[]*MatchResult
			var t *Target
			var ok bool
			var hTarget, h uint64
			var prevQuery string
			var floatMsSize float64
			var match *MatchResult
			var first bool

			onlyTopNScore := topNScore > 0
			var nScore int
			var pScore float64
			var processThisMatch bool

			taxids := make([]uint32, 0, 128)
			var taxid1, taxid2 uint32
			var theSameSpecies bool

			reader, err := breader.NewBufferedReader(file, opt.NumCPUs, chunkSize, fn)
			checkError(err)
			var data interface{}

			matches = make(map[uint64]*[]*MatchResult)
			pScore = 1024
			nScore = 0
			processThisMatch = true
			for chunk := range reader.Ch {
				checkError(chunk.Err)

				for _, data = range chunk.Data {
					match = data.(*MatchResult)

					if prevQuery != match.Query { // new query
						nReads++

						if len(matches) > 0 { // not the first query
							if levelSpecies {
								taxids = taxids[:0]
								for h, ms = range matches {
									taxid1, ok = taxidMap[(*ms)[0].Target]
									if !ok {
										checkError(fmt.Errorf("unknown taxid for %s, please check taxid mapping file(s)", (*ms)[0].Target))
									}
									taxids = append(taxids, taxid1)
								}
								// LCA
								theSameSpecies = false
								taxid1 = taxids[0]
								for _, taxid2 = range taxids[1:] {
									taxid1 = taxdb.LCA(taxid1, taxid2)
								}
								if taxdb.AtOrBelowRank(taxid1, "species") {
									theSameSpecies = true
								}
							}

							for h, ms = range matches {
								floatMsSize = float64(len(*ms))
								first = true
								for _, m = range *ms { // multiple matches in different fragments
									if t, ok = profile[h]; !ok {
										t0 := Target{
											Name:         m.Target,
											GenomeSize:   m.GSize,
											Match:        make([]float64, m.IdxNum),
											UniqMatch:    make([]float64, m.IdxNum),
											UniqMatchHic: make([]float64, m.IdxNum),
											// QLen:         make([]float64, m.IdxNum),
											// RelDepth:  make([]float64, m.IdxNum),
										}
										profile[h] = &t0
										t = &t0
									}

									if first { // count once
										if len(matches) == 1 || theSameSpecies {
											t.UniqMatch[m.FragIdx]++
											if m.QCov >= hicUreadsMinQcov {
												t.UniqMatchHic[m.FragIdx]++
											}
										}
										// t.QLen[m.FragIdx] += float64(m.QLen)

										first = false
									}

									// for a read matching multiple regions of a reference, distribute count to multiple regions,
									// the sum is still the one.
									t.Match[m.FragIdx] += floatOne / floatMsSize
								}
								poolMatchResults.Put(ms)
							}
						}

						matches = make(map[uint64]*[]*MatchResult)
						pScore = 1024
						nScore = 0
						processThisMatch = true
					} else if keepFullMatch { // not the first match
						if !processThisMatch {
							prevQuery = match.Query
							continue
						}

						if pScore == 1 && match.QCov < 1 {
							processThisMatch = false

							prevQuery = match.Query
							continue
						}
					} else if keepMainMatch && pScore <= 1 {
						if !processThisMatch {
							prevQuery = match.Query
							continue
						}

						if pScore-match.QCov > maxScoreGap {
							processThisMatch = false

							prevQuery = match.Query
							continue
						}
					}

					if onlyTopNScore {
						if !processThisMatch {
							prevQuery = match.Query
							continue
						}

						if match.QCov < pScore { // match with a smaller score
							nScore++
							if nScore > topNScore {
								processThisMatch = false

								prevQuery = match.Query
								continue
							}
						}
					}

					hTarget = wyhash.HashString(match.Target, 1)
					if ms, ok = matches[hTarget]; !ok {
						// tmp := []*MatchResult{match}
						tmp := poolMatchResults.Get().(*[]*MatchResult)
						*tmp = (*tmp)[:0]
						*tmp = append(*tmp, match)
						matches[hTarget] = tmp
					} else {
						*ms = append(*ms, match)
					}

					prevQuery = match.Query
					pScore = match.QCov
				}
			}

			if len(matches) > 0 {
				nReads++

				if levelSpecies {
					taxids = taxids[:0]
					for h, ms = range matches {
						taxid1, ok = taxidMap[(*ms)[0].Target]
						if !ok {
							checkError(fmt.Errorf("unknown taxid for %s, please check taxid mapping file(s)", (*ms)[0].Target))
						}
						taxids = append(taxids, taxid1)
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

				for h, ms = range matches {
					floatMsSize = float64(len(*ms))
					first = true
					for _, m = range *ms { // multiple matches in different fragments
						if t, ok = profile[h]; !ok {
							t0 := Target{
								Name:         m.Target,
								GenomeSize:   m.GSize,
								Match:        make([]float64, m.IdxNum),
								UniqMatch:    make([]float64, m.IdxNum),
								UniqMatchHic: make([]float64, m.IdxNum),
								// QLen:         make([]float64, m.IdxNum),
								// RelDepth:  make([]float64, m.IdxNum),
							}
							profile[h] = &t0
							t = &t0
						}

						if first { // count once
							if len(matches) == 1 || theSameSpecies {
								t.UniqMatch[m.FragIdx]++
								if m.QCov >= hicUreadsMinQcov {
									t.UniqMatchHic[m.FragIdx]++
								}
							}
							// t.QLen[m.FragIdx] += float64(m.QLen)

							first = false
						}

						// for a read matching multiple regions of a reference, distribute count to multiple regions,
						// the sum is still the one.
						t.Match[m.FragIdx] += floatOne / floatMsSize
					}
					poolMatchResults.Put(ms)
				}
			}
		}

		// --------------------
		// sum up #1

		if debug {
			outfhD.WriteString("#------------------ round 1 ------------------\n")
		}

		if opt.Verbose || opt.Log2File {
			log.Infof("  number of references in search result: %d", len(profile))
		}

		var c float64
		var c1 float64
		var c2 float64
		var hs []uint64
		hs = make([]uint64, 0, 10240) // list to delete
		for h, t := range profile {
			for _, c1 = range t.UniqMatch {
				t.SumUniqMatch += c1
			}
			if t.SumUniqMatch < 1 { // no enough unique match
				hs = append(hs, h)
				if debug {
					fmt.Fprintf(outfhD, "failed1: %s, %s: %.0f\n",
						t.Name, "no enough unique match", t.SumUniqMatch)
				}
				continue
			}

			for _, c1 = range t.UniqMatchHic {
				t.SumUniqMatchHic += c1
			}
			if t.SumUniqMatchHic < 1 { // no enough high-confidence unique match
				hs = append(hs, h)
				if debug {
					fmt.Fprintf(outfhD, "failed1: %s, %s: %.0f\n",
						t.Name, "no enough high-confidence unique match", t.SumUniqMatchHic)
				}
				continue
			}

			// the SumUniqMatchHic may increase later

			// if t.SumUniqMatchHic < t.SumUniqMatch*HicUreadsMinProp {
			// 	hs = append(hs, h)
			// 	continue
			// }

			// ---------------

			for _, c = range t.Match {
				if c > 0 {
					t.FragsProp++
				}
				t.SumMatch += c
			}
			t.FragsProp = t.FragsProp / float64(len(t.Match))
			if t.FragsProp < minFragsProp { // low coverage
				hs = append(hs, h)
				if debug {
					fmt.Fprintf(outfhD, "failed1: %s, %s: %.1f %v\n",
						t.Name, "low fragments fraction", t.FragsProp, t.Match)
				}
				continue
			}
		}

		for _, h := range hs {
			delete(profile, h)
		}

		if opt.Verbose || opt.Log2File {
			log.Infof("  number of estimated references: %d", len(profile))
			log.Infof("  elapsed time: %s", time.Since(timeStart1))
			log.Info()
		}

		// ---------------------------------------------------------------
		// stage 2/4, counting ambiguous reads/matches
		if opt.Verbose || opt.Log2File {
			log.Infof("stage 2/4: counting ambiguous matches for correcting matches")
		}
		timeStart1 = time.Now()

		// hashA -> hashB -> count
		ambMatch := make(map[uint64]map[uint64]float64, len(profile))

		for _, file := range files {
			if opt.Verbose || opt.Log2File {
				log.Infof("  parsing file: %s", file)
			}

			var matches map[uint64]*[]*MatchResult // target -> match result
			var ok bool
			// var ms *[]*MatchResult
			var hTarget, h, h1, h2 uint64
			var prevQuery string
			hs := make([]uint64, 0, 256)
			var match *MatchResult
			var amb map[uint64]float64
			var i, j int
			var n, np1 int

			onlyTopNScore := topNScore > 0
			var nScore int
			var pScore float64
			var processThisMatch bool

			reader, err := breader.NewBufferedReader(file, opt.NumCPUs, chunkSize, fn)
			checkError(err)
			var data interface{}

			matches = make(map[uint64]*[]*MatchResult)
			pScore = 1024
			nScore = 0
			processThisMatch = true
			for chunk := range reader.Ch {
				checkError(chunk.Err)

				for _, data = range chunk.Data {
					match = data.(*MatchResult)

					hTarget = wyhash.HashString(match.Target, 1)
					if _, ok = profile[hTarget]; !ok { // skip matches of unwanted targets
						continue
					}

					if prevQuery != match.Query { // new query
						if len(matches) > 1 { // skip uniq match
							hs = hs[:0]
							for h = range matches {
								hs = append(hs, h)
							}

							sorts.Quicksort(Uint64Slice(hs))

							n = len(hs)
							np1 = len(hs) - 1
							for i = 0; i < np1; i++ {
								for j = i + 1; j < n; j++ {
									h1, h2 = hs[i], hs[j]
									if amb, ok = ambMatch[h1]; !ok {
										tmp := make(map[uint64]float64, 128)
										tmp[h2]++ // count of cooccurence
										ambMatch[h1] = tmp
									} else {
										amb[h2]++
									}
								}
							}
						}

						matches = make(map[uint64]*[]*MatchResult)
						pScore = 1024
						nScore = 0
						processThisMatch = true
					} else if keepFullMatch {
						if !processThisMatch {
							prevQuery = match.Query
							continue
						}

						if pScore == 1 && match.QCov < 1 {
							processThisMatch = false

							prevQuery = match.Query
							continue
						}
					} else if keepMainMatch && pScore <= 1 {
						if !processThisMatch {
							prevQuery = match.Query
							continue
						}

						if pScore-match.QCov > maxScoreGap {
							processThisMatch = false

							prevQuery = match.Query
							continue
						}
					}

					if onlyTopNScore {
						if !processThisMatch {
							prevQuery = match.Query
							continue
						}

						if match.QCov < pScore { // match with a smaller score
							nScore++
							if nScore > topNScore {
								processThisMatch = false

								prevQuery = match.Query
								continue
							}
						}
					}

					// just need to collect keys
					matches[hTarget] = nil

					prevQuery = match.Query
					pScore = match.QCov
				}
			}

			if len(matches) > 1 {
				hs = hs[:0]
				for h = range matches {
					hs = append(hs, h)
				}

				sorts.Quicksort(Uint64Slice(hs))

				n = len(hs)
				np1 = len(hs) - 1
				for i = 0; i < np1; i++ {
					for j = i + 1; j < n; j++ {
						h1, h2 = hs[i], hs[j]
						if amb, ok = ambMatch[h1]; !ok {
							tmp := make(map[uint64]float64, 128)
							tmp[h2]++ // count of cooccurence
							ambMatch[h1] = tmp
						} else {
							amb[h2]++
						}
					}
				}
			}
		}

		if opt.Verbose || opt.Log2File {
			log.Infof("  elapsed time: %s", time.Since(timeStart1))
			log.Info()
		}

		// ---------------------------------------------------------------
		// stage 3/4
		if opt.Verbose || opt.Log2File {
			log.Infof("stage 3/4: recounting matches and unique matches")
		}
		timeStart1 = time.Now()

		profile2 := make(map[uint64]*Target, len(profile))

		for _, file := range files {
			if opt.Verbose || opt.Log2File {
				log.Infof("  parsing file: %s", file)
			}

			var matches map[uint64]*[]*MatchResult // target -> match result
			var m *MatchResult
			var ms *[]*MatchResult
			var t *Target
			var ok bool
			var hTarget, h, h1, h2 uint64
			var prevQuery string
			var floatMsSize float64
			var uniqMatch bool
			var first bool
			hss := make([]uint64, 0, 256) // for sorting hash value of reference
			hsm := make([]bool, 0, 256)   // marking hash values to delete
			var n, np1, i, j int
			var match *MatchResult

			onlyTopNScore := topNScore > 0
			var nScore int
			var pScore float64
			var processThisMatch bool

			taxids := make([]uint32, 0, 128)
			var taxid1, taxid2 uint32
			var theSameSpecies bool

			reader, err := breader.NewBufferedReader(file, opt.NumCPUs, chunkSize, fn)
			checkError(err)
			var data interface{}

			matches = make(map[uint64]*[]*MatchResult)
			pScore = 1024
			nScore = 0
			processThisMatch = true
			for chunk := range reader.Ch {
				checkError(chunk.Err)

				for _, data = range chunk.Data {
					match = data.(*MatchResult)

					hTarget = wyhash.HashString(match.Target, 1)
					if _, ok = profile[hTarget]; !ok { // skip matches of unwanted targets
						continue
					}

					if prevQuery != match.Query {
						uniqMatch = false
						if len(matches) > 1 {
							hss = hss[:0]
							hsm = hsm[:0]
							for h = range matches {
								hss = append(hss, h)
								hsm = append(hsm, false)
							}
							sort.Slice(hss, func(i, j int) bool {
								return (*matches[hss[i]])[0].QCov > (*matches[hss[j]])[0].QCov
							})

							n = len(hss)
							np1 = len(hss) - 1
							for i = 0; i < np1; i++ {
								if hsm[i] { // deleted
									continue
								}
								for j = i + 1; j < n; j++ {
									if hsm[j] { // deleted
										continue
									}

									h1, h2 = sortTwoUint64s(hss[i], hss[j]) // sort to extract data from ambMatch

									if profile[hss[i]].SumMatch*(1-minDReadsProp) >= ambMatch[h1][h2] &&
										profile[hss[j]].SumUniqMatch < profile[hss[i]].SumUniqMatch*maxMismatchErr {
										// remove hss[j]
										hsm[j] = true
										// fmt.Println(matches[hss[i]], matches[hss[j]])
									} else if profile[hss[j]].SumMatch*(1-minDReadsProp) >= ambMatch[h1][h2] &&
										profile[hss[i]].SumUniqMatch < profile[hss[j]].SumUniqMatch*maxMismatchErr {
										// remove hss[i]
										hsm[i] = true
										// fmt.Println(matches[hss[j]], matches[hss[i]])
									}
								}
							}
							for i, h = range hss {
								if hsm[i] {
									delete(matches, h)
								}
							}

							if len(matches) > 1 { // redistribute matches
								taxids = taxids[:0]
								if levelSpecies {
									for h, ms = range matches {
										taxid1, ok = taxidMap[(*ms)[0].Target]
										if !ok {
											checkError(fmt.Errorf("unknown taxid for %s, please check taxid mapping file(s)", (*ms)[0].Target))
										}
										taxids = append(taxids, taxid1)
									}
									// LCA
									theSameSpecies = false
									taxid1 = taxids[0]
									for _, taxid2 = range taxids[1:] {
										taxid1 = taxdb.LCA(taxid1, taxid2)
									}
									if taxdb.AtOrBelowRank(taxid1, "species") {
										theSameSpecies = true
									}
								}

								for h, ms = range matches {
									floatMsSize = float64(len(*ms))
									first = true
									for _, m = range *ms {
										if t, ok = profile2[h]; !ok {
											t0 := Target{
												Name:         m.Target,
												GenomeSize:   m.GSize,
												Match:        make([]float64, m.IdxNum),
												UniqMatch:    make([]float64, m.IdxNum),
												UniqMatchHic: make([]float64, m.IdxNum),
												QLen:         make([]float64, m.IdxNum),
												RelDepth:     make([]float64, m.IdxNum),
												// Stats:        stats.NewQuantiler(),
											}
											profile2[h] = &t0
											t = &t0
										}

										if first { // count once
											if levelSpecies && theSameSpecies {
												t.UniqMatch[m.FragIdx] += floatOne / floatMsSize
												if m.QCov >= hicUreadsMinQcov {
													t.UniqMatchHic[m.FragIdx] += floatOne / floatMsSize
												}

												// t.Stats.Add(m.QCov) // the best match on a subject
											}
											t.QLen[m.FragIdx] += float64(m.QLen)

											first = false
										}

										// for a read matching multiple regions of a reference, distribute count to multiple regions,
										// the sum is still the one.
										t.Match[m.FragIdx] += floatOne / floatMsSize
									}
									poolMatchResults.Put(ms)
								}
							} else { // len(matches) == 1
								uniqMatch = true
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
											Name:         m.Target,
											GenomeSize:   m.GSize,
											Match:        make([]float64, m.IdxNum),
											UniqMatch:    make([]float64, m.IdxNum),
											UniqMatchHic: make([]float64, m.IdxNum),
											QLen:         make([]float64, m.IdxNum),
											RelDepth:     make([]float64, m.IdxNum),
											// Stats:        stats.NewQuantiler(),
										}
										profile2[h] = &t0
										t = &t0
									}

									if first { // count once
										if len(matches) == 1 {
											t.UniqMatch[m.FragIdx]++
											if m.QCov >= hicUreadsMinQcov {
												t.UniqMatchHic[m.FragIdx]++
											}

											//	t.Stats.Add(m.QCov) // the best match on a subject
										}
										t.QLen[m.FragIdx] += float64(m.QLen)

										first = false
									}

									// for a read matching multiple regions of a reference, distribute count to multiple regions,
									// the sum is still the one.
									t.Match[m.FragIdx] += floatOne / floatMsSize
								}
								poolMatchResults.Put(ms)
							}
						}

						matches = make(map[uint64]*[]*MatchResult)
						pScore = 1024
						nScore = 0
						processThisMatch = true
					} else if keepFullMatch {
						if !processThisMatch {
							prevQuery = match.Query
							continue
						}

						if pScore == 1 && match.QCov < 1 {
							processThisMatch = false

							prevQuery = match.Query
							continue
						}
					} else if keepMainMatch && pScore <= 1 {
						if !processThisMatch {
							prevQuery = match.Query
							continue
						}

						if pScore-match.QCov > maxScoreGap {
							processThisMatch = false

							prevQuery = match.Query
							continue
						}
					}

					if onlyTopNScore {
						if !processThisMatch {
							prevQuery = match.Query
							continue
						}

						if match.QCov < pScore { // match with a smaller score
							nScore++
							if nScore > topNScore {
								processThisMatch = false

								prevQuery = match.Query
								continue
							}
						}
					}

					if ms, ok = matches[hTarget]; !ok {
						// tmp := []*MatchResult{match}
						tmp := poolMatchResults.Get().(*[]*MatchResult)
						*tmp = (*tmp)[:0]
						*tmp = append(*tmp, match)
						matches[hTarget] = tmp
					} else {
						*ms = append(*ms, match)
					}

					prevQuery = match.Query
					pScore = match.QCov
				}
			}

			uniqMatch = false
			if len(matches) > 1 {
				hss = hss[:0]
				hsm = hsm[:0]
				for h = range matches {
					hss = append(hss, h)
					hsm = append(hsm, false)
				}
				sort.Slice(hss, func(i, j int) bool {
					return (*matches[hss[i]])[0].QCov > (*matches[hss[j]])[0].QCov
				})

				n = len(hss)
				np1 = len(hss) - 1
				for i = 0; i < np1; i++ {
					if hsm[i] { // deleted
						continue
					}
					for j = i + 1; j < n; j++ {
						if hsm[j] { // deleted
							continue
						}

						h1, h2 = sortTwoUint64s(hss[i], hss[j]) // sort to extract data from ambMatch

						if profile[hss[i]].SumMatch*(1-minDReadsProp) >= ambMatch[h1][h2] &&
							profile[hss[j]].SumUniqMatch < profile[hss[i]].SumUniqMatch*maxMismatchErr {
							// remove hss[j]
							hsm[j] = true
							// fmt.Println(matches[hss[i]], matches[hss[j]])
						} else if profile[hss[j]].SumMatch*(1-minDReadsProp) >= ambMatch[h1][h2] &&
							profile[hss[i]].SumUniqMatch < profile[hss[j]].SumUniqMatch*maxMismatchErr {
							// remove hss[i]
							hsm[i] = true
							// fmt.Println(matches[hss[j]], matches[hss[i]])
						}
					}
				}
				for i, h = range hss {
					if hsm[i] {
						delete(matches, h)
					}
				}

				if len(matches) > 1 { // redistribute matches
					taxids = taxids[:0]
					if levelSpecies {
						for h, ms = range matches {
							taxid1, ok = taxidMap[(*ms)[0].Target]
							if !ok {
								checkError(fmt.Errorf("unknown taxid for %s, please check taxid mapping file(s)", (*ms)[0].Target))
							}
							taxids = append(taxids, taxid1)
						}
						// LCA
						theSameSpecies = false
						taxid1 = taxids[0]
						for _, taxid2 = range taxids[1:] {
							taxid1 = taxdb.LCA(taxid1, taxid2)
						}
						if taxdb.AtOrBelowRank(taxid1, "species") {
							theSameSpecies = true
						}
					}

					for h, ms = range matches {
						floatMsSize = float64(len(*ms))
						first = true
						for _, m = range *ms {
							if t, ok = profile2[h]; !ok {
								t0 := Target{
									Name:         m.Target,
									GenomeSize:   m.GSize,
									Match:        make([]float64, m.IdxNum),
									UniqMatch:    make([]float64, m.IdxNum),
									UniqMatchHic: make([]float64, m.IdxNum),
									QLen:         make([]float64, m.IdxNum),
									RelDepth:     make([]float64, m.IdxNum),
									// Stats:        stats.NewQuantiler(),
								}
								profile2[h] = &t0
								t = &t0
							}

							if first { // count once
								if levelSpecies && theSameSpecies {
									t.UniqMatch[m.FragIdx] += floatOne / floatMsSize
									if m.QCov >= hicUreadsMinQcov {
										t.UniqMatchHic[m.FragIdx] += floatOne / floatMsSize
									}

									// t.Stats.Add(m.QCov) // the best match on a subject
								}
								t.QLen[m.FragIdx] += float64(m.QLen)

								first = false
							}

							// for a read matching multiple regions of a reference, distribute count to multiple regions,
							// the sum is still the one.
							t.Match[m.FragIdx] += floatOne / floatMsSize
						}
						poolMatchResults.Put(ms)
					}
				} else { // len(matches) == 1
					uniqMatch = true
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
								Name:         m.Target,
								GenomeSize:   m.GSize,
								Match:        make([]float64, m.IdxNum),
								UniqMatch:    make([]float64, m.IdxNum),
								UniqMatchHic: make([]float64, m.IdxNum),
								QLen:         make([]float64, m.IdxNum),
								RelDepth:     make([]float64, m.IdxNum),
								// Stats:        stats.NewQuantiler(),
							}
							profile2[h] = &t0
							t = &t0
						}

						if first { // count once
							if len(matches) == 1 {
								t.UniqMatch[m.FragIdx]++
								if m.QCov >= hicUreadsMinQcov {
									t.UniqMatchHic[m.FragIdx]++
								}

								// t.Stats.Add(m.QCov) // the best match on a subject
							}
							t.QLen[m.FragIdx] += float64(m.QLen)

							first = false
						}

						// for a read matching multiple regions of a reference, distribute count to multiple regions,
						// the sum is still the one.
						t.Match[m.FragIdx] += floatOne / floatMsSize
					}
					poolMatchResults.Put(ms)
				}
			}
		}

		// --------------------
		// sum up #3

		if debug {
			outfhD.WriteString("\n\n")
			outfhD.WriteString("#------------------ round 2 ------------------\n")
		}

		hs = make([]uint64, 0, len(profile)) // list to delete
		for h, t := range profile2 {
			for _, c1 = range t.UniqMatch {
				t.SumUniqMatch += c1
			}
			if t.SumUniqMatch < minUReads {
				hs = append(hs, h)
				if debug {
					fmt.Fprintf(outfhD, "failed2: %s, %s: %.0f\n",
						t.Name, "no enough unique match", t.SumUniqMatch)
				}
				continue
			}

			for _, c1 = range t.UniqMatchHic {
				t.SumUniqMatchHic += c1
			}
			if t.SumUniqMatchHic < minHicUreads {
				hs = append(hs, h)
				if debug {
					fmt.Fprintf(outfhD, "failed2: %s, %s: %.0f\n",
						t.Name, "no enough high-confidence unique match", t.SumUniqMatchHic)
				}
				continue
			}

			if t.SumUniqMatchHic < HicUreadsMinProp*t.SumUniqMatch {
				hs = append(hs, h)
				if debug {
					fmt.Fprintf(outfhD, "failed2: %s, %s: %.4f (%.0f/%.0f)\n",
						t.Name, "no enough high-confidence unique match proportion", t.SumUniqMatchHic/t.SumUniqMatch, t.SumUniqMatchHic, t.SumUniqMatch)
				}
				continue
			}

			// ----------------------

			for _, c = range t.Match {
				if c >= minReads {
					t.FragsProp++
				}
				t.SumMatch += c
			}
			t.FragsProp = t.FragsProp / float64(len(t.Match))
			if t.FragsProp < minFragsProp {
				hs = append(hs, h)
				if debug {
					fmt.Fprintf(outfhD, "failed2: %s, %s: %.1f %v\n",
						t.Name, "low fragments fraction", t.FragsProp, t.Match)
				}
				continue
			}

			t.Qlens = 0
			for _, c2 = range t.QLen {
				t.Qlens += c2
			}
			for i, c2 := range t.QLen {
				t.RelDepth[i] = c2 / t.Qlens * float64(len(t.QLen))
			}

			_, t.RelDepthStd = MeanStdev(t.RelDepth)
			if t.RelDepthStd > maxFragsDepthStdev {
				hs = append(hs, h)
				if debug {
					fmt.Fprintf(outfhD, "failed2: %s, %s: %f\n",
						t.Name, "high FragsDepthStdev", t.RelDepthStd)
				}
				continue
			}

		}

		for _, h := range hs {
			delete(profile2, h)
		}

		if opt.Verbose || opt.Log2File {
			log.Infof("  number of estimated references: %d", len(profile2))
			log.Infof("  elapsed time: %s", time.Since(timeStart1))
			log.Info()
		}

		// ---------------------------------------------------------------
		// stage 4/4
		if opt.Verbose || opt.Log2File {
			log.Infof("stage 4/4: computing profile")
		}
		timeStart1 = time.Now()

		var outfhB *bufio.Writer
		var gwB io.WriteCloser
		var wB *os.File
		var nB uint64
		if outputBinningResult {
			outfhB, gwB, wB, err = outStream(binningFile, strings.HasSuffix(strings.ToLower(binningFile), ".gz"), opt.CompressionLevel)
			checkError(err)
			defer func() {
				outfhB.Flush()
				if gwB != nil {
					gwB.Close()
				}
				wB.Close()
			}()

			outfhB.WriteString("# This is the bioboxes.org binning output format at\n")
			outfhB.WriteString("# https://github.com/bioboxes/rfc/tree/master/data-format\n")
			outfhB.WriteString("@Version:0.10.0\n")
			outfhB.WriteString(fmt.Sprintf("@SampleID:%s\n", sampleID))
			outfhB.WriteString("@@SEQUENCEID	TAXID\n")
		}

		profile3 := make(map[uint64]*Target, len(profile2))

		var nAssignedReads float64

		for _, file := range files {
			if opt.Verbose || opt.Log2File {
				log.Infof("  parsing file: %s", file)
			}

			var matches map[uint64]*[]*MatchResult // target -> match result
			var m *MatchResult
			var ms *[]*MatchResult
			var t, t1 *Target
			var ok bool
			var hTarget, h uint64
			var prevQuery string
			var floatMsSize float64
			var uniqMatch bool
			var first bool
			var sumUReads, prop float64
			var uregionProp float64
			var match *MatchResult

			onlyTopNScore := topNScore > 0
			var nScore int
			var pScore float64
			var processThisMatch bool

			taxids := make([]uint32, 0, 128)
			var taxid1, taxid2 uint32
			var theSameSpecies bool

			reader, err := breader.NewBufferedReader(file, opt.NumCPUs, chunkSize, fn)
			checkError(err)
			var data interface{}

			matches = make(map[uint64]*[]*MatchResult)
			pScore = 1024
			nScore = 0
			processThisMatch = true
			for chunk := range reader.Ch {
				checkError(chunk.Err)

				for _, data = range chunk.Data {
					match = data.(*MatchResult)

					hTarget = wyhash.HashString(match.Target, 1)
					if _, ok = profile2[hTarget]; !ok { // skip matches of unwanted targets
						continue
					}

					if prevQuery != match.Query {
						nAssignedReads++
						uniqMatch = false
						if len(matches) > 1 { // redistribute matches
							sumUReads = 0

							taxids = taxids[:0]
							for h, ms = range matches {
								// consider unique sequence proportion of references.
								if considerUregionProp {
									if uregionProp, ok = uregionPropMap[profile2[h].Name]; ok {
										sumUReads += profile2[h].SumUniqMatch / uregionProp
									} else {
										sumUReads += profile2[h].SumUniqMatch
									}
								} else {
									sumUReads += profile2[h].SumUniqMatch
								}

								if mappingTaxids {
									taxid1, ok = taxidMap[(*ms)[0].Target]
									if !ok {
										checkError(fmt.Errorf("unknown taxid for %s, please check taxid mapping file(s)", (*ms)[0].Target))
									}
									taxids = append(taxids, taxid1)
								}
							}

							if mappingTaxids {
								// LCA
								theSameSpecies = false
								taxid1 = taxids[0]
								for _, taxid2 = range taxids[1:] {
									taxid1 = taxdb.LCA(taxid1, taxid2)
								}

								if levelSpecies && taxdb.AtOrBelowRank(taxid1, "species") {
									theSameSpecies = true
								}

								if outputBinningResult {
									outfhB.WriteString(fmt.Sprintf("%s\t%d\n", prevQuery, taxid1))
									nB++
								}
							}

							for h, ms = range matches {
								floatMsSize = float64(len(*ms))
								first = true
								t1 = profile2[h]

								// consider unique sequence proportion of references.
								if considerUregionProp {
									if uregionProp, ok = uregionPropMap[t1.Name]; ok {
										prop = t1.SumUniqMatch / uregionProp / sumUReads
									} else {
										prop = t1.SumUniqMatch / sumUReads
									}
								} else {
									prop = t1.SumUniqMatch / sumUReads
								}

								for _, m = range *ms {
									if t, ok = profile3[h]; !ok {
										t0 := Target{
											Name:         m.Target,
											GenomeSize:   m.GSize,
											Match:        make([]float64, m.IdxNum),
											UniqMatch:    make([]float64, m.IdxNum),
											UniqMatchHic: make([]float64, m.IdxNum),
											QLen:         make([]float64, m.IdxNum),
											RelDepth:     make([]float64, m.IdxNum),
											Stats:        stats.NewQuantiler(),
										}
										profile3[h] = &t0
										t = &t0
									}

									if first { // count once
										t.QLen[m.FragIdx] += float64(m.QLen) * prop
										if levelSpecies && theSameSpecies {
											t.Stats.Add(m.QCov) // the best match on a subject
										}
										first = false
									}

									t.Match[m.FragIdx] += prop / floatMsSize

									if levelSpecies && theSameSpecies {
										t.UniqMatch[m.FragIdx] += prop / floatMsSize
										if m.QCov >= hicUreadsMinQcov {
											t.UniqMatchHic[m.FragIdx] += prop / floatMsSize
										}
									}
								}
								poolMatchResults.Put(ms)
							}

							uniqMatch = false
						} else if len(matches) == 1 {
							uniqMatch = true
						} else {
							// should not happen here, but it may happen out the main loop
						}

						if uniqMatch {
							for h, ms = range matches {
								floatMsSize = float64(len(*ms))
								first = true
								for _, m = range *ms {
									if t, ok = profile3[h]; !ok {
										t0 := Target{
											Name:         m.Target,
											GenomeSize:   m.GSize,
											Match:        make([]float64, m.IdxNum),
											UniqMatch:    make([]float64, m.IdxNum),
											UniqMatchHic: make([]float64, m.IdxNum),
											QLen:         make([]float64, m.IdxNum),
											RelDepth:     make([]float64, m.IdxNum),
											Stats:        stats.NewQuantiler(),
										}
										profile3[h] = &t0
										t = &t0
									}

									if first { // count once
										if len(matches) == 1 {
											t.UniqMatch[m.FragIdx]++
											if m.QCov >= hicUreadsMinQcov {
												t.UniqMatchHic[m.FragIdx]++
											}
											t.Stats.Add(m.QCov) // the best match on a subject
										}

										t.QLen[m.FragIdx] += float64(m.QLen)
										first = false

										if outputBinningResult {
											outfhB.WriteString(fmt.Sprintf("%s\t%d\n", prevQuery, taxidMap[m.Target]))
											nB++
										}
									}

									t.Match[m.FragIdx] += floatOne / floatMsSize
								}
								poolMatchResults.Put(ms)
							}
						}

						matches = make(map[uint64]*[]*MatchResult)
						pScore = 1024
						nScore = 0
						processThisMatch = true
					} else if keepFullMatch {
						if !processThisMatch {
							prevQuery = match.Query
							continue
						}

						if pScore == 1 && match.QCov < 1 {
							processThisMatch = false

							prevQuery = match.Query
							continue
						}
					} else if keepMainMatch && pScore <= 1 {
						if !processThisMatch {
							prevQuery = match.Query
							continue
						}

						if pScore-match.QCov > maxScoreGap {
							processThisMatch = false

							prevQuery = match.Query
							continue
						}
					}

					if onlyTopNScore {
						if !processThisMatch {
							prevQuery = match.Query
							continue
						}

						if match.QCov < pScore { // match with a smaller score
							nScore++
							if nScore > topNScore {
								processThisMatch = false

								prevQuery = match.Query
								continue
							}
						}
					}

					if ms, ok = matches[hTarget]; !ok {
						// tmp := []*MatchResult{match}
						tmp := poolMatchResults.Get().(*[]*MatchResult)
						*tmp = (*tmp)[:0]
						*tmp = append(*tmp, match)
						matches[hTarget] = tmp
					} else {
						*ms = append(*ms, match)
					}

					prevQuery = match.Query
					pScore = match.QCov
				}
			}

			nAssignedReads++
			uniqMatch = false
			if len(matches) > 1 { // redistribute matches
				sumUReads = 0

				taxids = taxids[:0]
				for h, ms = range matches {
					// consider unique sequence proportion of references.
					if considerUregionProp {
						if uregionProp, ok = uregionPropMap[profile2[h].Name]; ok {
							sumUReads += profile2[h].SumUniqMatch / uregionProp
						} else {
							sumUReads += profile2[h].SumUniqMatch
						}
					} else {
						sumUReads += profile2[h].SumUniqMatch
					}

					if mappingTaxids {
						taxid1, ok = taxidMap[(*ms)[0].Target]
						if !ok {
							checkError(fmt.Errorf("unknown taxid for %s, please check taxid mapping file(s)", (*ms)[0].Target))
						}
						taxids = append(taxids, taxid1)
					}
				}

				if mappingTaxids {
					// LCA
					theSameSpecies = false
					taxid1 = taxids[0]
					for _, taxid2 = range taxids[1:] {
						taxid1 = taxdb.LCA(taxid1, taxid2)
					}

					if levelSpecies && taxdb.AtOrBelowRank(taxid1, "species") {
						theSameSpecies = true
					}

					if outputBinningResult {
						outfhB.WriteString(fmt.Sprintf("%s\t%d\n", prevQuery, taxid1))
						nB++
					}
				}

				for h, ms = range matches {
					floatMsSize = float64(len(*ms))
					first = true
					t1 = profile2[h]

					// consider unique sequence proportion of references.
					if considerUregionProp {
						if uregionProp, ok = uregionPropMap[t1.Name]; ok {
							prop = t1.SumUniqMatch / uregionProp / sumUReads
						} else {
							prop = t1.SumUniqMatch / sumUReads
						}
					} else {
						prop = t1.SumUniqMatch / sumUReads
					}

					for _, m = range *ms {
						if t, ok = profile3[h]; !ok {
							t0 := Target{
								Name:         m.Target,
								GenomeSize:   m.GSize,
								Match:        make([]float64, m.IdxNum),
								UniqMatch:    make([]float64, m.IdxNum),
								UniqMatchHic: make([]float64, m.IdxNum),
								QLen:         make([]float64, m.IdxNum),
								RelDepth:     make([]float64, m.IdxNum),
								Stats:        stats.NewQuantiler(),
							}
							profile3[h] = &t0
							t = &t0
						}

						if first { // count once
							t.QLen[m.FragIdx] += float64(m.QLen) * prop
							if levelSpecies && theSameSpecies {
								t.Stats.Add(m.QCov) // the best match on a subject
							}
							first = false
						}

						t.Match[m.FragIdx] += prop / floatMsSize

						if levelSpecies && theSameSpecies {
							t.UniqMatch[m.FragIdx] += prop / floatMsSize
							if m.QCov >= hicUreadsMinQcov {
								t.UniqMatchHic[m.FragIdx] += prop / floatMsSize
							}
						}
					}
					poolMatchResults.Put(ms)
				}

				uniqMatch = false
			} else if len(matches) == 1 {
				uniqMatch = true
			} else {
				// should not happen here, but it may happen out the main loop
			}

			if uniqMatch {
				for h, ms = range matches {
					floatMsSize = float64(len(*ms))
					first = true
					for _, m = range *ms {
						if t, ok = profile3[h]; !ok {
							t0 := Target{
								Name:         m.Target,
								GenomeSize:   m.GSize,
								Match:        make([]float64, m.IdxNum),
								UniqMatch:    make([]float64, m.IdxNum),
								UniqMatchHic: make([]float64, m.IdxNum),
								QLen:         make([]float64, m.IdxNum),
								RelDepth:     make([]float64, m.IdxNum),
								Stats:        stats.NewQuantiler(),
							}
							profile3[h] = &t0
							t = &t0
						}

						if first { // count once
							if len(matches) == 1 {
								t.UniqMatch[m.FragIdx]++
								if m.QCov >= hicUreadsMinQcov {
									t.UniqMatchHic[m.FragIdx]++
								}
								t.Stats.Add(m.QCov) // the best match on a subject
							}

							t.QLen[m.FragIdx] += float64(m.QLen)
							first = false

							if outputBinningResult {
								outfhB.WriteString(fmt.Sprintf("%s\t%d\n", prevQuery, taxidMap[m.Target]))
								nB++
							}
						}

						t.Match[m.FragIdx] += floatOne / floatMsSize
					}
					poolMatchResults.Put(ms)
				}
			}

		}

		// --------------------
		// sum up #4

		if debug {
			outfhD.WriteString("\n\n")
			outfhD.WriteString("#------------------ round 3 ------------------\n")
		}

		targets := make([]*Target, 0, 256)

		for _, t := range profile3 {
			for _, c1 = range t.UniqMatch {
				t.SumUniqMatch += c1
			}
			if t.SumUniqMatch < minUReads {
				if debug {
					fmt.Fprintf(outfhD, "failed3: %s, %s: %.0f\n",
						t.Name, "no enough unique match", t.SumUniqMatch)
				}
				continue
			}

			for _, c1 = range t.UniqMatchHic {
				t.SumUniqMatchHic += c1
			}
			if t.SumUniqMatchHic < minHicUreads {
				if debug {
					fmt.Fprintf(outfhD, "failed3: %s, %s: %.0f\n",
						t.Name, "no enough high-confidence unique match", t.SumUniqMatchHic)
				}
				continue
			}

			if t.SumUniqMatchHic < HicUreadsMinProp*t.SumUniqMatch {
				if debug {
					fmt.Fprintf(outfhD, "failed3: %s, %s: %.4f (%.0f/%.0f)\n",
						t.Name, "no enough high-confidence unique match proportion", t.SumUniqMatchHic/t.SumUniqMatch, t.SumUniqMatchHic, t.SumUniqMatch)
				}
				continue
			}

			// ----------------------

			for _, c = range t.Match {
				if c >= minReads {
					t.FragsProp++
				}
				t.SumMatch += c
			}
			t.FragsProp = t.FragsProp / float64(len(t.Match))
			if t.FragsProp < minFragsProp {
				if debug {
					fmt.Fprintf(outfhD, "failed3: %s, %s: %.1f %v\n",
						t.Name, "low fragments fraction", t.FragsProp, t.Match)
				}
				continue
			}

			t.Qlens = 0
			for _, c2 = range t.QLen {
				t.Qlens += c2
			}
			for i, c2 := range t.QLen {
				t.RelDepth[i] = c2 / t.Qlens * float64(len(t.QLen))
			}

			_, t.RelDepthStd = MeanStdev(t.RelDepth)
			if t.RelDepthStd > maxFragsDepthStdev {
				if debug {
					fmt.Fprintf(outfhD, "failed3: %s, %s: %f\n",
						t.Name, "high FragsDepthStdev", t.RelDepthStd)
				}
				continue
			}

			// ----------------------

			switch normAbund {
			case "mean":
				t.Coverage = t.Qlens / float64(t.GenomeSize)
			case "min":
				tmp := math.MaxFloat64
				for _, c2 = range t.QLen {
					if c2 == 0 {
						continue
					}
					if c2 < tmp {
						tmp = c2
					}
				}
				t.Coverage = tmp * float64(len(t.QLen)) / float64(t.GenomeSize)
			case "max":
				var tmp float64
				for _, c2 = range t.QLen {
					if c2 == 0 {
						continue
					}
					if c2 > tmp {
						tmp = c2
					}
				}
				t.Coverage = tmp * float64(len(t.QLen)) / float64(t.GenomeSize)
			}

			// t.Score = similarity(t.Stats.Percentile(90))
			t.Score = t.Stats.Percentile(90) * 100

			targets = append(targets, t)
		}

		if opt.Verbose || opt.Log2File {
			log.Infof("  number of estimated references: %d", len(targets))
			log.Infof("  elapsed time: %s", time.Since(timeStart1))
			log.Info()
			if outputBinningResult {
				log.Infof("%d binning results are save to %s", nB, binningFile)
			}
			log.Info()
			log.Infof("#input matched reads: %.0f, #reads belonging to references in profile: %0.f, proportion: %.6f%%",
				nReads, nAssignedReads, nAssignedReads/nReads*100)
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

		sorts.Quicksort(Targets(targets))

		var totalCoverage float64
		for _, t := range targets {
			totalCoverage += t.Coverage
		}

		for _, t := range targets {
			t.Percentage = t.Coverage / totalCoverage * 100
		}

		if fileterLowAbc && len(targets) > 1 {
			if opt.Verbose || opt.Log2File {
				log.Infof("filtering out predictions with the smallest relative abundances summing up %v%%", lowAbcPct)
			}
			var accPct float64
			var t *Target
			var i, n int

			for i = len(targets) - 1; i >= 0; i-- { // reverse order
				t = targets[i]
				accPct += t.Percentage

				if accPct > lowAbcPct {
					break
				}
				n++
			}

			if n > 0 {
				if opt.Verbose || opt.Log2File {
					log.Infof("  %d targets being filtered out", n)
				}
				targets = targets[:len(targets)-n]

				totalCoverage = 0
				for _, t := range targets {
					totalCoverage += t.Coverage
				}

				for _, t := range targets {
					t.Percentage = t.Coverage / totalCoverage * 100
				}
			} else if opt.Verbose || opt.Log2File {
				log.Infof("no targets being filtered out", n)
			}

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

		outfh.WriteString("ref\tpercentage\tcoverage\tscore\tfragsFrac\tfragsRelDepth\tfragsRelDepthStd\treads\tureads\thicureads\trefsize\trefname\ttaxid\trank\ttaxname\ttaxpath\ttaxpathsn\n")

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
			covs := make([]string, len(t.QLen))
			for i, v := range t.RelDepth {
				covs[i] = fmt.Sprintf("%.2f", v)
			}

			outfh.WriteString(fmt.Sprintf("%s\t%.6f\t%.2f\t%.2f\t%.2f\t%s\t%.2f\t%.0f\t%.0f\t%.0f\t%d\t%s\t%d\t%s\t%s\t%s\t%s\n",
				t.Name, t.Percentage, t.Coverage, t.Score,
				t.FragsProp, strings.Join(covs, ";"), t.RelDepthStd,
				t.SumMatch, t.SumUniqMatch, t.SumUniqMatchHic, t.GenomeSize,
				t.RefName,
				t.Taxid, t.Rank, t.TaxonName,
				strings.Join(t.LineageNames, separator),
				strings.Join(t.LineageTaxids, separator)))
		}

		// ---------------------------------------------------------------
		// more output

		var profile4 map[uint32]*ProfileNode
		var nodes []*ProfileNode

		if outputCamiReport || outputMetaphlanReport {
			profile4 = generateProfile(taxdb, targets)

			nodes = make([]*ProfileNode, 0, len(profile4))
			for _, node := range profile4 {
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

			outfh2.WriteString(fmt.Sprintf("#SampleID\t%s\n", sampleID))

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

			outfh3.WriteString(fmt.Sprintf("@SampleID:%s\n", sampleID))
			outfh3.WriteString("@Version:0.10.0\n")
			outfh3.WriteString("@Ranks:superkingdom|phylum|class|order|family|genus|species|strain\n")
			outfh3.WriteString(fmt.Sprintf("@TaxonomyID:%s\n", taxonomyID))
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
	profileCmd.Flags().Float64P("max-fpr", "f", 0.05, `maximal false positive rate of a read in search result`)
	profileCmd.Flags().Float64P("min-query-cov", "t", 0.55, `minimal query coverage of a read in search result`)
	profileCmd.Flags().IntP("keep-top-qcovs", "n", 0, `keep matches with the top N qcovs for a query, 0 for all`)
	profileCmd.Flags().BoolP("keep-perfect-match", "", false, `only keep the perfect matches (qcov == 1) if there are`)
	profileCmd.Flags().BoolP("keep-main-match", "", false, `only keep main matches, abandon matches with sharply decreased qcov (> --max-qcov-gap)`)
	profileCmd.Flags().Float64P("max-qcov-gap", "", 0.2, `max qcov gap between adjacent matches`)

	// for matches against a reference
	profileCmd.Flags().IntP("min-frags-reads", "r", 50, `minimal number of reads for a reference fragment`)
	profileCmd.Flags().IntP("min-uniq-reads", "u", 20, `minimal number of uniquely matched reads for a reference`)
	profileCmd.Flags().Float64P("min-frags-fraction", "p", 0.8, `minimal fraction of matched reference fragments with reads >= -r/--min-frags-reads`)
	profileCmd.Flags().Float64P("max-frags-depth-stdev", "d", 2, `maximal standard deviation of relative depths of all fragments`)

	profileCmd.Flags().IntP("min-hic-ureads", "U", 5, `minimal number of high-confidence uniquely matched reads for a reference`)
	profileCmd.Flags().Float64P("min-hic-ureads-qcov", "H", 0.75, `minimal query coverage of high-confidence uniquely matched reads`)
	profileCmd.Flags().Float64P("min-hic-ureads-prop", "P", 0.1, `minimal proportion of high-confidence uniquely matched reads`)

	// for the two-stage taxonomy assignment algorithm in MagaPath
	profileCmd.Flags().Float64P("min-dreads-prop", "D", 0.05, `minimal proportion of distinct reads, for determing the right reference for ambiguous reads. Range: (0, 1)`)
	profileCmd.Flags().Float64P("max-mismatch-err", "R", 0.05, `maximal error rate of a read being matched to a wrong reference, for determing the right reference for ambiguous reads. Range: (0, 1)`)

	// name mapping
	profileCmd.Flags().StringSliceP("name-map", "N", []string{}, `tabular two-column file(s) mapping reference IDs to reference names`)

	// taxonomy
	profileCmd.Flags().StringSliceP("taxid-map", "T", []string{}, `tabular two-column file(s) mapping reference IDs to TaxIds`)
	profileCmd.Flags().StringP("taxdump", "X", "", `directory of NCBI taxonomy dump files: names.dmp, nodes.dmp, optional with merged.dmp and delnodes.dmp`)
	profileCmd.Flags().StringP("separator", "S", ";", `separator of TaxIds and taxonomy names`)
	profileCmd.Flags().StringSliceP("show-rank", "", []string{"superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain"}, "only show TaxIds and names of these ranks")
	profileCmd.Flags().StringSliceP("rank-prefix", "", []string{"k__", "p__", "c__", "o__", "f__", "g__", "s__", "t__"}, "prefixes of taxon name in certain ranks, used with --metaphlan-report ")

	// other output formats
	profileCmd.Flags().StringP("sample-id", "s", "", `sample ID in result file`)
	profileCmd.Flags().StringP("taxonomy-id", "", "", `taxonomy ID in result file`)
	profileCmd.Flags().StringP("metaphlan-report", "M", "", `save extra metaphlan-like report`)
	profileCmd.Flags().StringP("cami-report", "C", "", `save extra CAMI-like report`)
	profileCmd.Flags().StringP("binning-result", "B", "", `save extra binning result in CAMI report`)

	profileCmd.Flags().Float64P("filter-low-pct", "F", 0, `filter out predictions with the smallest relative abundances summing up N%. Range: [0,100)`)

	// abundance
	profileCmd.Flags().StringP("norm-abund", "", "mean", `method for normalize abundance of a reference by the mean/min/max abundance in all fragments, available values: mean, min, max`)

	profileCmd.Flags().StringP("level", "", "species", `level to estimate abundance at. available values: species, strain/assembly`)

	// profileCmd.Flags().StringSliceP("uregion-prop-map", "Y", []string{}, `tabular two-column file(s) mapping reference IDs to unique region proportion (experimental)`)

	// profileCmd.Flags().StringP("debug", "", "", `debug output file`)

	profileCmd.Flags().IntP("mode", "m", 3, `profiling mode, type "kmcp profile -h" for details. available values: 0 (for pathogen detection), 1 (highest recall), 2 (high recall), 3 (default), 4 (higher precision), 5 (highest precision)`)

}

// s = lambda qcov: 87.456 + 26.410*qcov - 22.008*qcov*qcov + 7.325*qcov*qcov*qcov
func similarity(qcov float64) float64 {
	square := qcov * qcov
	return 87.456 + 26.410*qcov - -22.008*square + 7.325*square*qcov
}

var poolMatchResults = &sync.Pool{New: func() interface{} {
	tmp := make([]*MatchResult, 0, 128)
	return &tmp
}}
