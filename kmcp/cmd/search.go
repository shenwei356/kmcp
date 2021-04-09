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
	"io"
	"io/ioutil"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/pkg/errors"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/util/cliutil"
	"github.com/shenwei356/util/pathutil"
	"github.com/spf13/cobra"
	"github.com/twotwotwo/sorts/sortutil"
)

var searchCmd = &cobra.Command{
	Use:   "search",
	Short: "Search sequence from database",
	Long: `Search sequence from database

Attentions:
  1. Input format should be (gzipped) FASTA or FASTQ from files or stdin.
  2. Increase value of -j/--threads for acceleratation.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		seq.ValidateSeq = false

		timeStart := time.Now()
		defer func() {
			if opt.Verbose {
				log.Info()
				log.Infof("elapsed time: %s", time.Since(timeStart))
			}
		}()

		var err error

		// ---------------------------------------------------------------

		dbDir := getFlagString(cmd, "db-dir")
		if dbDir == "" {
			checkError(fmt.Errorf("flag -d/--db-dir needed"))
		}
		outFile := getFlagString(cmd, "out-file")
		queryCov := getFlagFloat64(cmd, "query-cov")
		targetCov := getFlagFloat64(cmd, "target-cov")
		minCount := getFlagNonNegativeInt(cmd, "min-count")
		useMmap := !getFlagBool(cmd, "low-mem")
		nameMappingFiles := getFlagStringSlice(cmd, "name-map")
		loadDefaultNameMap := getFlagBool(cmd, "default-name-map")
		keepUnmatched := getFlagBool(cmd, "keep-unmatched")
		// topN := getFlagNonNegativeInt(cmd, "keep-top")
		topN := 0
		topNScore := getFlagNonNegativeInt(cmd, "keep-top-scores")
		noHeaderRow := getFlagBool(cmd, "no-header-row")
		sortBy := getFlagString(cmd, "sort-by")
		// keepOrder := getFlagBool(cmd, "keep-order")
		keepOrder := true

		switch sortBy {
		case "qcov", "jacc", "tcov":
			break
		default:
			checkError(fmt.Errorf("invalid value for flag -s/--sort-by: %s. Available: qcov/tsov/jacc", sortBy))
		}

		if queryCov < 0 || queryCov > 1 {
			checkError(fmt.Errorf("value of -t/--query-cov should be in range [0, 1]"))
		}
		if targetCov < 0 || targetCov > 1 {
			checkError(fmt.Errorf("value of -T/-target-cov should be in range [0, 1]"))
		}

		// ---------------------------------------------------------------
		// check Database

		subFiles, err := ioutil.ReadDir(dbDir)
		if err != nil {
			checkError(fmt.Errorf("read database error: %s", err))
		}

		dbDirs := make([]string, 0, 8)
		for _, file := range subFiles {
			if file.Name() == "." || file.Name() == ".." {
				continue
			}
			path := filepath.Join(dbDir, file.Name())

			if !file.IsDir() {
				continue
			}
			existed, err := pathutil.Exists(filepath.Join(path, dbInfoFile))
			if err != nil {
				checkError(fmt.Errorf("read database error: %s", err))
			}
			if existed {
				dbDirs = append(dbDirs, path)
			}
		}
		if len(dbDirs) == 0 {
			checkError(fmt.Errorf("invalid kmcp database: %s", dbDir))
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
		// input files

		if opt.Verbose {
			log.Info("checking input files ...")
		}
		files := getFileListFromArgsAndFile(cmd, args, true, "infile-list", true)
		if opt.Verbose {
			if len(files) == 1 && isStdin(files[0]) {
				log.Info("  no files given, reading from stdin")
			} else {
				log.Infof("  %d input file(s) given", len(files))
			}
		}

		outFileClean := filepath.Clean(outFile)
		for _, file := range files {
			if !isStdin(file) && filepath.Clean(file) == outFileClean {
				checkError(fmt.Errorf("out file should not be one of the input file"))
			}
		}

		// ---------------------------------------------------------------
		// load db

		if opt.Verbose {
			if useMmap {
				log.Info("loading database with mmap enabled ...")
			} else {
				log.Info("loading database ...")
			}
		}
		searchOpt := SearchOptions{
			UseMMap: useMmap,
			Threads: opt.NumCPUs,

			TopN:       topN,
			TopNScores: topNScore,
			SortBy:     sortBy,

			MinMatched:   minCount,
			MinQueryCov:  queryCov,
			MinTargetCov: targetCov,

			LoadDefaultNameMap: loadDefaultNameMap,
			NameMap:            namesMap,
		}
		sg, err := NewUnikIndexDBSearchEngine(searchOpt, dbDirs...)
		if err != nil {
			checkError(err)
		}

		for _, db := range sg.DBs {
			if queryCov <= db.Info.FPR {
				checkError(fmt.Errorf("query coverage threshold (%f) should not be smaller than FPR of single bloom filter of index database (%f)", queryCov, db.Info.FPR))
			}
		}
		if opt.Verbose {
			log.Infof("database loaded")
			log.Info()
			log.Infof("-------------------- [important parameters] --------------------")
			log.Infof("  minimum  matched k-mers: %d", minCount)
			log.Infof("  minimum  query coverage: %f", queryCov)
			log.Infof("  minimum target coverage: %f", targetCov)
			log.Infof("  minimum target coverage: %f", targetCov)
			log.Infof("-------------------- [important parameters] --------------------")
			log.Info()
			log.Info("searching ...")
		}

		timeStart1 := time.Now()

		outfh, gw, w, err := outStream(outFile, strings.HasSuffix(outFile, ".gz"), opt.CompressionLevel)
		checkError(err)
		defer func() {
			outfh.Flush()
			if gw != nil {
				gw.Close()
			}
			w.Close()
		}()

		if !noHeaderRow {
			outfh.WriteString("#query\tqLen\tqKmers\tFPR\thits\ttarget\tfragIdx\tfrags\ttLen\tmKmers\tqCov\ttCov\tjacc\n")
		}

		var fastxReader *fastx.Reader
		var record *fastx.Record

		var prefix2 string

		// ---------------------------------------------------------------
		// receive result and output

		var total, matched uint64
		var speed float64 // k reads/second

		output := func(result QueryResult) {
			if opt.Verbose {
				total++
			}

			if len(result.Matches) == 0 && !keepUnmatched {
				return
			}

			if opt.Verbose {
				matched++

				if total&8191 == 0 {
					speed = float64(total) / 1000000 / time.Since(timeStart1).Minutes()
					fmt.Fprintf(os.Stderr, "processed queries: %d, speed: %.2f million queries per minute\r", total, speed)
				}
			}

			// query, len_query, num_kmers, fpr, num_matches,
			prefix2 = fmt.Sprintf("%s\t%d\t%d\t%e\t%d",
				result.QueryID, result.QueryLen,
				result.NumKmers, result.FPR, len(result.Matches))

			if keepUnmatched && len(result.Matches) == 0 {
				outfh.WriteString(fmt.Sprintf("%s\t%s\t%d\t%d\t%d\t%0.4f\t%0.4f\n",
					prefix2, "", -1, 0, 0, float64(0), float64(0)))
				return
			}

			for _, match := range result.Matches {
				// query, len_query, num_kmers, fpr, num_matches
				// target, fragIdx, idxNum, tlength, num_matched_kmers, qcov, tcov, jacc
				outfh.WriteString(fmt.Sprintf("%s\t%s\t%d\t%d\t%d\t%d\t%0.4f\t%0.4f\t%0.4f\n",
					prefix2, match.Target[0], uint16(match.TargetIdx[0]), match.TargetIdx[0]>>16, match.GenomeSize[0], match.NumKmers, match.QCov, match.TCov, match.JaccardIndex))
			}

			outfh.Flush()

			result.Recycle()
		}

		done := make(chan int)
		go func() {
			if !keepOrder {
				for result := range sg.OutCh {
					output(result)
				}
			} else {
				m := make(map[uint64]*[]QueryResult, opt.NumCPUs)
				var id, _id uint64
				var ok bool
				var _result QueryResult
				var _results *[]QueryResult
				ndb := len(sg.DBs)
				for result := range sg.OutCh {
					_id = result.QueryIdx
					if _results, ok = m[_id]; ok {
						*_results = append(*_results, result)
					} else {
						m[_id] = &[]QueryResult{result}
					}

					if _results, ok = m[id]; ok && len(*_results) == ndb {
						for _, _result = range *_results {
							output(_result)
						}
						delete(m, id)
						id++
					}
				}

				if len(m) > 0 {
					ids := make([]uint64, len(m))
					i := 0
					for id = range m {
						ids[i] = id
						i++
					}
					sortutil.Uint64s(ids)
					for _, id = range ids {
						_results = m[id]
						for _, _result = range *_results {
							output(_result)
						}
					}
				}
			}
			done <- 1
		}()

		// ---------------------------------------------------------------
		// send query

		var id uint64
		for _, file := range files {
			if opt.Verbose {
				log.Infof("reading sequence file: %s", file)
			}
			fastxReader, err = fastx.NewDefaultReader(file)
			checkError(errors.Wrap(err, file))
			for {
				record, err = fastxReader.Read()
				if err != nil {
					if err == io.EOF {
						break
					}
					errors.Wrap(err, file)
					break
				}

				recordID := make([]byte, len(record.ID))
				copy(recordID, record.ID)

				// may block due to search engine' control
				sg.InCh <- Query{
					Idx: id,
					ID:  recordID,
					Seq: record.Seq.Clone2(),
				}

				id++
			}
		}
		close(sg.InCh) // close Inch

		sg.Wait() // wait all searching finished
		<-done    // all result returned and outputed

		checkError(sg.Close()) // cleanup
		if opt.Verbose {
			fmt.Fprintf(os.Stderr, "\n")

			speed = float64(total) / 1000000 / time.Since(timeStart1).Minutes()
			log.Infof("")
			log.Infof("processed queries: %d, speed: %.2f million queries per minute\n", total, speed)
			log.Infof("%.4f%% (%d/%d) queries matched", float64(matched)/float64(total)*100, matched, total)
			log.Infof("done searching")
		}

	},
}

func init() {
	RootCmd.AddCommand(searchCmd)

	// database option
	searchCmd.Flags().StringP("db-dir", "d", "", `database directory created by "kmcp index"`)
	searchCmd.Flags().BoolP("low-mem", "m", false, `do not load all index files into memory, searching would be very slow`)

	// query option
	searchCmd.Flags().IntP("min-count", "c", 3, `minimal number of matched k-mers (sketch)`)
	searchCmd.Flags().Float64P("query-cov", "t", 0.7, `minimal query coverage, i.e., proportion of matched k-mers and unique k-mers of a query`)
	searchCmd.Flags().Float64P("target-cov", "T", 0, `minimal target coverage, i.e., proportion of matched k-mers and unique k-mers of a target`)

	// output
	searchCmd.Flags().StringP("out-file", "o", "-", `out file, supports and recommends a ".gz" suffix ("-" for stdout)`)
	searchCmd.Flags().StringSliceP("name-map", "M", []string{}, `tabular two-column file(s) mapping names to user-defined values`)
	searchCmd.Flags().BoolP("default-name-map", "D", false, `load ${db}/__name_mapping.tsv for mapping name first`)
	searchCmd.Flags().BoolP("keep-unmatched", "K", false, `keep unmatched query sequence information`)
	// searchCmd.Flags().BoolP("keep-order", "k", false, `keep results in order of input sequences`)
	// searchCmd.Flags().IntP("keep-top", "n", 0, `keep top N hits, 0 for all`)
	searchCmd.Flags().IntP("keep-top-scores", "n", 0, `keep matches with the top N score, 0 for all`)
	searchCmd.Flags().BoolP("no-header-row", "H", false, `do not print header row`)
	searchCmd.Flags().StringP("sort-by", "s", "qcov", `sort hits by "qcov" (Containment Index), "tcov" or "jacc" (Jaccard Index)`)

}
