// Copyright © 2018-2020 Wei Shen <shenwei356@gmail.com>
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
	"runtime"
	"sync"
	"time"

	"github.com/pkg/errors"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/util/cliutil"
	"github.com/spf13/cobra"
)

var searchCmd = &cobra.Command{
	Use:   "search",
	Short: "Search sequence from database",
	Long: `Search sequence from database

Attentions:
  0. Input format should be (gzipped) FASTA or FASTQ from file or stdin.
  1. Increase value of -j/--threads for acceleratation.
  2. Switch on -m/--use-mmap to load index files into memory to 
     accelerate searching, memory usage is roughly equal to size of index files.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)
		seq.ValidateSeq = false

		timeStart := time.Now()
		defer func() {
			if opt.Verbose {
				log.Infof("elapsed time: %s", time.Since(timeStart))
			}
		}()

		var err error

		// ---------------------------------------------------------------

		dbDirs := getFlagStringSlice(cmd, "db-dir")
		if len(dbDirs) == 0 {
			checkError(fmt.Errorf("flag -d/--db-dir needed"))
		}
		outFile := getFlagString(cmd, "out-prefix")
		queryCov := getFlagFloat64(cmd, "query-cov")
		targetCov := getFlagFloat64(cmd, "target-cov")
		useMmap := getFlagBool(cmd, "use-mmap")
		nameMappingFile := getFlagString(cmd, "name-map")
		// keepUnmatched := getFlagBool(cmd, "keep-unmatched")
		topN := getFlagNonNegativeInt(cmd, "keep-top")
		noHeaderRow := getFlagBool(cmd, "no-header-row")
		sortBy := getFlagString(cmd, "sort-by")
		// circular := getFlagBool(cmd, "circular")
		if !(sortBy == "qcov" || sortBy == "tcov" || sortBy == "sum") {
			checkError(fmt.Errorf("invalid value for flag -s/--sort-by: %s. Available: qcov/tsov/sum", sortBy))
		}

		if queryCov < 0 || queryCov > 1 {
			checkError(fmt.Errorf("value of -t/--query-cov should be in range [0, 1]"))
		}
		if targetCov < 0 || targetCov > 1 {
			checkError(fmt.Errorf("value of -T/-target-cov should be in range [0, 1]"))
		}

		var namesMap map[string]string
		mappingNames := nameMappingFile != ""
		if mappingNames {
			namesMap, err = cliutil.ReadKVs(nameMappingFile, false)
			checkError(errors.Wrap(err, nameMappingFile))
			if opt.Verbose {
				log.Infof("%d pairs of name mapping values loaded", len(namesMap))
			}
		}

		// ---------------------------------------------------------------
		// input files

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

			Align: false,

			TopN:   topN,
			SortBy: sortBy,

			MinQueryCov:  queryCov,
			MinTargetCov: targetCov,
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
			log.Infof("%d databases loaded", len(sg.DBs))

			log.Infof("-------------------- [important parameters] --------------------")
			log.Infof("query coverage threshold: %f", queryCov)
			log.Infof("target coverage threshold: %f", targetCov)
			log.Infof("-------------------- [important parameters] --------------------")
		}
		// if mappingNames {
		// 	var ok bool
		// 	var _n int
		// 	for _, name := range db.Info.Names {
		// 		if _, ok = namesMap[name]; !ok {
		// 			_n++
		// 		}
		// 	}
		// 	if _n > 0 {
		// 		log.Warningf("%d names are not defined in name mapping file: %s", _n, nameMappingFile)
		// 	}
		// }

		// if !isStdout(outFile) {
		// 	outFile += ".txt"
		// }
		outfh, gw, w, err := outStream(outFile, false, opt.CompressionLevel)
		checkError(err)
		defer func() {
			outfh.Flush()
			if gw != nil {
				gw.Close()
			}
			w.Close()
		}()

		if !noHeaderRow {
			outfh.WriteString("query\tqlength\tdb\tqKmers\tFPR\thits\ttarget\ttKmers\tqCov\ttCov\n")
		}

		var fastxReader *fastx.Reader
		var record *fastx.Record
		var ok bool

		var prefix2 string
		var t, target string

		done := make(chan int)
		go func() {
			for result := range sg.OutCh {
				// query, len_query,
				// db, num_kmers, fpr, num_matches,
				prefix2 = fmt.Sprintf("%s\t%d\t%s\t%d\t%e\t%d",
					result.Query.ID, len(result.Query.Seq.Seq),
					result.DBName, result.NumKmers, result.FPR, len(result.Matches))

				for _, match := range result.Matches {

					// if keepUnmatched && len(matched) == 0 {
					// 	outfh.WriteString(fmt.Sprintf("%s\t%s\t%d\t%d\t%d\n",
					// 		prefix2, t, len(kmers), 0, 0))
					// }
					target = match.Target
					if mappingNames {
						if t, ok = namesMap[target]; ok {
							t = target
						}
					} else {
						t = target
					}

					// query, len_query,
					// db, num_kmers, fpr,
					// target, num_target_kmers, qcov, tcov
					outfh.WriteString(fmt.Sprintf("%s\t%s\t%d\t%0.4f\t%0.4f\n",
						prefix2, t, match.Kmers, match.QCov, match.TCov))
				}

			}
			done <- 1
		}()

		var wg sync.WaitGroup
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

				wg.Add(1)
				go func(record *fastx.Record) {
					sg.InCh <- Query{
						ID:  record.ID,
						Seq: record.Seq,
					}
					wg.Done()
				}(record.Clone())
			}
		}
		wg.Wait()      // all sequences being sent
		close(sg.InCh) // close Inch

		sg.Wait() // wait all searching finished
		<-done    // all result returned and outputed

		checkError(sg.Close()) // cleanup
		if opt.Verbose {
			log.Infof("done searching")
		}

	},
}

func init() {
	RootCmd.AddCommand(searchCmd)

	searchCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)

	searchCmd.Flags().StringSliceP("db-dir", "d", []string{}, `database directorys created by "kmcp index"`)
	searchCmd.Flags().StringP("name-map", "M", "", `tabular two-column file mapping names to user-defined values`)
	searchCmd.Flags().BoolP("use-mmap", "m", false, `load index files into memory to accelerate searching (recommended)`)

	// searchCmd.Flags().BoolP("circular", "", false, `query sequence is circular genome`)

	searchCmd.Flags().Float64P("query-cov", "t", 0.6, `query coverage threshold, i.e., proportion of matched k-mers and unique k-mers of a query`)
	searchCmd.Flags().Float64P("target-cov", "T", 0, `target coverage threshold, i.e., proportion of matched k-mers and unique k-mers of a target`)

	searchCmd.Flags().BoolP("keep-unmatched", "k", false, `keep unmatched query sequence information`)
	searchCmd.Flags().IntP("keep-top", "n", 0, `keep top N hits, 0 for all`)
	searchCmd.Flags().BoolP("no-header-row", "H", false, `do not print header row`)
	searchCmd.Flags().StringP("sort-by", "s", "qcov", `sort hits by qcov, tcov or sum (qcov+tcov)`)

}
