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
	"container/heap"
	"fmt"
	"os"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/shenwei356/bio/seq"
	"github.com/spf13/cobra"
	"github.com/twotwotwo/sorts"
)

var mergeCmd = &cobra.Command{
	Use:   "merge",
	Short: "Merge search results from multiple databases",
	Long: `Merge search results from multiple databases

Attentions
  1. Referene IDs should be distinct accross all databases.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		seq.ValidateSeq = false

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

		outFile := getFlagString(cmd, "out-file")
		queryIdxField := getFlagPositiveInt(cmd, "field-queryIdx")
		hitsField := getFlagPositiveInt(cmd, "field-hits")
		noHeaderRow := getFlagBool(cmd, "no-header-row")

		sortBy := getFlagString(cmd, "sort-by")
		switch sortBy {
		case "qcov", "jacc", "tcov":
			break
		default:
			checkError(fmt.Errorf("invalid value for flag -s/--sort-by: %s. Available: qcov/tsov/jacc", sortBy))
		}

		var scoreField int
		switch sortBy {
		case "qcov":
			scoreField = queryIdxField - 3
		case "tcov":
			scoreField = queryIdxField - 2
		case "jacc":
			scoreField = queryIdxField - 1
		default:
			checkError(fmt.Errorf("invalid value for flag -s/--sort-by: %s. Available: qcov/tsov/jacc", sortBy))
		}

		// ---------------------------------------------------------------
		// input files

		if opt.Verbose || opt.Log2File {
			log.Info("checking input files ...")
		}
		files := getFileListFromArgsAndFile(cmd, args, true, "infile-list", true)
		if opt.Verbose || opt.Log2File {
			if len(files) == 1 && isStdin(files[0]) {
				log.Info("  no files given, reading from stdin")
			} else {
				log.Infof("  %d input file(s) given", len(files))
			}
		}

		if len(files) < 2 {
			checkError(fmt.Errorf(">= 2 files needed"))
		}

		// checking duplication
		fileMap := make(map[string]interface{})
		for _, file := range files {
			if _, ok := fileMap[file]; ok {
				checkError(fmt.Errorf("duplicated file: %s", file))
			}
			fileMap[file] = struct{}{}
		}

		// ---------------------------------------------------------------

		if opt.Verbose || opt.Log2File {
			log.Info("merging ...")
		}

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
			outfh.WriteString("#query\tqLen\tqKmers\tFPR\thits\ttarget\tfragIdx\tfrags\ttLen\tkSize\tmKmers\tqCov\ttCov\tjacc\tqueryIdx\n")
		}

		poolMatches := &sync.Pool{New: func() interface{} {
			tmp := make([]*_Match, 0, 32)
			return &tmp
		}}

		poolStrings := &sync.Pool{New: func() interface{} {
			tmp := make([]string, queryIdxField)
			return &tmp
		}}

		chOut := make(chan *SearchResult, 256)
		done := make(chan int)
		go func() {
			var preId uint64
			first := true
			var match *_Match
			buf0 := make([]*_Match, 0, 1024)
			var hits string
			hitsFieldM1 := hitsField - 1
			var item string
			for r := range chOut {
				if first {
					buf0 = append(buf0, *r.Matches...)
					// poolMatches.Put(r.Matches)

					preId = r.QueryIdx
					first = false
					continue
				}

				if r.QueryIdx != preId {
					sorts.Quicksort(_Matches(buf0))
					hits = strconv.Itoa(len(buf0))
					for _, match = range buf0 {
						(*match.Data)[hitsFieldM1] = hits

						for _, item = range *match.Data {
							outfh.WriteString(item)
							outfh.WriteByte('\t')
						}
						outfh.WriteByte('\n')

						// poolStrings.Put(match.Data)
						// pool_Match.Put(match)
					}

					buf0 = buf0[:0]
					buf0 = append(buf0, *r.Matches...)
					// poolMatches.Put(r.Matches)

					preId = r.QueryIdx
					continue
				}

				buf0 = append(buf0, *r.Matches...)
				// poolMatches.Put(r.Matches)
			}

			// last one
			sorts.Quicksort(_Matches(buf0))
			hits = strconv.Itoa(len(buf0))
			for _, match = range buf0 {
				(*match.Data)[hitsFieldM1] = hits

				for _, item = range *match.Data {
					outfh.WriteString(item)
					outfh.WriteByte('\t')
				}
				outfh.WriteByte('\n')

				// poolStrings.Put(match.Data)
				// pool_Match.Put(match)
			}

			done <- 1
		}()

		m := make([]chan *SearchResult, len(files))
		for i, file := range files {
			ch, err := NewSearchResultParser(file, poolStrings, poolMatches, hitsField, scoreField, queryIdxField, 128)
			checkError(err)
			m[i] = ch
		}

		var idx uint64

		var ch chan *SearchResult
		var found bool
		var allClosed bool
		var imin uint64
		var ids *Uint64Slice
		buf := make(map[uint64]*[]*SearchResult, 256)
		var results *[]*SearchResult
		for {
			found = false
			allClosed = true
			for _, ch = range m {
				r, ok := <-ch
				if !ok { // closed channel
					continue
				} else {
					allClosed = false
				}

				if r.QueryIdx != idx { // not found, save to buffer
					if results, ok = buf[r.QueryIdx]; !ok {
						buf[r.QueryIdx] = &[]*SearchResult{r}

						if ids == nil {
							ids = &Uint64Slice{r.QueryIdx}
							heap.Init(ids)
						} else {
							heap.Push(ids, r.QueryIdx)
						}
					} else {
						*results = append(*results, r)
					}
					continue
				}

				found = true
				chOut <- r
			}

			if allClosed {
				break
			}

			if found {
				idx++
				continue
			}

			// not found.
			// find the minimum idx in buffer
			imin = heap.Pop(ids).(uint64)

			for _, r := range *buf[imin] {
				chOut <- r
			}
			delete(buf, imin)

			idx = imin + 1
		}

		// left data in buffer, no need to sort
		for _, results = range buf {
			for _, r := range *results {
				chOut <- r
			}
		}

		close(chOut)
		<-done
	},
}

func init() {
	RootCmd.AddCommand(mergeCmd)

	mergeCmd.Flags().StringP("out-file", "o", "-", `out file, supports and recommends a ".gz" suffix ("-" for stdout)`)
	mergeCmd.Flags().IntP("field-queryIdx", "f", 15, `field of queryIdx`)
	mergeCmd.Flags().IntP("field-hits", "n", 5, `field of hits`)
	mergeCmd.Flags().StringP("sort-by", "s", "qcov", `sort hits by "qcov" (Containment Index), "tcov" or "jacc" (Jaccard Index)`)
	mergeCmd.Flags().BoolP("no-header-row", "H", false, `do not print header row`)
}

// do not use, it's slower
var pool_Match = &sync.Pool{New: func() interface{} {
	return &_Match{}
}}

type _Match struct {
	Data  *[]string
	Score float64

	Hits int
}

type _Matches []*_Match

func (s _Matches) Len() int           { return len(s) }
func (s _Matches) Less(i, j int) bool { return s[i].Score > s[j].Score }
func (s _Matches) Swap(i, j int)      { s[i], s[j] = s[j], s[i] }

type SearchResult struct {
	QueryIdx uint64
	Matches  *[]*_Match
}

func NewSearchResultParser(file string, poolStrings *sync.Pool, poolMatches *sync.Pool, matchesField int, scoreField int, numFields int, bufferSize int) (chan *SearchResult, error) {
	ch := make(chan *SearchResult, bufferSize)

	go func() {
		infh, r, _, err := inStream(file)
		checkError(err)

		field := numFields - 1
		fieldNum := matchesField - 1
		fieldScore := scoreField - 1
		scanner := bufio.NewScanner(infh)

		tmp := make([]*_Match, 0, 32)
		matches := &tmp
		// matches := poolMatches.Get().(*[]*_Match)
		// *matches = (*matches)[:0]

		var idx, i uint64
		var score float64
		var nmatches int
		var line string
		first := true
		for scanner.Scan() {
			line = scanner.Text()
			if line == "" || line[0] == '#' {
				continue
			}

			tmp := make([]string, numFields)
			items := &tmp

			// Strangely, it's much slower (nearly 1/2 speed) in server with Intel CPUs.
			// items := poolStrings.Get().(*[]string)

			stringSplitNByByte(line, '\t', numFields, items)
			if len(*items) < numFields {
				checkError(fmt.Errorf("number of fields (%d) < query index field (%d)", len(*items), numFields))
			}

			i, err = strconv.ParseUint((*items)[field], 10, 64)
			if err != nil {
				checkError(fmt.Errorf("invalid query index at field %d: %s", numFields, (*items)[field]))
			}

			nmatches, err = strconv.Atoi((*items)[fieldNum])
			if err != nil {
				checkError(fmt.Errorf("invalid matches number at field %d: %s", fieldNum, (*items)[fieldNum]))
			}

			score, err = strconv.ParseFloat((*items)[fieldScore], 64)
			if err != nil {
				checkError(fmt.Errorf("failed to parse score: %s", (*items)[fieldScore]))
			}

			if first {
				idx = i
				first = false
				*matches = append(*matches, &_Match{Data: items, Score: score, Hits: nmatches})

				// do not use, it's slower
				// _m := pool_Match.Get().(*_Match)
				// _m.Data = items
				// _m.Score = score
				// _m.Hits = nmatches
				// *matches = append(*matches, _m)
				continue
			}

			if i != idx { // new
				ch <- &SearchResult{
					QueryIdx: idx,
					Matches:  matches,
				}
				tmp := make([]*_Match, 0, 32)
				matches = &tmp
				// matches = poolMatches.Get().(*[]*_Match)
				// *matches = (*matches)[:0]
				idx = i
			}
			*matches = append(*matches, &_Match{Data: items, Score: score, Hits: nmatches})

			// do not use, it's slower
			// _m := pool_Match.Get().(*_Match)
			// _m.Data = items
			// _m.Score = score
			// _m.Hits = nmatches
			// *matches = append(*matches, _m)
		}

		checkError(scanner.Err())
		r.Close()

		ch <- &SearchResult{
			QueryIdx: idx,
			Matches:  matches,
		}
		close(ch)
	}()

	return ch, nil
}
