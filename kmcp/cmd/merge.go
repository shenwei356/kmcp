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
	"math"
	"os"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/shenwei356/bio/seq"
	"github.com/spf13/cobra"
)

var mergeCmd = &cobra.Command{
	Use:   "merge",
	Short: "Merge search results from multiple databases",
	Long: `Merge search results from multiple databases

Attentions
  0. Input files should contain queryIdx field.
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
		queryIdxField := getFlagPositiveInt(cmd, "queryIdx-field")
		noHeaderRow := getFlagBool(cmd, "no-header-row")

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

		// ---------------------------------------------------------------

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
			return make([]string, 0, 32)
		}}

		m := make([]chan *SearchResult, len(files))
		for i, file := range files {
			ch, err := NewSearchResultParser(file, poolMatches, queryIdxField, 128)
			checkError(err)
			m[i] = ch
		}

		var idx uint64

		var ch chan *SearchResult
		var line string
		var found bool
		var allClosed bool
		var i, imin uint64
		buf := make(map[uint64][]*SearchResult, 256)
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

				if r.QueryIdx != idx { // not found
					if _, ok = buf[r.QueryIdx]; !ok {
						buf[r.QueryIdx] = []*SearchResult{r}
					} else {
						buf[r.QueryIdx] = append(buf[r.QueryIdx], r)
					}
					continue
				}

				found = true
				for _, line = range r.Matches {
					outfh.WriteString(line + "\n")
				}
				poolMatches.Put(r.Matches)
			}

			if allClosed {
				break
			}

			if found {
				// still check the buffer
				if rs, ok := buf[idx]; ok {
					for _, r := range rs {
						for _, line = range r.Matches {
							outfh.WriteString(line + "\n")
						}
						poolMatches.Put(r.Matches)
					}

					delete(buf, idx)
				}

				idx++
				continue
			}

			// not found.
			// find the minimum idx in buffer
			imin = math.MaxUint64
			for i = range buf {
				if i < imin {
					imin = i
				}
			}

			for _, r := range buf[imin] {
				for _, line = range r.Matches {
					outfh.WriteString(line + "\n")
				}
				poolMatches.Put(r.Matches)
			}
			delete(buf, imin)

			idx = imin + 1
		}
		// left data in buffer
		for _, rs := range buf {
			for _, r := range rs {
				for _, line = range r.Matches {
					outfh.WriteString(line + "\n")
				}
				poolMatches.Put(r.Matches)
			}
		}
	},
}

func init() {
	RootCmd.AddCommand(mergeCmd)

	mergeCmd.Flags().StringP("out-file", "o", "-", `out file, supports and recommends a ".gz" suffix ("-" for stdout)`)
	mergeCmd.Flags().IntP("queryIdx-field", "f", 15, `field of queryIdx`)
	mergeCmd.Flags().BoolP("no-header-row", "H", false, `do not print header row`)
}

type SearchResult struct {
	QueryIdx uint64
	Matches  []string
}

func NewSearchResultParser(file string, poolMatches *sync.Pool, numFields int, bufferSize int) (chan *SearchResult, error) {
	ch := make(chan *SearchResult, bufferSize)

	go func() {
		infh, r, _, err := inStream(file)
		checkError(err)

		field := numFields - 1
		items := make([]string, numFields)
		scanner := bufio.NewScanner(infh)

		matches := poolMatches.Get().([]string)

		var idx, i uint64
		var line string
		first := true
		for scanner.Scan() {
			line = scanner.Text()
			if line == "" || line[0] == '#' {
				continue
			}

			stringSplitN(line, "\t", numFields, &items)
			if len(items) < numFields {
				checkError(fmt.Errorf("number of fields (%d) < query index field (%d)", len(items), numFields))
			}

			i, err = strconv.ParseUint(items[field], 10, 64)
			if err != nil {
				checkError(fmt.Errorf("invalid query index at field %d: %s", numFields, items[field]))
			}

			if first {
				idx = i
				first = false
			}

			if i != idx { // new
				ch <- &SearchResult{
					QueryIdx: idx,
					Matches:  matches,
				}
				matches = poolMatches.Get().([]string)
				matches = matches[:0]
				idx = i
			}
			matches = append(matches, line)
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
