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
	"math"
	"os"
	"path/filepath"
	"sort"
	"sync"

	"github.com/clausecker/pospop"
	"github.com/edsrzf/mmap-go"
	"github.com/pkg/errors"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/unikmer"
	"github.com/shenwei356/unikmer/index"
)

// ---------------------------------------------------------------
// messenging between search command and search engine

// Query strands for a query sequence
type Query struct {
	ID  []byte
	Seq *seq.Seq

	Ch chan QueryResult // result chanel
}

// QueryResult is search result of a query sequence
type QueryResult struct {
	Query Query

	DBName string
	FPR    float64 // fpr, p is related to database

	NumKmers int      // number of k-mers
	Kmers    []uint64 // hashes of k-mers (sketch), for alignment vs target

	Matches []Match // all matches
}

// Match is the type of matching detail
type Match struct {
	Target string  // target name
	Kmers  int     // matched k-mers
	QCov   float64 // coverage of query
	TCov   float64 // coverage of target
}

// ---------------------------------------------------------------
// messenging between database and indice

// IndexQuery is a query sent to multiple indice of a database
type IndexQuery struct {
	Hashes [][]uint64 // related to database

	Ch chan []Match // result chanel
}

// ---------------------------------------------------------------

// SearchOptions defines options for searching
type SearchOptions struct {
	UseMMap bool
	Threads int

	Align bool

	TopN   int
	SortBy string

	MinQueryCov  float64
	MinTargetCov float64
}

// UnikIndexDBSearchEngine search sequence on multiple database
type UnikIndexDBSearchEngine struct {
	Options SearchOptions

	DBs     []*UnikIndexDB
	DBNames []string

	wg   sync.WaitGroup
	done chan int

	InCh  chan Query // queries
	OutCh chan QueryResult
}

// NewUnikIndexDBSearchEngine returns a search engine based on multiple engines
func NewUnikIndexDBSearchEngine(opt SearchOptions, dbPaths ...string) (*UnikIndexDBSearchEngine, error) {
	dbs := make([]*UnikIndexDB, 0, len(dbPaths))
	names := make([]string, 0, len(dbPaths))
	for _, path := range dbPaths {
		db, err := NewUnikIndexDB(path, opt)
		if err != nil {
			return nil, errors.Wrapf(err, "open kmcp db: %s", path)
		}
		dbs = append(dbs, db)
		names = append(names, filepath.Base(path))
	}

	sg := &UnikIndexDBSearchEngine{Options: opt, DBs: dbs, DBNames: names}
	sg.done = make(chan int)
	sg.InCh = make(chan Query, opt.Threads)
	sg.OutCh = make(chan QueryResult, opt.Threads)

	go func() {
		for query := range sg.InCh {
			sg.wg.Add(1)
			go func(query Query) {
				query.Ch = make(chan QueryResult, len(sg.DBs))

				// send to all databases
				for _, db := range sg.DBs {
					db.InCh <- query
				}

				// get matches from all databases
				var _queryResult QueryResult
				for i := 0; i < len(sg.DBs); i++ {
					// block to read
					_queryResult = <-query.Ch
					if _queryResult.Matches == nil {
						continue
					}

					sg.OutCh <- _queryResult
				}

				sg.wg.Done()
			}(query)
		}
		sg.done <- 1
	}()

	return sg, nil
}

// Wait waits
func (sg *UnikIndexDBSearchEngine) Wait() {
	<-sg.done       //  confirm inCh being closed
	sg.wg.Wait()    // wait all results being sent
	close(sg.OutCh) // close OutCh
}

// Close closes the search engine.
func (sg *UnikIndexDBSearchEngine) Close() error {
	var err0 error
	for _, db := range sg.DBs {
		err := db.Close()
		if err != nil && err0 == nil {
			err0 = err
		}
	}
	return err0
}

// UnikIndexDB is database for multiple .unik indices.
type UnikIndexDB struct {
	Options SearchOptions
	path    string

	InCh chan Query

	wg0        sync.WaitGroup
	stop, done chan int

	Info   UnikIndexDBInfo
	Header index.Header

	Indices []*UnikIndex
}

func (db *UnikIndexDB) String() string {
	return fmt.Sprintf("unikmer index db v%d: #blocksize: %d, #blocks: %d, #%d-mers: %d, #hashes: %d",
		db.Info.Version, db.Info.BlockSize, len(db.Info.Files), db.Header.K, db.Info.Kmers, db.Header.NumHashes)
}

// NewUnikIndexDB opens and read from database directory.
func NewUnikIndexDB(path string, opt SearchOptions) (*UnikIndexDB, error) {
	info, err := UnikIndexDBInfoFromFile(filepath.Join(path, dbInfoFile))
	if err != nil {
		return nil, err
	}

	if len(info.Files) == 0 {
		checkError(fmt.Errorf("no index files"))
	}

	err = info.Check()
	if err != nil {
		return nil, err
	}

	indices := make([]*UnikIndex, 0, len(info.Files))

	// first idx
	idx1, err := NewUnixIndex(filepath.Join(path, info.Files[0]), opt)
	checkError(errors.Wrap(err, filepath.Join(path, info.Files[0])))

	if info.IndexVersion == idx1.Header.Version &&
		info.K == idx1.Header.K &&
		info.Canonical == idx1.Header.Canonical &&
		info.NumHashes == int(idx1.Header.NumHashes) {
	} else {
		checkError(fmt.Errorf("index files not compatible"))
	}

	indices = append(indices, idx1)

	db := &UnikIndexDB{Info: info, Header: idx1.Header, path: path}

	db.InCh = make(chan Query, opt.Threads)

	db.stop = make(chan int)
	db.done = make(chan int)

	if len(info.Files) == 1 {
		db.Indices = indices
		return db, nil
	}

	ch := make(chan *UnikIndex, len(info.Files)-1)
	done := make(chan int)
	go func() {
		for idx := range ch {
			indices = append(indices, idx)
		}
		done <- 1
	}()

	var wg sync.WaitGroup
	for _, f := range info.Files[1:] {
		f = filepath.Join(path, f)

		wg.Add(1)
		go func(f string) {
			defer wg.Done()

			idx, err := NewUnixIndex(f, opt)
			checkError(errors.Wrap(err, f))

			if !idx.Header.Compatible(idx1.Header) {
				checkError(fmt.Errorf("index files not compatible"))
			}

			ch <- idx
		}(f)
	}
	wg.Wait()
	close(ch)
	<-done

	db.Indices = indices

	numHashes := db.Info.NumHashes
	go func() {
	LOOP:
		for {
			select {
			case <-db.stop:
				break LOOP
			case query := <-db.InCh:
				// db.wg0.Add(1)
				go func(query Query) {
					// defer db.wg0.Done()

					// compute kmers
					kmers, err := db.generateKmers(query.Seq)

					if err != nil {
						checkError(err)
					}

					result := QueryResult{
						Query:    query,
						NumKmers: len(kmers),
						Kmers:    kmers,
						Matches:  nil,
					}

					if kmers == nil { // sequence shorter than k
						query.Ch <- result
						return
					}

					// compute hashes
					hashes := make([][]uint64, len(kmers))
					for i, kmer := range kmers {
						hashes[i] = hashValues(kmer, numHashes)
					}

					// send queries
					chMatches := make(chan []Match, len(db.Indices))
					for i := len(db.Indices) - 1; i >= 0; i-- { // start from bigger files
						db.Indices[i].InCh <- IndexQuery{
							Hashes: hashes,
							Ch:     chMatches,
						}
					}

					// get matches from all indice
					matches := make([]Match, 0, 8)
					var _matches []Match
					var _match Match
					for i := 0; i < len(db.Indices); i++ {
						// block to read
						_matches = <-chMatches
						for _, _match = range _matches {
							matches = append(matches, _match)
						}
					}
					// close(chMatches)

					db.filterMatches(matches)

					// send result
					result.FPR = maxFPR(db.Info.FPR, opt.MinQueryCov, len(kmers))
					result.DBName = filepath.Base(path)
					result.Kmers = kmers
					result.Matches = matches
					query.Ch <- result
				}(query)
			}

		}
		db.done <- 1
	}()

	return db, nil
}

func (db *UnikIndexDB) filterMatches(matches []Match) {
	switch db.Options.SortBy {
	case "qcov":
		sort.Slice(matches,
			func(i, j int) bool {
				return matches[i].QCov > matches[j].QCov
			})
	case "tcov":
		sort.Slice(matches,
			func(i, j int) bool {
				return matches[i].TCov > matches[j].TCov
			})
	case "sum":
		sort.Slice(matches,
			func(i, j int) bool {
				return matches[i].QCov+matches[i].TCov > matches[j].QCov+matches[j].TCov
			})
	}

	if db.Options.TopN > 0 && db.Options.TopN < len(matches) {
		matches = matches[0:db.Options.TopN]
	}
}

func (db *UnikIndexDB) generateKmers(sequence *seq.Seq) ([]uint64, error) {
	scaled := db.Info.Scaled
	scale := db.Info.Scale
	maxHash := ^uint64(0)
	if scaled {
		maxHash = uint64(float64(^uint64(0)) / float64(scale))
	}

	var err error
	var iter *unikmer.Iterator
	var sketch *unikmer.Sketch
	var code uint64
	var ok bool

	kmers := make([]uint64, 0, 256)

	// using ntHash
	if db.Info.Syncmer {
		sketch, err = unikmer.NewSyncmerSketch(sequence, db.Header.K, int(db.Info.SyncmerS), false)
	} else if db.Info.Minimizer {
		sketch, err = unikmer.NewMinimizerSketch(sequence, db.Header.K, int(db.Info.MinimizerW), false)
	} else {
		iter, err = unikmer.NewHashIterator(sequence, db.Header.K, db.Header.Canonical, false)
	}
	if err != nil {
		if err == unikmer.ErrShortSeq {
			return nil, nil
		}
		return nil, err
	}

	if db.Info.Syncmer {
		for {
			code, ok = sketch.NextSyncmer()
			if !ok {
				break
			}
			if scaled && code > maxHash {
				continue
			}
			kmers = append(kmers, code)
		}
	} else if db.Info.Minimizer {
		for {
			code, ok = sketch.NextMinimizer()
			if !ok {
				break
			}
			if scaled && code > maxHash {
				continue
			}
			kmers = append(kmers, code)
		}
	} else {
		for {
			code, ok = iter.NextHash()
			if !ok {
				break
			}
			if scaled && code > maxHash {
				continue
			}
			kmers = append(kmers, code)
		}
	}
	return kmers, nil
}

// CompatibleWith has loose restric tions for enabling searching from database of different perameters.
func (db *UnikIndexDB) CompatibleWith(db2 *UnikIndexDB) bool {
	if db.Info.Version == db2.Info.Version &&
		db.Info.IndexVersion == db2.Info.IndexVersion {
		return true
	}
	return false
}

// Close closes database.
func (db *UnikIndexDB) Close() error {
	close(db.stop)
	<-db.done
	close(db.InCh)

	var err0 error
	for _, idx := range db.Indices {
		err := idx.Close()
		if err != nil && err0 == nil {
			err0 = err
		}
	}
	return err0
}

// ------------------------------------------------------------------

// PosPopCountBufSize defines the buffer size of byte slice feeding to
// pospopcount (github.com/clausecker/pospop).
//
// Theoretically, size >240 is better, but in this scenario,
// we need firstly transposing the signature matrix, which is the performance
// bottleneck. Column size of the matrix is fixed, therefore we must control
// the row size to balance time of matrix transposing and popopcount.
//
// 128 is the best value for my machine (AMD ryzen 2700X).
const PosPopCountBufSize = 128

// UnikIndex defines a unik index struct.
type UnikIndex struct {
	Options SearchOptions

	stop, done chan int
	InCh       chan IndexQuery

	Path   string
	Header index.Header

	fh      *os.File
	reader  *index.Reader
	offset0 int64

	// -------------------------------

	moreThanOneHash bool

	useMmap bool
	sigs    mmap.MMap // mapped sigatures
	sigsB   []byte

	_data  [][]uint8
	buffs  [][]byte
	buffsT [][PosPopCountBufSize]byte // cache line size 64
}

func (idx *UnikIndex) String() string {
	return fmt.Sprintf("%s: %s", idx.Path, idx.Header.String())
}

// NewUnixIndex create a index from file.
func NewUnixIndex(file string, opt SearchOptions) (*UnikIndex, error) {
	fh, err := os.Open(file)
	if err != nil {
		return nil, fmt.Errorf("failed to open unikmer index file: %s", file)
	}

	reader, err := index.NewReader(fh)
	if err != nil {
		return nil, fmt.Errorf("failed to read unikmer index file: %s", file)
	}

	offset, err := fh.Seek(0, 1)
	if err != nil {
		return nil, fmt.Errorf("error on seek unikmer index file: %s", file)
	}

	h := index.Header{}
	h.Version = reader.Version
	h.K = reader.K
	h.Canonical = reader.Canonical
	h.NumHashes = reader.NumHashes
	h.Names = reader.Names
	h.Sizes = reader.Sizes
	h.NumRowBytes = reader.NumRowBytes
	h.NumSigs = reader.NumSigs
	idx := &UnikIndex{Options: opt, Path: file, Header: h, fh: fh, reader: reader, offset0: offset}
	idx.useMmap = opt.UseMMap

	idx.stop = make(chan int)
	idx.done = make(chan int)
	idx.InCh = make(chan IndexQuery, opt.Threads)

	if idx.useMmap {
		idx.sigs, err = mmap.Map(fh, mmap.RDONLY, 0)
		if err != nil {
			return nil, err
		}
		idx.sigsB = []byte(idx.sigs)
	}
	idx._data = make([][]uint8, reader.NumHashes)
	for i := 0; i < int(reader.NumHashes); i++ {
		idx._data[i] = make([]byte, reader.NumRowBytes)
	}

	// byte matrix for counting
	buffs := make([][]byte, PosPopCountBufSize)
	for i := 0; i < PosPopCountBufSize; i++ {
		buffs[i] = make([]byte, reader.NumRowBytes)
	}

	// transpose of buffs
	buffsT := make([][PosPopCountBufSize]byte, reader.NumRowBytes)
	for i := 0; i < reader.NumRowBytes; i++ {
		buffsT[i] = [PosPopCountBufSize]byte{}
	}
	idx.buffs = buffs
	idx.buffsT = buffsT

	idx.moreThanOneHash = reader.NumHashes > 1

	go func() {
	LOOP:
		for {
			select {
			case <-idx.stop:
				break LOOP
			case query := <-idx.InCh:
				query.Ch <- idx.Search(query.Hashes)
			}
		}
		idx.done <- 1
	}()
	return idx, nil
}

// Search searches with hashes.
func (idx *UnikIndex) Search(hashes [][]uint64) []Match {
	numNames := len(idx.Header.Names)
	numRowBytes := idx.Header.NumRowBytes
	numSigs := idx.Header.NumSigs
	numSigsInt := uint64(numSigs)
	offset0 := idx.offset0
	fh := idx.fh
	sizes := idx.Header.Sizes
	sigs := idx.sigsB
	useMmap := idx.useMmap
	data := idx._data
	moreThanOneHash := idx.moreThanOneHash
	queryCov := idx.Options.MinQueryCov
	targetCov := idx.Options.MinTargetCov

	counts := make([][8]int, numRowBytes)

	var offset int
	var offset2 int64
	var loc int
	var i, j int
	var hs []uint64
	var row []byte
	var b byte
	var h uint64

	buffs := idx.buffs
	buffsT := idx.buffsT
	bufIdx := 0
	var buf *[PosPopCountBufSize]byte

	for _, hs = range hashes {
		if useMmap {
			for i, h = range hs {
				loc = int(h % numSigsInt)
				offset = int(offset0 + int64(loc*numRowBytes))

				data[i] = sigs[offset : offset+numRowBytes]
			}
		} else {
			for i, h = range hs {
				loc = int(h % numSigsInt)
				offset2 = offset0 + int64(loc*numRowBytes)

				fh.Seek(offset2, 0)
				io.ReadFull(fh, data[i])
			}
		}

		// AND
		var and []byte // must creat a new local variable
		if moreThanOneHash {
			and = make([]byte, numRowBytes) // create new slice to avoid edit original data source
			copy(and, data[0])
			for _, row = range data[1:] {
				for i, b = range row {
					and[i] &= b
				}
			}
		} else if useMmap { // just point to the orginial data (mmaped)
			and = data[0]
		} else { // ！useMmap, where io.ReadFull(fh, data[i])
			and = make([]byte, numRowBytes) // create new slice because we don't want data in buffs point to the same data[0]
			copy(and, data[0])
		}

		// add to buffer for counting
		buffs[bufIdx] = and
		bufIdx++

		if bufIdx == PosPopCountBufSize {
			// transpose
			for i = 0; i < numRowBytes; i++ { // every column in matrix
				buf = &buffsT[i]
				for j = 0; j < PosPopCountBufSize; j++ {
					(*buf)[j] = buffs[j][i]
				}
			}

			// count
			for i = 0; i < numRowBytes; i++ { // every row in transposed matrix
				pospop.Count8(&counts[i], buffsT[i][:])
			}

			bufIdx = 0
		}
	}
	// left data in buffer
	if bufIdx > 0 {
		// transpose
		for i = 0; i < numRowBytes; i++ { // every column in matrix
			buf = &buffsT[i]
			for j = 0; j < bufIdx; j++ {
				(*buf)[j] = buffs[j][i]
			}
		}

		// count
		for i = 0; i < numRowBytes; i++ { // every row in transposed matrix
			pospop.Count8(&counts[i], buffsT[i][0:bufIdx])
		}
	}

	var _counts [8]int
	var count int
	var ix8 int
	var k int

	results := make([]Match, 0, 8)

	var c, t, T float64
	iLast := numRowBytes - 1
	for i, _counts = range counts {
		ix8 = i << 3
		for j = 0; j < 8; j++ {
			count = _counts[7-j] // because count in package pospop is in reversed order
			k = ix8 + j
			if i == iLast && k == numNames {
				break
			}
			if count == 0 {
				continue
			}

			c = float64(count)

			t = c / float64(len(hashes))
			if t < queryCov {
				continue
			}
			T = c / float64(sizes[k])
			if T < targetCov {
				continue
			}

			results = append(results, Match{
				Target: idx.Header.Names[k],
				Kmers:  count,
				QCov:   t,
				TCov:   t,
			})
		}
	}

	return results
}

// Close closes the index.
func (idx *UnikIndex) Close() error {
	close(idx.stop) // send stop signal
	<-idx.done      // confirm stopped
	close(idx.InCh) // close InCh

	if idx.useMmap {
		err := idx.sigs.Unmap()
		if err != nil {
			return err
		}
	}

	return idx.fh.Close()
}

func maxFPR(p float64, k float64, l int) float64 {
	return math.Exp(-float64(l) * (k - p) * (k - p) / 2 / (1 - p))
}
