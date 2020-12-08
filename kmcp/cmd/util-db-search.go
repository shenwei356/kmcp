// Copyright © 2020 Wei Shen <shenwei356@gmail.com>
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
	"sync"

	"github.com/edsrzf/mmap-go"
	"github.com/fuzxxl/pospop"
	"github.com/pkg/errors"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/unikmer"
	"github.com/shenwei356/unikmer/index"
	"github.com/twotwotwo/sorts"
)

// ---------------------------------------------------------------
// messenging between search command and search engine

// Query strands for a query sequence.
type Query struct {
	Idx uint64 // id for keep output in order
	ID  []byte
	Seq *seq.Seq

	Ch chan QueryResult // result chanel
}

// QueryResult is the search result of a query sequence.
type QueryResult struct {
	QueryIdx uint64 // id for keep output in order
	QueryID  []byte
	QueryLen int

	DBId int // id of database, for getting database name with few space

	FPR float64 // fpr, p is related to database

	NumKmers int      // number of k-mers
	Kmers    []uint64 // hashes of k-mers (sketch), for alignment vs target

	Matches []Match // all matches
}

// Clone returns a clone of QueryResult
func (r QueryResult) Clone() QueryResult {
	var matches []Match
	if len(r.Matches) > 0 {
		matches = make([]Match, len(r.Matches))
		for i := 0; i < len(r.Matches); i++ {
			matches[i] = r.Matches[i].Clone()
		}
	}

	return QueryResult{
		QueryIdx: r.QueryIdx,
		QueryID:  r.QueryID,
		QueryLen: r.QueryLen,
		DBId:     r.DBId,
		FPR:      r.FPR,
		NumKmers: r.NumKmers,
		Kmers:    r.Kmers,
		Matches:  matches,
	}
}

// Match is the struct of matching detail.
type Match struct {
	IndexID   int // index file id
	TargetIdx int // target name index, for saving space
	NumKmers  int // mat	ched k-mers

	QCov float64 // coverage of query
	TCov float64 // coverage of target

	PIdt float64 // percent of ident linearly matched kmers
	Loc  int
}

// Clone returns a clone of Match
func (m Match) Clone() Match {
	return Match{
		IndexID:   m.IndexID,
		TargetIdx: m.TargetIdx,
		NumKmers:  m.NumKmers,
		QCov:      m.QCov,
		TCov:      m.TCov,
		PIdt:      m.PIdt,
		Loc:       m.Loc,
	}
}

// Matches is list of Matches, for sorting.
type Matches []Match

// Len returns length of Matches.
func (ms Matches) Len() int { return len(ms) }

// Swap swaps two elements.
func (ms Matches) Swap(i int, j int) { ms[i], ms[j] = ms[j], ms[i] }

// Less judges if element is i is less than element in j.
func (ms Matches) Less(i int, j int) bool {
	return ms[i].QCov > ms[j].QCov && ms[i].NumKmers > ms[j].NumKmers
}

// SortByQCov is used to sort matches by qcov.
type SortByQCov struct{ Matches }

// SortByPIdt is used to sort matches by pindent.
type SortByPIdt struct{ Matches }

// Less judges if element is i is less than element in j.
func (ms SortByPIdt) Less(i int, j int) bool {
	return ms.Matches[i].PIdt > ms.Matches[j].PIdt && ms.Matches[i].NumKmers > ms.Matches[j].NumKmers
}

// SortByTCov is used to sort matches by tcov.
type SortByTCov struct{ Matches }

// Less judges if element is i is less than element in j.
func (ms SortByTCov) Less(i int, j int) bool {
	return ms.Matches[i].TCov > ms.Matches[j].TCov && ms.Matches[i].NumKmers > ms.Matches[j].NumKmers
}

// SortBySum12 is used to sort matches by qcov + pidt.
type SortBySum12 struct{ Matches }

// Less judges if element is i is less than element in j.
func (ms SortBySum12) Less(i int, j int) bool {
	return ms.Matches[i].QCov+ms.Matches[i].PIdt > ms.Matches[j].QCov+ms.Matches[j].PIdt && ms.Matches[i].NumKmers > ms.Matches[j].NumKmers
}

// SortBySum13 is used to sort matches by qcov + tcov.
type SortBySum13 struct{ Matches }

// Less judges if element is i is less than element in j.
func (ms SortBySum13) Less(i int, j int) bool {
	return ms.Matches[i].QCov+ms.Matches[i].TCov > ms.Matches[j].QCov+ms.Matches[j].TCov && ms.Matches[i].NumKmers > ms.Matches[j].NumKmers
}

// SortBySum123 is used to sort matches by qcov + pidt + tcov.
type SortBySum123 struct{ Matches }

// Less judges if element is i is less than element in j.
func (ms SortBySum123) Less(i int, j int) bool {
	return ms.Matches[i].QCov+ms.Matches[i].TCov+ms.Matches[i].TCov > ms.Matches[j].QCov+ms.Matches[j].TCov+ms.Matches[j].TCov && ms.Matches[i].NumKmers > ms.Matches[j].NumKmers
}

// ---------------------------------------------------------------
// messenging between database and indice

// IndexQuery is a query sent to multiple indice of a database.
type IndexQuery struct {
	// Kmers  []uint64
	Hashes [][]uint64 // related to database

	Ch chan []Match // result chanel
}

// ---------------------------------------------------------------

// SearchOptions defines options for searching
type SearchOptions struct {
	UseMMap bool
	Threads int

	// for align
	Align       bool
	SkipAlign   int
	RefuseAlign int
	MinIdentPct float64

	KeepUnmatched bool
	TopN          int
	SortBy        string

	MinMatched   int
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
	for i, path := range dbPaths {
		db, err := NewUnikIndexDB(path, opt, i)
		if err != nil {
			return nil, errors.Wrapf(err, "open kmcp db: %s", path)
		} // for returning target name by (DBId, nameIdx)
		dbs = append(dbs, db)
		names = append(names, filepath.Base(path))
	}

	sg := &UnikIndexDBSearchEngine{Options: opt, DBs: dbs, DBNames: names}
	sg.done = make(chan int)
	sg.InCh = make(chan Query, opt.Threads)
	sg.OutCh = make(chan QueryResult, opt.Threads)

	go func() {
		// have to control maximum concurrence number to prevent memory (goroutine) leak.
		tokens := make(chan int, sg.Options.Threads)
		wg := &sg.wg
		for query := range sg.InCh {
			wg.Add(1)
			tokens <- 1
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
					if _queryResult.Matches == nil && !opt.KeepUnmatched {
						continue
					}

					sg.OutCh <- _queryResult
				}

				wg.Done()
				<-tokens
			}(query)
		}
		sg.done <- 1
	}()

	return sg, nil
}

// Wait waits
func (sg *UnikIndexDBSearchEngine) Wait() {
	<-sg.done       // confirm inCh being closed
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

	DBId          int // id for current database
	CumBlockSizes []int

	InCh chan Query

	wg0        sync.WaitGroup
	stop, done chan int

	Info   UnikIndexDBInfo
	Header index.Header

	Indices []*UnikIndex
}

func (db *UnikIndexDB) String() string {
	return fmt.Sprintf("kmcp database v%d: %s,#blocksize: %d, #blocks: %d, #%d-mers: %d, #hashes: %d",
		db.Info.Version, db.Info.Alias, db.Info.BlockSize, len(db.Info.Files), db.Header.K, db.Info.Kmers, db.Header.NumHashes)
}

// NewUnikIndexDB opens and read from database directory.
func NewUnikIndexDB(path string, opt SearchOptions, dbID int) (*UnikIndexDB, error) {
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

	cumBlockSizes := make([]int, 0, len(info.Files))
	i := 0
	for _, nfile := range info.BlockSizes {
		cumBlockSizes = append(cumBlockSizes, i)
		i += nfile
	}

	unikPaths := make([]string, len(info.Paths))
	for j, p := range info.Paths {
		unikPaths[j] = filepath.Join(path, p)
	}

	indices := make([]*UnikIndex, 0, len(info.Files))

	// first idx
	idx1, err := NewUnixIndex(filepath.Join(path, info.Files[0]), opt, 0)
	checkError(errors.Wrap(err, filepath.Join(path, info.Files[0])))

	if info.IndexVersion == idx1.Header.Version &&
		info.K == idx1.Header.K &&
		info.Canonical == idx1.Header.Canonical &&
		info.NumHashes == int(idx1.Header.NumHashes) {
	} else {
		checkError(fmt.Errorf("index files not compatible"))
	}

	indices = append(indices, idx1)

	db := &UnikIndexDB{Options: opt, Info: info, Header: idx1.Header, path: path}

	db.InCh = make(chan Query, opt.Threads)

	db.CumBlockSizes = cumBlockSizes
	db.DBId = dbID

	db.stop = make(chan int)
	db.done = make(chan int)

	if len(info.Files) > 1 {

		ch := make(chan *UnikIndex, len(info.Files)-1)
		done := make(chan int)
		go func() {
			for idx := range ch {
				indices = append(indices, idx)
			}
			done <- 1
		}()

		var wg sync.WaitGroup
		for i, f := range info.Files[1:] {
			f = filepath.Join(path, f)

			wg.Add(1)
			go func(f string, indexID int) {
				defer wg.Done()

				idx, err := NewUnixIndex(f, opt, indexID)
				checkError(errors.Wrap(err, f))

				if !idx.Header.Compatible(idx1.Header) {
					checkError(fmt.Errorf("index files not compatible"))
				}

				ch <- idx
			}(f, i+1)
		}
		wg.Wait()
		close(ch)
		<-done
	}

	db.Indices = indices

	go func() {

		var kmerLocLists [][]uint64
		if opt.Align {
			kmerLocLists = make([][]uint64, len(db.Info.Names))
		}

		tokens := make(chan int, db.Options.Threads)
		for query := range db.InCh {
			tokens <- 1
			go func(query Query) {
				defer func() { <-tokens }()

				numHashes := db.Info.NumHashes
				indices := db.Indices
				align := opt.Align
				skipAlign := opt.SkipAlign
				refuseAlign := opt.RefuseAlign
				numIndices := len(indices)
				cumBlockSizes = db.CumBlockSizes
				pident := db.Options.MinIdentPct
				threads := db.Options.Threads

				// compute kmers
				// reuse []uint64 object, to reduce GC
				kmers := poolKmers.Get().([]uint64)
				kmers, err = db.generateKmers(query.Seq, kmers)
				// buf := poolIdxValues.Get().([]unikmer.IdxValue)
				// kmers, err = db.generateKmers(query.Seq, kmers, buf)
				if err != nil {
					checkError(err)
				}
				nKmers := len(kmers)
				nKmersF := float64(len(kmers))

				// recycle objects
				// buf = buf[:0]
				// poolIdxValues.Put(buf)

				result := QueryResult{
					QueryIdx: query.Idx,
					QueryID:  query.ID,
					QueryLen: query.Seq.Length(),
					NumKmers: nKmers,
					Kmers:    kmers,
					Matches:  nil,
				}

				// sequence shorter than k, or too few k-mer sketchs.
				if kmers == nil || nKmers < db.Options.MinMatched {
					query.Ch <- result // still send result!
					return
				}

				// compute hashes
				// reuse [][]uint64 object, to reduce GC
				hashes := poolHashes.Get().([][]uint64)
				for _, kmer := range kmers {
					hashes = append(hashes, hashValues(kmer, numHashes))
				}

				// send queries
				// reuse chan []Match object, to reduce GC
				chMatches := poolChanMatches.Get().(chan []Match)
				for i := numIndices - 1; i >= 0; i-- { // start from bigger files
					indices[i].InCh <- IndexQuery{
						// Kmers:  kmers,
						Hashes: hashes,
						Ch:     chMatches,
					}
				}

				// get matches from all indice
				// reuse []Match object
				matches := poolMatches.Get().([]Match)
				var _matches []Match
				var _match Match
				for i := 0; i < numIndices; i++ {
					// block to read
					_matches = <-chMatches
					for _, _match = range _matches {
						matches = append(matches, _match)
					}
				}

				// alignment
				if align {
					if len(matches) >= refuseAlign { // too much matches
						// recycle objects
						poolChanMatches.Put(chMatches)

						hashes = hashes[:0]
						poolHashes.Put(hashes)

						matches = matches[:0]
						poolMatches.Put(matches)

						query.Ch <- result // still send result!
						return
					}

					var wg sync.WaitGroup
					tokens2 := make(chan int, threads)
					done := make(chan int)

					chMatches2 := make(chan Match, threads)
					matches2 := poolMatches.Get().([]Match)

					go func() {
						for _match := range chMatches2 {
							matches2 = append(matches2, _match)
						}
						done <- 1
					}()

					for _, _match = range matches {
						if _match.NumKmers >= skipAlign {
							_match.PIdt = 100
							_match.Loc = -1
							chMatches2 <- _match
							continue
						}

						wg.Add(1)
						tokens2 <- 1
						go func(_match Match) {
							defer func() {
								wg.Done()
								<-tokens2
							}()

							k := cumBlockSizes[_match.IndexID] + _match.TargetIdx

							if kmerLocLists[k] == nil {
								kmerLocLists[k], err = readKmers(unikPaths[k])
								if err != nil {
									checkError(err)
								}
							}
							_m, loc := linearMatched(kmers, kmerLocLists[k])
							m := float64(_m) / nKmersF
							if m < pident {
								kmerLocLists[k] = nil
								return
							}

							_match.PIdt = m
							_match.Loc = loc
							chMatches2 <- _match
						}(_match)

					}
					wg.Wait()
					close(chMatches2)
					<-done

					matches = matches[:0]
					poolMatches.Put(matches)

					matches = matches2
				}

				// recycle objects
				poolChanMatches.Put(chMatches)

				hashes = hashes[:0]
				poolHashes.Put(hashes)

				// sort and filter
				if len(matches) > 0 {
					db.sortAndFilterMatches(matches)

					// send result
					result.FPR = maxFPR(db.Info.FPR, opt.MinQueryCov, nKmers)
					result.DBId = db.DBId
					result.Kmers = kmers
					result.Matches = matches
				}

				query.Ch <- result
			}(query)
		}
		db.done <- 1
	}()

	return db, nil
}

// TargetName return name of a target
func (db *UnikIndexDB) TargetName(indexID int, targetIdx int) string {
	return db.Info.Names[db.CumBlockSizes[indexID]+targetIdx]
}

func (db *UnikIndexDB) sortAndFilterMatches(matches []Match) {
	switch db.Options.SortBy {
	case "qcov":
		sorts.Quicksort(Matches(matches))
	case "pidt":
		sorts.Quicksort(SortByPIdt{Matches(matches)})
	case "tcov":
		sorts.Quicksort(SortByTCov{Matches(matches)})
	case "sum12":
		sorts.Quicksort(SortBySum12{Matches(matches)})
	case "sum13":
		sorts.Quicksort(SortBySum13{Matches(matches)})
	case "sum123":
		sorts.Quicksort(SortBySum123{Matches(matches)})
	}

	if db.Options.TopN > 0 && db.Options.TopN < len(matches) {
		matches = matches[0:db.Options.TopN]
	}
}

func (db *UnikIndexDB) generateKmers(sequence *seq.Seq, kmers []uint64) ([]uint64, error) {
	// func (db *UnikIndexDB) generateKmers(sequence *seq.Seq, kmers []uint64, buf []unikmer.IdxValue) ([]uint64, error) {
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

	// using ntHash
	if db.Info.Syncmer {
		sketch, err = unikmer.NewSyncmerSketch(sequence, db.Header.K, int(db.Info.SyncmerS), false)
		// sketch, err = unikmer.NewSyncmerSketchWithBuffer(sequence, db.Header.K, int(db.Info.SyncmerS), false, buf)
	} else if db.Info.Minimizer {
		sketch, err = unikmer.NewMinimizerSketch(sequence, db.Header.K, int(db.Info.MinimizerW), false)
		// sketch, err = unikmer.NewMinimizerSketchWithBuffer(sequence, db.Header.K, int(db.Info.MinimizerW), false, buf)
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
	close(db.InCh)
	<-db.done

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
	IndexID int // id
	Options SearchOptions

	done chan int
	InCh chan IndexQuery

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
func NewUnixIndex(file string, opt SearchOptions, indexID int) (*UnikIndex, error) {
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
	idx.IndexID = indexID

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

	// -------------------------------------------------------

	// receive query and execute
	go func() {

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
		minMatched := idx.Options.MinMatched
		buffs := idx.buffs
		buffsT := idx.buffsT

		iLast := numRowBytes - 1
		indexID := idx.IndexID

		var _count *[8]int

		var _offset int
		var offset2 int64
		var loc int
		var i, j int
		var hs []uint64
		var row []byte
		var b byte
		var _h uint64
		var hashes [][]uint64
		var nHashes float64
		var bufIdx int
		var buf *[PosPopCountBufSize]byte

		var _counts [8]int
		var count int
		var ix8 int
		var k int
		var c, t, T, m float64
		var lastRound bool

		counts := make([][8]int, numRowBytes)

		for query := range idx.InCh {

			m = -1

			hashes = query.Hashes
			nHashes = float64(len(hashes))

			// reset counts
			bufIdx = 0
			for i := 0; i < numRowBytes; i++ {
				_count = &counts[i]
				_counts = *_count
				(*_count)[0] = 0
				(*_count)[1] = 0
				(*_count)[2] = 0
				(*_count)[3] = 0
				(*_count)[4] = 0
				(*_count)[5] = 0
				(*_count)[6] = 0
				(*_count)[7] = 0
			}

			for _, hs = range hashes {
				if useMmap {
					for i, _h = range hs {
						loc = int(_h % numSigsInt)
						_offset = int(offset0 + int64(loc*numRowBytes))

						data[i] = sigs[_offset : _offset+numRowBytes]
					}
				} else {
					for i, _h = range hs {
						loc = int(_h % numSigsInt)
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

			results := make([]Match, 0, 1)

			for i, _counts = range counts {
				ix8 = i << 3
				lastRound = i == iLast

				for j = 0; j < 8; j++ {
					k = ix8 + j

					if lastRound && k == numNames {
						break
					}

					count = _counts[7-j] // because count in package pospop is in reversed order

					if count < minMatched {
						continue
					}

					c = float64(count)

					t = c / nHashes
					if t < queryCov {
						continue
					}
					T = c / float64(sizes[k])
					if T < targetCov {
						continue
					}

					results = append(results, Match{
						IndexID:   indexID,
						TargetIdx: k,
						NumKmers:  count,
						QCov:      t,
						TCov:      T,
						PIdt:      m,
						Loc:       -2,
					})
				}
			}

			query.Ch <- results
		}

		idx.done <- 1
	}()
	return idx, nil
}

func readKmers(file string) ([]uint64, error) {
	infh, r, _, err := inStream(file)
	if err != nil {
		return nil, err
	}
	defer r.Close()

	reader, err := unikmer.NewReader(infh)
	if err != nil {
		return nil, err
	}
	if reader.IsSorted() {
		return nil, fmt.Errorf(".unik file not supposed to be sorted")
	}

	kmers := make([]uint64, 0, reader.Number)

	var code uint64
	for {
		code, _, err = reader.ReadCodeWithTaxid()
		if err != nil {
			if err == io.EOF {
				break
			}
			checkError(errors.Wrap(err, file))
		}
		kmers = append(kmers, code)
	}
	return kmers, nil
}

type loc2vals [][2]uint64

func (l2v loc2vals) Len() int               { return len(l2v) }
func (l2v loc2vals) Less(i int, j int) bool { return l2v[i][1] < l2v[j][1] }
func (l2v loc2vals) Swap(i int, j int)      { l2v[i], l2v[j] = l2v[j], l2v[i] }

func linearMatched(kmers []uint64, kmers2 []uint64) (int, int) {
	// find a kmer that's in genome
	var kmer, kmer2 uint64
	var i0 int                // location in query kmers
	var j0 int                // location in target kmers2
	locs := make([]int, 0, 2) // kmer at i in kmers matchs kmers at locs in kmers2
	for i0, kmer = range kmers {
		for j0, kmer2 = range kmers2 {
			if kmer == kmer2 {
				locs = append(locs, j0)
			}
		}
		if len(locs) > 0 {
			break
		}
	}

	// find the longest match
	nKmers := len(kmers)
	nKmers2 := len(kmers2)
	max := 0
	maxLoc := -1
	var i int                // location in query kmers
	var j int                // location in target kmers2
	var k int                // a tmp variable
	for _, j0 = range locs { // _loc is location of the kmer in target sequence
		matched := -1

		k = 0
		for i = i0; i >= 0; i-- { // from start to 0 in kmer
			j = j0 - k // locatin in kmer2
			if j < 0 {
				break
			}
			if kmers[i] == kmers2[j] {
				matched++
			}
			k++
		}

		k = 0
		for i = i0; i < nKmers; i++ { // from start to end
			j = j0 + k
			if j >= nKmers2 {
				break
			}
			if kmers[i] == kmers2[j] {
				matched++
			}
			k++
		}

		if matched > max {
			max = matched
			maxLoc = j0
		}
	}

	return max, maxLoc
}

// Close closes the index.
func (idx *UnikIndex) Close() error {
	close(idx.InCh) // close InCh
	<-idx.done      // confirm stopped

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

var poolKmers = &sync.Pool{New: func() interface{} {
	return make([]uint64, 0, 8)
}}

var poolHashes = &sync.Pool{New: func() interface{} {
	return make([][]uint64, 0, 8)
}}

var poolMatches = &sync.Pool{New: func() interface{} {
	return make([]Match, 0, 1)
}}

var poolChanMatch = &sync.Pool{New: func() interface{} {
	return make(chan Match, 8)
}}

var poolChanMatches = &sync.Pool{New: func() interface{} {
	return make(chan []Match, 8)
}}

// Recycle put pooled objects back.
func (r QueryResult) Recycle() {
	// recycle kmer-sketch ([]uint64) object
	r.Kmers = r.Kmers[:0]
	poolKmers.Put(r.Kmers)

	r.Matches = r.Matches[:0]
	poolMatches.Put(r.Matches)
}

var poolIdxValues = &sync.Pool{New: func() interface{} {
	return make([]unikmer.IdxValue, 0, 128)
}}