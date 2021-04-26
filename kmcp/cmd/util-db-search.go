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
	"math"
	"os"
	"path/filepath"
	"sync"

	"github.com/clausecker/pospop"
	"github.com/pkg/errors"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/kmcp/kmcp/cmd/index"
	"github.com/shenwei356/mmap-go"
	"github.com/shenwei356/unikmer"
	"github.com/shenwei356/util/cliutil"
	"github.com/shenwei356/util/pathutil"
	"github.com/twotwotwo/sorts"
	"github.com/twotwotwo/sorts/sortutil"
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

	NumKmers int // number of k-mers
	// Kmers    []uint64 // hashes of k-mers (sketch), for alignment vs target

	Matches []Match // all matches
}

// Name2Idx is a struct of name and index
type Name2Idx struct {
	Name  string
	Index uint32
}

// Match is the struct of matching detail.
type Match struct {
	Target     []string // target name
	TargetIdx  []uint32
	GenomeSize []uint64
	NumKmers   int // matched k-mers

	QCov         float64 // |A∩B|/|A|, coverage of query. i.e., Containment Index
	TCov         float64 // |A∩B|/|B|, coverage of target
	JaccardIndex float64 // |A∩B|/|A∪B|, i.e., JaccardIndex
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

// SortByTCov is used to sort matches by tcov.
type SortByTCov struct{ Matches }

// Less judges if element is i is less than element in j.
func (ms SortByTCov) Less(i int, j int) bool {
	return ms.Matches[i].TCov > ms.Matches[j].TCov && ms.Matches[i].NumKmers > ms.Matches[j].NumKmers
}

// SortByJacc is used to sort matches by jaccard index.
type SortByJacc struct{ Matches }

// Less judges if element is i is less than element in j.
func (ms SortByJacc) Less(i int, j int) bool {
	return ms.Matches[i].JaccardIndex > ms.Matches[j].JaccardIndex && ms.Matches[i].NumKmers > ms.Matches[j].NumKmers
}

// ---------------------------------------------------------------
// messenging between databases and indices

// IndexQuery is a query sent to multiple indices of a database.
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

	DeduplicateThreshold int // deduplicate k-mers only number of kmers > this threshold

	Translated bool
	Frames     int // translate frames

	KeepUnmatched bool
	TopN          int
	TopNScores    int
	SortBy        string

	MinMatched   int
	MinQueryCov  float64
	MinTargetCov float64

	LoadDefaultNameMap bool
	NameMap            map[string]string
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
	multipleDBs := len(dbs) > 1
	mappingName := len(opt.NameMap) > 0

	go func() {
		// have to control maximum concurrence number to prevent memory (goroutine) leak.
		tokens := make(chan int, sg.Options.Threads)
		wg := &sg.wg
		nDBs := len(sg.DBs)
		sortBy := opt.SortBy
		nameMap := opt.NameMap

		topN := opt.TopN
		onlyTopN := topN > 0
		topNScore := opt.TopNScores
		onlyTopNScore := topNScore > 0

		if !multipleDBs {
			handleQuerySingleDB := func(query Query) {
				query.Ch = make(chan QueryResult, nDBs)

				sg.DBs[0].InCh <- query

				_queryResult := <-query.Ch

				if _queryResult.Matches != nil {
					switch sortBy {
					case "qcov":
						sorts.Quicksort(Matches(_queryResult.Matches))
					case "tcov":
						sorts.Quicksort(SortByTCov{Matches(_queryResult.Matches)})
					case "jacc":
						sorts.Quicksort(SortByJacc{Matches(_queryResult.Matches)})
					}

					// filter by scores
					if onlyTopNScore {
						var n, i int
						var score, pScore float64
						pScore = 1024
						var m Match
						for i, m = range _queryResult.Matches {
							switch sortBy {
							case "qcov":
								score = m.QCov
							case "tcov":
								score = m.TCov
							case "jacc":
								score = m.JaccardIndex
							}

							if score < pScore {
								n++
								if n > topNScore {
									break
								}

								pScore = score
							}

						}
						_queryResult.Matches = _queryResult.Matches[:i+1]
					}
				}

				if onlyTopN && len(_queryResult.Matches) > topN {
					_queryResult.Matches = _queryResult.Matches[:topN]
				}

				if mappingName {
					var _m *Match
					var ok bool
					var t string
					_dbInfo := dbs[_queryResult.DBId].Info
					for _, _match := range _queryResult.Matches {
						_m = &_match
						if t, ok = nameMap[_match.Target[0]]; ok {
							_m.Target[0] = t
						} else if opt.LoadDefaultNameMap {
							if t, ok = _dbInfo.NameMapping[_match.Target[0]]; ok {
								_m.Target[0] = t
							}
						}
					}
				}

				sg.OutCh <- _queryResult

				wg.Done()
				<-tokens
			}

			for query := range sg.InCh {
				wg.Add(1)
				tokens <- 1
				go handleQuerySingleDB(query)
			}

			sg.done <- 1

			return
		}

		handleQueryMultiDBs := func(query Query) {
			query.Ch = make(chan QueryResult, nDBs)

			// send to all databases
			for _, db := range sg.DBs {
				db.InCh <- query
			}

			// get matches from all databases
			var queryResult *QueryResult
			var m map[Name2Idx]*Match
			var m2 map[Name2Idx]interface{} // mark shared keys
			// var _match Match
			var _name string
			var key Name2Idx
			var j int
			var ok bool
			var firstDB bool
			var _match0 *Match
			var toDelete []Name2Idx
			toDelete = make([]Name2Idx, 1024)
			for i := 0; i < nDBs; i++ {
				// block to read
				_queryResult := <-query.Ch

				if queryResult == nil { // assign one
					queryResult = &_queryResult
				}

				if _queryResult.Matches == nil {
					continue
				}

				firstDB = i == 0

				if firstDB {
					m = make(map[Name2Idx]*Match, len(_queryResult.Matches))
				} else {
					m2 = make(map[Name2Idx]interface{}, len(_queryResult.Matches))
				}

				for _i := range _queryResult.Matches {
					_match := _queryResult.Matches[_i]
					for j, _name = range _match.Target {
						key = Name2Idx{Name: _name, Index: _match.TargetIdx[j]}
						if firstDB {
							m[key] = &_match
							continue
						}

						if _match0, ok = m[key]; ok { // shared
							if _match.NumKmers < _match0.NumKmers { // update numkmers with smaller value
								m[key] = &_match
							}
							m2[key] = struct{}{} // mark shared keys
						}
					}
				}
				if firstDB {
					continue
				}

				// delete targets not found in previous result
				toDelete = toDelete[:0]
				for key = range m {
					if _, ok = m2[key]; !ok {
						toDelete = append(toDelete, key)
					}
				}
				for _, key = range toDelete {
					delete(m, key)
				}

				// recycle matches
				_queryResult.Matches = _queryResult.Matches[:0]
				poolMatches.Put(_queryResult.Matches)
			}

			_matches2 := poolMatches.Get().([]Match)
			for key, _match := range m {
				_match.Target = []string{key.Name}
				_match.TargetIdx = []uint32{key.Index}
				if multipleDBs {
					_match.TCov = 0
				}
				_matches2 = append(_matches2, *_match)
			}

			switch sortBy {
			case "qcov":
				sorts.Quicksort(Matches(_matches2))
			case "tcov":
				sorts.Quicksort(SortByTCov{Matches(_matches2)})
			case "jacc":
				sorts.Quicksort(SortByJacc{Matches(_matches2)})
			}

			queryResult.Matches = _matches2

			// filter by scores
			if onlyTopNScore {
				var n, i int
				var score, pScore float64
				pScore = 1024
				var m Match
				for i, m = range queryResult.Matches {
					switch sortBy {
					case "qcov":
						score = m.QCov
					case "tcov":
						score = m.TCov
					case "jacc":
						score = m.JaccardIndex
					}

					if score < pScore {
						n++
						if n > topNScore {
							break
						}

						pScore = score
					}

				}
				queryResult.Matches = queryResult.Matches[:i+1]
			}

			if onlyTopN && len(queryResult.Matches) > topN {
				queryResult.Matches = queryResult.Matches[:topN]
			}

			if mappingName {
				var _m *Match
				var ok bool
				var t string
				_dbInfo := dbs[queryResult.DBId].Info
				for _, _match := range queryResult.Matches {
					_m = &_match
					if t, ok = nameMap[_match.Target[0]]; ok {
						_m.Target[0] = t
					} else if opt.LoadDefaultNameMap {
						if t, ok = _dbInfo.NameMapping[_match.Target[0]]; ok {
							_m.Target[0] = t
						}
					}
				}
			}

			sg.OutCh <- *queryResult

			wg.Done()
			<-tokens
		}

		for query := range sg.InCh {
			wg.Add(1)
			tokens <- 1
			go handleQueryMultiDBs(query)
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

	ch := make(chan error)
	done := make(chan int)
	go func() {
		for err := range ch {
			if err != nil && err0 == nil {
				err0 = err
			}
		}
		done <- 1
	}()

	var wg sync.WaitGroup

	for _, db := range sg.DBs {
		wg.Add(1)
		go func(db *UnikIndexDB) {
			ch <- db.Close()
			wg.Done()
		}(db)
	}

	wg.Wait()
	close(ch)
	<-done

	return err0
}

// UnikIndexDB is database for multiple .unik indices.
type UnikIndexDB struct {
	Options SearchOptions
	path    string

	DBId int // id for current database

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

	if opt.LoadDefaultNameMap {
		fileNameMapping := filepath.Join(path, dbNameMappingFile)
		var existed bool
		existed, err = pathutil.Exists(fileNameMapping)
		if existed {
			info.NameMapping, err = cliutil.ReadKVs(fileNameMapping, false)
			checkError(err)
		}
		info.MappingNames = len(info.NameMapping) > 0
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

	db := &UnikIndexDB{Options: opt, Info: info, Header: idx1.Header, path: path}

	db.InCh = make(chan Query, opt.Threads)

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
	}

	db.Indices = indices

	go func() {
		tokens := make(chan int, db.Options.Threads)

		numHashes := db.Info.NumHashes
		indices := db.Indices
		numIndices := len(indices)

		handleQuery := func(query Query) {
			defer func() { <-tokens }()

			// compute kmers
			// reuse []uint64 object, to reduce GC
			kmers := poolKmers.Get().([]uint64)
			kmers, err = db.generateKmers(query.Seq, kmers)
			// buf := poolIdxValues.Get().([]unikmer.IdxValue)
			// kmers, err = db.generateKmers(query.Seq, kmers, buf)
			// kmers, err := db.generateKmers(query.Seq)
			if err != nil {
				checkError(err)
			}
			nKmers := len(kmers)

			// recycle objects
			// buf = buf[:0]
			// poolIdxValues.Put(buf)

			result := QueryResult{
				QueryIdx: query.Idx,
				QueryID:  query.ID,
				QueryLen: query.Seq.Length(),
				NumKmers: nKmers,
				// Kmers:    kmers,
				Matches: nil,
			}

			// sequence shorter than k, or too few k-mer sketchs.
			if kmers == nil || nKmers < db.Options.MinMatched {
				query.Ch <- result // still send result!
				return
			}

			if nKmers > opt.DeduplicateThreshold {
				// map is slower than sorting

				// m := make(map[uint64]interface{}, nKmers)
				// i := 0
				// var ok bool
				// for _, v := range kmers {
				// 	if _, ok = m[v]; !ok {
				// 		m[v] = struct{}{}
				// 		kmers[i] = v
				// 		i++
				// 	}
				// }
				// kmers = kmers[:i]

				sortutil.Uint64s(kmers)
				_kmers := poolKmers.Get().([]uint64)
				var p uint64 = math.MaxUint64
				i := 0
				for _, v := range kmers {
					if p != v {
						_kmers = append(_kmers, v)
						i++
						p = v
					}
				}
				kmers = kmers[:0]
				poolKmers.Put(kmers)
				kmers = _kmers
			}

			// update nKmers
			nKmers = len(kmers)

			result.NumKmers = nKmers

			// compute hashes
			// reuse [][]uint64 object, to reduce GC
			hashes := poolHashes.Get().([][]uint64)
			for _, kmer := range kmers {
				// for kmer := range kmers {
				hashes = append(hashes, hashValues(kmer, numHashes))
			}

			// recycle kmer-sketch ([]uint64) object
			kmers = kmers[:0]
			poolKmers.Put(kmers)

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

			// get matches from all indices
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

			// recycle objects
			poolChanMatches.Put(chMatches)

			hashes = hashes[:0]
			poolHashes.Put(hashes)

			// filter
			if len(matches) > 0 {
				// send result
				result.FPR = maxFPR(db.Info.FPR, opt.MinQueryCov, nKmers)
				result.DBId = db.DBId
				result.Matches = matches
			}

			query.Ch <- result
		}

		for query := range db.InCh {
			tokens <- 1
			go handleQuery(query)
		}
		db.done <- 1
	}()

	return db, nil
}

func (db *UnikIndexDB) generateKmers(sequence *seq.Seq, kmers []uint64) ([]uint64, error) {
	// func (db *UnikIndexDB) generateKmers(sequence *seq.Seq, kmers []uint64, buf []unikmer.IdxValue) ([]uint64, error) {
	// func (db *UnikIndexDB) generateKmers(sequence *seq.Seq) (map[uint64]interface{}, error) {
	// kmers := make(map[uint64]interface{}, 128)
	scaled := db.Info.Scaled
	scale := db.Info.Scale
	maxHash := ^uint64(0)
	if scaled {
		maxHash = uint64(float64(^uint64(0)) / float64(scale))
	}

	var err error
	var iter *unikmer.Iterator
	var sketch *unikmer.Sketch
	var sketchProt *unikmer.ProteinMinimizerSketch
	var flagContinue bool
	var frame int
	var code uint64
	var ok bool

	if db.Info.Protein {
		if db.Options.Translated {
			sequence.Alphabet = seq.Protein
		}
		flagContinue = false
		for _, frame = range transFrames[db.Options.Frames] {
			sketchProt, err = unikmer.NewProteinMinimizerSketch(
				sequence, db.Header.K, db.Info.CodonTable, frame, int(db.Info.MinimizerW))
			if err != nil {
				if err == unikmer.ErrShortSeq {
					flagContinue = true
					break
				} else {
					return nil, err
				}
			}

			for {
				code, ok = sketchProt.Next()
				if !ok {
					break
				}
				if scaled && code > maxHash {
					continue
				}
				kmers = append(kmers, code)
			}

			if flagContinue {
				continue
			}
		}

		return kmers, nil
	}

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
			// kmers[code] = struct{}{}
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
			// kmers[code] = struct{}{}
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
			// kmers[code] = struct{}{}
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

	ch := make(chan error)
	done := make(chan int)
	go func() {
		for err := range ch {
			if err != nil && err0 == nil {
				err0 = err
			}
		}
		done <- 1
	}()

	var wg sync.WaitGroup

	for _, idx := range db.Indices {
		wg.Add(1)
		go func(idx *UnikIndex) {
			ch <- idx.Close()
			wg.Done()
		}(idx)
	}

	wg.Wait()
	close(ch)
	<-done

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
// 64 is the best value for my machine (AMD ryzen 2700X).
const PosPopCountBufSize = 64

// UnikIndex defines a unik index struct.
type UnikIndex struct {
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

	_data [][]uint8
	buffs [][]byte
}

func (idx *UnikIndex) String() string {
	return fmt.Sprintf("%s: %s", idx.Path, idx.Header.String())
}

// NewUnixIndex create a index from file.
func NewUnixIndex(file string, opt SearchOptions) (*UnikIndex, error) {
	fh, err := os.Open(file)
	if err != nil {
		return nil, err
	}

	reader, err := index.NewReader(fh)
	if err != nil {
		return nil, err
	}

	offset, err := fh.Seek(0, 1)
	if err != nil {
		return nil, err
	}

	h := index.Header{}
	h.Version = reader.Version
	h.K = reader.K
	h.Canonical = reader.Canonical
	h.NumHashes = reader.NumHashes
	h.Names = reader.Names
	h.GSizes = reader.GSizes
	h.Indices = reader.Indices
	h.Sizes = reader.Sizes
	h.NumRowBytes = reader.NumRowBytes
	h.NumSigs = reader.NumSigs

	idx := &UnikIndex{Options: opt, Path: file, Header: h, fh: fh, reader: reader, offset0: offset}
	idx.useMmap = opt.UseMMap

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

	idx.buffs = buffs

	idx.moreThanOneHash = reader.NumHashes > 1

	// -------------------------------------------------------

	// receive query and execute
	go func() {

		names := idx.Header.Names
		gsizes := idx.Header.GSizes
		indices := idx.Header.Indices
		numNames := len(idx.Header.Names)
		numRowBytes := idx.Header.NumRowBytes
		numSigs := idx.Header.NumSigs
		numSigsUint := uint64(numSigs)
		numSigsUintM1 := numSigsUint - 1
		offset0 := idx.offset0
		fh := idx.fh
		sizes := idx.Header.Sizes
		sigs := idx.sigsB
		// useMmap := idx.useMmap
		useMmap := opt.UseMMap
		data := idx._data
		moreThanOneHash := idx.moreThanOneHash
		queryCov := idx.Options.MinQueryCov
		targetCov := idx.Options.MinTargetCov
		minMatched := idx.Options.MinMatched
		buffs := idx.buffs

		iLast := numRowBytes - 1
		sizesFloat := make([]float64, len(sizes))
		for i, s := range sizes {
			sizesFloat[i] = float64(s)
		}

		var offset int
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

		var _counts [8]int
		var count int
		var ix8 int
		var k int
		var nHashesTarget, c, t, T float64
		var lastRound bool
		var tmp int

		counts0 := make([][8]int, numRowBytes)
		counts := make([][8]int, numRowBytes)

		buf := make([]byte, PosPopCountBufSize)

		for query := range idx.InCh {
			hashes = query.Hashes
			nHashes = float64(len(hashes))

			// reset counts
			bufIdx = 0
			copy(counts, counts0)

			for _, hs = range hashes {
				if useMmap {
					for i, _h = range hs {
						// loc = int(_h % numSigsUint)
						loc = int(_h & numSigsUintM1) // & X is faster than % X when X is power of 2
						offset = int(offset0 + int64(loc*numRowBytes))

						data[i] = sigs[offset : offset+numRowBytes]
					}
				} else {
					for i, _h = range hs {
						// loc = int(_h % numSigsUint)
						loc = int(_h & numSigsUintM1) // & X is faster than % X when X is power of 2
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

						// for j = 0; j < bufIdx; j++ {
						// 	buf[j] = buffs[j][i]
						// }

						// unroll loop
						buf[0] = buffs[0][i]
						buf[1] = buffs[1][i]
						buf[2] = buffs[2][i]
						buf[3] = buffs[3][i]
						buf[4] = buffs[4][i]
						buf[5] = buffs[5][i]
						buf[6] = buffs[6][i]
						buf[7] = buffs[7][i]
						buf[8] = buffs[8][i]
						buf[9] = buffs[9][i]
						buf[10] = buffs[10][i]
						buf[11] = buffs[11][i]
						buf[12] = buffs[12][i]
						buf[13] = buffs[13][i]
						buf[14] = buffs[14][i]
						buf[15] = buffs[15][i]
						buf[16] = buffs[16][i]
						buf[17] = buffs[17][i]
						buf[18] = buffs[18][i]
						buf[19] = buffs[19][i]
						buf[20] = buffs[20][i]
						buf[21] = buffs[21][i]
						buf[22] = buffs[22][i]
						buf[23] = buffs[23][i]
						buf[24] = buffs[24][i]
						buf[25] = buffs[25][i]
						buf[26] = buffs[26][i]
						buf[27] = buffs[27][i]
						buf[28] = buffs[28][i]
						buf[29] = buffs[29][i]
						buf[30] = buffs[30][i]
						buf[31] = buffs[31][i]
						buf[32] = buffs[32][i]
						buf[33] = buffs[33][i]
						buf[34] = buffs[34][i]
						buf[35] = buffs[35][i]
						buf[36] = buffs[36][i]
						buf[37] = buffs[37][i]
						buf[38] = buffs[38][i]
						buf[39] = buffs[39][i]
						buf[40] = buffs[40][i]
						buf[41] = buffs[41][i]
						buf[42] = buffs[42][i]
						buf[43] = buffs[43][i]
						buf[44] = buffs[44][i]
						buf[45] = buffs[45][i]
						buf[46] = buffs[46][i]
						buf[47] = buffs[47][i]
						buf[48] = buffs[48][i]
						buf[49] = buffs[49][i]
						buf[50] = buffs[50][i]
						buf[51] = buffs[51][i]
						buf[52] = buffs[52][i]
						buf[53] = buffs[53][i]
						buf[54] = buffs[54][i]
						buf[55] = buffs[55][i]
						buf[56] = buffs[56][i]
						buf[57] = buffs[57][i]
						buf[58] = buffs[58][i]
						buf[59] = buffs[59][i]
						buf[60] = buffs[60][i]
						buf[61] = buffs[61][i]
						buf[62] = buffs[62][i]
						buf[63] = buffs[63][i]

						// count
						pospop.Count8(&counts[i], buf)
					}

					bufIdx = 0
				}
			}

			// left data in buffer
			if bufIdx > 0 {
				// transpose
				for i = 0; i < numRowBytes; i++ { // every column in matrix
					// for j = 0; j < bufIdx; j++ {
					// 	buf[j] = buffs[j][i]
					// }

					// unroll loop
					tmp = bufIdx - 8
					for j = 0; j < tmp; j++ {
						buf[j] = buffs[j][i]
						j++
						buf[j] = buffs[j][i]
						j++
						buf[j] = buffs[j][i]
						j++
						buf[j] = buffs[j][i]
						j++
						buf[j] = buffs[j][i]
						j++
						buf[j] = buffs[j][i]
						j++
						buf[j] = buffs[j][i]
						j++
						buf[j] = buffs[j][i]
					}
					for ; j < bufIdx; j++ {
						buf[j] = buffs[j][i]
					}

					// count
					pospop.Count8(&counts[i], buf[:bufIdx])
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

					t = c / nHashes // Containment index
					if t < queryCov {
						continue
					}
					nHashesTarget = sizesFloat[k]
					T = c / nHashesTarget
					if T < targetCov {
						continue
					}

					results = append(results, Match{
						Target:     names[k],
						GenomeSize: gsizes[k],
						TargetIdx:  indices[k],
						NumKmers:   count,
						QCov:       t,
						TCov:       T,

						JaccardIndex: c / (nHashes + nHashesTarget - c), // Jaccard Index
					})
				}
			}

			query.Ch <- results
		}

		idx.done <- 1
	}()
	return idx, nil
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

var poolKmers = &sync.Pool{New: func() interface{} {
	return make([]uint64, 0, 256)
}}

var poolHashes = &sync.Pool{New: func() interface{} {
	return make([][]uint64, 0, 256)
}}

var poolMatches = &sync.Pool{New: func() interface{} {
	return make([]Match, 0, 8)
}}

var poolChanMatch = &sync.Pool{New: func() interface{} {
	return make(chan Match, 8)
}}

var poolChanMatches = &sync.Pool{New: func() interface{} {
	return make(chan []Match, 8)
}}

// Recycle put pooled objects back.
func (r QueryResult) Recycle() {
	r.Matches = r.Matches[:0]
	poolMatches.Put(r.Matches)
}

var poolIdxValues = &sync.Pool{New: func() interface{} {
	return make([]unikmer.IdxValue, 0, 128)
}}

var poolBytes = &sync.Pool{New: func() interface{} {
	return make([]byte, 0, 10<<20)
}}
