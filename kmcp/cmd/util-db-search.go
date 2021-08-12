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
	"github.com/shenwei356/pand"
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

	K        int
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
	Hashes  [][]uint64 // related to database
	Hashes1 []uint64

	Ch chan []Match // result chanel
}

// ---------------------------------------------------------------

// SearchOptions defines options for searching
type SearchOptions struct {
	UseMMap bool
	Threads int

	DeduplicateThreshold int // deduplicate k-mers only number of kmers > this threshold

	KeepUnmatched bool
	TopN          int
	TopNScores    int
	SortBy        string
	DoNotSort     bool

	MinQLen      int
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

func channelBuffSize(v int) int {
	return int(float64(v) * 1)
}
func tokenNum(v int) int {
	return int(float64(v) * 1)
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
	sg.InCh = make(chan Query, channelBuffSize(opt.Threads))
	sg.OutCh = make(chan QueryResult, channelBuffSize(opt.Threads))
	multipleDBs := len(dbs) > 1
	mappingName := len(opt.NameMap) > 0

	go func() {
		// have to control maximum concurrence number to prevent memory (goroutine) leak.
		tokens := make(chan int, tokenNum(sg.Options.Threads))
		wg := &sg.wg
		nDBs := len(sg.DBs)
		sortBy := opt.SortBy
		doNotSort := opt.DoNotSort
		nameMap := opt.NameMap

		topN := opt.TopN
		onlyTopN := topN > 0
		topNScore := opt.TopNScores
		onlyTopNScore := topNScore > 0 && !doNotSort

		if !multipleDBs {
			handleQuerySingleDB := func(query Query) {
				query.Ch = make(chan QueryResult, nDBs)

				// send to DB
				sg.DBs[0].InCh <- query

				// wait and receive result from it
				_queryResult := <-query.Ch

				if _queryResult.Matches != nil {
					if len(_queryResult.Matches) > 1 && !doNotSort {
						switch sortBy {
						case "qcov":
							sorts.Quicksort(Matches(_queryResult.Matches))
						case "tcov":
							sorts.Quicksort(SortByTCov{Matches(_queryResult.Matches)})
						case "jacc":
							sorts.Quicksort(SortByJacc{Matches(_queryResult.Matches)})
						}
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
			var queryResult QueryResult
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
			firstDB = true
			var noInter bool
			for i := 0; i < nDBs; i++ {
				// block to read
				_queryResult := <-query.Ch

				if noInter { // skip
					continue
				}

				if firstDB { // assign once
					queryResult = QueryResult{
						QueryIdx: _queryResult.QueryIdx,
						QueryID:  _queryResult.QueryID,
						QueryLen: _queryResult.QueryLen,
						DBId:     _queryResult.DBId,
						FPR:      _queryResult.FPR,
						K:        _queryResult.K,
						NumKmers: _queryResult.NumKmers,
					}
				}

				if _queryResult.Matches == nil {
					continue
				}

				if firstDB {
					m = make(map[Name2Idx]*Match, len(_queryResult.Matches)*nDBs)
				} else {
					m2 = make(map[Name2Idx]interface{}, len(_queryResult.Matches))
				}

				for _, _match := range _queryResult.Matches {
					for j, _name = range _match.Target {
						key = Name2Idx{Name: _name, Index: _match.TargetIdx[j] & 65535}
						if firstDB {
							m[key] = &Match{
								Target:     []string{_match.Target[j]},
								TargetIdx:  []uint32{_match.TargetIdx[j]},
								GenomeSize: []uint64{_match.GenomeSize[j]},
								NumKmers:   _match.NumKmers,

								QCov:         _match.QCov,
								TCov:         _match.TCov,
								JaccardIndex: _match.JaccardIndex,
							}
							continue
						}

						if _match0, ok = m[key]; ok { // shared
							if _match.NumKmers < _match0.NumKmers { // update numkmers with smaller value
								_match0.QCov = _match.QCov
								_match0.TCov = _match.TCov
								_match0.JaccardIndex = _match.JaccardIndex
							}
							m2[key] = struct{}{} // mark shared keys
						}
					}
				}

				// recycle matches
				_queryResult.Matches = _queryResult.Matches[:0]
				poolMatches.Put(_queryResult.Matches)

				if firstDB {
					firstDB = false
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

				if len(m) == 0 {
					noInter = true
				}
			}

			if noInter {
				queryResult.Matches = nil
				sg.OutCh <- queryResult

				wg.Done()
				<-tokens
				return
			}

			_matches2 := poolMatches.Get().([]Match)
			for _, _match := range m {
				_matches2 = append(_matches2, *_match)
			}

			if len(_matches2) > 1 {
				switch sortBy {
				case "qcov":
					sorts.Quicksort(Matches(_matches2))
				case "tcov":
					sorts.Quicksort(SortByTCov{Matches(_matches2)})
				case "jacc":
					sorts.Quicksort(SortByJacc{Matches(_matches2)})
				}
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

			sg.OutCh <- queryResult

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

	// wg0        sync.WaitGroup
	stop, done chan int

	Info   UnikIndexDBInfo
	Header index.Header

	Indices []*UnikIndex
}

func (db *UnikIndexDB) String() string {
	return fmt.Sprintf("kmcp database v%d: name: %s, path: %s, #blocksize: %d, #blocks: %d, #%d-mers: %d, #hashes: %d",
		db.Info.Version, db.Info.Alias, db.path, db.Info.BlockSize, len(db.Info.Files), db.Header.K, db.Info.Kmers, db.Header.NumHashes)
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

	// the first idx
	idx1, err := NewUnixIndex(filepath.Join(path, info.Files[0]), opt)
	checkError(errors.Wrap(err, filepath.Join(path, info.Files[0])))

	if info.IndexVersion == idx1.Header.Version &&
		info.Ks[len(info.Ks)-1] == idx1.Header.K &&
		info.Canonical == idx1.Header.Canonical &&
		info.NumHashes == int(idx1.Header.NumHashes) {
	} else {
		checkError(fmt.Errorf("index files not compatible"))
	}

	indices = append(indices, idx1)

	db := &UnikIndexDB{Options: opt, Info: info, Header: idx1.Header, path: path}

	db.InCh = make(chan Query, channelBuffSize(opt.Threads))

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
		tokens := make(chan int, tokenNum(db.Options.Threads))

		numHashes := db.Info.NumHashes
		singleHash := numHashes == 1
		indices := db.Indices
		numIndices := len(indices)
		ks := db.Info.Ks
		sortutil.Ints(ks)
		minLen := db.Options.MinQLen
		// sort ks in descent order
		for i, j := 0, len(ks)-1; i < j; i, j = i+1, j-1 {
			ks[i], ks[j] = ks[j], ks[i]
		}
		lastIk := len(ks) - 1

		handleQuery := func(query Query) {
			for _ik, k := range ks {
				if len(query.Seq.Seq) < minLen { // skip short query
					query.Ch <- QueryResult{
						QueryIdx: query.Idx,
						QueryID:  query.ID,
						QueryLen: query.Seq.Length(),
						K:        k,
						NumKmers: 0,
						// Kmers:    kmers,
						Matches: nil,
					} // still send result!
					<-tokens
					return
				}

				// compute kmers
				// reuse []uint64 object, to reduce GC
				kmers := poolKmers.Get().([]uint64)
				kmers, err = db.generateKmers(query.Seq, k, kmers)
				if err != nil {
					checkError(err)
				}
				nKmers := len(kmers)

				result := QueryResult{
					QueryIdx: query.Idx,
					QueryID:  query.ID,
					QueryLen: query.Seq.Length(),
					K:        k,
					NumKmers: nKmers,
					// Kmers:    kmers,
					Matches: nil,
				}

				// sequence shorter than k, or too few k-mer sketchs.
				if kmers == nil || nKmers < db.Options.MinMatched {
					query.Ch <- result // still send result!
					<-tokens
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
				var hashes [][]uint64
				if !singleHash {
					hashes = poolHashes.Get().([][]uint64)
					for _, kmer := range kmers {
						// for kmer := range kmers {
						hashes = append(hashes, hashValues(kmer, numHashes))
					}

					// recycle kmer-sketch ([]uint64) object
					kmers = kmers[:0]
					poolKmers.Put(kmers)
				}

				// send queries
				// reuse chan []Match object, to reduce GC
				chMatches := poolChanMatches.Get().(chan []Match)
				if !singleHash {
					for i := numIndices - 1; i >= 0; i-- { // start from bigger files
						indices[i].InCh <- IndexQuery{
							// Kmers:  kmers,
							Hashes: hashes,
							Ch:     chMatches,
						}
					}
				} else {
					for i := numIndices - 1; i >= 0; i-- { // start from bigger files
						indices[i].InCh <- IndexQuery{
							// Kmers:  kmers,
							Hashes1: kmers,
							Ch:      chMatches,
						}
					}
				}

				// get matches from all indices
				// reuse []Match object
				matches := poolMatches.Get().([]Match)
				// var _matches []Match
				for i := 0; i < numIndices; i++ {
					// block to read
					_matches := <-chMatches
					if len(_matches) > 0 {
						matches = append(matches, _matches...)
					}
					_matches = _matches[:0]
					poolMatches.Put(_matches)

				}

				// recycle objects
				poolChanMatches.Put(chMatches)

				if !singleHash {
					hashes = hashes[:0]
					poolHashes.Put(hashes)
				} else {
					kmers = kmers[:0]
					poolKmers.Put(kmers)
				}

				// found
				if len(matches) > 0 {
					// send result
					result.FPR = maxFPR(db.Info.FPR, opt.MinQueryCov, nKmers)
					result.DBId = db.DBId
					result.Matches = matches

					query.Ch <- result
					<-tokens
					return
				}

				// not found, try smaller k
				if _ik == lastIk { // already the smallest k, give up
					query.Ch <- result
					<-tokens
					return
				}
			}
		}

		for query := range db.InCh {
			tokens <- 1
			go handleQuery(query)
		}
		db.done <- 1
	}()

	return db, nil
}

func (db *UnikIndexDB) generateKmers(sequence *seq.Seq, k int, kmers []uint64) ([]uint64, error) {
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
		sketch, err = unikmer.NewSyncmerSketch(sequence, k, int(db.Info.SyncmerS), false)
	} else if db.Info.Minimizer {
		sketch, err = unikmer.NewMinimizerSketch(sequence, k, int(db.Info.MinimizerW), false)
	} else {
		iter, err = unikmer.NewHashIterator(sequence, k, db.Header.Canonical, false)
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
	h.Compact = reader.Compact
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
	idx.InCh = make(chan IndexQuery, channelBuffSize(opt.Threads))

	if idx.useMmap {
		idx.sigs, err = mmap.Map(fh, mmap.RDONLY, 0)
		if err != nil {
			return nil, err
		}
		idx.sigsB = []byte(idx.sigs)
	}

	idx.moreThanOneHash = reader.NumHashes > 1

	// -------------------------------------------------------

	// receive query and execute
	fn := func() {
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
		moreThanOneHash := idx.moreThanOneHash
		moreThanTwoHashes := idx.reader.NumHashes > 2
		queryCov := idx.Options.MinQueryCov
		targetCov := idx.Options.MinTargetCov
		minMatched := idx.Options.MinMatched
		compactSize := idx.Header.Compact

		// bit matrix
		data := make([][]byte, reader.NumHashes)
		for i := 0; i < int(reader.NumHashes); i++ {
			data[i] = make([]byte, numRowBytes)
		}

		// byte matrix for counting
		buffs := make([][]byte, PosPopCountBufSize)
		for i := 0; i < PosPopCountBufSize; i++ {
			buffs[i] = make([]byte, numRowBytes)
		}

		b0 := &buffs[0]
		b1 := &buffs[1]
		b2 := &buffs[2]
		b3 := &buffs[3]
		b4 := &buffs[4]
		b5 := &buffs[5]
		b6 := &buffs[6]
		b7 := &buffs[7]
		b8 := &buffs[8]
		b9 := &buffs[9]
		b10 := &buffs[10]
		b11 := &buffs[11]
		b12 := &buffs[12]
		b13 := &buffs[13]
		b14 := &buffs[14]
		b15 := &buffs[15]
		b16 := &buffs[16]
		b17 := &buffs[17]
		b18 := &buffs[18]
		b19 := &buffs[19]
		b20 := &buffs[20]
		b21 := &buffs[21]
		b22 := &buffs[22]
		b23 := &buffs[23]
		b24 := &buffs[24]
		b25 := &buffs[25]
		b26 := &buffs[26]
		b27 := &buffs[27]
		b28 := &buffs[28]
		b29 := &buffs[29]
		b30 := &buffs[30]
		b31 := &buffs[31]
		b32 := &buffs[32]
		b33 := &buffs[33]
		b34 := &buffs[34]
		b35 := &buffs[35]
		b36 := &buffs[36]
		b37 := &buffs[37]
		b38 := &buffs[38]
		b39 := &buffs[39]
		b40 := &buffs[40]
		b41 := &buffs[41]
		b42 := &buffs[42]
		b43 := &buffs[43]
		b44 := &buffs[44]
		b45 := &buffs[45]
		b46 := &buffs[46]
		b47 := &buffs[47]
		b48 := &buffs[48]
		b49 := &buffs[49]
		b50 := &buffs[50]
		b51 := &buffs[51]
		b52 := &buffs[52]
		b53 := &buffs[53]
		b54 := &buffs[54]
		b55 := &buffs[55]
		b56 := &buffs[56]
		b57 := &buffs[57]
		b58 := &buffs[58]
		b59 := &buffs[59]
		b60 := &buffs[60]
		b61 := &buffs[61]
		b62 := &buffs[62]
		b63 := &buffs[63]

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
		// var b byte
		var _h uint64
		var hashes [][]uint64
		var hashes1 []uint64
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
			if moreThanOneHash {
				hashes = query.Hashes
				nHashes = float64(len(hashes))
			} else {
				hashes1 = query.Hashes1
				nHashes = float64(len(hashes1))
			}

			// reset counts
			bufIdx = 0
			copy(counts, counts0)

			// copy code block to reduce branch testing in loop

			// -------------------------------------------------------------------------
			// useMmap

			if useMmap {
				if compactSize {
					if moreThanOneHash {
						for _, hs = range hashes {
							for i, _h = range hs {
								loc = int(_h % numSigsUint)
								// loc = int(_h & numSigsUintM1) // & X is faster than % X when X is power of 2
								offset = int(offset0 + int64(loc*numRowBytes))

								data[i] = sigs[offset : offset+numRowBytes]
							}

							// first two rows
							pand.AndUnsafe(buffs[bufIdx], data[0], data[1])

							// more rows
							if moreThanTwoHashes {
								for _, row = range data[2:] {
									pand.AndUnsafeInplace(buffs[bufIdx], row)
								}
							}

							// add to buffer for counting
							bufIdx++

							if bufIdx == PosPopCountBufSize {
								// transpose
								for i = 0; i < numRowBytes; i++ { // every column in matrix

									// for j = 0; j < bufIdx; j++ {
									// 	buf[j] = buffs[j][i]
									// }

									// unroll loop
									buf[0] = (*b0)[i]
									buf[1] = (*b1)[i]
									buf[2] = (*b2)[i]
									buf[3] = (*b3)[i]
									buf[4] = (*b4)[i]
									buf[5] = (*b5)[i]
									buf[6] = (*b6)[i]
									buf[7] = (*b7)[i]
									buf[8] = (*b8)[i]
									buf[9] = (*b9)[i]
									buf[10] = (*b10)[i]
									buf[11] = (*b11)[i]
									buf[12] = (*b12)[i]
									buf[13] = (*b13)[i]
									buf[14] = (*b14)[i]
									buf[15] = (*b15)[i]
									buf[16] = (*b16)[i]
									buf[17] = (*b17)[i]
									buf[18] = (*b18)[i]
									buf[19] = (*b19)[i]
									buf[20] = (*b20)[i]
									buf[21] = (*b21)[i]
									buf[22] = (*b22)[i]
									buf[23] = (*b23)[i]
									buf[24] = (*b24)[i]
									buf[25] = (*b25)[i]
									buf[26] = (*b26)[i]
									buf[27] = (*b27)[i]
									buf[28] = (*b28)[i]
									buf[29] = (*b29)[i]
									buf[30] = (*b30)[i]
									buf[31] = (*b31)[i]
									buf[32] = (*b32)[i]
									buf[33] = (*b33)[i]
									buf[34] = (*b34)[i]
									buf[35] = (*b35)[i]
									buf[36] = (*b36)[i]
									buf[37] = (*b37)[i]
									buf[38] = (*b38)[i]
									buf[39] = (*b39)[i]
									buf[40] = (*b40)[i]
									buf[41] = (*b41)[i]
									buf[42] = (*b42)[i]
									buf[43] = (*b43)[i]
									buf[44] = (*b44)[i]
									buf[45] = (*b45)[i]
									buf[46] = (*b46)[i]
									buf[47] = (*b47)[i]
									buf[48] = (*b48)[i]
									buf[49] = (*b49)[i]
									buf[50] = (*b50)[i]
									buf[51] = (*b51)[i]
									buf[52] = (*b52)[i]
									buf[53] = (*b53)[i]
									buf[54] = (*b54)[i]
									buf[55] = (*b55)[i]
									buf[56] = (*b56)[i]
									buf[57] = (*b57)[i]
									buf[58] = (*b58)[i]
									buf[59] = (*b59)[i]
									buf[60] = (*b60)[i]
									buf[61] = (*b61)[i]
									buf[62] = (*b62)[i]
									buf[63] = (*b63)[i]

									// count
									pospop.Count8(&counts[i], buf)
								}

								bufIdx = 0
							}
						}
					} else {
						for _, _h = range hashes1 {
							loc = int(_h % numSigsUint)
							// loc = int(_h & numSigsUintM1) // & X is faster than % X when X is power of 2
							offset = int(offset0 + int64(loc*numRowBytes))

							buffs[bufIdx] = sigs[offset : offset+numRowBytes] // just point to the orginial data (mmaped)

							// add to buffer for counting
							bufIdx++

							if bufIdx == PosPopCountBufSize {
								// transpose
								for i = 0; i < numRowBytes; i++ { // every column in matrix

									// for j = 0; j < bufIdx; j++ {
									// 	buf[j] = buffs[j][i]
									// }

									// unroll loop
									buf[0] = (*b0)[i]
									buf[1] = (*b1)[i]
									buf[2] = (*b2)[i]
									buf[3] = (*b3)[i]
									buf[4] = (*b4)[i]
									buf[5] = (*b5)[i]
									buf[6] = (*b6)[i]
									buf[7] = (*b7)[i]
									buf[8] = (*b8)[i]
									buf[9] = (*b9)[i]
									buf[10] = (*b10)[i]
									buf[11] = (*b11)[i]
									buf[12] = (*b12)[i]
									buf[13] = (*b13)[i]
									buf[14] = (*b14)[i]
									buf[15] = (*b15)[i]
									buf[16] = (*b16)[i]
									buf[17] = (*b17)[i]
									buf[18] = (*b18)[i]
									buf[19] = (*b19)[i]
									buf[20] = (*b20)[i]
									buf[21] = (*b21)[i]
									buf[22] = (*b22)[i]
									buf[23] = (*b23)[i]
									buf[24] = (*b24)[i]
									buf[25] = (*b25)[i]
									buf[26] = (*b26)[i]
									buf[27] = (*b27)[i]
									buf[28] = (*b28)[i]
									buf[29] = (*b29)[i]
									buf[30] = (*b30)[i]
									buf[31] = (*b31)[i]
									buf[32] = (*b32)[i]
									buf[33] = (*b33)[i]
									buf[34] = (*b34)[i]
									buf[35] = (*b35)[i]
									buf[36] = (*b36)[i]
									buf[37] = (*b37)[i]
									buf[38] = (*b38)[i]
									buf[39] = (*b39)[i]
									buf[40] = (*b40)[i]
									buf[41] = (*b41)[i]
									buf[42] = (*b42)[i]
									buf[43] = (*b43)[i]
									buf[44] = (*b44)[i]
									buf[45] = (*b45)[i]
									buf[46] = (*b46)[i]
									buf[47] = (*b47)[i]
									buf[48] = (*b48)[i]
									buf[49] = (*b49)[i]
									buf[50] = (*b50)[i]
									buf[51] = (*b51)[i]
									buf[52] = (*b52)[i]
									buf[53] = (*b53)[i]
									buf[54] = (*b54)[i]
									buf[55] = (*b55)[i]
									buf[56] = (*b56)[i]
									buf[57] = (*b57)[i]
									buf[58] = (*b58)[i]
									buf[59] = (*b59)[i]
									buf[60] = (*b60)[i]
									buf[61] = (*b61)[i]
									buf[62] = (*b62)[i]
									buf[63] = (*b63)[i]

									// count
									pospop.Count8(&counts[i], buf)
								}

								bufIdx = 0
							}
						}
					}
				} else {
					if moreThanOneHash {
						for _, hs = range hashes {
							for i, _h = range hs {
								// loc = int(_h % numSigsUint)
								loc = int(_h & numSigsUintM1) // & X is faster than % X when X is power of 2
								offset = int(offset0 + int64(loc*numRowBytes))

								data[i] = sigs[offset : offset+numRowBytes]
							}

							// first two rows
							pand.AndUnsafe(buffs[bufIdx], data[0], data[1])

							// more rows
							if moreThanTwoHashes {
								for _, row = range data[2:] {
									pand.AndUnsafeInplace(buffs[bufIdx], row)
								}
							}

							// add to buffer for counting
							bufIdx++

							if bufIdx == PosPopCountBufSize {
								// transpose
								for i = 0; i < numRowBytes; i++ { // every column in matrix

									// for j = 0; j < bufIdx; j++ {
									// 	buf[j] = buffs[j][i]
									// }

									// unroll loop
									buf[0] = (*b0)[i]
									buf[1] = (*b1)[i]
									buf[2] = (*b2)[i]
									buf[3] = (*b3)[i]
									buf[4] = (*b4)[i]
									buf[5] = (*b5)[i]
									buf[6] = (*b6)[i]
									buf[7] = (*b7)[i]
									buf[8] = (*b8)[i]
									buf[9] = (*b9)[i]
									buf[10] = (*b10)[i]
									buf[11] = (*b11)[i]
									buf[12] = (*b12)[i]
									buf[13] = (*b13)[i]
									buf[14] = (*b14)[i]
									buf[15] = (*b15)[i]
									buf[16] = (*b16)[i]
									buf[17] = (*b17)[i]
									buf[18] = (*b18)[i]
									buf[19] = (*b19)[i]
									buf[20] = (*b20)[i]
									buf[21] = (*b21)[i]
									buf[22] = (*b22)[i]
									buf[23] = (*b23)[i]
									buf[24] = (*b24)[i]
									buf[25] = (*b25)[i]
									buf[26] = (*b26)[i]
									buf[27] = (*b27)[i]
									buf[28] = (*b28)[i]
									buf[29] = (*b29)[i]
									buf[30] = (*b30)[i]
									buf[31] = (*b31)[i]
									buf[32] = (*b32)[i]
									buf[33] = (*b33)[i]
									buf[34] = (*b34)[i]
									buf[35] = (*b35)[i]
									buf[36] = (*b36)[i]
									buf[37] = (*b37)[i]
									buf[38] = (*b38)[i]
									buf[39] = (*b39)[i]
									buf[40] = (*b40)[i]
									buf[41] = (*b41)[i]
									buf[42] = (*b42)[i]
									buf[43] = (*b43)[i]
									buf[44] = (*b44)[i]
									buf[45] = (*b45)[i]
									buf[46] = (*b46)[i]
									buf[47] = (*b47)[i]
									buf[48] = (*b48)[i]
									buf[49] = (*b49)[i]
									buf[50] = (*b50)[i]
									buf[51] = (*b51)[i]
									buf[52] = (*b52)[i]
									buf[53] = (*b53)[i]
									buf[54] = (*b54)[i]
									buf[55] = (*b55)[i]
									buf[56] = (*b56)[i]
									buf[57] = (*b57)[i]
									buf[58] = (*b58)[i]
									buf[59] = (*b59)[i]
									buf[60] = (*b60)[i]
									buf[61] = (*b61)[i]
									buf[62] = (*b62)[i]
									buf[63] = (*b63)[i]

									// count
									pospop.Count8(&counts[i], buf)
								}

								bufIdx = 0
							}
						}
					} else {
						for _, _h = range hashes1 {
							// loc = int(_h % numSigsUint)
							loc = int(_h & numSigsUintM1) // & X is faster than % X when X is power of 2
							offset = int(offset0 + int64(loc*numRowBytes))

							buffs[bufIdx] = sigs[offset : offset+numRowBytes] // just point to the orginial data (mmaped)

							// add to buffer for counting
							bufIdx++

							if bufIdx == PosPopCountBufSize {
								// transpose
								for i = 0; i < numRowBytes; i++ { // every column in matrix

									// for j = 0; j < bufIdx; j++ {
									// 	buf[j] = buffs[j][i]
									// }

									// unroll loop
									buf[0] = (*b0)[i]
									buf[1] = (*b1)[i]
									buf[2] = (*b2)[i]
									buf[3] = (*b3)[i]
									buf[4] = (*b4)[i]
									buf[5] = (*b5)[i]
									buf[6] = (*b6)[i]
									buf[7] = (*b7)[i]
									buf[8] = (*b8)[i]
									buf[9] = (*b9)[i]
									buf[10] = (*b10)[i]
									buf[11] = (*b11)[i]
									buf[12] = (*b12)[i]
									buf[13] = (*b13)[i]
									buf[14] = (*b14)[i]
									buf[15] = (*b15)[i]
									buf[16] = (*b16)[i]
									buf[17] = (*b17)[i]
									buf[18] = (*b18)[i]
									buf[19] = (*b19)[i]
									buf[20] = (*b20)[i]
									buf[21] = (*b21)[i]
									buf[22] = (*b22)[i]
									buf[23] = (*b23)[i]
									buf[24] = (*b24)[i]
									buf[25] = (*b25)[i]
									buf[26] = (*b26)[i]
									buf[27] = (*b27)[i]
									buf[28] = (*b28)[i]
									buf[29] = (*b29)[i]
									buf[30] = (*b30)[i]
									buf[31] = (*b31)[i]
									buf[32] = (*b32)[i]
									buf[33] = (*b33)[i]
									buf[34] = (*b34)[i]
									buf[35] = (*b35)[i]
									buf[36] = (*b36)[i]
									buf[37] = (*b37)[i]
									buf[38] = (*b38)[i]
									buf[39] = (*b39)[i]
									buf[40] = (*b40)[i]
									buf[41] = (*b41)[i]
									buf[42] = (*b42)[i]
									buf[43] = (*b43)[i]
									buf[44] = (*b44)[i]
									buf[45] = (*b45)[i]
									buf[46] = (*b46)[i]
									buf[47] = (*b47)[i]
									buf[48] = (*b48)[i]
									buf[49] = (*b49)[i]
									buf[50] = (*b50)[i]
									buf[51] = (*b51)[i]
									buf[52] = (*b52)[i]
									buf[53] = (*b53)[i]
									buf[54] = (*b54)[i]
									buf[55] = (*b55)[i]
									buf[56] = (*b56)[i]
									buf[57] = (*b57)[i]
									buf[58] = (*b58)[i]
									buf[59] = (*b59)[i]
									buf[60] = (*b60)[i]
									buf[61] = (*b61)[i]
									buf[62] = (*b62)[i]
									buf[63] = (*b63)[i]

									// count
									pospop.Count8(&counts[i], buf)
								}

								bufIdx = 0
							}
						}
					}
				}
			} else { // !useMmap {
				if compactSize {
					if moreThanOneHash {
						for _, hs = range hashes {
							for i, _h = range hs {
								loc = int(_h % numSigsUint)
								// loc = int(_h & numSigsUintM1) // & X is faster than % X when X is power of 2
								// offset = int(offset0 + int64(loc*numRowBytes))

								// data[i] = sigs[offset : offset+numRowBytes]

								offset2 = offset0 + int64(loc*numRowBytes)
								fh.Seek(offset2, 0)
								io.ReadFull(fh, data[i])
							}

							// first two rows
							pand.AndUnsafe(buffs[bufIdx], data[0], data[1])

							// more rows
							if moreThanTwoHashes {
								for _, row = range data[2:] {
									pand.AndUnsafeInplace(buffs[bufIdx], row)
								}
							}

							// add to buffer for counting
							bufIdx++

							if bufIdx == PosPopCountBufSize {
								// transpose
								for i = 0; i < numRowBytes; i++ { // every column in matrix

									// for j = 0; j < bufIdx; j++ {
									// 	buf[j] = buffs[j][i]
									// }

									// unroll loop
									buf[0] = (*b0)[i]
									buf[1] = (*b1)[i]
									buf[2] = (*b2)[i]
									buf[3] = (*b3)[i]
									buf[4] = (*b4)[i]
									buf[5] = (*b5)[i]
									buf[6] = (*b6)[i]
									buf[7] = (*b7)[i]
									buf[8] = (*b8)[i]
									buf[9] = (*b9)[i]
									buf[10] = (*b10)[i]
									buf[11] = (*b11)[i]
									buf[12] = (*b12)[i]
									buf[13] = (*b13)[i]
									buf[14] = (*b14)[i]
									buf[15] = (*b15)[i]
									buf[16] = (*b16)[i]
									buf[17] = (*b17)[i]
									buf[18] = (*b18)[i]
									buf[19] = (*b19)[i]
									buf[20] = (*b20)[i]
									buf[21] = (*b21)[i]
									buf[22] = (*b22)[i]
									buf[23] = (*b23)[i]
									buf[24] = (*b24)[i]
									buf[25] = (*b25)[i]
									buf[26] = (*b26)[i]
									buf[27] = (*b27)[i]
									buf[28] = (*b28)[i]
									buf[29] = (*b29)[i]
									buf[30] = (*b30)[i]
									buf[31] = (*b31)[i]
									buf[32] = (*b32)[i]
									buf[33] = (*b33)[i]
									buf[34] = (*b34)[i]
									buf[35] = (*b35)[i]
									buf[36] = (*b36)[i]
									buf[37] = (*b37)[i]
									buf[38] = (*b38)[i]
									buf[39] = (*b39)[i]
									buf[40] = (*b40)[i]
									buf[41] = (*b41)[i]
									buf[42] = (*b42)[i]
									buf[43] = (*b43)[i]
									buf[44] = (*b44)[i]
									buf[45] = (*b45)[i]
									buf[46] = (*b46)[i]
									buf[47] = (*b47)[i]
									buf[48] = (*b48)[i]
									buf[49] = (*b49)[i]
									buf[50] = (*b50)[i]
									buf[51] = (*b51)[i]
									buf[52] = (*b52)[i]
									buf[53] = (*b53)[i]
									buf[54] = (*b54)[i]
									buf[55] = (*b55)[i]
									buf[56] = (*b56)[i]
									buf[57] = (*b57)[i]
									buf[58] = (*b58)[i]
									buf[59] = (*b59)[i]
									buf[60] = (*b60)[i]
									buf[61] = (*b61)[i]
									buf[62] = (*b62)[i]
									buf[63] = (*b63)[i]

									// count
									pospop.Count8(&counts[i], buf)
								}

								bufIdx = 0
							}
						}
					} else {
						for _, _h = range hashes1 {
							loc = int(_h % numSigsUint)
							// loc = int(_h & numSigsUintM1) // & X is faster than % X when X is power of 2
							// offset = int(offset0 + int64(loc*numRowBytes))

							offset2 = offset0 + int64(loc*numRowBytes)
							fh.Seek(offset2, 0)
							io.ReadFull(fh, buffs[bufIdx])

							// add to buffer for counting
							bufIdx++

							if bufIdx == PosPopCountBufSize {
								// transpose
								for i = 0; i < numRowBytes; i++ { // every column in matrix

									// for j = 0; j < bufIdx; j++ {
									// 	buf[j] = buffs[j][i]
									// }

									// unroll loop
									buf[0] = (*b0)[i]
									buf[1] = (*b1)[i]
									buf[2] = (*b2)[i]
									buf[3] = (*b3)[i]
									buf[4] = (*b4)[i]
									buf[5] = (*b5)[i]
									buf[6] = (*b6)[i]
									buf[7] = (*b7)[i]
									buf[8] = (*b8)[i]
									buf[9] = (*b9)[i]
									buf[10] = (*b10)[i]
									buf[11] = (*b11)[i]
									buf[12] = (*b12)[i]
									buf[13] = (*b13)[i]
									buf[14] = (*b14)[i]
									buf[15] = (*b15)[i]
									buf[16] = (*b16)[i]
									buf[17] = (*b17)[i]
									buf[18] = (*b18)[i]
									buf[19] = (*b19)[i]
									buf[20] = (*b20)[i]
									buf[21] = (*b21)[i]
									buf[22] = (*b22)[i]
									buf[23] = (*b23)[i]
									buf[24] = (*b24)[i]
									buf[25] = (*b25)[i]
									buf[26] = (*b26)[i]
									buf[27] = (*b27)[i]
									buf[28] = (*b28)[i]
									buf[29] = (*b29)[i]
									buf[30] = (*b30)[i]
									buf[31] = (*b31)[i]
									buf[32] = (*b32)[i]
									buf[33] = (*b33)[i]
									buf[34] = (*b34)[i]
									buf[35] = (*b35)[i]
									buf[36] = (*b36)[i]
									buf[37] = (*b37)[i]
									buf[38] = (*b38)[i]
									buf[39] = (*b39)[i]
									buf[40] = (*b40)[i]
									buf[41] = (*b41)[i]
									buf[42] = (*b42)[i]
									buf[43] = (*b43)[i]
									buf[44] = (*b44)[i]
									buf[45] = (*b45)[i]
									buf[46] = (*b46)[i]
									buf[47] = (*b47)[i]
									buf[48] = (*b48)[i]
									buf[49] = (*b49)[i]
									buf[50] = (*b50)[i]
									buf[51] = (*b51)[i]
									buf[52] = (*b52)[i]
									buf[53] = (*b53)[i]
									buf[54] = (*b54)[i]
									buf[55] = (*b55)[i]
									buf[56] = (*b56)[i]
									buf[57] = (*b57)[i]
									buf[58] = (*b58)[i]
									buf[59] = (*b59)[i]
									buf[60] = (*b60)[i]
									buf[61] = (*b61)[i]
									buf[62] = (*b62)[i]
									buf[63] = (*b63)[i]

									// count
									pospop.Count8(&counts[i], buf)
								}

								bufIdx = 0
							}
						}
					}
				} else {
					if moreThanOneHash {
						for _, hs = range hashes {
							for i, _h = range hs {
								// loc = int(_h % numSigsUint)
								loc = int(_h & numSigsUintM1) // & X is faster than % X when X is power of 2
								// offset = int(offset0 + int64(loc*numRowBytes))

								// data[i] = sigs[offset : offset+numRowBytes]

								offset2 = offset0 + int64(loc*numRowBytes)
								fh.Seek(offset2, 0)
								io.ReadFull(fh, data[i])
							}

							// first two rows
							pand.AndUnsafe(buffs[bufIdx], data[0], data[1])

							// more rows
							if moreThanTwoHashes {
								for _, row = range data[2:] {
									pand.AndUnsafeInplace(buffs[bufIdx], row)
								}
							}

							// add to buffer for counting
							bufIdx++

							if bufIdx == PosPopCountBufSize {
								// transpose
								for i = 0; i < numRowBytes; i++ { // every column in matrix

									// for j = 0; j < bufIdx; j++ {
									// 	buf[j] = buffs[j][i]
									// }

									// unroll loop
									buf[0] = (*b0)[i]
									buf[1] = (*b1)[i]
									buf[2] = (*b2)[i]
									buf[3] = (*b3)[i]
									buf[4] = (*b4)[i]
									buf[5] = (*b5)[i]
									buf[6] = (*b6)[i]
									buf[7] = (*b7)[i]
									buf[8] = (*b8)[i]
									buf[9] = (*b9)[i]
									buf[10] = (*b10)[i]
									buf[11] = (*b11)[i]
									buf[12] = (*b12)[i]
									buf[13] = (*b13)[i]
									buf[14] = (*b14)[i]
									buf[15] = (*b15)[i]
									buf[16] = (*b16)[i]
									buf[17] = (*b17)[i]
									buf[18] = (*b18)[i]
									buf[19] = (*b19)[i]
									buf[20] = (*b20)[i]
									buf[21] = (*b21)[i]
									buf[22] = (*b22)[i]
									buf[23] = (*b23)[i]
									buf[24] = (*b24)[i]
									buf[25] = (*b25)[i]
									buf[26] = (*b26)[i]
									buf[27] = (*b27)[i]
									buf[28] = (*b28)[i]
									buf[29] = (*b29)[i]
									buf[30] = (*b30)[i]
									buf[31] = (*b31)[i]
									buf[32] = (*b32)[i]
									buf[33] = (*b33)[i]
									buf[34] = (*b34)[i]
									buf[35] = (*b35)[i]
									buf[36] = (*b36)[i]
									buf[37] = (*b37)[i]
									buf[38] = (*b38)[i]
									buf[39] = (*b39)[i]
									buf[40] = (*b40)[i]
									buf[41] = (*b41)[i]
									buf[42] = (*b42)[i]
									buf[43] = (*b43)[i]
									buf[44] = (*b44)[i]
									buf[45] = (*b45)[i]
									buf[46] = (*b46)[i]
									buf[47] = (*b47)[i]
									buf[48] = (*b48)[i]
									buf[49] = (*b49)[i]
									buf[50] = (*b50)[i]
									buf[51] = (*b51)[i]
									buf[52] = (*b52)[i]
									buf[53] = (*b53)[i]
									buf[54] = (*b54)[i]
									buf[55] = (*b55)[i]
									buf[56] = (*b56)[i]
									buf[57] = (*b57)[i]
									buf[58] = (*b58)[i]
									buf[59] = (*b59)[i]
									buf[60] = (*b60)[i]
									buf[61] = (*b61)[i]
									buf[62] = (*b62)[i]
									buf[63] = (*b63)[i]

									// count
									pospop.Count8(&counts[i], buf)
								}

								bufIdx = 0
							}
						}
					} else {
						for _, _h = range hashes1 {
							// loc = int(_h % numSigsUint)
							loc = int(_h & numSigsUintM1) // & X is faster than % X when X is power of 2
							// offset = int(offset0 + int64(loc*numRowBytes))

							offset2 = offset0 + int64(loc*numRowBytes)
							fh.Seek(offset2, 0)
							io.ReadFull(fh, buffs[bufIdx])

							// add to buffer for counting
							bufIdx++

							if bufIdx == PosPopCountBufSize {
								// transpose
								for i = 0; i < numRowBytes; i++ { // every column in matrix

									// for j = 0; j < bufIdx; j++ {
									// 	buf[j] = buffs[j][i]
									// }

									// unroll loop
									buf[0] = (*b0)[i]
									buf[1] = (*b1)[i]
									buf[2] = (*b2)[i]
									buf[3] = (*b3)[i]
									buf[4] = (*b4)[i]
									buf[5] = (*b5)[i]
									buf[6] = (*b6)[i]
									buf[7] = (*b7)[i]
									buf[8] = (*b8)[i]
									buf[9] = (*b9)[i]
									buf[10] = (*b10)[i]
									buf[11] = (*b11)[i]
									buf[12] = (*b12)[i]
									buf[13] = (*b13)[i]
									buf[14] = (*b14)[i]
									buf[15] = (*b15)[i]
									buf[16] = (*b16)[i]
									buf[17] = (*b17)[i]
									buf[18] = (*b18)[i]
									buf[19] = (*b19)[i]
									buf[20] = (*b20)[i]
									buf[21] = (*b21)[i]
									buf[22] = (*b22)[i]
									buf[23] = (*b23)[i]
									buf[24] = (*b24)[i]
									buf[25] = (*b25)[i]
									buf[26] = (*b26)[i]
									buf[27] = (*b27)[i]
									buf[28] = (*b28)[i]
									buf[29] = (*b29)[i]
									buf[30] = (*b30)[i]
									buf[31] = (*b31)[i]
									buf[32] = (*b32)[i]
									buf[33] = (*b33)[i]
									buf[34] = (*b34)[i]
									buf[35] = (*b35)[i]
									buf[36] = (*b36)[i]
									buf[37] = (*b37)[i]
									buf[38] = (*b38)[i]
									buf[39] = (*b39)[i]
									buf[40] = (*b40)[i]
									buf[41] = (*b41)[i]
									buf[42] = (*b42)[i]
									buf[43] = (*b43)[i]
									buf[44] = (*b44)[i]
									buf[45] = (*b45)[i]
									buf[46] = (*b46)[i]
									buf[47] = (*b47)[i]
									buf[48] = (*b48)[i]
									buf[49] = (*b49)[i]
									buf[50] = (*b50)[i]
									buf[51] = (*b51)[i]
									buf[52] = (*b52)[i]
									buf[53] = (*b53)[i]
									buf[54] = (*b54)[i]
									buf[55] = (*b55)[i]
									buf[56] = (*b56)[i]
									buf[57] = (*b57)[i]
									buf[58] = (*b58)[i]
									buf[59] = (*b59)[i]
									buf[60] = (*b60)[i]
									buf[61] = (*b61)[i]
									buf[62] = (*b62)[i]
									buf[63] = (*b63)[i]

									// count
									pospop.Count8(&counts[i], buf)
								}

								bufIdx = 0
							}
						}
					}
				}
			}

			// -------------------------------------------------------------------------
			// -------------------------------------------------------------------------

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

			// results := make([]Match, 0, 8)
			results := poolMatches.Get().([]Match)

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
	}

	// start two goroutines
	go fn()

	if opt.UseMMap {
		go fn()
	}

	return idx, nil
}

// Close closes the index.
func (idx *UnikIndex) Close() error {
	close(idx.InCh) // close InCh
	<-idx.done      // confirm stopped
	// another one
	if idx.useMmap {
		<-idx.done
	}

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
	return make([]Match, 0, 256)
}}

var poolChanMatches = &sync.Pool{New: func() interface{} {
	return make(chan []Match, 256)
}}

// not used
var poolIdxValues = &sync.Pool{New: func() interface{} {
	return make([]unikmer.IdxValue, 0, 128)
}}

// not used
var poolBytes = &sync.Pool{New: func() interface{} {
	return make([]byte, 0, 10<<20)
}}
