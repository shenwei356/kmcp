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
	"io/ioutil"
	"math"
	"os"
	"path/filepath"
	"sync"

	"github.com/clausecker/pospop"
	"github.com/edsrzf/mmap-go"
	"github.com/pkg/errors"
	"github.com/shenwei356/unikmer/index"
	"github.com/shenwei356/util/pathutil"
	"github.com/smallnest/ringbuffer"
	"gopkg.in/yaml.v2"
)

const dbInfoFile = "_db.yml"

var ErrVersionMismatch = errors.New("kmcp/index: version mismatch")

const UnikIndexDBVersion uint8 = 3

type UnikIndexDBInfo struct {
	Version      uint8 `yaml:"version"`
	IndexVersion uint8 `yaml:"unikiVersion"`
	K            int   `yaml:"k"`
	Hashed       bool  `yaml:"hashed"`
	Canonical    bool  `yaml:"canonical"`

	Scaled     bool   `yaml:"scaled"`
	Scale      uint32 `yaml:"scale"`
	Minizimer  bool   `yaml:"minimizer"`
	MinizimerW uint32 `yaml:"minimizer-w"`
	Syncmer    bool   `yaml:"syncmer"`
	SyncmerS   uint32 `yaml:"syncmer-s"`

	NumHashes int      `yaml:"hashes"`
	FPR       float64  `yaml:"fpr"`
	BlockSize int      `yaml:"blocksize"`
	Kmers     int      `yaml:"totalKmers"`
	Files     []string `yaml:"files"`
	NumNames  int      `yaml:"numNames"`
	Names     []string `yaml:"names"`
	Sizes     []uint64 `yaml:"kmers"`

	path string
}

func (i UnikIndexDBInfo) String() string {
	return fmt.Sprintf("kmcp database v%d: k: %d, hashed: %v,  canonical: %v, #hashes: %d, fpr:%f, #blocksize: %d, #blocks: %d, #%d-mers: %d",
		i.Version, i.K, i.Hashed, i.Canonical, i.NumHashes, i.FPR, i.BlockSize, len(i.Files), i.K, i.Kmers)
}

func NewUnikIndexDBInfo(files []string) UnikIndexDBInfo {
	return UnikIndexDBInfo{Version: UnikIndexDBVersion, IndexVersion: index.Version, Files: files}
}

func UnikIndexDBInfoFromFile(file string) (UnikIndexDBInfo, error) {
	info := UnikIndexDBInfo{}

	r, err := os.Open(file)
	if err != nil {
		return info, fmt.Errorf("fail to read unikmer index db info file: %s", file)
	}

	data, err := ioutil.ReadAll(r)
	if err != nil {
		return info, fmt.Errorf("fail to read unikmer index db info file: %s", file)
	}

	err = yaml.Unmarshal(data, &info)
	if err != nil {
		return info, fmt.Errorf("fail to unmarshal unikmer index db info")
	}

	r.Close()

	if info.Version != UnikIndexDBVersion {
		return info, ErrVersionMismatch
	}

	p, _ := filepath.Abs(file)
	info.path = filepath.Dir(p)
	return info, nil
}

func (i UnikIndexDBInfo) WriteTo(file string) error {
	data, err := yaml.Marshal(i)
	if err != nil {
		return fmt.Errorf("fail to marshal uniker index db info")
	}

	w, err := os.Create(file)
	if err != nil {
		return fmt.Errorf("fail to write unikmer index db info file: %s", file)
	}
	_, err = w.Write(data)
	if err != nil {
		return fmt.Errorf("fail to write unikmer index db info file: %s", file)
	}

	w.Close()
	return nil
}

func (i UnikIndexDBInfo) Check() error {
	for _, file := range i.Files {
		file = filepath.Join(i.path, file)
		ok, err := pathutil.Exists(file)
		if err != nil {
			return fmt.Errorf("error on checking unikmer index file: %s: %s", file, err)
		}
		if !ok {
			return fmt.Errorf("unikmer index file missing: %s", file)
		}
	}
	return nil
}

// ------------------------------------------------------------------

type UnikIndexDB struct {
	path string

	Info   UnikIndexDBInfo
	Header index.Header

	Indices []*UnikIndex
}

func (db *UnikIndexDB) String() string {
	return fmt.Sprintf("unikmer index db v%d: #blocksize: %d, #blocks: %d, #%d-mers: %d, #hashes: %d",
		db.Info.Version, db.Info.BlockSize, len(db.Info.Files), db.Header.K, db.Info.Kmers, db.Header.NumHashes)
}

func NewUnikIndexDB(path string, useMmap bool) (*UnikIndexDB, error) {
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
	idx1, err := NewUnixIndex(filepath.Join(path, info.Files[0]), useMmap)
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

			idx, err := NewUnixIndex(f, useMmap)
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
	return db, nil
}

func (db *UnikIndexDB) Close() error {
	var err0 error
	for _, idx := range db.Indices {
		err := idx.Close()
		if err != nil && err0 == nil {
			err0 = err
		}
	}
	return err0
}

func (db *UnikIndexDB) SearchMap(kmers map[uint64]interface{}, threads int, queryCov float64, targetCov float64) map[string][3]float64 {
	numHashes := db.Info.NumHashes
	hashes := make([][]uint64, len(kmers))
	i := 0

	if db.Info.Hashed {
		for kmer := range kmers {
			hashes[i] = hashValues(kmer, numHashes)
			i++
		}
	} else {
		for kmer := range kmers {
			hashes[i] = hashValues(hash64(kmer), numHashes)
			i++
		}
	}

	return db.searchHashes(hashes, threads, queryCov, targetCov)
}

func (db *UnikIndexDB) Search(kmers []uint64, threads int, queryCov float64, targetCov float64) map[string][3]float64 {
	numHashes := db.Info.NumHashes
	hashes := make([][]uint64, len(kmers))
	if db.Info.Hashed {
		for i, kmer := range kmers {
			hashes[i] = hashValues(kmer, numHashes)
		}
	} else {
		for i, kmer := range kmers {
			hashes[i] = hashValues(hash64(kmer), numHashes)
		}
	}

	return db.searchHashes(hashes, threads, queryCov, targetCov)
}

func (db *UnikIndexDB) searchHashes(hashes [][]uint64, threads int, queryCov float64, targetCov float64) map[string][3]float64 {

	m := make(map[string][3]float64, 8)

	if threads <= 0 {
		threads = len(db.Indices)
	}

	ch := make([]*map[string][3]float64, len(db.Indices))

	var wg sync.WaitGroup
	// tokens := make(chan int, threads)
	tokens := ringbuffer.New(threads) // ringbufer is faster than channel
	// for i, idx := range db.Indices {
	for i := len(db.Indices) - 1; i >= 0; i-- { // start from bigger files
		idx := db.Indices[i]

		wg.Add(1)
		// tokens <- 1
		tokens.WriteByte(0)
		go func(idx *UnikIndex, ch []*map[string][3]float64, i int) {
			_m := idx.Search(hashes, queryCov, targetCov)
			ch[i] = &_m

			wg.Done()
			// <-tokens
			tokens.ReadByte()
		}(idx, ch, i)
	}
	wg.Wait()

	var k string
	var v [3]float64
	for _, _m := range ch {
		for k, v = range *_m {
			m[k] = v
		}
	}

	return m
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

type UnikIndex struct {
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

func NewUnixIndex(file string, useMmap bool) (*UnikIndex, error) {
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
	idx := &UnikIndex{Path: file, Header: h, fh: fh, reader: reader, offset0: offset}
	idx.useMmap = useMmap

	if useMmap {
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
	return idx, nil
}

func (idx *UnikIndex) Search(hashes [][]uint64, queryCov float64, targetCov float64) map[string][3]float64 {
	numNames := len(idx.Header.Names)
	numRowBytes := idx.Header.NumRowBytes
	numSigs := idx.Header.NumSigs
	numSigsInt := uint64(numSigs)
	offset0 := idx.offset0
	fh := idx.fh
	names := idx.Header.Names
	sizes := idx.Header.Sizes
	sigs := idx.sigsB
	useMmap := idx.useMmap
	data := idx._data
	moreThanOneHash := idx.moreThanOneHash

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

	result := make(map[string][3]float64, 8)

	var c, t, T float64
	iLast := numRowBytes - 1
	for i, _counts = range counts {
		ix8 = i << 3
		// for j, count = range _counts {
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

			result[names[k]] = [3]float64{c, t, T}
		}
	}

	return result
}

func (idx *UnikIndex) Close() error {
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
