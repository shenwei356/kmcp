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

package index

import (
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"strings"
)

// Version is the version of index format
const Version uint8 = 4

// Magic number of index file.
var Magic = [8]byte{'.', 'k', 'm', 'c', 'p', 'i', 'd', 'x'}

// ErrInvalidIndexFileFormat means invalid index format.
var ErrInvalidIndexFileFormat = errors.New("kmcp: invalid index format")

// ErrUnfishedWrite means writing not finished
var ErrUnfishedWrite = errors.New("kmcp: index not fished writing")

// ErrTruncateIndexFile means the file is truncated
var ErrTruncateIndexFile = errors.New("kmcp: truncated index file")

// ErrWrongWriteDataSize means the size of data to write is invalid
var ErrWrongWriteDataSize = errors.New("kmcp: write data with wrong size")

// ErrVersionMismatch means version mismatch between files and program
var ErrVersionMismatch = errors.New("kmcp: version mismatch")

// ErrNameAndSizeMismatch means size of names and sizes are not equal.
var ErrNameAndSizeMismatch = errors.New("kmcp: size of names and sizes unequal")

// ErrNameAndIndexMismatch means size of names and sizes are not equal.
var ErrNameAndIndexMismatch = errors.New("kmcp: size of names and indices unequal")

var be = binary.BigEndian

const (
	CANONICAL = 1 << iota
	COMPACT
)

// Header contains metadata
type Header struct {
	Version   uint8 // uint8
	K         int   // uint8
	flag      uint8 //uint8
	Canonical bool
	Compact   bool
	NumHashes uint8 // uint8
	NumSigs   uint64

	Names  [][]string // one bloom filter contains union of multiple sets
	GSizes [][]uint64 // genome sizes
	// Kmers   [][]uint64 // kmer numbers
	Indices [][]uint32 // coresponding fragment indices of all sets.
	Sizes   []uint64

	NumRowBytes int // length of bytes for storing one row of signiture for n names
}

func (h Header) String() string {
	return fmt.Sprintf("kmcp index file v%d: k: %d, canonical: %v, #hashes: %d, #signatures: %d, #names: %d",
		h.Version, h.K, h.Canonical, h.NumHashes, h.NumSigs, len(h.Names))
}

// Compatible checks compatibility
func (h Header) Compatible(b Header) bool {
	if h.Version == b.Version &&
		h.K == b.K &&
		h.Canonical == b.Canonical &&
		h.NumHashes == b.NumHashes {

		return true
	}
	return false
}

// ------------------------------------------------------------------------

// Writer writes KmerCode.
type Writer struct {
	Header
	w           io.Writer
	wroteHeader bool

	count uint64
}

// NewWriter creates a Writer.
// func NewWriter(w io.Writer, k int, canonical bool, compact bool, numHashes uint8, numSigs uint64,
// 	names [][]string, gsizes [][]uint64, kmers [][]uint64, indices [][]uint32, sizes []uint64) (*Writer, error) {
func NewWriter(w io.Writer, k int, canonical bool, compact bool, numHashes uint8, numSigs uint64,
	names [][]string, gsizes [][]uint64, indices [][]uint32, sizes []uint64) (*Writer, error) {
	if len(names) != len(sizes) {
		return nil, ErrNameAndSizeMismatch
	}
	if len(names) != len(indices) {
		return nil, ErrNameAndIndexMismatch
	}

	writer := &Writer{
		Header: Header{
			Version:   Version,
			K:         k,
			Canonical: canonical,
			Compact:   compact,
			NumHashes: numHashes,
			NumSigs:   numSigs,
			Names:     names,
			GSizes:    gsizes,
			// Kmers:     kmers,
			Indices: indices,
			Sizes:   sizes,
		},
		w: w,
	}
	writer.NumRowBytes = int((len(names) + 7) / 8)
	writer.flag = 0
	if canonical {
		writer.flag |= CANONICAL
	}
	if compact {
		writer.flag |= COMPACT
	}

	return writer, nil
}

// WriteHeader writes file header
func (writer *Writer) WriteHeader() (err error) {
	if writer.wroteHeader {
		return nil
	}
	w := writer.w

	// 8 bytes magic number
	err = binary.Write(w, be, Magic)
	if err != nil {
		return err
	}

	// 4 bytes meta info

	err = binary.Write(w, be, [4]uint8{writer.Version, uint8(writer.K), writer.flag, writer.NumHashes})
	if err != nil {
		return err
	}

	// 8 bytes signature size
	err = binary.Write(w, be, writer.NumSigs)
	if err != nil {
		return err
	}

	// ----------------------------------------------------------
	// names

	// 4 bytes length of name groups
	err = binary.Write(w, be, uint32(len(writer.Names)))
	if err != nil {
		return err
	}

	// N := 24

	var name string
	for _, names := range writer.Names {
		n := 0

		// 4 bytes length of Names
		for _, name := range names {
			n += len(name) + 1
		}

		err = binary.Write(w, be, uint32(n))
		if err != nil {
			return err
		}
		// N += 4

		// Names
		for _, name = range names {
			err = binary.Write(w, be, []byte(name+"\n"))
			if err != nil {
				return err
			}
			// N += len(name + "\n")
		}
	}

	// ----------------------------------------------------------
	// GSizes

	// 4 bytes length of GSizes groups
	err = binary.Write(w, be, uint32(len(writer.GSizes)))
	if err != nil {
		return err
	}
	// N += 4

	for _, gsizes := range writer.GSizes {
		err = binary.Write(w, be, uint32(len(gsizes)))
		if err != nil {
			return err
		}
		// N += 4
		err = binary.Write(w, be, gsizes)
		if err != nil {
			return err
		}
		// N += 8 * len(gsizes)
	}

	// // ----------------------------------------------------------
	// // Kmers

	// // 4 bytes length of Kmers groups
	// err = binary.Write(w, be, uint32(len(writer.Kmers)))
	// if err != nil {
	// 	return err
	// }

	// for _, kmers := range writer.Kmers {
	// 	err = binary.Write(w, be, uint32(len(kmers)))
	// 	if err != nil {
	// 		return err
	// 	}
	// 	err = binary.Write(w, be, kmers)
	// 	if err != nil {
	// 		return err
	// 	}
	// }

	// ----------------------------------------------------------
	// Indices

	// 4 bytes length of Indices groups
	err = binary.Write(w, be, uint32(len(writer.Indices)))
	if err != nil {
		return err
	}
	// N += 4

	for _, indices := range writer.Indices {
		err = binary.Write(w, be, uint32(len(indices)))
		if err != nil {
			return err
		}
		// N += 4

		err = binary.Write(w, be, indices)
		if err != nil {
			return err
		}
		// N += 4 * len(indices)
	}

	// Sizes
	err = binary.Write(w, be, writer.Sizes)
	if err != nil {
		return err
	}
	// N += 8 * len(writer.Sizes)

	// // padding to 64X
	// nPadding := 56 - (N % 64)
	// err = binary.Write(w, be, uint64(nPadding))
	// if err != nil {
	// 	return err
	// }
	// // N += 8
	// padding := make([]byte, nPadding)
	// err = binary.Write(w, be, padding)
	// if err != nil {
	// 	return err
	// }
	// // N += len(padding)

	writer.wroteHeader = true
	return nil
}

// Write writes one row of signitures
func (writer *Writer) Write(data []byte) (err error) {
	if len(data) != writer.NumRowBytes {
		return ErrWrongWriteDataSize
	}

	// lazily write header
	if !writer.wroteHeader {
		err = writer.WriteHeader()
		if err != nil {
			return err
		}
		writer.wroteHeader = true
	}

	_, err = writer.w.Write(data)
	if err != nil {
		return err
	}

	writer.count++
	return nil
}

// WriteBatch writes a batch of data
func (writer *Writer) WriteBatch(data []byte, n int) (err error) {
	// lazily write header
	if !writer.wroteHeader {
		err = writer.WriteHeader()
		if err != nil {
			return err
		}
		writer.wroteHeader = true
	}

	_, err = writer.w.Write(data)
	if err != nil {
		return err
	}

	writer.count += uint64(n)
	return nil
}

// Flush check completeness
func (writer *Writer) Flush() (err error) {
	if !writer.wroteHeader {
		writer.WriteHeader()
	}
	if writer.count != writer.NumSigs {
		return ErrUnfishedWrite
	}
	return nil
}

// ------------------------------------------------------------------------

// Reader is for reading signatures.
type Reader struct {
	Header
	r io.Reader

	count uint64
}

// NewReader returns a Reader.
func NewReader(r io.Reader) (reader *Reader, err error) {
	reader = &Reader{r: r}
	err = reader.readHeader()
	if err != nil {
		return nil, err
	}

	reader.NumRowBytes = int((len(reader.Names) + 7) / 8)
	return reader, nil
}

func (reader *Reader) readHeader() (err error) {
	buf := make([]byte, 56)
	r := reader.r

	// check Magic number (8 byte)
	_, err = io.ReadFull(r, buf[:8])
	if err != nil {
		return err
	}
	same := true
	for i := 0; i < 8; i++ {
		if Magic[i] != buf[i] {
			same = false
			break
		}
	}
	if !same {
		return ErrInvalidIndexFileFormat
	}

	// 4 bytes meta info
	_, err = io.ReadFull(r, buf[:4])
	if err != nil {
		return err
	}
	// check compatibility
	if Version != buf[0] {
		return ErrVersionMismatch
	}
	reader.Version = buf[0]
	reader.K = int(buf[1])

	if buf[2]&CANONICAL > 0 {
		reader.Canonical = true
	}
	if buf[2]&COMPACT > 0 {
		reader.Compact = true
	}

	reader.NumHashes = buf[3]

	// 8 bytes signature size
	_, err = io.ReadFull(r, buf[:8])
	if err != nil {
		return err
	}
	reader.NumSigs = be.Uint64(buf[:8])

	// ----------------------------------------------------------
	// names

	// 4 bytes length of Names groups
	var n uint32
	_, err = io.ReadFull(r, buf[:4])
	if err != nil {
		return err
	}
	n = be.Uint32(buf[:4])

	names := make([][]string, n)

	for i := 0; i < int(n); i++ {
		// 4 bytes length of Names
		_, err = io.ReadFull(r, buf[:4])
		if err != nil {
			return err
		}

		// Names
		namesData := make([]byte, be.Uint32(buf[:4]))
		_, err = io.ReadFull(r, namesData)
		if err != nil {
			return err
		}

		_names := strings.Split(string(namesData), "\n")
		names[i] = _names[0 : len(_names)-1]
	}
	reader.Names = names

	// ----------------------------------------------------------
	// GenomeSizes

	// 4 bytes length of Gsizes groups
	_, err = io.ReadFull(r, buf[:4])
	if err != nil {
		return err
	}
	n = be.Uint32(buf[:4])
	gsizes := make([][]uint64, n)

	var _n int
	var j int
	var buf2 []byte
	var k int

	for i := 0; i < int(n); i++ {
		_, err = io.ReadFull(r, buf[:4])
		if err != nil {
			return err
		}
		_n = int(be.Uint32(buf[:4]))
		buf2 = make([]byte, _n<<3)
		_, err = io.ReadFull(r, buf2)
		if err != nil {
			return err
		}

		gsizesData := make([]uint64, _n)
		for j = 0; j < _n; j++ {
			k = j << 3
			gsizesData[j] = be.Uint64(buf2[k : k+8])
		}
		gsizes[i] = gsizesData
	}
	reader.GSizes = gsizes

	// // ----------------------------------------------------------
	// // Kmers

	// // 4 bytes length of Kmers groups
	// _, err = io.ReadFull(r, buf[:4])
	// if err != nil {
	// 	return err
	// }
	// n = be.Uint32(buf[:4])
	// kmers := make([][]uint64, n)

	// for i := 0; i < int(n); i++ {
	// 	_, err = io.ReadFull(r, buf[:4])
	// 	if err != nil {
	// 		return err
	// 	}
	// 	_n = int(be.Uint32(buf[:4]))
	// 	buf2 = make([]byte, _n<<3)
	// 	_, err = io.ReadFull(r, buf2)
	// 	if err != nil {
	// 		return err
	// 	}

	// 	kmersData := make([]uint64, _n)
	// 	for j = 0; j < _n; j++ {
	// 		k = j << 3
	// 		kmersData[j] = be.Uint64(buf2[k : k+8])
	// 	}
	// 	kmers[i] = kmersData
	// }
	// reader.Kmers = kmers

	// ----------------------------------------------------------
	// Indices

	// 4 bytes length of Indices groups
	_, err = io.ReadFull(r, buf[:4])
	if err != nil {
		return err
	}
	n = be.Uint32(buf[:4])
	indices := make([][]uint32, n)

	for i := 0; i < int(n); i++ {
		_, err = io.ReadFull(r, buf[:4])
		if err != nil {
			return err
		}
		_n = int(be.Uint32(buf[:4]))
		buf2 = make([]byte, _n<<2)
		_, err = io.ReadFull(r, buf2)
		if err != nil {
			return err
		}

		indicesData := make([]uint32, _n)
		for j = 0; j < _n; j++ {
			k = j << 2
			indicesData[j] = be.Uint32(buf2[k : k+4])
		}
		indices[i] = indicesData
	}
	reader.Indices = indices

	// Sizes
	buf2 = make([]byte, len(names)<<3)
	_, err = io.ReadFull(r, buf2)
	if err != nil {
		return err
	}

	sizesData := make([]uint64, len(names))
	for j = 0; j < len(names); j++ {
		k = j << 3
		sizesData[j] = be.Uint64(buf2[k : k+8])
	}

	reader.Sizes = sizesData

	// // 8 bytes padding size
	// _, err = io.ReadFull(r, buf[:8])
	// if err != nil {
	// 	return err
	// }
	// nPadding := be.Uint64(buf[:8])

	// buf2 = make([]byte, nPadding)
	// _, err = io.ReadFull(r, buf2)
	// if err != nil {
	// 	return err
	// }

	return nil
}

// Read reads one code.
func (reader *Reader) Read() ([]byte, error) {
	data := make([]byte, reader.NumRowBytes)
	nReaded, err := io.ReadFull(reader.r, data)
	if err != nil {
		if err == io.EOF {
			if reader.count != reader.NumSigs {
				return nil, ErrTruncateIndexFile
			}
		}
		return nil, err
	}
	if nReaded < reader.NumRowBytes {
		return nil, ErrTruncateIndexFile
	}
	reader.count++
	return data, nil
}
