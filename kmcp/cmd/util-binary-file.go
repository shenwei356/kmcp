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

	"github.com/pkg/errors"
	"github.com/shenwei356/unik/v5"
)

const extDataFile = ".unik"

func checkCompatibility(reader0 *unik.Reader, reader *unik.Reader, file string, meta0 *Meta, meta *Meta) {
	if reader0.K != reader.K {
		checkError(fmt.Errorf(`k-mer length not consistent (%d != %d), please check with "unikmer info": %s`, reader0.K, reader.K, file))
	}
	if reader0.IsCanonical() != reader.IsCanonical() {
		checkError(fmt.Errorf(`'canonical' flags not consistent, please check with "unikmer info": %s`, file))
	}
	if reader0.IsHashed() != reader.IsHashed() {
		checkError(fmt.Errorf(`'hashed' flags not consistent, please check with "unikmer info": %s`, file))
	}
	if reader0.IsScaled() != reader.IsScaled() {
		checkError(fmt.Errorf(`'scaled' flags not consistent, please check with "unikmer info": %s`, file))
	}

	if meta0.MinimizerW == meta.MinimizerW &&
		meta0.SyncmerS == meta.SyncmerS &&
		meta0.SplitSize == meta.SplitSize &&
		meta0.SplitOverlap == meta.SplitOverlap {
		return
	}
	checkError(fmt.Errorf(`sketch information (description) not consistent, please check with "unikmer info -a ": %s. file1: %s, file: %s`,
		file, meta0, meta))

}

// not used
func readKmers(file string) ([]uint64, error) {
	infh, r, _, err := inStream(file)
	if err != nil {
		return nil, err
	}
	defer r.Close()

	reader, err := unik.NewReader(infh)
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
