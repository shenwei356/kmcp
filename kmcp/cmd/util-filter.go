// Copyright © 2020-2022 Wei Shen <shenwei356@gmail.com>
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
	"strconv"
	"sync"
)

type MatchResult2 struct {
	Query string
	// QLen    int
	// QKmers  int
	FPR float64
	// Hits    int
	Target string
	// FragIdx int
	// IdxNum  int
	// GSize   uint64
	// K       int
	// MKmers  int
	QCov float64

	Line *string
}

func parseMatchResult2(line string, numFields int, items *[]string, maxPFR float64, minQcov float64) (*MatchResult2, bool) {
	stringSplitNByByte(line, '\t', numFields, items)
	if len(*items) < numFields {
		checkError(fmt.Errorf("invalid kmcp search result format"))
	}

	m := &MatchResult2{} // do not use sync.Pool, which is slower for this case

	var err error

	// too slow
	m.FPR, err = strconv.ParseFloat((*items)[3], 64)
	if err != nil {
		checkError(fmt.Errorf("failed to parse FPR: %s", (*items)[3]))
	}
	if m.FPR > maxPFR {
		return m, false
	}

	// too slow
	m.QCov, err = strconv.ParseFloat((*items)[11], 64)
	if err != nil {
		checkError(fmt.Errorf("failed to parse qCov: %s", (*items)[11]))
	}
	if m.QCov < minQcov {
		return m, false
	}

	// -----------

	m.Query = (*items)[0]
	m.Target = (*items)[5]

	m.Line = &line

	return m, true
}

var poolMatchResults2 = &sync.Pool{New: func() interface{} {
	tmp := make([]*MatchResult2, 0, 128)
	return &tmp
}}
