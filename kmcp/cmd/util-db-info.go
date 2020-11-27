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
	"io/ioutil"
	"os"
	"path/filepath"

	"github.com/pkg/errors"
	"github.com/shenwei356/unikmer/index"
	"github.com/shenwei356/util/pathutil"
	"gopkg.in/yaml.v2"
)

const dbInfoFile = "__db.yml"

// ErrVersionMismatch indicates mismatched version
var ErrVersionMismatch = errors.New("kmcp/index: version mismatch")

// UnikIndexDBVersion is the version of database.
const UnikIndexDBVersion uint8 = 3

// UnikIndexDBInfo is the meta data of a database.
type UnikIndexDBInfo struct {
	Version      uint8 `yaml:"version"`
	IndexVersion uint8 `yaml:"unikiVersion"`
	K            int   `yaml:"k"`
	Hashed       bool  `yaml:"hashed"`
	Canonical    bool  `yaml:"canonical"`

	Scaled     bool   `yaml:"scaled"`
	Scale      uint32 `yaml:"scale"`
	Minimizer  bool   `yaml:"minimizer"`
	MinimizerW uint32 `yaml:"minimizer-w"`
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
	Paths     []string `yaml:"unik-path"`

	path string
}

func (i UnikIndexDBInfo) String() string {
	return fmt.Sprintf("kmcp database v%d: k: %d, hashed: %v,  canonical: %v, #hashes: %d, fpr:%f, #blocksize: %d, #blocks: %d, #%d-mers: %d",
		i.Version, i.K, i.Hashed, i.Canonical, i.NumHashes, i.FPR, i.BlockSize, len(i.Files), i.K, i.Kmers)
}

// NewUnikIndexDBInfo creates UnikIndexDBInfo from index files, but you have to manually assign other values.
func NewUnikIndexDBInfo(files []string) UnikIndexDBInfo {
	return UnikIndexDBInfo{Version: UnikIndexDBVersion, IndexVersion: index.Version, Files: files}
}

// UnikIndexDBInfoFromFile creates UnikIndexDBInfo from files.
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

// WriteTo dumps UnikIndexDBInfo to file.
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

// Check check if all index files exist.
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
