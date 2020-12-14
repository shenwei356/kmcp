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
	"io/ioutil"
	"os"
	"path/filepath"

	"github.com/pkg/errors"
	"github.com/shenwei356/unikmer/index"
	"github.com/shenwei356/util/cliutil"
	"github.com/shenwei356/util/pathutil"
	"gopkg.in/yaml.v2"
)

const dbInfoFile = "__db.yml"
const dbNameMappingFile = "__name_mapping.tsv"

// ErrVersionMismatch indicates mismatched version
var ErrVersionMismatch = errors.New("kmcp/index: version mismatch")

// UnikIndexDBVersion is the version of database.
const UnikIndexDBVersion uint8 = 3

// UnikIndexDBInfo is the meta data of a database.
type UnikIndexDBInfo struct {
	Version      uint8  `yaml:"version"`
	IndexVersion uint8  `yaml:"unikiVersion"`
	Alias        string `yaml:"alias"`
	K            int    `yaml:"k"`
	Hashed       bool   `yaml:"hashed"`
	Canonical    bool   `yaml:"canonical"`

	Scaled     bool   `yaml:"scaled"`
	Scale      uint32 `yaml:"scale"`
	Minimizer  bool   `yaml:"minimizer"`
	MinimizerW uint32 `yaml:"minimizer-w"`
	Syncmer    bool   `yaml:"syncmer"`
	SyncmerS   uint32 `yaml:"syncmer-s"`

	SplitSeq     bool `yaml:"split-seq"`
	SplitSize    int  `yaml:"split-size"`
	SplitOverlap int  `yaml:"split-overlap"`

	NumHashes int      `yaml:"hashes"`
	FPR       float64  `yaml:"fpr"`
	BlockSize int      `yaml:"blocksize"`
	Kmers     int      `yaml:"totalKmers"`
	Files     []string `yaml:"files"`
	NumNames  int      `yaml:"numNames"`

	path         string            `yaml:"path,omitempty"`
	NameMapping  map[string]string `yaml:"name-mapping,omitempty"`
	MappingNames bool              `yaml:"mapping-names,omitempty"`
}

func (i UnikIndexDBInfo) String() string {
	return fmt.Sprintf("kmcp database (v%d): %s, k: %d, hashed: %v, canonical: %v, #hashes: %d, fpr:%f, #blocksize: %d, #blocks: %d, #%d-mers: %d",
		i.Version, i.Alias, i.K, i.Hashed, i.Canonical, i.NumHashes, i.FPR, i.BlockSize, len(i.Files), i.K, i.Kmers)
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
		return info, fmt.Errorf("fail to open kmcp database info file: %s", file)
	}

	data, err := ioutil.ReadAll(r)
	if err != nil {
		return info, fmt.Errorf("fail to read kmcp database info file: %s", file)
	}

	err = yaml.Unmarshal(data, &info)
	if err != nil {
		return info, fmt.Errorf("fail to unmarshal kmcp database info")
	}

	r.Close()

	if info.Version != UnikIndexDBVersion {
		return info, ErrVersionMismatch
	}

	p, _ := filepath.Abs(file)
	info.path = filepath.Dir(p)

	fileNameMapping := filepath.Join(filepath.Dir(file), dbNameMappingFile)
	var existed bool
	existed, err = pathutil.Exists(fileNameMapping)
	if existed {
		info.NameMapping, err = cliutil.ReadKVs(fileNameMapping, false)
		checkError(err)
	}
	info.MappingNames = len(info.NameMapping) > 0
	return info, nil
}

// WriteTo dumps UnikIndexDBInfo to file.
func (i UnikIndexDBInfo) WriteTo(file string) (int, error) {
	data, err := yaml.Marshal(i)
	if err != nil {
		return 0, fmt.Errorf("fail to marshal database info")
	}

	w, err := os.Create(file)
	if err != nil {
		return 0, fmt.Errorf("fail to write kmcp database info file: %s", file)
	}
	var n int
	n, err = w.Write(data)
	if err != nil {
		return 0, fmt.Errorf("fail to write kmcp database info file: %s", file)
	}

	w.Close()
	return n, nil
}

// CompatibleWith checks whether two databases have the same parameters.
func (i UnikIndexDBInfo) CompatibleWith(j UnikIndexDBInfo) bool {
	if i.Version == j.Version &&
		i.IndexVersion == j.IndexVersion &&
		i.K == j.K &&
		i.Hashed == j.Hashed &&
		i.Canonical == j.Canonical &&
		i.Scale == j.Scale &&
		i.Scaled == j.Scaled &&
		i.Minimizer == j.Minimizer &&
		i.MinimizerW == j.MinimizerW &&
		i.Syncmer == j.Syncmer &&
		i.SyncmerS == j.SyncmerS {
		return true
	}

	return false
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
