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
	"strings"

	"github.com/shenwei356/breader"
)

const extIndex = ".uniki"

// UnikFileInfo store basic info of .unik file.
type UnikFileInfo struct {
	Path       string
	Name       string
	GenomeSize uint64
	Index      uint32
	Indexes    uint32
	Kmers      uint64
}

func (i UnikFileInfo) String() string {
	return fmt.Sprintf("UnikFile{Kmers: %d, Path: %s, Name: %s}", i.Kmers, i.Path, i.Name)
}

// UnikFileInfos is list of UnikFileInfo.
type UnikFileInfos []UnikFileInfo

func (l UnikFileInfos) Len() int               { return len(l) }
func (l UnikFileInfos) Less(i int, j int) bool { return l[i].Kmers < l[j].Kmers }
func (l UnikFileInfos) Swap(i int, j int)      { l[i], l[j] = l[j], l[i] }

// UnikFileInfosByName is used to sort infos by name and indices
type UnikFileInfosByName []UnikFileInfo

func (l UnikFileInfosByName) Len() int { return len(l) }
func (l UnikFileInfosByName) Less(i int, j int) bool {
	v := strings.Compare(l[i].Name, l[j].Name)
	switch {
	case v < 0:
		return true
	case v > 0:
		return false
	}
	return l[i].Index < l[j].Index
}
func (l UnikFileInfosByName) Swap(i int, j int) { l[i], l[j] = l[j], l[i] }

// UnikFileInfoGroup represents a slice of UnikFileInfos
type UnikFileInfoGroup struct {
	Infos []UnikFileInfo
	Kmers uint64
}

func (i UnikFileInfoGroup) String() string {
	return fmt.Sprintf("UnikFileGroups{Files: %d, Kmers: %d}", len(i.Infos), i.Kmers)
}

// UnikFileInfoGroups is just a slice of UnikFileInfoGroup
type UnikFileInfoGroups []UnikFileInfoGroup

func (l UnikFileInfoGroups) Len() int               { return len(l) }
func (l UnikFileInfoGroups) Less(i int, j int) bool { return l[i].Kmers < l[j].Kmers }
func (l UnikFileInfoGroups) Swap(i int, j int)      { l[i], l[j] = l[j], l[i] }

var fnParseUnikInfoFile = func(line string) (interface{}, bool, error) {
	if len(line) > 0 && line[len(line)-1] == '\n' {
		line = line[:len(line)-1]
	}
	if len(line) == 0 {
		return nil, false, nil
	}
	if line[0] == '#' {
		return nil, false, nil
	}

	items := strings.Split(line, "\t")
	if len(items) < 5 {
		return nil, false, nil
	}
	idx, err := strconv.Atoi(items[2])
	if err != nil || idx < 0 {
		return nil, false, err
	}
	idxNum, err := strconv.Atoi(items[3])
	if err != nil || idxNum < 0 {
		return nil, false, err
	}
	gSize, err := strconv.ParseUint(items[4], 10, 64)
	if err != nil {
		return nil, false, err
	}
	kmers, err := strconv.Atoi(items[5])
	if err != nil || kmers < 0 {
		return nil, false, err
	}
	return UnikFileInfo{
		Path:       items[0],
		Name:       items[1],
		Index:      uint32(idx),
		Indexes:    uint32(idxNum),
		GenomeSize: gSize,
		Kmers:      uint64(kmers),
	}, true, nil
}

func readUnikFileInfos(file string) ([]UnikFileInfo, error) {
	infos := make([]UnikFileInfo, 0, mapInitSize)

	reader, err := breader.NewBufferedReader(file, 4, 100, fnParseUnikInfoFile)
	if err != nil {
		return nil, err
	}
	var data interface{}
	for chunk := range reader.Ch {
		if chunk.Err != nil {
			return nil, err
		}
		for _, data = range chunk.Data {
			infos = append(infos, data.(UnikFileInfo))
		}
	}
	return infos, nil
}

func dumpUnikFileInfos(fileInfos []UnikFileInfo, file string) {
	outfh, gw, w, err := outStream(file, false, -1)
	checkError(err)
	defer func() {
		outfh.Flush()
		if gw != nil {
			gw.Close()
		}
		w.Close()
	}()

	outfh.WriteString("#path\tname\tfragIdx\tidxNum\tgSize\tkmers\n")
	for _, info := range fileInfos {
		outfh.WriteString(fmt.Sprintf("%s\t%s\t%d\t%d\t%d\t%d\n", info.Path, info.Name, info.Index, info.Indexes, info.GenomeSize, info.Kmers))
	}
}

// Meta contains some meta information
type Meta struct {
	SeqID      string `json:"id"`   // sequence ID
	FragIdx    uint32 `json:"idx"`  // sequence location index
	GenomeSize uint64 `json:"gn-s"` // genome length

	Ks []int `json:"ks"` // ks

	Syncmer  bool `json:"sm"` // syncmer
	SyncmerS int  `json:"sm-s"`

	Minimizer  bool `json:"mm"` // minimizer
	MinimizerW int  `json:"mm-w"`

	SplitSeq     bool `json:"sp"` // split sequence
	SplitSize    int  `json:"sp-s"`
	SplitNum     int  `json:"sp-n"`
	SplitOverlap int  `json:"sp-o"`
}

func (m Meta) String() string {
	return fmt.Sprintf("seqID: %s; fragIdx: %d; genome-len: %d, syncmer: %v, size %d; minimizer: %v, window: %d; split-seq: %v, number:%d / size: %d, overlap: %d",
		m.SeqID, m.FragIdx, m.GenomeSize, m.Syncmer, m.SyncmerS, m.Minimizer, m.MinimizerW, m.SplitSeq, m.SplitNum, m.SplitSize, m.SplitOverlap)
}
