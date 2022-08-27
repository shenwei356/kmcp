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
	"bytes"
	"fmt"
	"math"
	"strconv"

	"github.com/shenwei356/bio/taxdump"
	"github.com/shenwei356/util/stats"
)

type MatchResult struct {
	Query   string
	QLen    int
	QKmers  int
	FPR     float64
	Hits    int
	Target  string
	FragIdx int
	IdxNum  int
	GSize   uint64
	K       int
	MKmers  int
	QCov    float64
}

var float64powm10 = []float64{
	1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9,
	1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 1e-16, 1e-17, 1e-18, 1e-19,
	1e-20, 1e-21, 1e-22,
}

func parseQcov(value string) float64 {
	if value[0] == '1' {
		return 1
	}

	var v float64
	for i, b := range value[2:] {
		v += float64powm10[i] * float64(b-'0')
	}
	return v
}

var float64pow10 = []float64{
	1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9,
	1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19,
	1e20, 1e21, 1e22,
}

func parseFPR(value string) float64 {
	var v float64
	v += float64(value[0] - '0')

	var ie int
	for i, b := range value[2:] {
		if b == 'e' {
			ie = i
			break
		}
		v += float64powm10[i] * float64(b-'0')
	}
	ie += 2
	var e float64
	lm1 := len(value) - 1
	end := len(value) - ie - 2
	for i := 0; i < end; i++ {
		e += float64pow10[i] * float64(value[lm1-i]-'0')
	}
	return v * math.Pow10(-int(e))
}

func parseMatchResult(line string, numFields int, items *[]string, maxPFR float64, minQcov float64) (*MatchResult, bool) {
	stringSplitNByByte(line, '\t', numFields, items)
	if len(*items) < numFields {
		checkError(fmt.Errorf("invalid kmcp search result format"))
	}

	m := &MatchResult{} // do not use sync.Pool, which is slower for this case

	var err error

	// slow
	m.FPR, err = strconv.ParseFloat((*items)[3], 64)
	if err != nil {
		checkError(fmt.Errorf("failed to parse FPR: %s", (*items)[3]))
	}
	// did not bring final speedup
	// if len((*items)[3]) < 2 {
	// 	checkError(fmt.Errorf("failed to parse qCov: %s", (*items)[3]))
	// }
	// m.FPR = parseFPR((*items)[3])

	if m.FPR > maxPFR {
		return m, false
	}

	// slow
	m.QCov, err = strconv.ParseFloat((*items)[11], 64)
	if err != nil {
		checkError(fmt.Errorf("failed to parse qCov: %s", (*items)[11]))
	}

	// did not bring final speedup
	// if len((*items)[11]) < 2 {
	// 	checkError(fmt.Errorf("failed to parse qCov: %s", (*items)[11]))
	// }
	// m.QCov = parseQcov((*items)[11])

	if m.QCov < minQcov {
		return m, false
	}

	// -----------

	m.Query = (*items)[0]

	m.QLen, err = strconv.Atoi((*items)[1])
	if err != nil {
		checkError(fmt.Errorf("failed to parse qLen: %s", (*items)[1]))
	}

	m.QKmers, err = strconv.Atoi((*items)[2])
	if err != nil {
		checkError(fmt.Errorf("failed to parse qKmers: %s", (*items)[2]))
	}

	m.Hits, err = strconv.Atoi((*items)[4])
	if err != nil {
		checkError(fmt.Errorf("failed to parse hits: %s", (*items)[4]))
	}

	m.Target = (*items)[5]

	m.FragIdx, err = strconv.Atoi((*items)[6])
	if err != nil {
		checkError(fmt.Errorf("failed to parse chunkIdx: %s", (*items)[6]))
	}

	m.IdxNum, err = strconv.Atoi((*items)[7])
	if err != nil {
		checkError(fmt.Errorf("failed to parse IdxNum: %s", (*items)[7]))
	}

	m.GSize, err = strconv.ParseUint((*items)[8], 10, 64)
	if err != nil {
		checkError(fmt.Errorf("failed to parse genomeSize: %s", (*items)[8]))
	}

	m.K, err = strconv.Atoi((*items)[9])
	if err != nil {
		checkError(fmt.Errorf("failed to parse K: %s", (*items)[9]))
	}

	m.MKmers, err = strconv.Atoi((*items)[10])
	if err != nil {
		checkError(fmt.Errorf("failed to parse mKmers: %s", (*items)[10]))
	}

	return m, true
}

type Target struct {
	Name string

	GenomeSize uint64

	// Counting matches in all chunks
	// some reads match multiple sites in the same genome,
	// the count should be divided by number of sites.
	Match []float64

	// sum of read (query) length
	QLen []float64

	// unique match
	UniqMatch []float64

	// unique match with high confidence
	UniqMatchHic []float64

	SumMatch        float64 // depth
	SumUniqMatch    float64
	SumUniqMatchHic float64

	FragsProp   float64 // coverage
	Coverage    float64
	Qlens       float64
	RelDepth    []float64
	RelDepthStd float64

	//
	RefName string

	// Taxonomy information
	Taxid         uint32
	Rank          string
	TaxonName     string
	LineageNames  []string
	LineageTaxids []string

	CompleteLineageNames  []string
	CompleteLineageTaxids []uint32

	Percentage float64 // relative abundance

	Stats  *stats.Quantiler // for computing percentil of qcov of unique matches
	StatsA *stats.Quantiler // for computing percentil of qcov of all matches

	Score float64
}

func (t *Target) AddTaxonomy(taxdb *taxdump.Taxonomy, showRanksMap map[string]interface{}, taxid uint32) {
	t.Taxid, _ = taxdb.TaxId(taxid)
	t.Rank = taxdb.Rank(taxid)
	t.TaxonName = taxdb.Name(taxid)

	_taxids := taxdb.LineageTaxIds(taxid)

	t.CompleteLineageTaxids = _taxids
	t.CompleteLineageNames = taxdb.LineageNames(taxid)

	var _taxids2 []uint32
	var ok bool
	if len(showRanksMap) > 0 {
		_taxids2 = make([]uint32, 0, len(_taxids))
		for _, _taxid := range _taxids {
			if _, ok = showRanksMap[taxdb.Rank(_taxid)]; ok {
				_taxids2 = append(_taxids2, _taxid)
			}
		}
		_taxids = _taxids2
	}

	t.LineageTaxids = make([]string, len(_taxids))
	for i, _taxid := range _taxids {
		t.LineageTaxids[i] = strconv.Itoa(int(_taxid))
	}

	t.LineageNames = make([]string, len(_taxids))
	for i, _taxid := range _taxids {
		t.LineageNames[i] = taxdb.Names[_taxid]
	}
}

func (t Target) String() string {
	var buf bytes.Buffer
	buf.WriteString(t.Name)
	for i := range t.Match {
		buf.WriteString(fmt.Sprintf(", %d: %.0f(%.0f)", i, t.Match[i], t.UniqMatch[i]))
	}
	buf.WriteString("\n")
	return buf.String()
}

type Targets []*Target

func (t Targets) Len() int { return len(t) }
func (t Targets) Less(i, j int) bool {
	d := t[i].Coverage - t[j].Coverage
	if d == 0 {
		return t[i].FragsProp > t[j].FragsProp
	}
	return t[i].Coverage > t[j].Coverage
}
func (t Targets) Swap(i, j int) {
	t[i], t[j] = t[j], t[i]
}

type ProfileNode struct {
	Taxid         uint32
	Rank          string
	TaxonName     string
	LineageNames  []string // complete lineage
	LineageTaxids []uint32
	Percentage    float64
}

func generateProfile(taxdb *taxdump.Taxonomy, targets []*Target) map[uint32]*ProfileNode {

	profile := make(map[uint32]*ProfileNode, len(targets))

	for _, target := range targets {
		for _, taxid := range target.CompleteLineageTaxids {
			if node, ok := profile[taxid]; !ok {
				profile[taxid] = &ProfileNode{
					Taxid:         taxid,
					Rank:          taxdb.Rank(taxid),
					TaxonName:     taxdb.Names[taxid],
					LineageNames:  taxdb.LineageNames(taxid),
					LineageTaxids: taxdb.LineageTaxIds(taxid),

					Percentage: target.Percentage,
				}
			} else {
				node.Percentage += target.Percentage
			}
		}
	}

	return profile
}
