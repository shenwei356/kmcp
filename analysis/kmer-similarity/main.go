// -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore mismatch gapopen"
package main

import (
	"fmt"
	"os"
	"regexp"
	"runtime"
	"strconv"
	"strings"
	"sync"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/breader"
	"github.com/shenwei356/unikmer"
	"github.com/shenwei356/xopen"
)

func main() {
	if len(os.Args) < 4 || len(os.Args) > 5 {
		checkError(fmt.Errorf(`usage: %s <sequence file> <blastn result> <k-mer size> [out file]

blastn output format:
    -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore mismatch gapopen gaps"

`, os.Args[0]))
	}

	sfile := os.Args[1]

	sequences, err := getSeqsMap(sfile, seq.Unlimit, 4, 100, "")
	checkError(err)

	file := os.Args[2]

	k, err := strconv.Atoi(os.Args[3])
	checkError(err)

	var outFile string
	if len(os.Args) > 4 {
		outFile = os.Args[4]
	} else {
		outFile = "-"
	}

	// -----------------------------------------------

	// kmers := make(map[string]map[uint64]interface{}, len(sequences))
	// var ok bool
	// for id, record := range sequences {
	// 	kmers[id], ok, err = getKmers(record.Seq, k)
	// 	if err != nil || !ok {
	// 		checkError(fmt.Errorf("fail to compute kmers for %s", id))
	// 	}
	// }

	// -----------------------------------------------

	numFields := 15
	pool := &sync.Pool{New: func() interface{} {
		tmp := make([]string, numFields)
		return &tmp
	}}

	fn := func(line string) (interface{}, bool, error) {
		if line == "" || line[0] == '#' { // ignoring blank line and comment line
			return "", false, nil
		}
		if line[len(line)-1] == '\n' {
			line = line[:len(line)-1]
		}

		items := pool.Get().(*[]string)
		defer pool.Put(items)

		record, ok := parseBlastResult(line, numFields, items)
		if !ok {
			return nil, false, nil
		}
		// fmt.Println(record)

		// -----------------------------------------------------
		qseq, ok := sequences[record.Query]
		if !ok {
			return nil, false, nil
		}
		sseq, ok := sequences[record.Subject]
		if !ok {
			return nil, false, nil
		}

		_qseq := qseq.Seq.SubSeq(record.QueryStart, record.QueryEnd)
		qkmers, ok, err := getKmers(_qseq, k)
		if err != nil {
			return nil, false, err
		}
		if !ok {
			return nil, false, nil
		}

		var s, e int
		if record.SStart < record.SEnd {
			s, e = record.SStart, record.SEnd
		} else {
			e, s = record.SStart, record.SEnd
		}

		_sseq := sseq.Seq.SubSeq(s, e)
		skmers, ok, err := getKmers(_sseq, k)
		if err != nil {
			return nil, false, err
		}
		if !ok {
			return nil, false, nil
		}
		// skmers := kmers[record.Subject]

		var m int
		for kmer := range qkmers {
			if _, ok = skmers[kmer]; ok {
				m++
			}
		}

		match := Match{
			Query:   record.Query,
			Subject: record.Subject,
			QSeq:    _qseq.Seq,
			SSeq:    _sseq.Seq,

			ALen: record.ALen,
			AQov: float64(record.ALen) / float64(len(_qseq.Seq)),

			Pident:   record.Pident,
			KQ:       len(qkmers),
			KS:       len(skmers),
			KM:       m,
			KQcov:    float64(m) / float64(len(qkmers)),
			Mismatch: record.Mismatch,
			Gapopen:  record.Gapopen,
			Gaps:     record.Gaps,
		}

		// -----------------------------------------------------
		return match, true, nil
	}

	reader, err := breader.NewBufferedReader(file, runtime.NumCPU(), 100, fn)
	checkError(err)
	var data interface{}
	var m Match

	outfh, err := xopen.Wopen(outFile)
	checkError(err)
	defer outfh.Close()

	outfh.WriteString("alen\tacov\tpident\tqcov\tqkmers\tskmers\tmatches\tmismatches\tgapopen\tgaps\tquery\tsubject\tqlen\tslen\tqseq\tsseq\n")

	for chunk := range reader.Ch {
		checkError(chunk.Err)

		for _, data = range chunk.Data {
			m = data.(Match)

			// outfh.WriteString(m.String() + "\n")

			outfh.WriteString(fmt.Sprintf("%d\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\n",
				m.ALen, m.AQov, m.Pident, m.KQcov, m.KQ, m.KS, m.KM, m.Mismatch, m.Gapopen, m.Gaps,
				m.Query, m.Subject, len(m.QSeq), len(m.SSeq), m.QSeq, m.SSeq,
			))
		}
	}
}

func getKmers(sequence *seq.Seq, k int) (map[uint64]interface{}, bool, error) {
	iter, err := unikmer.NewHashIterator(sequence, k, true, false)
	if err != nil {
		if err == unikmer.ErrShortSeq {
			return nil, false, nil
		}
		return nil, false, err
	}

	var code uint64
	var ok bool
	kmers := make(map[uint64]interface{}, 256)
	for {
		code, ok = iter.NextHash()
		if !ok {
			break
		}
		kmers[code] = struct{}{}
	}
	return kmers, true, nil
}

type Match struct {
	Query, Subject string
	QSeq, SSeq     []byte

	ALen int
	AQov float64

	Pident     float64
	KQ, KS, KM int
	KQcov      float64

	Mismatch, Gapopen, Gaps int
}

// func (m Match) String() string {
// 	return fmt.Sprintf("%f\t%f\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s",
// 		m.Pident, m.KQcov, m.KQ, m.KS, m.KM, m.Mismatch, m.Gap,
// 		m.Query, m.Subject, m.QSeq, m.SSeq,
// 	)
// }

type BlastRecord struct {
	Query      string
	QueryStart int
	QueryEnd   int
	Subject    string
	Pident     float64

	ALen         int
	QStart, QEnd int
	SStart, SEnd int

	Mismatch int
	Gapopen  int
	Gaps     int
}

// type

var reQLoc = regexp.MustCompile(`^(.+)_sliding:(\d+)\-(\d+)$`)

// BAIM01000080.1_sliding:1-150   NZ_CP033092.2       100.000   150   150   4903501   1   150   1997567   1997716   6.07e-75   278   0    0
// BAIM01000080.1_sliding:1-150   NZ_CP033092.2       96.000    150   150   4903501   2   150   4153180   4153329   2.21e-64   243   5    1
// BAIM01000080.1_sliding:1-150   NZ_JMST01000036.1   100.000   150   150   1165      1   150   1073      924       6.07e-75   278   0    0
// BAIM01000080.1_sliding:1-150   NZ_KK583188.1       100.000   150   150   4906844   1   150   4743371   4743520   6.07e-75   278   0    0
func parseBlastResult(line string, numFields int, items *[]string) (BlastRecord, bool) {
	stringSplitN(line, "\t", numFields, items)
	if len(*items) < numFields {
		checkError(fmt.Errorf("invalid blastn outfmt6 format"))
	}

	var m BlastRecord
	var err error

	m.Query = (*items)[0]
	found := reQLoc.FindAllStringSubmatch((*items)[0], 2)
	if len(found) == 0 {
		checkError(fmt.Errorf("fail to parse query locations: %s", (*items)[0]))
	}
	m.Query = found[0][1]
	m.QueryStart, _ = strconv.Atoi(found[0][2])
	m.QueryEnd, _ = strconv.Atoi(found[0][3])

	m.Subject = (*items)[1]

	if m.Query == m.Subject { // skip self-compare
		return m, false
	}

	m.Pident, err = strconv.ParseFloat((*items)[2], 64)
	if err != nil {
		checkError(fmt.Errorf("fail to parse pident: %s", (*items)[2]))
	}
	m.ALen, err = strconv.Atoi((*items)[3])
	if err != nil {
		checkError(fmt.Errorf("fail to parse alignment length: %s", (*items)[3]))
	}

	m.QStart, err = strconv.Atoi((*items)[6])
	if err != nil {
		checkError(fmt.Errorf("fail to parse qstart: %s", (*items)[6]))
	}
	m.QEnd, err = strconv.Atoi((*items)[7])
	if err != nil {
		checkError(fmt.Errorf("fail to parse qend: %s", (*items)[7]))
	}

	m.SStart, err = strconv.Atoi((*items)[8])
	if err != nil {
		checkError(fmt.Errorf("fail to parse sstart: %s", (*items)[8]))
	}
	m.SEnd, err = strconv.Atoi((*items)[9])
	if err != nil {
		checkError(fmt.Errorf("fail to parse send: %s", (*items)[9]))
	}

	m.Mismatch, err = strconv.Atoi((*items)[12])
	if err != nil {
		checkError(fmt.Errorf("fail to parse mismatch: %s", (*items)[12]))
	}
	m.Gapopen, err = strconv.Atoi((*items)[13])
	if err != nil {
		checkError(fmt.Errorf("fail to parse gapopen: %s", (*items)[13]))
	}
	m.Gaps, err = strconv.Atoi((*items)[14])
	if err != nil {
		checkError(fmt.Errorf("fail to parse gaps: %s", (*items)[14]))
	}

	return m, true
}

func checkError(err error) {
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(-1)
	}
}

func stringSplitN(s string, sep string, n int, a *[]string) {
	if a == nil {
		tmp := make([]string, n)
		a = &tmp
	}

	n--
	i := 0
	for i < n {
		m := strings.Index(s, sep)
		if m < 0 {
			break
		}
		(*a)[i] = s[:m]
		s = s[m+len(sep):]
		i++
	}
	(*a)[i] = s

	(*a) = (*a)[:i+1]
}

func getSeqsMap(file string, alphabet *seq.Alphabet, bufferSize int, chunkSize int, idRegexp string) (map[string]*fastx.Record, error) {
	m := make(map[string]*fastx.Record)
	records, err := fastx.GetSeqs(file, alphabet, bufferSize, chunkSize, idRegexp)
	if err != nil {
		return m, err
	}
	for _, record := range records {
		m[string(record.ID)] = record
	}
	return m, nil
}
