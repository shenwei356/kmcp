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
	"math/bits"
	"os"
	"path/filepath"
	"regexp"
	"runtime"
	"strconv"
	"strings"

	"github.com/iafan/cwalk"
	"github.com/pkg/errors"
	"github.com/shenwei356/util/pathutil"
	"github.com/spf13/cobra"
	"github.com/twotwotwo/sorts"
)

var mapInitSize = 1 << 20 // 1M

// Options contains the global flags
type Options struct {
	NumCPUs int
	Verbose bool

	LogFile  string
	Log2File bool

	Compress         bool
	CompressionLevel int
}

func getOptions(cmd *cobra.Command) *Options {
	threads := getFlagNonNegativeInt(cmd, "threads")
	if threads == 0 {
		threads = runtime.NumCPU()
	}

	sorts.MaxProcs = threads
	runtime.GOMAXPROCS(threads)

	logfile := getFlagString(cmd, "log")
	return &Options{
		NumCPUs: threads,
		// Verbose: getFlagBool(cmd, "verbose"),
		Verbose: !getFlagBool(cmd, "quiet"),

		LogFile:  logfile,
		Log2File: logfile != "",

		Compress:         true,
		CompressionLevel: -1,
	}
}

func checkFileSuffix(opt *Options, suffix string, files ...string) {
	for _, file := range files {
		if isStdin(file) {
			continue
		}

		if suffix != "" && !strings.HasSuffix(file, suffix) {
			checkError(fmt.Errorf("input should be stdin or %s file: %s", suffix, file))
		}
	}
}

func makeOutDir(outDir string, force bool) {
	pwd, _ := os.Getwd()
	if outDir != "./" && outDir != "." && pwd != filepath.Clean(outDir) {
		existed, err := pathutil.DirExists(outDir)
		checkError(errors.Wrap(err, outDir))
		if existed {
			empty, err := pathutil.IsEmpty(outDir)
			checkError(errors.Wrap(err, outDir))
			if !empty {
				if force {
					log.Infof("removing old output directory: %s", outDir)
					checkError(os.RemoveAll(outDir))
				} else {
					checkError(fmt.Errorf("out-dir not empty: %s, use --force to overwrite", outDir))
				}
			} else {
				checkError(os.RemoveAll(outDir))
			}
		}
		checkError(os.MkdirAll(outDir, 0777))
	}
}

func getFileListFromDir(path string, pattern *regexp.Regexp, threads int) ([]string, error) {
	files := make([]string, 0, 512)
	ch := make(chan string, threads)
	done := make(chan int)
	go func() {
		for file := range ch {
			files = append(files, file)
		}
		done <- 1
	}()

	cwalk.NumWorkers = threads
	err := cwalk.WalkWithSymlinks(path, func(_path string, info os.FileInfo, err error) error {
		if err != nil {
			return err
		}
		if !info.IsDir() && pattern.MatchString(info.Name()) {
			ch <- filepath.Join(path, _path)
		}
		return nil
	})
	close(ch)
	<-done
	if err != nil {
		return nil, err
	}

	return files, err
}

func filepathTrimExtension(file string) (string, string) {
	unik := strings.HasSuffix(file, ".unik")
	if unik {
		file = file[0 : len(file)-5]
	}
	gz := strings.HasSuffix(file, ".gz") || strings.HasSuffix(file, ".GZ")
	if gz {
		file = file[0 : len(file)-3]
	}

	fasta := strings.HasSuffix(file, ".fasta") || strings.HasSuffix(file, ".FASTA")
	fastq := strings.HasSuffix(file, ".fastq") || strings.HasSuffix(file, ".FASTQ")
	var fa, fq bool
	if fasta || fastq {
		file = file[0 : len(file)-6]
	} else {
		fa = strings.HasSuffix(file, ".fa") || strings.HasSuffix(file, ".FA") || strings.HasSuffix(file, ".fna") || strings.HasSuffix(file, ".FNA")
		fq = strings.HasSuffix(file, ".fq") || strings.HasSuffix(file, ".FQ")
	}

	extension := filepath.Ext(file)
	name := file[0 : len(file)-len(extension)]
	switch {
	case fasta:
		extension += ".fasta"
	case fastq:
		extension += ".fastq"
	case fa:
		extension += ".fa"
	case fq:
		extension += ".fq"
	}
	if gz {
		extension += ".gz"
	}
	if unik {
		extension += ".unik"
	}
	return name, extension
}

func roundup32(x uint32) uint32 {
	if x == 0 {
		return 1
	}
	x--
	x |= x >> 1
	x |= x >> 2
	x |= x >> 4
	x |= x >> 8
	x |= x >> 16
	return (x | x>>32) + 1
}

func roundup64(x uint64) uint64 {
	if x == 0 {
		return 1
	}
	x--
	x |= x >> 1
	x |= x >> 2
	x |= x >> 4
	x |= x >> 8
	x |= x >> 16
	x |= x >> 32
	return (x | x>>64) + 1
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

func bytesSplitN(s []byte, sep []byte, n int, a *[][]byte) {
	if a == nil {
		tmp := make([][]byte, n)
		a = &tmp
	}

	n--
	i := 0
	for i < n {
		m := bytes.Index(s, sep)
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

func stringSplitNByByte(s string, sep byte, n int, a *[]string) {
	if a == nil {
		tmp := make([]string, n)
		a = &tmp
	}

	n--
	i := 0
	for i < n {
		m := strings.IndexByte(s, sep)
		if m < 0 {
			break
		}
		(*a)[i] = s[:m]
		s = s[m+1:]
		i++
	}
	(*a)[i] = s

	(*a) = (*a)[:i+1]
}

// modify from https://github.com/mxschmitt/golang-combinations/blob/master/combinations.go
// too slow for big n.
func Combinations(set []uint64, n int) (subsets [][]uint64) {
	length := uint(len(set))

	if n > len(set) {
		n = len(set)
	}

	// Go through all possible combinations of objects
	// from 1 (only first object in subset) to 2^length (all objects in subset)
	for subsetBits := 1; subsetBits < (1 << length); subsetBits++ {
		if n > 0 && bits.OnesCount(uint(subsetBits)) != n {
			continue
		}

		var subset []uint64

		for object := uint(0); object < length; object++ {
			// checks if object is contained in subset
			// by checking if bit 'object' is set in subsetBits
			if (subsetBits>>object)&1 == 1 {
				// add object to subset
				subset = append(subset, set[object])
			}
		}
		// add subset to subsets
		subsets = append(subsets, subset)
	}
	return subsets
}

type Uint64Slice []uint64

func (s Uint64Slice) Len() int           { return len(s) }
func (s Uint64Slice) Less(i, j int) bool { return s[i] < s[j] }
func (s Uint64Slice) Swap(i, j int)      { s[i], s[j] = s[j], s[i] }

func (s *Uint64Slice) Push(x interface{}) {
	*s = append(*s, x.(uint64))
}

func (s *Uint64Slice) Pop() interface{} {
	old := *s
	n := len(old)
	x := old[n-1]
	*s = old[0 : n-1]
	return x
}

// Note: set should not have duplicates
func Combinations2(set []uint64) [][2]uint64 {
	if len(set) < 2 {
		return nil
	}
	var i, j int
	var n = len(set)
	var np1 = n - 1
	comb := make([][2]uint64, 0, n*(n-1)/2)
	for i = 0; i < np1; i++ {
		for j = i + 1; j < n; j++ {
			comb = append(comb, [2]uint64{set[i], set[j]})
		}
	}
	return comb
}

func sortTwoUint64s(a, b uint64) (uint64, uint64) {
	if a < b {
		return a, b
	}
	return b, a
}

func minInt(a int, vals ...int) int {
	min := a
	for _, v := range vals {
		if v < min {
			min = v
		}
	}
	return min
}

func maxInt(a int, vals ...int) int {
	min := a
	for _, v := range vals {
		if v > min {
			min = v
		}
	}
	return min
}

func IntSlice2StringSlice(vals []int) []string {
	s := make([]string, len(vals))
	for i, v := range vals {
		s[i] = strconv.Itoa(v)
	}
	return s
}

func MeanStdev(values []float64) (float64, float64) {
	n := len(values)

	if n == 0 {
		return 0, 0
	}

	if n == 1 {
		return values[0], 0
	}

	var sum float64
	for _, v := range values {
		sum += v
	}

	mean := sum / float64(n)

	var variance float64
	for _, v := range values {
		variance += (v - mean) * (v - mean)
	}

	return mean, math.Sqrt(variance / float64(n))
}

func wrapByteSlice(s []byte, width int, buffer *bytes.Buffer) ([]byte, *bytes.Buffer) {
	if width < 1 {
		return s, buffer
	}
	l := len(s)
	if l == 0 {
		return s, buffer
	}

	var lines int
	if l%width == 0 {
		lines = l/width - 1
	} else {
		lines = int(l / width)
	}

	if buffer == nil {
		buffer = bytes.NewBuffer(make([]byte, 0, l+lines))
	} else {
		buffer.Reset()
	}

	var start, end int
	for i := 0; i <= lines; i++ {
		start = i * width
		end = (i + 1) * width
		if end > l {
			end = l
		}

		buffer.Write(s[start:end])
		if i < lines {
			buffer.Write(_mark_newline)
		}
	}
	return buffer.Bytes(), buffer
}

var _mark_fasta = []byte{'>'}
var _mark_fastq = []byte{'@'}
var _mark_plus_newline = []byte{'+', '\n'}
var _mark_newline = []byte{'\n'}
