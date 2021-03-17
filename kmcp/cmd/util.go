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
	"os"
	"path/filepath"
	"regexp"
	"runtime"
	"strings"

	"github.com/iafan/cwalk"
	"github.com/pkg/errors"
	"github.com/shenwei356/util/cliutil"
	"github.com/shenwei356/util/pathutil"
	"github.com/spf13/cobra"
	"github.com/twotwotwo/sorts"
)

var mapInitSize = 1 << 20 // 1M

// Options contains the global flags
type Options struct {
	NumCPUs int
	Verbose bool

	Compress         bool
	CompressionLevel int
}

func getOptions(cmd *cobra.Command) *Options {
	threads := cliutil.GetFlagNonNegativeInt(cmd, "threads")
	if threads == 0 {
		threads = runtime.NumCPU()
	}

	sorts.MaxProcs = threads
	runtime.GOMAXPROCS(threads)

	return &Options{
		NumCPUs: threads,
		// Verbose: cliutil.GetFlagBool(cmd, "verbose"),
		Verbose: !cliutil.GetFlagBool(cmd, "quiet"),

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
