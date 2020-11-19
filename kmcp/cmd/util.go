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
	"os"
	"path/filepath"
	"runtime"
	"strings"

	"github.com/pkg/errors"
	"github.com/shenwei356/util/cliutil"
	"github.com/shenwei356/util/pathutil"
	"github.com/spf13/cobra"
	"github.com/twotwotwo/sorts"
)

var mapInitSize = 1 << 20 // 1M

const (
	flagContinue = iota
	flagBreak
	flagReturn
)

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

	return &Options{
		NumCPUs: threads,
		Verbose: cliutil.GetFlagBool(cmd, "verbose"),

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
