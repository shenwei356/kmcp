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
	"os"

	"github.com/shenwei356/go-logging"
	"github.com/shenwei356/util/cliutil"
	"github.com/spf13/cobra"
)

var log = logging.MustGetLogger("kmcp")

func checkError(err error) {
	if err != nil {
		log.Error(err)
		os.Exit(-1)
	}
}

func isStdin(file string) bool {
	return file == "-"
}

func isStdout(file string) bool {
	return file == "-"
}

func getFileListFromArgsAndFile(cmd *cobra.Command, args []string, checkFileFromArgs bool, flag string, checkFileFromFile bool) []string {
	infileList := cliutil.GetFlagString(cmd, flag)
	files := cliutil.GetFileList(args, checkFileFromArgs)
	if infileList != "" {
		_files, err := cliutil.GetFileListFromFile(infileList, checkFileFromFile)
		checkError(err)
		if len(_files) == 0 {
			log.Warningf("no files found in file list: %s", infileList)
			return files
		}

		if len(files) == 1 && isStdin(files[0]) {
			return _files
		}
		files = append(files, _files...)
	}
	return files
}
