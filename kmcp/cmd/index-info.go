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
	"path/filepath"
	"strings"

	"github.com/pkg/errors"
	"github.com/shenwei356/kmcp/kmcp/cmd/index"
	"github.com/spf13/cobra"
)

var infoCmd = &cobra.Command{
	Use:   "info",
	Short: "Information of index file",
	Long: `Information of index file

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)

		var err error

		outFile := getFlagString(cmd, "out-prefix")
		all := getFlagBool(cmd, "all")
		basename := getFlagBool(cmd, "basename")

		if opt.Verbose {
			log.Info("checking input files ...")
		}
		files := getFileListFromArgsAndFile(cmd, args, true, "infile-list", true)
		if opt.Verbose {
			if len(files) == 1 && isStdin(files[0]) {
				log.Info("no files given, reading from stdin")
			} else {
				log.Infof("%d input file(s) given", len(files))
			}
		}

		checkFileSuffix(opt, extIndex, files...)

		outfh, gw, w, err := outStream(outFile, strings.HasSuffix(strings.ToLower(outFile), ".gz"), opt.CompressionLevel)
		checkError(err)
		defer func() {
			outfh.Flush()
			if gw != nil {
				gw.Close()
			}
			w.Close()
		}()

		if all {
			outfh.WriteString(fmt.Sprintf("file\tk\tcanonical\tnum-hashes\tnum-sigs\tnum-names\tnames\tindices\n"))
		} else {
			outfh.WriteString(fmt.Sprintf("file\tk\tcanonical\tnum-hashes\tnum-sigs\tnum-names\n"))
		}
		for _, file := range files {
			infh, r, _, err := inStream(file)
			checkError(err)

			reader, err := index.NewReader(infh)
			checkError(errors.Wrap(err, file))

			h := reader.Header

			if basename {
				file = filepath.Base(file)
			}

			if all {
				names := make([]string, len(reader.Names))
				for _i, _names := range reader.Names {
					names[_i] = strings.Join(_names, ",")
				}

				indices := make([]string, len(reader.Names))
				for _i, _indices := range reader.Indices {
					_indicesS := make([]string, len(_indices))
					for _j, _idx := range _indices {
						_indicesS[_j] = fmt.Sprintf("%d", _idx)
					}
					indices[_i] = strings.Join(_indicesS, ",")
				}

				outfh.WriteString(fmt.Sprintf("%s\t%d\t%v\t%d\t%d\t%d\t%s\t%s\n",
					file, h.K, h.Canonical, h.NumHashes, h.NumSigs, len(h.Names),
					strings.Join(names, "; "), strings.Join(indices, "; ")))
			} else {
				outfh.WriteString(fmt.Sprintf("%s\t%d\t%v\t%d\t%d\t%d\n", file, h.K, h.Canonical, h.NumHashes, h.NumSigs, len(h.Names)))
			}
			outfh.Flush()

			r.Close()
		}

	},
}

func init() {
	RootCmd.AddCommand(infoCmd)

	infoCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
	infoCmd.Flags().BoolP("all", "a", false, "all information")
	infoCmd.Flags().BoolP("basename", "b", false, "only output basenames of files")
}
