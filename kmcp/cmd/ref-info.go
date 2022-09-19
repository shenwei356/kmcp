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
	"os"
	"path/filepath"
	"strings"

	"github.com/pkg/errors"
	"github.com/shenwei356/kmcp/kmcp/cmd/index"
	"github.com/shenwei356/util/pathutil"
	"github.com/spf13/cobra"
)

var refInfoCmd = &cobra.Command{
	Use:   "ref-info",
	Short: "Print information of reference chunks in a database",
	Long: `Print information of reference chunks in a database

Columns:

    file,     the base name of index file
    i,        the idx of a reference chunk in the index file, 1-based
    target,   reference name
    chunkIdx, the idx of the chunk, 0-based
    chunks,   the number of chunks of the reference
    kmers,    the number of k-mers of the chunk
    fpr,      the actual false-positive reate of the chunk

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)

		dbDir := getFlagString(cmd, "db-dir")
		if dbDir == "" {
			checkError(fmt.Errorf("flag -d/--db-dir needed"))
		}

		outFile := getFlagString(cmd, "out-file")

		// basename := getFlagBool(cmd, "basename")

		noHeaderRow := getFlagBool(cmd, "no-header-row")

		// ---------------------------------------------------------------
		// check Database

		subFiles, err := os.ReadDir(dbDir)
		if err != nil {
			checkError(fmt.Errorf("read database error: %s", err))
		}

		dbDirs := make([]string, 0, 8)
		for _, file := range subFiles {
			if file.Name() == "." || file.Name() == ".." {
				continue
			}
			path := filepath.Join(dbDir, file.Name())

			if !file.IsDir() {
				continue
			}
			existed, err := pathutil.Exists(filepath.Join(path, dbInfoFile))
			if err != nil {
				checkError(fmt.Errorf("read database error: %s", err))
			}
			if existed {
				dbDirs = append(dbDirs, path)
			}
		}
		if len(dbDirs) == 0 {
			checkError(fmt.Errorf("invalid kmcp database: %s", dbDir))
		}

		// ---------------------------------------------------------------

		outfh, gw, w, err := outStream(outFile, strings.HasSuffix(strings.ToLower(outFile), ".gz"), opt.CompressionLevel)
		checkError(err)
		defer func() {
			outfh.Flush()
			if gw != nil {
				gw.Close()
			}
			w.Close()
		}()

		if !noHeaderRow {
			fmt.Fprintf(outfh, "file\ti\ttarget\tchunkIdx\tchunks\tkmers\tfpr\n")
		}

		for _, dbDir := range dbDirs {
			info, err := UnikIndexDBInfoFromFile(filepath.Join(dbDir, dbInfoFile))
			if err != nil {
				checkError(err)
			}

			if len(info.Files) == 0 {
				checkError(fmt.Errorf("no index files"))
			}

			err = info.Check()
			if err != nil {
				checkError(err)
			}

			// --------------

			for _, file := range info.Files {
				infh, r, _, err := inStream(filepath.Join(dbDir, file))
				checkError(err)

				reader, err := index.NewReader(infh)
				checkError(errors.Wrap(err, file))

				h := reader.Header

				// if basename {
				// 	file = filepath.Base(file)
				// }

				var fpr float64
				var idx uint32
				names := h.Names
				indices := h.Indices
				for i, nKmers := range h.Sizes {
					fpr = CalcFPR(nKmers, int(h.NumHashes), h.NumSigs)
					idx = indices[i][0]
					fmt.Fprintf(outfh, "%s\t%d\t%s\t%d\t%d\t%d\t%f\n", file, i+1, names[i][0], uint16(idx), idx>>16, nKmers, fpr)
				}

				r.Close()

			}
		}

	},
}

func init() {
	utilsCmd.AddCommand(refInfoCmd)

	refInfoCmd.Flags().StringP("db-dir", "d", "",
		formatFlagUsage(`Database directory created by "kmcp index".`))

	refInfoCmd.Flags().StringP("out-file", "o", "-", formatFlagUsage(`Out file, supports and recommends a ".gz" suffix ("-" for stdout).`))

	refInfoCmd.Flags().BoolP("no-header-row", "H", false,
		formatFlagUsage(`Do not print header row.`))

	// refInfoCmd.Flags().BoolP("basename", "b", false, formatFlagUsage("Only output basenames of files."))

}
