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
	"bufio"
	"fmt"
	"strconv"
	"strings"

	"github.com/spf13/cobra"
)

var cov2simiCmd = &cobra.Command{
	Use:   "cov2simi",
	Short: "Convert k-mer coverage to sequence similarity",
	Long: `Convert k-mer coverage to sequence similarity

The polynomial model of degree 3 is fitted based on simulated 150-bp E.coli reads.
Visit https://github.com/shenwei356/kmcp/tree/main/analysis/kmer-similarity

    similarity = 87.456 + 26.410*qcov - 22.008*qcov*qcov + 7.325*qcov*qcov*qcov

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)

		var err error

		outFile := getFlagString(cmd, "out-prefix")

		queryCov := getFlagFloat64(cmd, "query-cov")
		if queryCov < 0 || queryCov > 1 {
			checkError(fmt.Errorf("value of -t/--query-cov should be in range [0, 1]"))
		}

		outfh, gw, w, err := outStream(outFile, strings.HasSuffix(strings.ToLower(outFile), ".gz"), opt.CompressionLevel)
		checkError(err)
		defer func() {
			outfh.Flush()
			if gw != nil {
				gw.Close()
			}
			w.Close()
		}()

		if queryCov > 0 {
			fmt.Fprintf(outfh, "%.6f\t%.6f\n", queryCov, similarity(queryCov))
			return
		}

		if opt.Verbose {
			log.Info("the value of flag -t/--query-cov is 0, reading from stdin/file")
		}

		files := getFileListFromArgsAndFile(cmd, args, true, "infile-list", true)
		if opt.Verbose {
			if len(files) == 1 && isStdin(files[0]) {
				log.Info("no files given, reading from stdin")
			} else {
				log.Infof("%d input file(s) given", len(files))
			}
		}

		var line string
		var qcov float64
		for _, file := range files {
			infh, r, _, err := inStream(file)
			checkError(err)

			scanner := bufio.NewScanner(infh)

			for scanner.Scan() {
				line = scanner.Text()
				if line == "" {
					continue
				}

				qcov, err = strconv.ParseFloat(line, 64)
				if err != nil {
					checkError(fmt.Errorf("input should be a number: %s", line))
				}

				fmt.Fprintf(outfh, "%s\t%.6f\n", line, similarity(qcov))
			}

			if err := scanner.Err(); err != nil {
				checkError(err)
			}
			r.Close()
		}

	},
}

func init() {
	utilsCmd.AddCommand(cov2simiCmd)

	cov2simiCmd.Flags().StringP("out-prefix", "o", "-", formatFlagUsage(`Out file prefix ("-" for stdout).`))

	cov2simiCmd.Flags().Float64P("query-cov", "t", 0,
		formatFlagUsage(`K-mer query coverage, i.e., proportion of matched k-mers and unique k-mers of a query. range: [0, 1]`))
}
