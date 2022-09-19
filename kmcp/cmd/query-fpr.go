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

	"github.com/spf13/cobra"
)

var queryFPRCmd = &cobra.Command{
	Use:   "query-fpr",
	Short: "Compute the false positive rate of a query",
	Long: `Compute the false positive rate of a query

When the flag '-a/--all' is given, the Chernoff bound (column 'cbound')
is also output along with input parameters.

> Given K ≥ p, Solomon and Kingsford also apply a Chernoff bound and
show that the false positive probability for a query to be detected
in a document is ≤ exp(−l(K − p)^2 /(2(1 − p)))

Reference:
  1. Theorem 2 in https://doi.org/10.1038/nbt.3442
  2. Theorem 1 in https://arxiv.org/abs/1905.09624v2

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)

		var err error

		fpr := getFlagPositiveFloat64(cmd, "false-positive-rate")
		if fpr >= 1 {
			checkError(fmt.Errorf("value of -f/--false-positive-rate too big: %f", fpr))
		}

		mKmers := getFlagPositiveInt(cmd, "matched-kmers")

		nKmers := getFlagPositiveInt(cmd, "num-kmers")

		all := getFlagBool(cmd, "all")
		addHeader := getFlagBool(cmd, "add-header")

		outFile := getFlagString(cmd, "out-prefix")

		outfh, gw, w, err := outStream(outFile, strings.HasSuffix(strings.ToLower(outFile), ".gz"), opt.CompressionLevel)
		checkError(err)
		defer func() {
			outfh.Flush()
			if gw != nil {
				gw.Close()
			}
			w.Close()
		}()

		_fpr := QueryFPR(nKmers, mKmers, fpr)
		_fpr2 := maxFPR(fpr, float64(mKmers)/float64(nKmers), nKmers) // Chernoff bound

		if all {
			if addHeader {
				fmt.Fprintf(outfh, "fpr\tcbound\tfpr0\tnKmers\tmKmers\n")
			}
			fmt.Fprintf(outfh, "%s\t%f\t%f\t%d\t%d\n", strconv.FormatFloat(_fpr, 'e', 4, 64), _fpr2, fpr, nKmers, mKmers)
		} else {
			if addHeader {
				fmt.Fprintf(outfh, "fpr\n")
			}
			fmt.Fprintf(outfh, "%s\n", strconv.FormatFloat(_fpr, 'e', 4, 64))
		}

	},
}

func init() {
	utilsCmd.AddCommand(queryFPRCmd)

	queryFPRCmd.Flags().BoolP("add-header", "H", false, formatFlagUsage(`Add header line (column names`))
	queryFPRCmd.Flags().BoolP("all", "a", false, formatFlagUsage(`Also show the value of -f, -n, and -t`))

	queryFPRCmd.Flags().StringP("out-prefix", "o", "-", formatFlagUsage(`Out file prefix ("-" for stdout).`))

	queryFPRCmd.Flags().Float64P("false-positive-rate", "f", 0.3,
		formatFlagUsage(`False positive rate of a single k-mer, i.e., FPR of the bloom filters in the database. range: (0, 1)`))
	queryFPRCmd.Flags().IntP("matched-kmers", "m", 35,
		formatFlagUsage(`The number of matched k-mers of a query.`))
	queryFPRCmd.Flags().IntP("num-kmers", "n", 70, formatFlagUsage("Number of unique k-mers of the query."))
}
