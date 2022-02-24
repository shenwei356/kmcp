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
	Short: "Compute the maximal false positive rate of a query",
	Long: `Compute the maximal false positive rate of a query

Solomon and Kingsford apply a Chernoff bound and show that the 
false positive probability for a query is:

    fpr ≤ exp( -n(t-f)^2 / (2(1-f)) )

Where:

    f,  the false positive rate of the bloom filters
    t,  the minimal proportion of matched k-mers and unique k-mers of a query
    n,  the number of unique k-mers of the query 

Reference:
  1. SBT: https://doi.org/10.1038/nbt.3442
  2. COBS: https://arxiv.org/abs/1905.09624v2

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)

		var err error

		fpr := getFlagPositiveFloat64(cmd, "false-positive-rate")
		if fpr >= 1 {
			checkError(fmt.Errorf("value of -f/--false-positive-rate too big: %f", fpr))
		}

		queryCov := getFlagFloat64(cmd, "min-query-cov")
		if queryCov < 0 || queryCov > 1 {
			checkError(fmt.Errorf("value of -t/--min-query-cov should be in range [0, 1]"))
		}

		nKmers := getFlagPositiveInt(cmd, "num-kmers")

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

		/*
			p, fpr of single bloom filter.
			k, theshold of query coverage.
			l, number of k-mers.

			import math
			fpr = lambda p,k,l: math.exp(-l * (k - p) * (k - p) / 2 / (1 - p))

			fpr(0.3, 0.8, 60)
		*/

		fmt.Fprintf(outfh, "%s\n", strconv.FormatFloat(maxFPR(fpr, queryCov, nKmers), 'e', 4, 64))

	},
}

func init() {
	utilsCmd.AddCommand(queryFPRCmd)

	queryFPRCmd.Flags().StringP("out-prefix", "o", "-", formatFlagUsage(`Out file prefix ("-" for stdout).`))

	queryFPRCmd.Flags().Float64P("false-positive-rate", "f", 0.3,
		formatFlagUsage(`False positive rate of the bloom filters in the database. range: (0, 1)`))
	queryFPRCmd.Flags().Float64P("min-query-cov", "t", 0.55,
		formatFlagUsage(`Minimal query coverage, i.e., proportion of matched k-mers and unique k-mers of a query. range: [0, 1]`))
	queryFPRCmd.Flags().IntP("num-kmers", "n", 80, formatFlagUsage("Number of unique k-mers of the query."))
}
