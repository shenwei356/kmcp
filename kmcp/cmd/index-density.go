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
	"image"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/edsrzf/mmap-go"
	"github.com/shenwei356/kmcp/kmcp/cmd/index"
	"github.com/spf13/cobra"

	"image/color"
	"image/jpeg"
	_ "image/jpeg"
)

var densityCmd = &cobra.Command{
	Use:   "index-density",
	Short: "Plot the element density of bloom filters for an index file",
	Long: `Plot the element density of bloom filters for an index file

Purposes:
  1. Checking whether elements (Ones) in bloom filters are uniformly distributed
     via an intuitive grayscale image.

Outputs:
  1. default output (a TSV file), columns:
      1) target:   reference id
      2) chunkIdx: the index of genome chunk
      3) bins:     the number of bins in bloom filters for counting 1s
      4) binSize:  the size/width of a bin
      5) counts:   comma-seperated counts in each bin
  2. the density image (a grayscale JPEG image):
      - X: bins. The width is the number of bins
      - Y: bloom filters, with each representing a genome (chunk).
        The height is the number of names (genome or genome chunks)
      - greyscale/darkness of a pixel: the density of a bin, calculated as:
            255 - 255 * ${the number of 1s in the bin} / ${bin-size}

Examples:
  1. common use:
      kmcp utils index-density gtdb.kmcp/R001/_block001.uniki \
          --bins 1024  --out-file t.tsv --out-img t.jpg
  2. export every bit of each position, the image could fail to create:
      kmcp utils index-density gtdb.kmcp/R001/_block001.uniki \
          --bin-size 1 --out-file t.tsv
`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)

		var err error

		outFile := getFlagString(cmd, "out-file")
		binSize := getFlagNonNegativeInt(cmd, "bin-size")
		bins0 := getFlagNonNegativeInt(cmd, "bins")
		outImg := getFlagString(cmd, "out-img")

		timeStart := time.Now()
		defer func() {
			if opt.Verbose {
				log.Infof("elapsed time: %s", time.Since(timeStart))
			}
		}()

		if opt.Verbose {
			log.Info("checking input files ...")
		}
		files := getFileListFromArgsAndFile(cmd, args, true, "infile-list", true)
		if opt.Verbose {
			if len(files) == 1 && isStdin(files[0]) {
				log.Warningf("stdin not supported, a .uniki file is needed")
				return
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

		// -----------------------------------------------

		file := files[0]

		if filepath.Clean(outFile) == filepath.Clean(file) {
			checkError(fmt.Errorf("intput and output paths should not be the same: %s", file))
		}
		if filepath.Clean(outImg) == filepath.Clean(file) {
			checkError(fmt.Errorf("intput and output paths should not be the same: %s", file))
		}

		fh, err := os.Open(file)
		checkError(err)

		reader, err := index.NewReader(fh)
		checkError(err)

		_offset0, err := fh.Seek(0, 1)
		checkError(err)
		offset0 := int(_offset0)

		sigsM, err := mmap.Map(fh, mmap.RDONLY, 0)
		checkError(err)

		sigs := []byte(sigsM)

		numRowBytes := reader.NumRowBytes
		numSig := int(reader.NumSigs)

		var offset0P1 int
		var ix8 int
		var b int
		var c0, c1, c2, c3, c4, c5, c6, c7 int
		numNames := len(reader.Names)

		var bins int

		if binSize > 0 {
			bins = numSig/binSize + 1
		} else {
			bins = bins0
			binSize = numSig / bins
		}

		if opt.Verbose {
			log.Infof("#names: %d, #sigs: %d", len(reader.Names), reader.NumSigs)
			log.Infof("#bins: %d, bin size: %d", bins, binSize)
		}

		if outImg != "" && bins >= 65536 {
			log.Warningf("the number of bins is too large for ploting: %d.", bins)
		}

		counts := make([][]int, numRowBytes<<3)
		for i := 0; i < len(counts); i++ {
			counts[i] = make([]int, 0, bins)
		}
		var c int
		for i := 0; i < numRowBytes; i++ { // every column in matrix
			offset0P1 = offset0 + i

			c0, c1, c2, c3, c4, c5, c6, c7 = 0, 0, 0, 0, 0, 0, 0, 0
			c = 0
			for loc := 0; loc < numSig; loc++ {
				b = int(sigs[offset0P1+loc*numRowBytes])
				c0 += b >> 7
				c1 += b >> 6 & 1
				c2 += b >> 5 & 1
				c3 += b >> 4 & 1
				c4 += b >> 3 & 1
				c5 += b >> 2 & 1
				c6 += b >> 1 & 1
				c7 += b & 1

				c++
				ix8 = i << 3
				if c == binSize {
					counts[ix8] = append(counts[ix8], c0)
					counts[ix8+1] = append(counts[ix8+1], c1)
					counts[ix8+2] = append(counts[ix8+2], c2)
					counts[ix8+3] = append(counts[ix8+3], c3)
					counts[ix8+4] = append(counts[ix8+4], c4)
					counts[ix8+5] = append(counts[ix8+5], c5)
					counts[ix8+6] = append(counts[ix8+6], c6)
					counts[ix8+7] = append(counts[ix8+7], c7)

					c0, c1, c2, c3, c4, c5, c6, c7 = 0, 0, 0, 0, 0, 0, 0, 0
					c = 0
				}
			}

			ix8 = i << 3
			counts[ix8] = append(counts[ix8], c0)
			counts[ix8+1] = append(counts[ix8+1], c1)
			counts[ix8+2] = append(counts[ix8+2], c2)
			counts[ix8+3] = append(counts[ix8+3], c3)
			counts[ix8+4] = append(counts[ix8+4], c4)
			counts[ix8+5] = append(counts[ix8+5], c5)
			counts[ix8+6] = append(counts[ix8+6], c6)
			counts[ix8+7] = append(counts[ix8+7], c7)
		}

		counts = counts[:numNames]

		var m, M int // min and max count
		m = binSize << 1
		var count []int

		fmt.Fprintf(outfh, "target\tchunkIdx\tbins\tbinSize\tcounts\n")
		for i := 0; i < len(counts); i++ {
			count = counts[i]
			fmt.Fprintf(outfh, "%s\t%d\t%d\t%d\t%d",
				reader.Names[i][0], reader.Indices[i][0]&65535, bins, binSize, count[0])
			if count[0] < m {
				m = count[0]
			}
			if count[0] > M {
				M = count[0]
			}
			for _, c = range count[1 : len(count)-1] {
				fmt.Fprintf(outfh, ",%d", c)
				if c < m {
					m = c
				}
				if c > M {
					M = c
				}
			}
			fmt.Fprint(outfh, "\n")
		}

		if opt.Verbose {
			log.Infof("minimum count in bins of %d: %d (%f)", binSize, m, float64(m)/float64(binSize))
			log.Infof("maximum count in bins of %d: %d (%f)", binSize, M, float64(M)/float64(binSize))
		}

		outfh.Flush()
		checkError(fh.Close())

		// -----------------------------------------------

		if outImg == "" {
			return
		}

		if bins >= 65536 {
			log.Warningf("generating image skipped for a number of bins >= 65536")
			return
		}

		r := float64(255) / float64(binSize)
		img := image.NewGray(image.Rectangle{Max: image.Point{X: bins, Y: numNames}})
		for i := 0; i < len(counts); i++ {
			for j, c := range counts[i] {
				img.SetGray(j, i, color.Gray{255 - uint8(float64(c)*r)})
			}
		}

		f2, err := os.Create(outImg)
		checkError(err)

		checkError(jpeg.Encode(f2, img, &jpeg.Options{Quality: 100}))

		checkError(f2.Close())

		if opt.Verbose {
			log.Infof("out image saved to: %s", outImg)
		}
	},
}

func init() {
	utilsCmd.AddCommand(densityCmd)

	densityCmd.Flags().StringP("out-file", "o", "-", formatFlagUsage(`Out file, supports and recommends a ".gz" suffix ("-" for stdout).`))
	densityCmd.Flags().IntP("bins", "b", 1024, formatFlagUsage(`number of bins for counting the number of 1s.`))
	densityCmd.Flags().IntP("bin-size", "s", 0, formatFlagUsage(`bin size/width`))

	densityCmd.Flags().StringP("out-img", "", "", formatFlagUsage(`Out density image, in format of jpeg`))

}
