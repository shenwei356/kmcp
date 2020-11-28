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
	"io"
	"os"
	"path/filepath"
	"runtime"
	"sync"
	"time"

	"github.com/pkg/errors"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/unikmer"
	"github.com/spf13/cobra"
	"github.com/twotwotwo/sorts/sortutil"
	"github.com/vbauerster/mpb"
	"github.com/vbauerster/mpb/decor"
)

var computeCmd = &cobra.Command{
	Use:   "compute",
	Short: "Generate k-mers (sketch) from FASTA/Q sequences",
	Long: `Generate k-mers (sketchs) from FASTA/Q sequences

Attentions:
  1. K-mers (sketchs) are not sorted and duplicates remained.

K-mer:
  1. ntHash of k-mer (-k)

K-mer sketchs:
  1. Scaled MinHash (-k -D)
  2. Minimizer      (-k -W), optionally scaling/down-sampling (-D)
  3. Syncmer        (-k -S), optionally scaling/down-sampling (-D)

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)
		seq.ValidateSeq = false

		timeStart := time.Now()
		defer func() {
			if opt.Verbose {
				log.Infof("elapsed time: %s", time.Since(timeStart))
			}
		}()

		// ---------------------------------------------------------------
		// basic flags

		var err error

		k := getFlagPositiveInt(cmd, "kmer-len")
		circular := getFlagBool(cmd, "circular")

		outFile := getFlagString(cmd, "out-prefix")
		outDir := getFlagString(cmd, "out-dir")
		force := getFlagBool(cmd, "force")
		exactNumber := getFlagBool(cmd, "exact-number")

		if outDir == "" && outFile == "" {
			checkError(fmt.Errorf("flag -o/--out-prefix OR -O/--out-dir is needed"))
		} else if outDir != "" && outFile != "" {
			checkError(fmt.Errorf("either flag -o/--out-prefix OR -O/--out-dir is allowed"))
		}

		// ---------------------------------------------------------------
		// flags of sketch

		scale := getFlagPositiveInt(cmd, "scale")
		if scale > 1<<32-1 {
			checkError(fmt.Errorf("value of flag --scale is too big"))
		}
		scaled := scale > 1
		maxHash := uint64(float64(^uint64(0)) / float64(scale))

		minimizerW := getFlagNonNegativeInt(cmd, "minimizer-w")
		if minimizerW > 1<<32-1 {
			checkError(fmt.Errorf("value of flag --minimizer-w is too big"))
		}
		minimizer := minimizerW > 0

		syncmerS := getFlagNonNegativeInt(cmd, "syncmer-s")
		if syncmerS > 1<<32-1 {
			checkError(fmt.Errorf("value of flag --syncmer-s is too big"))
		}
		syncmer := syncmerS > 0

		if minimizer && syncmer {
			checkError(fmt.Errorf("flag --minimizer-w and --syncmer-s can not be given simultaneously"))
		}

		// ---------------------------------------------------------------
		// out dir

		outputDir := outDir != ""
		if outputDir {
			makeOutDir(outDir, force)
		}

		// ---------------------------------------------------------------
		// input files

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

		// ---------------------------------------------------------------
		// log

		if opt.Verbose {
			log.Infof("-------------------- [main parameters] --------------------")
			log.Infof("k: %d", k)
			log.Infof("circular genome: %v", circular)
			if scaled {
				log.Infof("down-sampling scale: %d", scale)
			}
			if minimizer {
				log.Infof("minizimer window: %d", minimizerW)
			}
			if syncmer {
				log.Infof("syncmer s: %d", syncmerS)
			}
			log.Infof("-------------------- [main parameters] --------------------")
		}

		// ---------------------------------------------------------------
		// compute for every file and output to a directory
		if outputDir {

			// process bar
			var pbs *mpb.Progress
			var bar *mpb.Bar
			var chDuration chan time.Duration
			var doneDuration chan int

			if opt.Verbose {
				pbs = mpb.New(mpb.WithWidth(60), mpb.WithOutput(os.Stderr))
				bar = pbs.AddBar(int64(len(files)),
					mpb.BarStyle("[=>-]<+"),
					mpb.PrependDecorators(
						decor.Name("processing file: ", decor.WC{W: len("processing file: "), C: decor.DidentRight}),
						decor.Name("", decor.WCSyncSpaceR),
						decor.CountersNoUnit("%d / %d", decor.WCSyncWidth),
					),
					mpb.AppendDecorators(
						decor.Name("ETA: ", decor.WC{W: len("ETA: ")}),
						decor.EwmaETA(decor.ET_STYLE_GO, 60),
						decor.OnComplete(decor.Name(""), ". done"),
					),
				)

				chDuration = make(chan time.Duration, opt.NumCPUs)
				doneDuration = make(chan int)
				go func() {
					for t := range chDuration {
						bar.Increment()
						bar.DecoratorEwmaUpdate(t)
					}
					doneDuration <- 1
				}()
			}

			// wait group
			var wg sync.WaitGroup
			tokens := make(chan int, opt.NumCPUs)

			for _, file := range files {
				tokens <- 1
				wg.Add(1)

				go func(file string) {
					startTime := time.Now()

					outFile := filepath.Join(outDir, fmt.Sprintf("%s%s", filepath.Base(file), extDataFile))
					outfh, gw, w, err := outStream(outFile, opt.Compress, opt.CompressionLevel)
					checkError(err)
					defer func() {
						checkError(outfh.Flush())
						if gw != nil {
							checkError(gw.Close())
						}
						checkError(w.Close())

						wg.Done()
						<-tokens

						if opt.Verbose {
							chDuration <- time.Duration(float64(time.Since(startTime)) / float64(opt.NumCPUs))
						}
					}()

					var writer *unikmer.Writer
					var mode uint32
					var n int

					mode |= unikmer.UnikCanonical
					mode |= unikmer.UnikHashed

					writer, err = unikmer.NewWriter(outfh, k, mode)
					checkError(errors.Wrap(err, outFile))
					if scaled {
						writer.SetScale(uint32(scale))
					}

					var record *fastx.Record
					var fastxReader *fastx.Reader
					var ok bool
					var code uint64
					var iter *unikmer.Iterator
					var sketch *unikmer.Sketch
					codes := make([]uint64, 0, mapInitSize)

					fastxReader, err = fastx.NewDefaultReader(file)
					checkError(errors.Wrap(err, file))

					for {
						record, err = fastxReader.Read()
						if err != nil {
							if err == io.EOF {
								break
							}
							checkError(errors.Wrap(err, file))
							break
						}

						if syncmer {
							sketch, err = unikmer.NewSyncmerSketch(record.Seq, k, syncmerS, circular)
						} else if minimizer {
							sketch, err = unikmer.NewMinimizerSketch(record.Seq, k, minimizerW, circular)
						} else {
							iter, err = unikmer.NewHashIterator(record.Seq, k, true, circular)
						}
						if err != nil {
							if err == unikmer.ErrShortSeq {
								continue
							} else {
								checkError(errors.Wrapf(err, "seq: %s", record.Name))
							}
						}

						if syncmer {
							for {
								code, ok = sketch.NextSyncmer()
								if !ok {
									break
								}
								if scaled && code > maxHash {
									continue
								}
								codes = append(codes, code)
							}
							continue
						}

						if minimizer {
							for {
								code, ok = sketch.NextMinimizer()
								if !ok {
									break
								}
								if scaled && code > maxHash {
									continue
								}
								codes = append(codes, code)
							}
							continue
						}

						for {
							code, ok = iter.NextHash()
							if !ok {
								break
							}
							if scaled && code > maxHash {
								continue
							}
							codes = append(codes, code)
						}
					}

					n = len(codes)

					if exactNumber {
						// compute exact number of unique k-mers
						codes2 := make([]uint64, n)
						copy(codes2, codes)
						sortutil.Uint64s(codes2)
						var pre uint64 = 0
						n = 0
						for _, code = range codes2 {
							if code != pre {
								n++
								pre = code
							}
						}
					}

					writer.Number = int64(n)

					for _, code = range codes {
						writer.WriteCode(code)
					}

					checkError(writer.Flush())
				}(file)
			}

			wg.Wait()

			if opt.Verbose {
				close(chDuration)
				<-doneDuration
				pbs.Wait()
			}

			return
		}

		// ---------------------------------------------------------------
		// compute from all input files and output to a single file

		if !isStdout(outFile) {
			outFile += extDataFile
		}
		outfh, gw, w, err := outStream(outFile, opt.Compress, opt.CompressionLevel)
		checkError(err)
		defer func() {
			checkError(outfh.Flush())
			if gw != nil {
				checkError(gw.Close())
			}
			checkError(w.Close())
		}()

		var writer *unikmer.Writer
		var mode uint32
		var n int

		mode |= unikmer.UnikCanonical
		mode |= unikmer.UnikHashed

		writer, err = unikmer.NewWriter(outfh, k, mode)
		checkError(errors.Wrap(err, outFile))
		if scaled {
			writer.SetScale(uint32(scale))
		}

		var record *fastx.Record
		var fastxReader *fastx.Reader
		var ok bool
		var code uint64
		var iter *unikmer.Iterator
		var sketch *unikmer.Sketch
		codes := make([]uint64, 0, mapInitSize)

		var pbs *mpb.Progress
		var bar *mpb.Bar
		var startTime time.Time
		if opt.Verbose {
			pbs = mpb.New(mpb.WithWidth(60), mpb.WithOutput(os.Stderr))
			bar = pbs.AddBar(int64(len(files)),
				mpb.BarStyle("[=>-]<+"),
				mpb.PrependDecorators(
					decor.Name("processing file: ", decor.WC{W: len("processing file: "), C: decor.DidentRight}),
					decor.Name("", decor.WCSyncSpaceR),
					decor.CountersNoUnit("%d / %d", decor.WCSyncWidth),
				),
				mpb.AppendDecorators(
					decor.Name("ETA: ", decor.WC{W: len("ETA: ")}),
					decor.EwmaETA(decor.ET_STYLE_GO, 60),
					decor.OnComplete(decor.Name(""), ". done"),
				),
			)
		}

		for _, file := range files {
			startTime = time.Now()

			fastxReader, err = fastx.NewDefaultReader(file)
			checkError(errors.Wrap(err, file))

			for {
				record, err = fastxReader.Read()
				if err != nil {
					if err == io.EOF {
						break
					}
					checkError(errors.Wrap(err, file))
					break
				}

				if syncmer {
					sketch, err = unikmer.NewSyncmerSketch(record.Seq, k, syncmerS, circular)
				} else if minimizer {
					sketch, err = unikmer.NewMinimizerSketch(record.Seq, k, minimizerW, circular)
				} else {
					iter, err = unikmer.NewHashIterator(record.Seq, k, true, circular)
				}
				if err != nil {
					if err == unikmer.ErrShortSeq {
						continue
					} else {
						checkError(errors.Wrapf(err, "seq: %s", record.Name))
					}
				}

				if syncmer {
					for {
						code, ok = sketch.NextSyncmer()
						if !ok {
							break
						}
						if scaled && code > maxHash {
							continue
						}
						codes = append(codes, code)
					}
					continue
				}

				if minimizer {
					for {
						code, ok = sketch.NextMinimizer()
						if !ok {
							break
						}
						if scaled && code > maxHash {
							continue
						}
						codes = append(codes, code)
					}
					continue
				}

				for {
					code, ok = iter.NextHash()
					if !ok {
						break
					}
					if scaled && code > maxHash {
						continue
					}
					codes = append(codes, code)
				}

				if opt.Verbose {
					bar.Increment()
					bar.DecoratorEwmaUpdate(time.Since(startTime))
				}
			}
		}

		if opt.Verbose {
			pbs.Wait()
		}

		n = len(codes)

		if exactNumber {
			// compute exact number of unique k-mers
			if opt.Verbose {
				log.Infof("sorting and counting unique k-mers ...")
			}
			codes2 := make([]uint64, n)
			copy(codes2, codes)
			sortutil.Uint64s(codes2)
			var pre uint64 = 0
			n = 0
			for _, code = range codes2 {
				if code != pre {
					n++
					pre = code
				}
			}
		}

		writer.Number = int64(n)

		if opt.Verbose {
			log.Infof("starting writing to %s ...", outFile)
		}
		for _, code = range codes {
			writer.WriteCode(code)
		}

		checkError(writer.Flush())
		if opt.Verbose {
			log.Infof("%d unique k-mers saved to %s", n, outFile)
		}
	},
}

func init() {
	RootCmd.AddCommand(computeCmd)

	computeCmd.Flags().StringP("out-dir", "O", "", `output directory, overide -o/--out-prefix`)
	computeCmd.Flags().StringP("out-prefix", "o", "", `out file prefix ("-" for stdout)`)
	computeCmd.Flags().BoolP("force", "", false, `overwrite output directory`)

	computeCmd.Flags().IntP("kmer-len", "k", 31, "k-mer length")
	computeCmd.Flags().BoolP("circular", "", false, "circular genome")

	computeCmd.Flags().IntP("scale", "D", 1, `scale/down-sample factor`)
	computeCmd.Flags().IntP("minimizer-w", "W", 0, `minimizer window size`)
	computeCmd.Flags().IntP("syncmer-s", "S", 0, `syncmer s`)

	computeCmd.Flags().BoolP("exact-number", "e", false, `save exact number of unique k-mers`)
}
