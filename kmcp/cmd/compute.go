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
	"path/filepath"
	"runtime"
	"sync"
	"time"

	"github.com/pkg/errors"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/unikmer"
	"github.com/smallnest/ringbuffer"
	"github.com/spf13/cobra"
	"github.com/vbauerster/mpb"
	"github.com/vbauerster/mpb/decor"
)

var computeCmd = &cobra.Command{
	Use:   "compute",
	Short: "Generate k-mers (sketch) from FASTA/Q sequences",
	Long: `Generate k-mers (sketch) from FASTA/Q sequences

K-mer:
  1. ntHash of k-mer (-k)

K-mer sketchs:
  1. Scaled MinHash (-D)
  2. Minimizer (-k -W), optionally scaling/down-sampling (-D)
  3. Syncmer   (-k -S), optionally scaling/down-sampling (-D)

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

		if outDir == "" && outFile == "" {
			checkError(fmt.Errorf("flag -o/--out-prefix OR -O/--out-dir is needed"))
		} else if outDir != "" && outFile != "-" {
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
			pbs := mpb.New(mpb.WithWidth(79))
			bar := pbs.AddBar(int64(len(files)),
				mpb.BarStyle("[=>-]<+"),
				mpb.PrependDecorators(
					decor.Name("processing file: ", decor.WC{W: len("compute") + 1, C: decor.DidentRight}),
					decor.Name("", decor.WCSyncSpaceR),
					decor.CountersNoUnit("%d / %d", decor.WCSyncWidth),
				),
				mpb.AppendDecorators(
					decor.EwmaETA(decor.ET_STYLE_GO, 60),
				),
			)

			chDuration := make(chan time.Duration, opt.NumCPUs)
			done := make(chan int)
			go func() {
				for t := range chDuration {
					bar.Increment()
					bar.DecoratorEwmaUpdate(t)
				}
				done <- 1
			}()

			// wait group
			var wg sync.WaitGroup
			tokens := ringbuffer.New(opt.NumCPUs)

			for _, file := range files {
				tokens.WriteByte(0)
				wg.Add(1)

				go func(file string) {
					startTime := time.Now()

					defer func() {
						wg.Done()
						tokens.ReadByte()

						chDuration <- time.Since(startTime)
					}()

					outFile := filepath.Join(outDir, fmt.Sprintf("%s%s", filepath.Base(file), extDataFile))
					outfh, gw, w, err := outStream(outFile, opt.Compress, opt.CompressionLevel)
					checkError(err)
					defer func() {
						outfh.Flush()
						if gw != nil {
							gw.Close()
						}
						w.Close()
					}()

					var writer *unikmer.Writer
					var mode uint32
					var n int

					mode |= unikmer.UnikCompact
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
					writer.Number = int64(n)
					for _, code = range codes {
						writer.WriteCode(code)
					}

					checkError(writer.Flush())
				}(file)
			}

			wg.Wait()

			close(chDuration)
			<-done
			pbs.Wait()

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
			outfh.Flush()
			if gw != nil {
				gw.Close()
			}
			w.Close()
		}()

		var writer *unikmer.Writer
		var mode uint32
		var n int

		mode |= unikmer.UnikCompact
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

		for _, file := range files {
			if opt.Verbose {
				log.Infof("reading sequence file: %s", file)
			}
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
		}

		n = len(codes)
		writer.Number = int64(n)
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
	computeCmd.Flags().StringP("out-prefix", "o", "-", `out file prefix ("-" for stdout)`)
	computeCmd.Flags().BoolP("force", "", false, `overwrite output directory`)

	computeCmd.Flags().IntP("kmer-len", "k", 31, "k-mer length")
	computeCmd.Flags().BoolP("circular", "", false, "circular genome")

	computeCmd.Flags().IntP("scale", "D", 1, `scale/down-sample factor`)
	computeCmd.Flags().IntP("minimizer-w", "W", 0, `minimizer window size`)
	computeCmd.Flags().IntP("syncmer-s", "S", 0, `syncmer s`)
}
