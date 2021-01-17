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
	"encoding/json"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"regexp"
	"runtime"
	"sync"
	"time"

	"github.com/pkg/errors"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/unikmer"
	"github.com/shenwei356/util/pathutil"
	"github.com/spf13/cobra"
	"github.com/twotwotwo/sorts/sortutil"
	"github.com/vbauerster/mpb/v5"
	"github.com/vbauerster/mpb/v5/decor"
	"github.com/zeebo/xxh3"
)

var computeCmd = &cobra.Command{
	Use:   "compute",
	Short: "Generate k-mers (sketch) from FASTA/Q sequences",
	Long: `Generate k-mers (sketchs) from FASTA/Q sequences

Attentions:
  1. Input files can be given as list of FASTA/Q files via
     positional arguments or a directory containing files
     via flag -I/--in-dir. Regular expression for matching
     sequencing file is given by -r/--file-regexp.
  2. K-mers (sketchs) are not sorted and duplicates remained.
  3. Sequence IDs in all input files should be unique.
     Multiple sequences belonging to a same group can be mapped
     to their group name via name mapping file after indexing.
  4. By default, we compute k-mers (sketches) of every file,
     you can also use --by-seq to compute for every file.
  5. It also supports splitting sequences into fragment via
     flag -s/--split-size and -l/--split-overlap.

Supported k-mer (sketches) types:
  1. K-mer:
     1. ntHash of k-mer (-k)
  2. K-mer sketchs (all using ntHash):
     1. Scaled MinHash (-k -D)
     2. Minimizer      (-k -W), optionally scaling/down-sampling (-D)
     3. Syncmer        (-k -S), optionally scaling/down-sampling (-D)

Splitting sequences:
  1. Sequences can be splitted into fragments (-s/--split-size)
     with overlap (-l/--split-overlap).
  2. Both sequence IDs and fragments indices are saved for later use,
     in form of meta/description data in .unik files.

Output:
  1. All outputted .unik files are saved in ${outdir}, with path
     ${outdir}/xxx/yyy/zzz/${infile}-id${seqID}.unik
     where dirctory tree '/xxx/yyy/zzz/' is built for > 1000 output files.
  2. For splitting sequence mode (--split-size > 0), output files are
     ${outdir}//xxx/yyy/zzz/${infile}/{seqID}-frag${fragIdx}.unik

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)
		seq.ValidateSeq = false

		timeStart := time.Now()
		defer func() {
			if opt.Verbose {
				log.Info()
				log.Infof("elapsed time: %s", time.Since(timeStart))
			}
		}()

		// ---------------------------------------------------------------
		// basic flags

		k := getFlagPositiveInt(cmd, "kmer-len")
		circular0 := getFlagBool(cmd, "circular")

		outDir := getFlagString(cmd, "out-dir")
		force := getFlagBool(cmd, "force")
		exactNumber := getFlagBool(cmd, "exact-number")
		compress := getFlagBool(cmd, "compress")
		bySeq := getFlagBool(cmd, "by-seq")

		if outDir == "" {
			checkError(fmt.Errorf("flag -O/--out-dir is needed"))
		}

		var err error

		inDir := getFlagString(cmd, "in-dir")

		readFromDir := inDir != ""
		if readFromDir {
			var isDir bool
			isDir, err = pathutil.IsDir(inDir)
			if err != nil {
				checkError(errors.Wrapf(err, "checking -I/--in-dir"))
			}
			if !isDir {
				checkError(fmt.Errorf("value of -I/--in-dir should be a directory: %s", inDir))
			}
		}

		reFileStr := getFlagString(cmd, "file-regexp")
		var reFile *regexp.Regexp
		if reFileStr != "" {
			if !reIgnoreCase.MatchString(reFileStr) {
				reFileStr = reIgnoreCaseStr + reFileStr
			}
			reFile, err = regexp.Compile(reFileStr)
			checkError(errors.Wrapf(err, "parsing regular expression for matching file: %s", reFileStr))
		}

		// ---------------------------------------------------------------
		// flags for splitting sequence

		splitSize0 := getFlagNonNegativeInt(cmd, "split-size")
		splitOverlap := getFlagNonNegativeInt(cmd, "split-overlap")
		splitSeq := splitSize0 > 0
		if splitSeq {
			if splitSize0 < k {
				checkError(fmt.Errorf("value of flag -s/--split-size should >= k"))
			}
			if splitSize0 <= splitOverlap {
				checkError(fmt.Errorf("value of flag -s/--split-size should > value of -l/--split-overlap"))
			}
			bySeq = true
		}
		step := splitSize0 - splitOverlap

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

		var files []string
		if readFromDir {
			files, err = getFileListFromDir(inDir, reFile, opt.NumCPUs)
			if err != nil {
				checkError(errors.Wrapf(err, "walking dir: %s", inDir))
			}
			if len(files) == 0 {
				log.Warningf("no files matching regular expression: %s", reFileStr)
			}
		} else {
			files = getFileListFromArgsAndFile(cmd, args, true, "infile-list", true)
			if opt.Verbose {
				if len(files) == 1 && isStdin(files[0]) {
					log.Info("no files given, reading from stdin")
				}
			}
		}
		if len(files) < 1 {
			checkError(fmt.Errorf("FASTA/Q files needed"))
		} else if opt.Verbose {
			log.Infof("%d input file(s) given", len(files))
		}

		// ---------------------------------------------------------------
		// log

		if opt.Verbose {
			log.Info()
			log.Infof("-------------------- [main parameters] --------------------")
			log.Infof("k: %d", k)
			log.Infof("circular genome: %v", circular0)
			if splitSeq {
				log.Infof("split seqequence size: %d, overlap: %d", splitSize0, splitOverlap)
			}
			if minimizer {
				log.Infof("minimizer window: %d", minimizerW)
			}
			if syncmer {
				log.Infof("bounded syncmer size: %d", syncmerS)
			}
			if scaled {
				log.Infof("down-sampling scale: %d", scale)
			}
			log.Infof("-------------------- [main parameters] --------------------")
			log.Info()
			log.Infof("computing ...")
		}

		// ---------------------------------------------------------------

		// file info
		outfh, gw, w, err := outStream(filepath.Join(outDir, fileUnikInfos), false, -1)
		checkError(err)
		defer func() {
			outfh.Flush()
			if gw != nil {
				gw.Close()
			}
			w.Close()
		}()

		ch := make(chan UnikFileInfo, opt.NumCPUs)
		done := make(chan int)
		go func() {
			outfh.WriteString(fmt.Sprintf("#path\tname\tfragIdx\tkmers\n"))
			for info := range ch {
				outfh.WriteString(fmt.Sprintf("%s\t%s\t%d\t%d\n", info.Path, info.Name, info.Index, info.Kmers))
			}
			done <- 1
		}()

		// process bar
		var pbs *mpb.Progress
		var bar *mpb.Bar
		var chDuration chan time.Duration
		var doneDuration chan int

		if opt.Verbose {
			pbs = mpb.New(mpb.WithWidth(40), mpb.WithOutput(os.Stderr))
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
		threadsFloat := float64(opt.NumCPUs)

		multiLevelFileTree := bySeq || len(files) > 1000

		for _, file := range files {
			tokens <- 1
			wg.Add(1)

			go func(file string, multiLevelFileTree bool) {
				startTime := time.Now()
				defer func() {
					wg.Done()
					<-tokens

					if opt.Verbose {
						chDuration <- time.Duration(float64(time.Since(startTime)) / threadsFloat)
					}
				}()

				var err error
				var record *fastx.Record
				var fastxReader *fastx.Reader
				var ok bool
				var code uint64
				var iter *unikmer.Iterator
				var sketch *unikmer.Sketch
				var n int

				var slider func() (*seq.Seq, bool)
				var _seq *seq.Seq
				var _ok bool
				var slidIdx uint32
				var greedy bool = true
				var seqID string
				var outFile string
				var baseFile = filepath.Base(file)
				var circular bool
				var splitSize int

				var codes []uint64

				// multiple level file tree for saving files
				var dir1, dir2 string
				var fileHash uint64
				if multiLevelFileTree {
					fileHash = xxh3.HashString(baseFile)
					dir1 = fmt.Sprintf("%03d", fileHash&1023)
					dir2 = fmt.Sprintf("%03d", (fileHash>>10)&1023)
				}
				var outFileBase, dir3 string

				if !splitSeq {
					codes = make([]uint64, 0, mapInitSize)
					circular = circular0
				} else {
					codes = make([]uint64, 0, splitSize0)
					circular = false // split seq is linear, isn't it?
				}

				fastxReader, err = fastx.NewDefaultReader(file)
				checkError(errors.Wrap(err, file))

				slidIdx = 0

				for {
					record, err = fastxReader.Read()
					if err != nil {
						if err == io.EOF {
							break
						}
						checkError(errors.Wrap(err, file))
						break
					}
					seqID = string(record.ID)

					splitSize = splitSize0
					if !splitSeq { // whole sequence
						splitSize = len(record.Seq.Seq)
						step = 1
						greedy = false
					}

					slider = record.Seq.Slider(splitSize, step, circular0, greedy)

					if bySeq {
						slidIdx = 0
					}

					for {
						_seq, _ok = slider()
						if !_ok {
							break
						}

						if bySeq {
							codes = codes[:0] // reset
						}

						if syncmer {
							sketch, err = unikmer.NewSyncmerSketch(_seq, k, syncmerS, circular)
						} else if minimizer {
							sketch, err = unikmer.NewMinimizerSketch(_seq, k, minimizerW, circular)
						} else {
							iter, err = unikmer.NewHashIterator(_seq, k, true, circular)
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
						} else if minimizer {
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
						} else {
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

						if !bySeq {
							break
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

						if splitSeq {
							outFileBase = fmt.Sprintf("%s/%s-frag%d%s", baseFile, seqID, slidIdx, extDataFile)
						} else {
							outFileBase = fmt.Sprintf("%s-id%s%s", baseFile, seqID, extDataFile)
						}
						// write to file
						if multiLevelFileTree {
							fileHash = xxh3.HashString(baseFile)
							dir3 = fmt.Sprintf("%03d", fileHash&1023)
						}
						outFile = filepath.Join(outDir, dir1, dir2, dir3, outFileBase)

						meta := Meta{
							SeqID:   seqID,
							FragIdx: slidIdx,

							Syncmer:      syncmer,
							SyncmerS:     syncmerS,
							Minimizer:    minimizer,
							MinimizerW:   minimizerW,
							SplitSeq:     splitSeq,
							SplitSize:    splitSize0,
							SplitOverlap: splitOverlap,
						}
						writeKmers(k, codes, uint64(n), outFile, compress, opt.CompressionLevel,
							scaled, scale, meta)

						ch <- UnikFileInfo{
							Path:  outFile,
							Name:  seqID,
							Index: slidIdx,
							Kmers: uint64(n),
						}

						slidIdx++
					}

				}

				if bySeq {
					return
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

				// write to file
				outFile = filepath.Join(outDir, dir1, dir2, fmt.Sprintf("%s%s", baseFile, extDataFile))

				fileprefix, _ := filepathTrimExtension(filepath.Base(file))
				meta := Meta{
					SeqID:   fileprefix,
					FragIdx: slidIdx,

					Syncmer:      syncmer,
					SyncmerS:     syncmerS,
					Minimizer:    minimizer,
					MinimizerW:   minimizerW,
					SplitSeq:     splitSeq,
					SplitSize:    splitSize0,
					SplitOverlap: splitOverlap,
				}
				writeKmers(k, codes, uint64(n), outFile, compress, opt.CompressionLevel,
					scaled, scale, meta)

				ch <- UnikFileInfo{
					Path:  outFile,
					Name:  fileprefix,
					Index: slidIdx,
					Kmers: uint64(n),
				}

				slidIdx++

			}(file, multiLevelFileTree)
		}

		wg.Wait()

		if opt.Verbose {
			close(chDuration)
			<-doneDuration
			pbs.Wait()
		}

		close(ch)
		<-done
	},
}

func writeKmers(k int, codes []uint64, n uint64,
	outFile string, compress bool, compressLevel int,
	scaled bool, scale int, meta Meta) {

	outfh, gw, w, err := outStream(outFile, compress, compressLevel)
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

	mode |= unikmer.UnikCanonical
	mode |= unikmer.UnikHashed

	writer, err = unikmer.NewWriter(outfh, k, mode)
	if err != nil {
		checkError(errors.Wrap(err, outFile))
	}

	if scaled {
		writer.SetScale(uint32(scale))
	}

	writer.Number = n

	metatext, err := json.Marshal(meta)
	if err != nil {
		checkError(errors.Wrap(err, outFile))
	}
	writer.Description = metatext

	for _, code := range codes {
		writer.WriteCode(code)
	}

	checkError(writer.Flush())
}

func init() {
	RootCmd.AddCommand(computeCmd)

	computeCmd.Flags().StringP("in-dir", "I", "", `directory containing FASTA/Q files. directory symlinks are followed`)
	computeCmd.Flags().StringP("file-regexp", "r", `\.(f[aq](st[aq])?|fna)(.gz)?$`, `regular expression for matching files in -I/--in-dir to compute, case ignored`)

	computeCmd.Flags().StringP("out-dir", "O", "", `output directory`)
	computeCmd.Flags().BoolP("force", "", false, `overwrite output directory`)

	computeCmd.Flags().IntP("kmer-len", "k", 31, `k-mer length`)
	computeCmd.Flags().BoolP("circular", "", false, `input sequence is curcular`)

	computeCmd.Flags().IntP("scale", "D", 1, `scale/down-sample factor`)
	computeCmd.Flags().IntP("minimizer-w", "W", 0, `minimizer window size`)
	computeCmd.Flags().IntP("syncmer-s", "S", 0, `bounded syncmer length`)

	computeCmd.Flags().BoolP("exact-number", "e", false, `save exact number of unique k-mers for indexing`)
	computeCmd.Flags().BoolP("compress", "c", false, `output gzipped .unik files (can save little space`)

	computeCmd.Flags().IntP("split-size", "s", 0, `fragment size for splitting sequences`)
	computeCmd.Flags().IntP("split-overlap", "l", 0, `fragment overlap for splitting sequences`)

	computeCmd.Flags().BoolP("by-seq", "", false, `compute k-mers (sketches) for every sequence, instead of whole file`)

}

var reIgnoreCaseStr = "(?i)"
var reIgnoreCase = regexp.MustCompile(`\(\?i\)`)

var fileUnikInfos = "_info.txt"
