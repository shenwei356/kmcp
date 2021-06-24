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
	"encoding/json"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"regexp"
	"strings"
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
     positional arguments or a directory containing sequence files
     via the flag -I/--in-dir. A regular expression for matching
     sequencing files is available by the flag -r/--file-regexp.
  2. Multiple sizes of k-mers are supported.
     K-mers (sketchs) are not sorted, and duplicates are kept
     unless the flag -e/--exact-number is on.
  3. By default, we compute k-mers (sketches) of every file,
     you can also use --by-seq to compute for every sequence,
     where sequence IDs in all input files better be distinct.
  4. Unwanted sequence like plasmid can be filtered out by
     regular expressions via -B/--seq-name-filter.
  5. It also supports splitting sequences into fragments, this
     could increase the specificity in profiling result in cost
     of searching speed.

Supported k-mer (sketches) types:
  1. K-mer:
     1). ntHash of k-mer (-k)
  2. K-mer sketchs (all using ntHash):
     1). Scaled MinHash (-k -D)
     2). Minimizer      (-k -W), optionally scaling/down-sampling (-D)
     3). Closed Syncmer (-k -S), optionally scaling/down-sampling (-D)

Splitting sequences:
  1. Sequences can be splitted into fragments by a fragment size 
     (-s/--split-size) or number of fragments (-n/--split-number)
     with overlap (-l/--split-overlap).
  2. When splitting by number of fragments, all sequences (except for
     these mathching any regular expression given by -B/--seq-name-filter)
     in a sequence file are concatenated before splitting.
  3. Both sequence IDs and fragments indices are saved for later use,
     in form of meta/description data in .unik files.

Meta data:
  1. Every outputted .unik file contains the sequence ID/reference name,
     fragment index, number of fragments, and genome size of reference.
  2. When parsing whole sequence files or splitting by number of fragments,
     the identifier of a reference is the basename of the input file
     by default. It can also be extracted from the input file name via
     -N/--ref-name-regexp, e.g., "^(\w{3}_\d{9}\.\d+)" for refseq records.

Output:
  1. All outputted .unik files are saved in ${outdir}, with path
     ${outdir}/xxx/yyy/zzz/${infile}-id_${seqID}.unik
     where dirctory tree '/xxx/yyy/zzz/' is built for > 1000 output files.
  2. For splitting sequence mode (--split-size > 0 or --split-number > 0),
     output files are:
     ${outdir}//xxx/yyy/zzz/${infile}/{seqID}-frag_${fragIdx}.unik

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		seq.ValidateSeq = false

		var fhLog *os.File
		if opt.LogFile != "" {
			fhLog = addLog(opt.LogFile)
		}
		timeStart := time.Now()
		defer func() {
			if opt.Verbose {
				log.Info()
				log.Infof("elapsed time: %s", time.Since(timeStart))
				log.Info()
			}
			if opt.LogFile != "" {
				fhLog.Close()
			}
		}()

		// ---------------------------------------------------------------
		// basic flags

		ks := getFlagIntSlice(cmd, "kmer")
		if len(ks) == 0 {
			checkError(fmt.Errorf("flag -k/--kmer needed"))
		}
		for _, k := range ks {
			if k < 1 {
				checkError(fmt.Errorf("invalid k: %d", k))
			}
		}
		sortutil.Ints(ks)
		kMax := ks[len(ks)-1]

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

		if filepath.Clean(inDir) == filepath.Clean(outDir) {
			checkError(fmt.Errorf("intput and output paths should not be the same"))
		}

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
			checkError(errors.Wrapf(err, "failed to parse regular expression for matching file: %s", reFileStr))
		}

		reRefNameStr := getFlagString(cmd, "ref-name-regexp")
		var reRefName *regexp.Regexp
		var extractRefName bool
		if reRefNameStr != "" {
			if !regexp.MustCompile(`\(.+\)`).MatchString(reRefNameStr) {
				checkError(fmt.Errorf(`value of --ref-name-regexp must contains "(" and ")" to capture the ref name from file name`))
			}
			if !reIgnoreCase.MatchString(reRefNameStr) {
				reRefNameStr = reIgnoreCaseStr + reRefNameStr
			}

			reRefName, err = regexp.Compile(reRefNameStr)
			if err != nil {
				checkError(errors.Wrapf(err, "failed to parse regular expression for matching sequence header: %s", reRefName))
			}
			extractRefName = true
		}

		reSeqNameStrs := getFlagStringSlice(cmd, "seq-name-filter")
		reSeqNames := make([]*regexp.Regexp, 0, len(reSeqNameStrs))
		for _, kw := range reSeqNameStrs {
			if !reIgnoreCase.MatchString(reRefNameStr) {
				kw = reIgnoreCaseStr + kw
			}
			re, err := regexp.Compile("(?i)" + kw)
			if err != nil {
				checkError(errors.Wrapf(err, "failed to parse regular expression for matching sequence header: %s", kw))
			}
			reSeqNames = append(reSeqNames, re)
		}
		filterNames := len(reSeqNames) > 0

		// ---------------------------------------------------------------
		// flags for splitting sequence

		splitNumber0 := getFlagNonNegativeInt(cmd, "split-number")
		splitSize0 := getFlagNonNegativeInt(cmd, "split-size")
		splitOverlap := getFlagNonNegativeInt(cmd, "split-overlap")
		splitMinRef := getFlagNonNegativeInt(cmd, "split-min-ref")
		if splitNumber0 == 0 {
			splitNumber0 = 1
		}
		if splitSize0 > 0 && splitNumber0 > 1 {
			checkError(fmt.Errorf("flag -s/--split-size and -n/--split-number are incompatible"))
		}
		splitSeq := splitSize0 > 0 || splitNumber0 > 1
		if splitSeq {
			if splitSize0 > 0 { // split by size
				if splitSize0 < kMax {
					checkError(fmt.Errorf("value of flag -s/--split-size should >= k"))
				}
				if splitSize0 <= splitOverlap {
					checkError(fmt.Errorf("value of flag -s/--split-size should > value of -l/--split-overlap"))
				}
			} else { // split by number
				if splitNumber0 > 65535 {
					checkError(fmt.Errorf(("value of flag -s/--split-number should not be greater than 65535")))
				}
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
				log.Warningf("  no files matching regular expression: %s", reFileStr)
			}
		} else {
			files = getFileListFromArgsAndFile(cmd, args, true, "infile-list", true)
			if opt.Verbose {
				if len(files) == 1 && isStdin(files[0]) {
					log.Info("  no files given, reading from stdin")
				}
			}
		}
		if len(files) < 1 {
			checkError(fmt.Errorf("FASTA/Q files needed"))
		} else if opt.Verbose {
			log.Infof("  %d input file(s) given", len(files))
		}

		// ---------------------------------------------------------------
		// log

		if opt.Verbose {
			log.Info()
			log.Infof("-------------------- [main parameters] --------------------")

			log.Info("input and output:")
			log.Infof("  input directory: %s", inDir)
			log.Infof("    regular expression of input files: %s", reFileStr)
			log.Infof("    *regular expression for extracting reference name from file name: %s", reRefNameStr)
			log.Infof("    *regular expressions for filtering out sequences: %s", reSeqNameStrs)
			log.Infof("  output directory: %s", outDir)
			log.Info()

			log.Infof("sequences splitting: %v", splitSeq)
			if splitSeq {
				if splitNumber0 > 1 {
					log.Infof("  split parts: %d, overlap: %d bp", splitNumber0, splitOverlap)
				} else {
					log.Infof("  split sequence size: %d bp, overlap: %d bp", splitSize0, splitOverlap)
				}
			}
			log.Info()

			log.Info("k-mer (sketches) computing:")

			if !splitSeq && bySeq {
				log.Infof("  computing k-mers (sketches) for every sequence: %v", bySeq)
			}

			log.Infof("  k-mer size(s): %s", strings.Join(IntSlice2StringSlice(ks), ", "))

			log.Infof("  circular genome: %v", circular0)
			if minimizer {
				log.Infof("  minimizer window: %d", minimizerW)
			}
			if syncmer {
				log.Infof("  bounded syncmer size: %d", syncmerS)
			}
			if scaled {
				log.Infof("  down-sampling scale: %d", scale)
			}
			log.Infof("  saving exact number of k-mers: %v", exactNumber)
			log.Info()

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
			outfh.WriteString("#path\tname\tfragIdx\tidxNum\tgenomeSize\tkmers\n")
			for info := range ch {
				outfh.WriteString(fmt.Sprintf("%s\t%s\t%d\t%d\t%d\t%d\n", info.Path, info.Name, info.Index, info.Indexes, info.GenomeSize, info.Kmers))
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
					decor.Name("processed files: ", decor.WC{W: len("processed files: "), C: decor.DidentRight}),
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

				var k int
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
				var seqLen int
				var splitSize int
				var splitNumber int
				splitByNumber := splitNumber0 > 1

				var genomeSize uint64

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

				var allSeqs [][]byte
				var bigSeq []byte
				var record1 *fastx.Record

				var ignoreSeq bool
				var re *regexp.Regexp
				if splitByNumber { // concatenate all seqs
					allSeqs = make([][]byte, 0, 8)
					lenSum := 0
					first := true
					for {
						record, err = fastxReader.Read()
						if err != nil {
							if err == io.EOF {
								break
							}
							checkError(errors.Wrap(err, file))
							break
						}

						if filterNames {
							ignoreSeq = false
							for _, re = range reSeqNames {
								if re.Match(record.Name) {
									ignoreSeq = true
									break
								}
							}
							if ignoreSeq {
								continue
							}
						}

						if first {
							record1 = record
							first = false
						}
						aseq := make([]byte, len(record.Seq.Seq))
						copy(aseq, record.Seq.Seq)
						allSeqs = append(allSeqs, aseq)
						lenSum += len(aseq)
					}
					if len(allSeqs) == 1 {
						bigSeq = allSeqs[0]
					} else {
						bigSeq = make([]byte, lenSum)
						i := 0
						for _, aseq := range allSeqs {
							copy(bigSeq[i:i+len(aseq)], aseq)
							i += len(aseq)
						}
					}
					if lenSum == 0 {
						log.Warningf("skipping %s: no invalid sequences", file)
						log.Info()
						return
					}
					record1.Seq.Seq = bigSeq
					record = record1

					genomeSize = uint64(len(record.Seq.Seq))
				}

				slidIdx = 0

				once := true
				for {
					if splitByNumber {
						if !once {
							break
						}
						once = false
					} else {
						record, err = fastxReader.Read()
						if err != nil {
							if err == io.EOF {
								break
							}
							checkError(errors.Wrap(err, file))
							break
						}

						if filterNames {
							ignoreSeq = false
							for _, re = range reSeqNames {
								if re.Match(record.Name) {
									ignoreSeq = true
									break
								}
							}
							if ignoreSeq {
								continue
							}
						}

						genomeSize += uint64(len(record.Seq.Seq))
					}

					seqID = string(record.ID)
					seqLen = len(record.Seq.Seq)

					splitSize = splitSize0
					splitNumber = splitNumber0
					if !splitSeq || seqLen < splitMinRef {
						splitSize = seqLen
						step = 1
						greedy = false
					} else if splitByNumber {
						if splitNumber == 1 {
							splitSize = seqLen
							step = 1
							greedy = false
						} else {
							splitSize = (seqLen + (splitNumber-1)*splitOverlap + splitNumber - 1) / splitNumber
							step = splitSize - splitOverlap
						}
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

						if len(_seq.Seq)-1 <= splitOverlap {
							continue
						}

						for _, k = range ks {
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
							if splitByNumber {
								if extractRefName {
									if reRefName.MatchString(baseFile) {
										seqID = reRefName.FindAllStringSubmatch(baseFile, 1)[0][1]
									} else {
										seqID, _ = filepathTrimExtension(baseFile)
									}
								} else {
									seqID, _ = filepathTrimExtension(baseFile)
								}
							}
							outFileBase = fmt.Sprintf("%s/%s-frag_%d%s", baseFile, seqID, slidIdx, extDataFile)
						} else {
							outFileBase = fmt.Sprintf("%s-id_%s%s", baseFile, seqID, extDataFile)
						}
						// write to file
						if multiLevelFileTree {
							fileHash = xxh3.HashString(baseFile)
							dir3 = fmt.Sprintf("%03d", fileHash&1023)
						}

						outFile = filepath.Join(outDir, dir1, dir2, dir3, outFileBase)

						meta := Meta{
							SeqID:      seqID,
							FragIdx:    slidIdx,
							GenomeSize: uint64(len(record.Seq.Seq)),

							Ks: ks,

							Syncmer:      syncmer,
							SyncmerS:     syncmerS,
							Minimizer:    minimizer,
							MinimizerW:   minimizerW,
							SplitSeq:     splitSeq,
							SplitNum:     splitNumber0,
							SplitSize:    splitSize0,
							SplitOverlap: splitOverlap,
						}
						writeKmers(kMax, codes, uint64(n), outFile, compress, opt.CompressionLevel,
							scaled, scale, meta)

						ch <- UnikFileInfo{
							Path:       outFile,
							Name:       seqID,
							GenomeSize: uint64(len(record.Seq.Seq)),
							Index:      slidIdx,
							Indexes:    uint32(splitNumber0),
							Kmers:      uint64(n),
						}

						slidIdx++
					}

				}

				if bySeq {
					return
				}

				n = len(codes)

				if n == 0 {
					log.Warningf("skipping %s: no invalid sequences", file)
					log.Info()
					return
				}

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

				var fileprefix string
				if extractRefName {
					if reRefName.MatchString(baseFile) {
						fileprefix = reRefName.FindAllStringSubmatch(baseFile, 1)[0][1]
					} else {
						fileprefix, _ = filepathTrimExtension(baseFile)
					}
				} else {
					fileprefix, _ = filepathTrimExtension(baseFile)
				}
				meta := Meta{
					SeqID:      fileprefix,
					FragIdx:    slidIdx,
					GenomeSize: genomeSize,

					Ks: ks,

					Syncmer:      syncmer,
					SyncmerS:     syncmerS,
					Minimizer:    minimizer,
					MinimizerW:   minimizerW,
					SplitSeq:     splitSeq,
					SplitNum:     splitNumber0,
					SplitSize:    splitSize0,
					SplitOverlap: splitOverlap,
				}
				writeKmers(kMax, codes, uint64(n), outFile, compress, opt.CompressionLevel,
					scaled, scale, meta)

				ch <- UnikFileInfo{
					Path:       outFile,
					Name:       fileprefix,
					GenomeSize: genomeSize,
					Index:      slidIdx,
					Indexes:    uint32(splitNumber0),
					Kmers:      uint64(n),
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

	computeCmd.Flags().IntSliceP("kmer", "k", []int{21}, `k-mer size(s)`)
	computeCmd.Flags().BoolP("circular", "", false, `input sequence is circular`)

	computeCmd.Flags().IntP("scale", "D", 1, `scale/down-sample factor`)
	computeCmd.Flags().IntP("minimizer-w", "W", 0, `minimizer window size`)
	computeCmd.Flags().IntP("syncmer-s", "S", 0, `bounded syncmer length`)

	computeCmd.Flags().BoolP("exact-number", "e", false, `save exact number of unique k-mers for indexing`)
	computeCmd.Flags().BoolP("compress", "c", false, `output gzipped .unik files, it's slower and can saves little space`)

	computeCmd.Flags().IntP("split-number", "n", 0, `fragment number, incompatible with -s/--split-size`)
	computeCmd.Flags().IntP("split-size", "s", 0, `fragment size for splitting sequences, incompatible with -n/--split-number`)
	computeCmd.Flags().IntP("split-overlap", "l", 0, `fragment overlap for splitting sequences`)
	computeCmd.Flags().IntP("split-min-ref", "m", 1000, `only splitting sequences >= M bp`)

	computeCmd.Flags().BoolP("by-seq", "", false, `compute k-mers (sketches) for every sequence, instead of whole file`)

	computeCmd.Flags().StringP("ref-name-regexp", "N", "", `regular expression (must contains "(" and ")") for extracting reference name from file name`)
	computeCmd.Flags().StringSliceP("seq-name-filter", "B", []string{}, `list of regular expressions for filtering out sequences by header/name, case ignored`)
}

var reIgnoreCaseStr = "(?i)"
var reIgnoreCase = regexp.MustCompile(`\(\?i\)`)

var fileUnikInfos = "_info.txt"
