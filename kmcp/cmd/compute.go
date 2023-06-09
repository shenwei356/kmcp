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
	"bytes"
	"encoding/json"
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"regexp"
	"strings"
	"sync"
	"time"

	"github.com/pkg/errors"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/bio/sketches"
	"github.com/shenwei356/unik/v5"
	"github.com/shenwei356/util/pathutil"
	"github.com/spf13/cobra"
	"github.com/twotwotwo/sorts/sortutil"
	"github.com/vbauerster/mpb/v5"
	"github.com/vbauerster/mpb/v5/decor"
	"github.com/zeebo/xxh3"
)

var computeCmd = &cobra.Command{
	Use:   "compute",
	Short: "Generate k-mers (sketches) from FASTA/Q sequences",
	Long: `Generate k-mers (sketches) from FASTA/Q sequences

Input:
  1. Input plain or gzipped FASTA/Q files can be given via positional
     arguments or the flag -i/--infile-list with the list of input files,
  2. Or a directory containing sequence files via the flag -I/--in-dir,
     with multiple-level sub-directories allowed. A regular expression
     for matching sequencing files is available via the flag -r/--file-regexp.
 *3. For taxonomic profiling, the sequences of each reference genome should be
     saved in a separate file, with the reference identifier in the file name.

  Attention:
    You may rename the sequence files for convenience because the 
  sequence/genome identifier in the index and search results would be:
    1). For the default mode (computing k-mers for the whole file):
          the basename of file with common FASTA/Q file extension removed,
          captured via the flag -N/--ref-name-regexp.  
    2). For splitting sequence mode (see details below):
          same to 1).
    3). For computing k-mers for each sequence:
          the sequence identifier.

Attentions:
  1. Unwanted sequences like plasmid can be filtered out by
     the name via regular expressions (-B/--seq-name-filter).
  2. By default, kmcp computes k-mers (sketches) of every file,
     you can also use --by-seq to compute for every sequence,
     where sequence IDs in all input files are better to be distinct.
  3. It also supports splitting sequences into chunks, this
     could increase the specificity in search results at the cost
     of a slower searching speed.
  4. Multiple sizes of k-mers are supported.

Supported k-mer (sketches) types:
  1. K-mer:
     1). ntHash of k-mer (-k)
  2. K-mer sketchs (all using ntHash):
     1). FracMinHash    (-k -D), previously named Scaled MinHash
     2). Minimizer      (-k -W), optionally scaling/down-sampling (-D)
     3). Closed Syncmer (-k -S), optionally scaling/down-sampling (-D)

Splitting sequences:
  1. Sequences can be splitted into chunks by a chunk size 
     (-s/--split-size) or number of chunks (-n/--split-number)
     with overlap (-l/--split-overlap).
     In this mode, the sequences of each genome should be saved in an
     individual file.
  2. When splitting by number of chunks, all sequences (except for
     these matching any regular expression given by -B/--seq-name-filter)
     in a sequence file are concatenated with k-1 N's before splitting.
  3. Both sequence IDs and chunks indices are saved for later use,
     in form of meta/description data in .unik files.

Metadata:
  1. Every outputted .unik file contains the sequence/reference ID,
     chunk index, number of chunks, and genome size of reference.
  2. When parsing whole sequence files or splitting by the number of chunks,
     the identifier of a reference is the basename of the input file
     by default. It can also be extracted from the input file name via
     -N/--ref-name-regexp, e.g., "^(\w{3}_\d{9}\.\d+)" for RefSeq records.

Output:
  1. All outputted .unik files are saved in ${outdir}, with path
     ${outdir}/xxx/yyy/zzz/${infile}-id_${seqID}.unik
     where dirctory tree '/xxx/yyy/zzz/' is built for > 1000 output files.
  2. For splitting sequence mode (--split-size > 0 or --split-number > 0),
     output files are:
     ${outdir}//xxx/yyy/zzz/${infile}/{seqID}-chunk_${chunkIdx}.unik
  3. A summary file ("${outdir}/_info.txt") is generated for later use.
     Users need to check if the reference IDs (column "name") are what
     supposed to be.

Performance tips:
  1. Decrease the value of -j/--threads for data in hard disk drives to
     reduce I/O pressure.

Next step:
  1. Check the summary file (${outdir}/_info.txt) to see if the reference
     IDs (column "name") are what supposed to be.
  2. Run "kmcp index" with the output directory.

Examples:
  1. From few sequence files:

        kmcp compute -k 21 -n 10 -l 150 -O tmp-k21-n10-l150 NC_045512.2.fna.gz

  2. From a list file:

        kmcp compute -k 21 -n 10 -l 150 -O tmp-k21-210-l150 -i list.txt

  3. From a directory containing many sequence files:

        kmcp compute -k 21 -n 10 -l 150 -B plasmid \
            -O gtdb-k21-n10-l150 -I gtdb-genomes/

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		seq.ValidateSeq = false

		var fhLog *os.File
		if opt.Log2File {
			fhLog = addLog(opt.LogFile, opt.Verbose)
		}
		timeStart := time.Now()
		defer func() {
			if opt.Verbose || opt.Log2File {
				log.Info()
				log.Infof("elapsed time: %s", time.Since(timeStart))
				log.Info()
			}
			if opt.Log2File {
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
			if k > 64 {
				checkError(fmt.Errorf("k-mer size (%d) should be <=64", k))
			}
		}
		sortutil.Ints(ks)
		kMax := ks[len(ks)-1]
		kMin := ks[0]

		circular0 := getFlagBool(cmd, "circular")

		outDir := getFlagString(cmd, "out-dir")
		force := getFlagBool(cmd, "force")
		// exactNumber := getFlagBool(cmd, "exact-number")
		exactNumber := true
		compress := getFlagBool(cmd, "compress")
		bySeq := getFlagBool(cmd, "by-seq")

		if outDir == "" {
			checkError(fmt.Errorf("flag -O/--out-dir is needed"))
		}

		var err error

		inDir := getFlagString(cmd, "in-dir")

		if filepath.Clean(inDir) == filepath.Clean(outDir) {
			checkError(fmt.Errorf("intput and output paths should not be the same: %s", outDir))
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
			if !reIgnoreCase.MatchString(kw) {
				kw = reIgnoreCaseStr + kw
			}
			re, err := regexp.Compile(kw)
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
		if !cmd.Flags().Lookup("split-overlap").Changed {
			splitOverlap = kMax - 1
		}
		splitMinRef := getFlagNonNegativeInt(cmd, "split-min-ref")
		if splitNumber0 == 0 {
			splitNumber0 = 1
		}
		if splitSize0 > 0 && splitNumber0 > 1 {
			checkError(fmt.Errorf("flag -s/--split-size and -n/--split-number are incompatible"))
		}
		if bySeq && splitNumber0 > 1 {
			checkError(fmt.Errorf("flag --by-seq and -n/--split-number are incompatible"))
		}
		if bySeq && splitSize0 > 0 {
			checkError(fmt.Errorf("flag --by-seq and -s/--split-size are incompatible"))
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

		var circular bool // for computting k-mers
		if !splitSeq {
			circular = circular0
		} else {
			circular = false // split seq is linear, isn't it?
		}

		// ---------------------------------------------------------------
		// flags of sketch

		scale := getFlagPositiveInt(cmd, "scale")
		if scale > math.MaxInt32 {
			checkError(fmt.Errorf("value of flag --scale is too big"))
		}
		scaled := scale > 1
		maxHash := uint64(float64(^uint64(0)) / float64(scale))

		minimizerW := getFlagNonNegativeInt(cmd, "minimizer-w")
		if minimizerW > math.MaxInt32 {
			checkError(fmt.Errorf("value of flag --minimizer-w is too big"))
		}
		minimizer := minimizerW > 0

		syncmerS := getFlagNonNegativeInt(cmd, "syncmer-s")
		if syncmerS > kMax {
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

		if opt.Verbose || opt.Log2File {
			log.Infof("kmcp v%s", VERSION)
			log.Info("  https://github.com/shenwei356/kmcp")
			log.Info()

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
			if opt.Verbose || opt.Log2File {
				if len(files) == 1 && isStdin(files[0]) {
					log.Info("  no files given, reading from stdin")
				}
			}
		}
		if len(files) < 1 {
			checkError(fmt.Errorf("FASTA/Q files needed"))
		} else if opt.Verbose || opt.Log2File {
			log.Infof("  %d input file(s) given", len(files))
		}

		// ---------------------------------------------------------------
		// log

		if opt.Verbose || opt.Log2File {
			log.Info()
			log.Infof("-------------------- [main parameters] --------------------")

			log.Info("input and output:")
			log.Infof("  input directory: %s", inDir)
			log.Infof("    regular expression of input files: %s", reFileStr)
			log.Infof("    *regular expression for extracting reference name from file name: %s", reRefNameStr)
			log.Infof("    *regular expressions for filtering out sequences: %s", reSeqNameStrs)
			log.Infof("  output directory: %s", outDir)
			log.Info()

			if !splitSeq && bySeq {
				log.Infof("compute k-mers (sketches) for each sequence: %v", bySeq)
			}

			log.Infof("sequences splitting: %v", splitSeq)
			if splitSeq {
				if splitNumber0 > 1 {
					log.Infof("  split number: %d, overlap: %d bp", splitNumber0, splitOverlap)
				} else {
					log.Infof("  split sequence size: %d bp, overlap: %d bp", splitSize0, splitOverlap)
				}
				if circular0 {
					log.Infof("  circular genome: %v (only applies to genomes with a single chromosome)", circular0) // circular2
				} else {
					log.Infof("  circular genome: %v", circular0)
				}
			}
			log.Info()

			log.Info("k-mer (sketches) computing:")

			if !splitSeq && bySeq {
				log.Infof("  computing k-mers (sketches) for every sequence: %v", bySeq)
			}

			log.Infof("  k-mer size(s): %s", strings.Join(IntSlice2StringSlice(ks), ", "))

			log.Infof("  circular genome: %v (only applies to non-split modes)", circular) // circular

			if minimizer {
				log.Infof("  minimizer window: %d", minimizerW)
			}
			if syncmer {
				log.Infof("  closed syncmer size: %d", syncmerS)
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
			outfh.WriteString("#path\tname\tchunkIdx\tidxNum\tgenomeSize\tkmers\n")
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
					decor.EwmaETA(decor.ET_STYLE_GO, 10),
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
		var wg sync.WaitGroup                 // ensure all jobs done
		tokens := make(chan int, opt.NumCPUs) // control the max concurrency number
		threadsFloat := float64(opt.NumCPUs)  // just avoid repeated type conversion

		multiLevelFileTree := bySeq || len(files) > 1000

		for _, file := range files {
			tokens <- 1
			wg.Add(1)

			go func(file string, multiLevelFileTree bool) {
				startTime := time.Now()
				defer func() {
					wg.Done()
					<-tokens

					if opt.Verbose || opt.Log2File {
						// update duration of a file
						chDuration <- time.Duration(float64(time.Since(startTime)) / threadsFloat)
					}
				}()

				var k int
				var err error
				var record *fastx.Record
				var fastxReader *fastx.Reader
				var ok bool
				var code uint64
				var iter *sketches.Iterator
				var sketch *sketches.Sketch
				var n int

				var slider func() (*seq.Seq, bool)
				var _seq *seq.Seq
				var _ok bool
				var slidIdx uint32
				var greedy bool = true
				var seqID string
				var outFile string
				var baseFile = filepath.Base(file)
				var seqLen int
				var splitSize int
				var splitNumber int
				var step int
				splitByNumber := splitNumber0 > 1

				var genomeSize uint64

				// var codes, codes2 []uint64
				codes := poolCodes.Get().(*[]uint64)
				codes2 := poolCodes.Get().(*[]uint64)
				*codes = (*codes)[:0]
				*codes2 = (*codes2)[:0]
				defer func() {
					*codes = (*codes)[:0]
					*codes2 = (*codes2)[:0]
					poolCodes.Put(codes)
					poolCodes.Put(codes2)
				}()

				// multiple level file tree for saving files
				var dir1, dir2 string
				var fileHash uint64
				if multiLevelFileTree {
					fileHash = xxh3.HashString(baseFile)
					dir1 = fmt.Sprintf("%03d", fileHash&1023)
					dir2 = fmt.Sprintf("%03d", (fileHash>>10)&1023)
				}
				var outFileBase, dir3 string

				fastxReader, err = fastx.NewDefaultReader(file)
				checkError(errors.Wrap(err, file))

				var allSeqs [][]byte
				var bigSeq []byte
				nnn := bytes.Repeat([]byte{'N'}, kMax-1)

				var ignoreSeq bool
				var re *regexp.Regexp
				if splitSeq { // concatenate all seqs. given split number of by split size
					allSeqs = make([][]byte, 0, 8)
					lenSum := 0
					for {
						record, err = fastxReader.Read()
						if err != nil {
							if err == io.EOF {
								break
							}
							checkError(errors.Wrap(err, file))
							break
						}

						// filter out sequences with names in the blast list
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

						aseq := make([]byte, len(record.Seq.Seq))
						copy(aseq, record.Seq.Seq)
						allSeqs = append(allSeqs, aseq)
						lenSum += len(aseq)
					}

					if lenSum == 0 {
						log.Warningf("skipping %s: no valid sequences", file)
						log.Info()
						return
					}

					if len(allSeqs) == 1 {
						bigSeq = allSeqs[0]
					} else {
						bigSeq = make([]byte, lenSum+(len(allSeqs)-1)*(kMax-1))
						i := 0
						for j, aseq := range allSeqs {
							copy(bigSeq[i:i+len(aseq)], aseq)
							if j < len(allSeqs)-1 {
								copy(bigSeq[i+len(aseq):i+len(aseq)+kMax-1], nnn)
							}
							i += len(aseq) + kMax - 1
						}
					}

					record, _ = fastx.NewRecordWithoutValidation(seq.Unlimit, nil, []byte{}, nil, bigSeq)

					genomeSize = uint64(len(record.Seq.Seq))
				}

				nSeqs := len(allSeqs)

				slidIdx = 0
				var _splitNumber int

				circular2 := circular0 // for computing sliding window and step

				once := true
				for {
					if splitSeq { // only one loop for split seq mode, the records has already loaded into record
						if !once {
							break
						}
						once = false
					} else { // regularly read one sequence
						record, err = fastxReader.Read()
						if err != nil {
							if err == io.EOF {
								break
							}
							checkError(errors.Wrap(err, file))
							break
						}

						// filter out sequences with names in the blast list
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
						step = seqLen
						greedy = false
						splitNumber = 1
					} else if splitSeq {
						if splitNumber > 1 {
							if circular0 && nSeqs == 1 {
								circular2 = true
								greedy = false
								splitSize = (seqLen + splitNumber*splitOverlap + splitNumber - 1) / splitNumber
								step = splitSize - splitOverlap
							} else {
								circular2 = false
								splitSize = (seqLen + (splitNumber-1)*splitOverlap + splitNumber - 1) / splitNumber
								step = splitSize - splitOverlap
							}
						} else if splitSize > 0 { // --split-size
							step = splitSize - splitOverlap
						} else { // splitNumber == 1 , no split
							splitSize = seqLen
							step = seqLen
							greedy = false
						}
					}

					// get the actual split number
					_splitNumber = 1
					if splitSeq {
						slider = record.Seq.Slider(splitSize, step, circular2, greedy)
						_splitNumber = 0
						for {
							_seq, _ok = slider()
							if !_ok {
								break
							}
							if len(_seq.Seq)-1 <= splitOverlap || len(_seq.Seq) < kMin {
								continue
							}

							_splitNumber++
						}

						if _splitNumber == 0 {
							log.Warningf("sequence is too short to split into %d chunks with an overlap of %d: %s", splitNumber, splitOverlap, file)
							return
						}
					}

					// a method to extract subsequence with given window size and step size
					slider = record.Seq.Slider(splitSize, step, circular2, greedy)

					if bySeq {
						slidIdx = 0
					}
					for {
						_seq, _ok = slider()
						if !_ok {
							break
						}

						if bySeq {
							*codes = (*codes)[:0] // reset
						}

						if len(_seq.Seq)-1 <= splitOverlap || len(_seq.Seq) < kMin {
							continue
						}

						for _, k = range ks {
							if syncmer {
								sketch, err = sketches.NewSyncmerSketch(_seq, k, syncmerS, circular)
							} else if minimizer {
								sketch, err = sketches.NewMinimizerSketch(_seq, k, minimizerW, circular)
							} else {
								iter, err = sketches.NewHashIterator(_seq, k, true, circular)
							}
							if err != nil {
								if err == sketches.ErrShortSeq {
									continue
								} else {
									checkError(errors.Wrapf(err, "seq: %s", record.Name))
								}
							}

							// place the if branch of deciding sketch type out of the for loop
							if syncmer {
								for {
									code, ok = sketch.NextSyncmer()
									if !ok {
										break
									}
									if scaled && code > maxHash {
										continue
									}
									if code > 0 {
										*codes = append(*codes, code)
									}
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
									if code > 0 {
										*codes = append(*codes, code)
									}
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
									if code > 0 {
										*codes = append(*codes, code)
									}
								}
							}
						}

						if !bySeq { // whole file
							break
						}

						n = len(*codes)

						*codes2 = (*codes2)[:0]
						// it's default now.
						// k-mer hashes are sorted and then duplicated values are removed.
						if exactNumber {
							sortutil.Uint64s(*codes)
							var pre uint64 = 0
							for _, code = range *codes {
								if code != pre {
									*codes2 = append(*codes2, code)
									pre = code
								}
							}
							n = len(*codes2)
						}

						// output file
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
							outFileBase = fmt.Sprintf("%s/%s-chunk_%d%s", baseFile, seqID, slidIdx, extDataFile)
						} else {
							outFileBase = fmt.Sprintf("%s-id_%s%s", baseFile, seqID, extDataFile)
						}
						// write to file
						if multiLevelFileTree {
							fileHash = xxh3.HashString(baseFile)
							dir3 = fmt.Sprintf("%03d", fileHash&1023)
						}

						outFile = filepath.Join(outDir, dir1, dir2, dir3, outFileBase)

						// meta data to save into the .unik files
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
							SplitNum:     _splitNumber,
							SplitSize:    splitSize0,
							SplitOverlap: splitOverlap,
						}
						if exactNumber {
							writeKmers(kMax, *codes2, uint64(n), outFile, compress, opt.CompressionLevel,
								scaled, scale, meta)
						} else {
							writeKmers(kMax, *codes, uint64(n), outFile, compress, opt.CompressionLevel,
								scaled, scale, meta)
						}

						ch <- UnikFileInfo{
							Path:       outFile,
							Name:       seqID,
							GenomeSize: uint64(len(record.Seq.Seq)),
							Index:      slidIdx,
							Indexes:    uint32(_splitNumber),
							Kmers:      uint64(n),
						}

						slidIdx++
					}

				}

				if bySeq {
					return
				}

				// whole file, no splitting

				n = len(*codes)

				if n == 0 {
					log.Warningf("skipping %s: no valid sequences", file)
					log.Info()
					return
				}

				*codes2 = (*codes2)[:0]
				if exactNumber {
					sortutil.Uint64s(*codes)
					var pre uint64 = 0
					for _, code = range *codes {
						if code != pre {
							*codes2 = append(*codes2, code)
							pre = code
						}
					}
					n = len(*codes2)
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
				slidIdx = 0
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
					SplitNum:     splitNumber, // 1
					SplitSize:    splitSize0,
					SplitOverlap: splitOverlap,
				}
				if exactNumber {
					writeKmers(kMax, *codes2, uint64(n), outFile, compress, opt.CompressionLevel,
						scaled, scale, meta)
				} else {
					writeKmers(kMax, *codes, uint64(n), outFile, compress, opt.CompressionLevel,
						scaled, scale, meta)
				}

				ch <- UnikFileInfo{
					Path:       outFile,
					Name:       fileprefix,
					GenomeSize: genomeSize,
					Index:      slidIdx,
					Indexes:    uint32(splitNumber),
					Kmers:      uint64(n),
				}
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

	var writer *unik.Writer
	var mode uint32

	mode |= unik.UnikCanonical // cononical k-mers
	mode |= unik.UnikHashed    // yes, we use ntHash
	mode |= unik.UnikSorted    // sorted k-mers could reduce the .unik file for k <= 23. https://github.com/shenwei356/unik#compression-rate-comparison

	writer, err = unik.NewWriter(outfh, k, mode)
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

	computeCmd.Flags().StringP("in-dir", "I", "",
		formatFlagUsage(`Directory containing FASTA/Q files. Directory symlinks are followed.`))

	computeCmd.Flags().StringP("file-regexp", "r", `\.(f[aq](st[aq])?|fna)(.gz)?$`,
		formatFlagUsage(`Regular expression for matching sequence files in -I/--in-dir, case ignored.`))

	computeCmd.Flags().StringP("out-dir", "O", "",
		formatFlagUsage(`Output directory.`))

	computeCmd.Flags().BoolP("force", "", false,
		formatFlagUsage(`Overwrite existed output directory.`))

	computeCmd.Flags().IntSliceP("kmer", "k", []int{21}, formatFlagUsage(`K-mer size(s). K needs to be <=64. Multiple values are supported, e.g., "-k 21,31" or "-k 21 -k 31"`))

	computeCmd.Flags().BoolP("circular", "", false,
		formatFlagUsage(`Input sequences are circular. Note that it only applies to genomes with a single chromosome.`))

	computeCmd.Flags().IntP("scale", "D", 1,
		formatFlagUsage(`Scale of the FracMinHash (Scaled MinHash), or down-sample factor for Syncmers and Minimizer.`))

	computeCmd.Flags().IntP("minimizer-w", "W", 0,
		formatFlagUsage(`Minimizer window size.`))

	computeCmd.Flags().IntP("syncmer-s", "S", 0,
		formatFlagUsage(`Length of the s-mer in Closed Syncmers.`))

	// computeCmd.Flags().BoolP("exact-number", "e", false, `save exact number of unique k-mers for indexing (recommended)`)

	computeCmd.Flags().BoolP("compress", "c", false,
		formatFlagUsage(`Output gzipped .unik files, it's slower and can save little space.`))

	computeCmd.Flags().IntP("split-number", "n", 0,
		formatFlagUsage(`Chunk number for splitting sequences, incompatible with -s/--split-size.`))

	computeCmd.Flags().IntP("split-size", "s", 0,
		formatFlagUsage(`Chunk size for splitting sequences, incompatible with -n/--split-number.`))

	computeCmd.Flags().IntP("split-overlap", "l", 0,
		formatFlagUsage(`Chunk overlap for splitting sequences. The default value will be set to k-1 unless you change it.`))

	computeCmd.Flags().IntP("split-min-ref", "m", 1000,
		formatFlagUsage(`Only splitting sequences >= X bp.`))

	computeCmd.Flags().BoolP("by-seq", "", false,
		formatFlagUsage(`Compute k-mers (sketches) for each sequence, instead of the whole file.`))

	computeCmd.Flags().StringP("ref-name-regexp", "N", `(?i)(.+)\.(f[aq](st[aq])?|fna)(.gz)?$`,
		formatFlagUsage(`Regular expression (must contains "(" and ")") for extracting reference name from filename.`))

	computeCmd.Flags().StringSliceP("seq-name-filter", "B", []string{},
		formatFlagUsage(`List of regular expressions for filtering out sequences by header/name, case ignored.`))

	computeCmd.SetUsageTemplate(usageTemplate("[-k <k>] [-n <chunks>] [-l <overlap>] {[-I <seqs dir>] | <seq files>} -O <out dir>"))

}

var reIgnoreCaseStr = "(?i)"
var reIgnoreCase = regexp.MustCompile(`\(\?i\)`)

var fileUnikInfos = "_info.txt"

var poolCodes = &sync.Pool{New: func() interface{} {
	tmp := make([]uint64, 0, mapInitSize)
	return &tmp
}}
