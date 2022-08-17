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
	"fmt"
	"io"
	"os"
	"path/filepath"
	"regexp"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/pkg/errors"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/bio/sketches"
	"github.com/shenwei356/util/pathutil"
	"github.com/shenwei356/xopen"
	"github.com/spf13/cobra"
)

var mergeGenomeCmd = &cobra.Command{
	Use:   "merge-genomes",
	Short: "merge genomes of the same species/strain and split them into chunks",
	Long: `merge genomes of the same species/strain and split them into chunks

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

		k := getFlagPositiveInt(cmd, "kmer")
		fragSize := getFlagPositiveInt(cmd, "frag-size")
		if fragSize < k {
			checkError(fmt.Errorf("value of flag -f/--frag-size should be >= that of -k/--kmer"))
		}

		circular0 := getFlagBool(cmd, "circular")

		outDir := getFlagString(cmd, "out-dir")
		outFile := getFlagString(cmd, "out-file")
		force := getFlagBool(cmd, "force")
		fileInfo := getFlagString(cmd, "info-file")
		outputInfo := fileInfo != ""

		if outDir == "" && outFile == "" {
			checkError(fmt.Errorf("flag -O/--out-dir or -o/--out-file is needed"))
		}

		var err error

		inDir := getFlagString(cmd, "in-dir")

		if outDir != "" && filepath.Clean(inDir) == filepath.Clean(outDir) {
			checkError(fmt.Errorf("input and output paths should not be the same"))
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

		// filterNames := len(reSeqNames) > 0

		// ---------------------------------------------------------------
		// flags for splitting sequence

		splitNumber0 := getFlagPositiveInt(cmd, "split-number")
		splitOverlap := getFlagNonNegativeInt(cmd, "split-overlap")
		splitMinRef := getFlagNonNegativeInt(cmd, "split-min-ref")
		if splitNumber0 <= 1 {
			checkError(fmt.Errorf("value of -n/--split-number should be > 1"))
		}
		if splitNumber0 > 65535 {
			checkError(fmt.Errorf(("value of flag -s/--split-number should not be greater than 65535")))
		}

		// ---------------------------------------------------------------
		// out dir
		var outdir string
		outputDir := outDir != ""

		if outputDir {
			makeOutDir(outDir, force)
			outdir = outDir
		} else {
			outdir = filepath.Join(os.TempDir(), fmt.Sprintf("kmcp.%d", os.Getpid()))
			makeOutDir(outdir, force)
			defer func() {
				checkError(os.RemoveAll(outdir))
				if opt.Verbose || opt.Log2File {
					log.Infof("remove tmp dir: %s", outdir)
				}
			}()
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
			checkError(fmt.Errorf("genome FASTA files needed"))
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
			log.Infof("    *regular expressions for filtering out sequences: %s", reSeqNameStrs)
			if outDir != "" {
				log.Infof("  output directory: %s", outDir)
			} else {
				log.Infof("       output file: %s", outFile)
			}
			if outputInfo {
				log.Infof("   extra info file: %s", fileInfo)
			}

			log.Info()

			log.Infof("sequences splitting: %v", true)
			log.Infof("  split parts: %d, overlap: %d bp", splitNumber0, splitOverlap)
			log.Info()

			log.Info("k-mer (sketches) computing:")

			log.Infof("  k-mer size: %d", k)

			log.Infof("  circular genome: %v", circular0)
			log.Info()

			log.Infof("-------------------- [main parameters] --------------------")
			log.Info()
		}

		// ---------------------------------------------------------------
		// choose the reference genome

		if opt.Verbose || opt.Log2File {
			log.Infof("choose the reference genome with the least contigs and the biggest genome size:")
		}
		refGenome := chooseRef(opt, files, reSeqNames)
		if opt.Verbose || opt.Log2File {
			log.Infof("  chosen ref genome: %s of %d bp with %d contig(s) ", refGenome.File, refGenome.Size, refGenome.Contigs)
		}

		// ---------------------------------------------------------------
		// splitting the reference genome

		if opt.Verbose || opt.Log2File {
			log.Info()
			log.Infof("splitting the reference genome ...")
		}

		seqs, hashes := splitGenome(opt, refGenome, reSeqNames, k, circular0,
			splitNumber0, splitOverlap, splitMinRef)

		outfhs := make([]*xopen.Writer, len(hashes))
		chunkFiles := make([]string, len(hashes))
		var chunkfile string
		for i := range hashes {
			chunkfile = filepath.Join(outdir, fmt.Sprintf("chunk%03d.fa.gz", i+1))
			outfhs[i], err = xopen.Wopen(chunkfile)
			if err != nil {
				checkError(err)
			}
			chunkFiles[i] = chunkfile
		}
		var infofh *xopen.Writer
		if outputInfo {
			infofh, err = xopen.Wopen(fileInfo)
			if err != nil {
				checkError(err)
			}
			fmt.Fprintf(infofh, "file\tseqId\tmKmers\tchunkId\tfragLoc\n")
		}
		defer func() {
			for _, outfh := range outfhs {
				outfh.Close()
			}
			if outputInfo {
				infofh.Close()
			}
		}()

		var buffer *bytes.Buffer
		var text []byte
		header := []byte(">seq\n")
		for i, _seq := range seqs {
			outfhs[i].Write(header)
			text, buffer = wrapByteSlice(_seq.Seq, 70, buffer)
			outfhs[i].Write(text)
			outfhs[i].Write(_mark_newline)
		}

		var nSizes []string
		var nHashes []string

		// ---------------------------------------------------------------
		// splitting other genomes

		if opt.Verbose || opt.Log2File {
			nSizes = make([]string, len(hashes))
			nHashes = make([]string, len(hashes))
			for i := 0; i < len(hashes); i++ {
				nHashes[i] = strconv.Itoa(len(hashes[i]))
				nSizes[i] = strconv.Itoa(len(seqs[i].Seq))
			}

			log.Infof("  seq length   of the %d chunks: %s", len(hashes), strings.Join(nSizes, ", "))
			log.Infof("  k-mer number of the %d chunks: %s", len(hashes), strings.Join(nHashes, ", "))
			log.Info()
			log.Infof("splitting other genomes and add subsequences to ref chunks ...")
		}

		step := fragSize - k + 1
		filterNames := len(reSeqNames) > 0
		for _, file := range files {
			if file == refGenome.File {
				continue
			}

			fastxReader, err := fastx.NewDefaultReader(file)
			checkError(errors.Wrap(err, file))

			var record *fastx.Record
			var ignoreSeq bool
			var re *regexp.Regexp

			var slider func() (*seq.Seq, bool)
			var _seq *seq.Seq
			var ok, _ok bool
			var code uint64
			var iter *sketches.Iterator

			codes := make([]uint64, 0, fragSize)
			hits := make([]int, len(hashes))
			var m map[uint64]interface{}
			var i, hit int
			var max int
			perfectN := fragSize - k + 1

			var loc int

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

				slider = record.Seq.Slider(fragSize, step, false, true)
				loc = 0
				for {
					_seq, _ok = slider()
					if !_ok {
						break
					}

					if len(_seq.Seq)-1 < k {
						loc += step
						continue
					}

					seqs = append(seqs, _seq)

					iter, err = sketches.NewHashIterator(_seq, k, true, false)
					if err != nil {
						if err == sketches.ErrShortSeq {
							continue
						} else {
							checkError(errors.Wrapf(err, "seq: %s", record.Name))
						}
					}

					for i = range hits {
						hits[i] = 0
					}
					codes = codes[:0]

					for {
						code, ok = iter.NextHash()
						if !ok {
							break
						}
						if code > 0 {
							for i, m = range hashes {
								if _, ok = m[code]; ok {
									hits[i]++
								}
							}
							codes = append(codes, code)
						}
					}

					max = 0
					for i, hit = range hits {
						if hit > max {
							max = hit
						}
					}
					if max == perfectN { // a fragment perfectly match one of the chunks
						loc += step
						continue
					}
					// the max may be 0, which means this fragment is new, so we just add it to all chunks
					for i, hit = range hits {
						if hit == max {
							outfhs[i].Write(header)
							text, buffer = wrapByteSlice(_seq.Seq, 70, buffer)
							outfhs[i].Write(text)
							outfhs[i].Write(_mark_newline)

							for _, code = range codes {
								hashes[i][code] = struct{}{}
							}

							if outputInfo {
								fmt.Fprintf(infofh, "%s\t%s\t%d\t%d\t%d\n",
									file, record.ID, hit, i+1, loc+1)
							}
						}
					}

					loc += step
				}

			}

		}

		if opt.Verbose || opt.Log2File {
			for i := 0; i < len(hashes); i++ {
				nHashes[i] = strconv.Itoa(len(hashes[i]))
			}
			log.Infof("  k-mer number of the %d chunks: %s", len(hashes), strings.Join(nHashes, ", "))
		}

		// ---------------------------------------------------
		// concatenate sequences in all chunks
		if !outputDir {
			for _, outfh := range outfhs {
				outfh.Close()
			}

			var outfh *xopen.Writer
			outfh, err = xopen.Wopen(outFile)
			if err != nil {
				checkError(err)
			}

			markSeq, _ := fastx.NewRecord(seq.Unlimit, []byte("marker-seq"), []byte("marker-seq"), nil, markerSeq)
			for _, file := range chunkFiles {
				fastxReader, err := fastx.NewDefaultReader(file)
				checkError(errors.Wrap(err, file))
				var record *fastx.Record
				for {
					record, err = fastxReader.Read()
					if err != nil {
						if err == io.EOF {
							break
						}
						checkError(errors.Wrap(err, file))
						break
					}

					record.FormatToWriter(outfh, 70)
				}

				markSeq.FormatToWriter(outfh, 70)
			}

			outfh.Close()

			if opt.Verbose || opt.Log2File {
				log.Infof("sequences saved to file: %s", outFile)
			}
		}
	},
}

func splitGenome(opt *Options, info *GenomeInfo, reSeqNames []*regexp.Regexp,
	k int, circular0 bool,
	splitNumber0 int, splitOverlap int, splitMinRef int) ([]*seq.Seq, []map[uint64]interface{}) {

	filterNames := len(reSeqNames) > 0

	// ------------------------------------------------------------------
	// concatenate all seqs

	file := info.File
	fastxReader, err := fastx.NewDefaultReader(file)
	checkError(errors.Wrap(err, file))

	var record *fastx.Record
	var ignoreSeq bool
	var re *regexp.Regexp

	var record1 *fastx.Record

	lenSum := 0
	first := true

	if info.Contigs == 1 { // just read the first valid sequence
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
		}
	} else { // concatenate all seqs
		var allSeqs [][]byte
		var bigSeq []byte

		allSeqs = make([][]byte, 0, 8)
		nnn := bytes.Repeat([]byte{'N'}, k-1)

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
			bigSeq = make([]byte, lenSum+(len(allSeqs)-1)*(k-1))
			i := 0
			for j, aseq := range allSeqs {
				copy(bigSeq[i:i+len(aseq)], aseq)
				if j < len(allSeqs)-1 {
					copy(bigSeq[i+len(aseq):i+len(aseq)+k-1], nnn)
				}
				i += len(aseq) + k - 1
			}
		}

		record1.Seq.Seq = bigSeq
	}

	record1.Seq.Qual = nil
	record = record1

	// ------------------------------------------------------------------
	// split

	var step int
	var greedy bool = true
	var seqLen int
	var splitSize int
	var splitNumber int
	splitByNumber := splitNumber0 > 1

	seqLen = len(record.Seq.Seq)

	splitNumber = splitNumber0
	if seqLen < splitMinRef {
		splitSize = seqLen
		step = seqLen
		greedy = false
		splitNumber = 1
	} else if splitByNumber {
		if splitNumber == 1 {
			splitSize = seqLen
			step = seqLen
			greedy = false
		} else {
			splitSize = (seqLen + (splitNumber-1)*(splitOverlap+k-1) + splitNumber - 1) / splitNumber
			step = splitSize - splitOverlap
		}
	}
	if opt.Verbose || opt.Log2File {
		log.Infof("  genome size: %d, splitNumber: %d, splitSize: %d, step: %d, greedy: %v\n", seqLen, splitNumber, splitSize, step, greedy)
	}
	var code uint64
	var iter *sketches.Iterator

	var slider func() (*seq.Seq, bool)
	var _seq *seq.Seq
	var ok, _ok bool
	var slidIdx uint32
	hashes := make([]map[uint64]interface{}, 0, splitNumber0)

	seqs := make([]*seq.Seq, 0, splitNumber0)

	slider = record.Seq.Slider(splitSize, step, circular0, greedy)
	for {
		_seq, _ok = slider()
		if !_ok {
			break
		}

		if len(_seq.Seq)-1 <= splitOverlap {
			continue
		}

		seqs = append(seqs, _seq)

		iter, err = sketches.NewHashIterator(_seq, k, true, false)
		if err != nil {
			if err == sketches.ErrShortSeq {
				continue
			} else {
				checkError(errors.Wrapf(err, "seq: %s", record.Name))
			}
		}

		codes := make(map[uint64]interface{}, mapInitSize)
		for {
			code, ok = iter.NextHash()
			if !ok {
				break
			}
			if code > 0 {
				codes[code] = struct{}{}
			}
		}

		hashes = append(hashes, codes)

		slidIdx++
	}

	return seqs, hashes
}

type GenomeInfo struct {
	File    string
	Contigs int
	Size    int
}

// choose the genome with minimal contigs and biggest genome size
func chooseRef(opt *Options, files []string, reSeqNames []*regexp.Regexp) *GenomeInfo {
	filterNames := len(reSeqNames) > 0

	infos := make([]*GenomeInfo, 0, len(files))

	var wg sync.WaitGroup
	tokens := make(chan int, opt.NumCPUs)
	ch := make(chan *GenomeInfo, len(files))
	done := make(chan int)

	go func() {
		for info := range ch {
			infos = append(infos, info)
		}
		done <- 1
	}()

	for _, file := range files {
		tokens <- 1
		wg.Add(1)

		go func(file string) {
			defer func() {
				wg.Done()
				<-tokens
			}()

			fastxReader, err := fastx.NewDefaultReader(file)
			checkError(errors.Wrap(err, file))

			var record *fastx.Record
			var contigs int
			var size int
			var ignoreSeq bool
			var re *regexp.Regexp
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

				if len(record.Seq.Seq) == 0 {
					log.Warningf("skipping %s: no valid sequences", file)
					continue
				}

				contigs++
				size += len(record.Seq.Seq)
			}

			ch <- &GenomeInfo{File: file, Contigs: contigs, Size: size}
		}(file)
	}

	wg.Wait()
	close(ch)
	<-done

	sort.Slice(infos, func(i, j int) bool { // sort by contigs and then size
		a := infos[i].Contigs - infos[j].Contigs
		if a == 0 {
			return infos[i].Size > infos[j].Size
		}
		return a < 0
	})

	return infos[0]
}

func init() {
	utilsCmd.AddCommand(mergeGenomeCmd)

	mergeGenomeCmd.Flags().StringP("in-dir", "I", "",
		formatFlagUsage(`Directory containing FASTA files. Directory symlinks are followed.`))

	mergeGenomeCmd.Flags().StringP("file-regexp", "r", `\.(f[aq](st[aq])?|fna)(.gz)?$`,
		formatFlagUsage(`Regular expression for matching sequence files in -I/--in-dir, case ignored.`))

	mergeGenomeCmd.Flags().StringP("out-file", "o", "",
		formatFlagUsage(`Out file.`))

	mergeGenomeCmd.Flags().StringP("out-dir", "O", "",
		formatFlagUsage(`Output directory.`))

	mergeGenomeCmd.Flags().StringP("info-file", "", "",
		formatFlagUsage(`An extra output file to show which chunk(s) are assigned to for each genome fragment.`))

	mergeGenomeCmd.Flags().BoolP("force", "", false,
		formatFlagUsage(`Overwrite existed output directory.`))

	mergeGenomeCmd.Flags().IntP("kmer", "k", 21, formatFlagUsage(`K-mer size.`))

	mergeGenomeCmd.Flags().BoolP("circular", "", false,
		formatFlagUsage(`Input sequences are circular.`))

	mergeGenomeCmd.Flags().IntP("split-number", "n", 0,
		formatFlagUsage(`Chunk number for splitting sequences, incompatible with -s/--split-size.`))

	mergeGenomeCmd.Flags().IntP("split-overlap", "l", 0,
		formatFlagUsage(`Chunk overlap for splitting sequences.`))

	mergeGenomeCmd.Flags().IntP("split-min-ref", "m", 1000,
		formatFlagUsage(`Only splitting sequences >= X bp.`))

	mergeGenomeCmd.Flags().StringSliceP("seq-name-filter", "B", []string{},
		formatFlagUsage(`List of regular expressions for filtering out sequences by header/name, case ignored.`))

	mergeGenomeCmd.Flags().IntP("frag-size", "f", 100,
		formatFlagUsage(`size of sequence fragments to be re-assigned to the reference genome chunks.`))

	mergeGenomeCmd.SetUsageTemplate(usageTemplate("[-k <k>] [-n <chunks>] [-l <overlap>] {[-I <seqs dir>] | <seq files>} {-O <out dir> | -o <out file>"))

}

var _mark_fasta = []byte{'>'}
var _mark_fastq = []byte{'@'}
var _mark_plus_newline = []byte{'+', '\n'}
var _mark_newline = []byte{'\n'}

var markerSeq = []byte("--kmcp--marker-seq--")
