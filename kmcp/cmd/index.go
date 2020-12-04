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
	"bufio"
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"regexp"
	"runtime"
	"sort"
	"sync"
	"time"

	"github.com/pkg/errors"
	"github.com/shenwei356/unikmer"
	"github.com/shenwei356/unikmer/index"
	"github.com/shenwei356/util/bytesize"
	"github.com/shenwei356/util/pathutil"
	"github.com/spf13/cobra"
	"github.com/twotwotwo/sorts"
	"github.com/twotwotwo/sorts/sortutil"
	"github.com/vbauerster/mpb"
	"github.com/vbauerster/mpb/decor"
)

var indexCmd = &cobra.Command{
	Use:   "index",
	Short: "Construct database from k-mer files",
	Long: `Construct database from k-mer files

Attentions:
  0. All input files should be generated by "kmcp compute".
  1. If input files contain minimizer/syncmer sketchs, you must set 
     the same parameter in this command, or manually edit "__db.yml".

Tips:
  1. Increase value of -j/--threads for acceleratation in cost of more
     memory occupation and I/O pressure, sqrt(#cpus) is recommended.
     #threads * #threads files are simultaneously opened, and max number
     of opened files is limited by flag -F/--max-open-files.
  2. Value of block size -b/--block-size better be multiple of 64.
  3. Use --dry-run to adjust parameters and see final #index files 
     and total file size.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		runtime.GOMAXPROCS(opt.NumCPUs)

		timeStart := time.Now()
		defer func() {
			if opt.Verbose {
				log.Infof("elapsed time: %s", time.Since(timeStart))
			}
		}()

		// ---------------------------------------------------------------
		// basic flags

		var err error

		dryRun := getFlagBool(cmd, "dry-run")
		if dryRun {
			opt.Verbose = true
		}

		outDir := getFlagString(cmd, "out-dir")

		inDir := getFlagString(cmd, "in-dir")

		if outDir == "" && inDir == "" {
			checkError(fmt.Errorf("flag -I/--in-dir OR -O/--out-dir is needed"))
		}

		readFromDir := inDir != ""
		if readFromDir {
			var isDir bool
			isDir, err = pathutil.IsDir(inDir)
			checkError(errors.Wrapf(err, "checking -I/--in-dir: %s", inDir))
			if !isDir {
				checkError(fmt.Errorf("value of -I/--in-dir should be a directory: %s", inDir))
			}
		}

		reFileStr := getFlagString(cmd, "file-regexp")
		var reFile *regexp.Regexp
		if reFileStr != "" {
			reFile, err = regexp.Compile(reFileStr)
			checkError(errors.Wrapf(err, "parsing regular expression for matching file: %s", reFileStr))
		}

		force := getFlagBool(cmd, "force")

		alias := getFlagString(cmd, "alias")

		// ---------------------------------------------------------------
		// index flags

		sBlock := getFlagInt(cmd, "block-size")

		fpr := getFlagPositiveFloat64(cmd, "false-positive-rate")
		numHashes := getFlagPositiveInt(cmd, "num-hash")
		if numHashes > 255 {
			checkError(fmt.Errorf("value of -n/--num-hash too big: %d", numHashes))
		}

		maxOpenFiles := getFlagPositiveInt(cmd, "max-open-files")

		// block-max-kmers-t1
		kmerThreshold8Str := getFlagString(cmd, "block-max-kmers-t1")
		kmerThreshold8Float, err := bytesize.ParseByteSize(kmerThreshold8Str)
		if err != nil {
			checkError(fmt.Errorf("invalid size: %s", kmerThreshold8Str))
		}
		kmerThreshold8 := int64(kmerThreshold8Float)
		if kmerThreshold8 <= 0 {
			checkError(fmt.Errorf("value of flag -m/--block-max-kmers-t1 should be positive: %d", kmerThreshold8))
		}

		// block-max-kmers-t2
		kmerThresholdSStr := getFlagString(cmd, "block-max-kmers-t2")
		kmerThresholdSFloat, err := bytesize.ParseByteSize(kmerThresholdSStr)
		if err != nil {
			checkError(fmt.Errorf("invalid size: %s", kmerThresholdSStr))
		}
		kmerThresholdS := int64(kmerThresholdSFloat)
		if kmerThresholdS <= 0 {
			checkError(fmt.Errorf("value of flag -M/--block-max-kmers-t2 should be positive: %d", kmerThresholdS))
		}

		if kmerThreshold8 >= kmerThresholdS {
			checkError(fmt.Errorf("value of flag -m/--block-max-kmers-t1 (%d) should be small than -M/--block-max-kmers-t2 (%d)", kmerThreshold8, kmerThresholdS))
		}

		// ---------------------------------------------------------------
		// flags of sketch

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

		// -------------------------------------------------------
		// flags for extracting names

		nameRegexp := getFlagString(cmd, "name-regexp")
		extractName := nameRegexp != ""
		var reName *regexp.Regexp
		if extractName {
			if !regexp.MustCompile(`\(.+\)`).MatchString(nameRegexp) {
				checkError(fmt.Errorf(`value of -r (--name-regexp) must contains "(" and ")" to extract name`))
			}
			reName, err = regexp.Compile(nameRegexp)
			checkError(errors.Wrapf(err, "parsing regular expression for extracting name: %s", nameRegexp))
		}

		// ---------------------------------------------------------------
		// out dir

		if outDir == "" {
			if inDir != "" {
				outDir = inDir
			}
		} else if !dryRun {
			makeOutDir(outDir, force)
		}
		if alias == "" {
			alias = filepath.Base(outDir)
		}

		// ---------------------------------------------------------------
		// input files

		if opt.Verbose {
			log.Info("checking input files ...")
		}

		var files []string
		if readFromDir {
			files, err = getFileListFromDir(inDir, reFile)
			checkError(errors.Wrapf(err, "err on walking dir: %s", inDir))
		} else {
			files = getFileListFromArgsAndFile(cmd, args, true, "infile-list", true)
			if opt.Verbose {
				if len(files) == 1 && isStdin(files[0]) {
					log.Info("no files given, reading from stdin")
				} else {
					log.Infof("%d input file(s) given", len(files))
				}
			}
		}
		if len(files) < 2 {
			checkError(fmt.Errorf("at least two .unik files needed"))
		} else if opt.Verbose {
			log.Infof("%d input file(s) given", len(files))
		}

		// ---------------------------------------------------------------
		// log

		if opt.Verbose {
			log.Infof("-------------------- [main parameters] --------------------")
			if minimizer {
				log.Infof("Minimizer window: %d", minimizerW)
			}
			if syncmer {
				log.Infof("syncmer s: %d", syncmerS)
			}

			log.Infof("number of hashes: %d", numHashes)
			log.Infof("false positive rate: %f", fpr)
			log.Infof("-------------------- [main parameters] --------------------")

		}

		// ---------------------------------------------------------------
		// check unik files and read k-mers numbers

		if opt.Verbose {
			log.Info("checking .unik files ...")
		}

		var pbs *mpb.Progress
		var bar *mpb.Bar
		var chDuration chan time.Duration
		var doneDuration chan int

		if opt.Verbose {
			pbs = mpb.New(mpb.WithWidth(50), mpb.WithOutput(os.Stderr))
			bar = pbs.AddBar(int64(len(files)),
				mpb.BarStyle("[=>-]<+"),
				mpb.PrependDecorators(
					decor.Name("checking .unik file: ", decor.WC{W: len("checking .unik file: "), C: decor.DidentRight}),
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

		fileInfos := make([]UnikFileInfo, 0, len(files))

		var k int = -1
		var hashed bool
		var canonical bool
		var scaled bool
		var scale uint32
		var n int64
		var nfiles = len(files)

		var reader0 *unikmer.Reader

		getInfo := func(file string, first bool) UnikFileInfo {
			infh, r, _, err := inStream(file)
			checkError(err)

			reader, err := unikmer.NewReader(infh)
			checkError(errors.Wrap(err, file))

			if first {
				reader0 = reader
				k = reader.K
				hashed = reader.IsHashed()
				if !hashed {
					if syncmer {
						checkError(fmt.Errorf(`flag --syncmer-s given but the k-mers in .unik file are not hashed: %s`, file))
					}
					if minimizer {
						checkError(fmt.Errorf(`flag --Minimizer-w given but the k-mers in .unik file are not hashed: %s`, file))
					}
				}

				canonical = reader.IsCanonical()
				if !canonical {
					checkError(fmt.Errorf(`files with 'canonical' flag needed: %s`, file))
				}
				scaled = reader.IsScaled()
				scale = reader.GetScale()
			} else {
				checkCompatibility(reader0, reader, file)
				if scaled && scale != reader.GetScale() {
					checkError(fmt.Errorf(`scales not consistent, please check with "unikmer stats": %s`, file))
				}
			}

			if reader.Number < 0 {
				checkError(fmt.Errorf("binary file not sorted or no k-mers number found: %s", file))
			}
			var name string
			if extractName {
				found := reName.FindAllStringSubmatch(filepath.Base(file), -1)
				if len(found) > 0 {
					name = found[0][1]
				} else {
					name = filepath.Base(file)
				}
			} else {
				name = filepath.Base(file)
			}

			var relPath string
			relPath, err = filepath.Rel(outDir, file)
			checkError(err)

			checkError(r.Close())
			return UnikFileInfo{Path: file, RelPath: relPath, Name: name, Kmers: reader.Number}
		}

		file := files[0]
		var t time.Time
		if opt.Verbose {
			t = time.Now()
		}

		info := getInfo(file, true)

		if opt.Verbose {
			bar.Increment()
			bar.DecoratorEwmaUpdate(time.Since(t))
		}

		fileInfos = append(fileInfos, info)
		namesMap := make(map[string]interface{}, nfiles)
		namesMap[info.Name] = struct{}{}

		// left files
		var wgGetInfo sync.WaitGroup
		chInfos := make(chan UnikFileInfo, opt.NumCPUs)
		tokensGetInfo := make(chan int, opt.NumCPUs)
		doneGetInfo := make(chan int)
		go func() {
			var ok bool
			for info := range chInfos {
				fileInfos = append(fileInfos, info)
				n += info.Kmers

				if _, ok = namesMap[info.Name]; ok {
					log.Warningf("duplicated name: %s", info.Name)
				} else {
					namesMap[info.Name] = struct{}{}
				}
			}
			doneGetInfo <- 1
		}()

		for _, file := range files[1:] {
			wgGetInfo.Add(1)
			tokensGetInfo <- 1
			go func(file string) {
				defer func() {
					wgGetInfo.Done()
					<-tokensGetInfo
				}()
				var t time.Time
				if opt.Verbose {
					t = time.Now()
				}

				chInfos <- getInfo(file, false)

				if opt.Verbose {
					chDuration <- time.Duration(float64(time.Since(t)) / float64(opt.NumCPUs))
				}
			}(file)
		}

		wgGetInfo.Wait()
		close(chInfos)
		<-doneGetInfo

		if opt.Verbose {
			close(chDuration)
			<-doneDuration
			pbs.Wait()
		}

		if opt.Verbose {
			log.Infof("finished checking %d .unik files", nfiles)
		}

		sorts.Quicksort(UnikFileInfos(fileInfos))

		names0 := make([]string, 0, len(files))
		sizes0 := make([]uint64, 0, len(files))
		paths0 := make([]string, 0, len(files))
		for _, info := range fileInfos {
			n += info.Kmers
			names0 = append(names0, info.Name)
			sizes0 = append(sizes0, uint64(info.Kmers))
			paths0 = append(paths0, info.RelPath)
		}

		if opt.Verbose {
			log.Infof("names (first and last 5 names):")
			for i, info := range fileInfos {
				if i < 5 || i >= len(fileInfos)-5 {
					log.Infof("name [%d/%d]: %s, #k-mers: %d, file: %s ", i+1, len(fileInfos), info.Name, info.Kmers, info.Path)
				}
			}
		}

		// ------------------------------------------------------------------------------------
		// begin creating index
		if opt.Verbose {
			log.Infof("starting indexing ...")
		}
		nFiles := len(fileInfos)
		if sBlock <= 0 {
			sBlock = int((math.Sqrt(float64(nFiles))+7)/8) * 8
		}
		if sBlock < 8 {
			sBlock = 8
		} else if sBlock > nFiles {
			sBlock = nFiles
		}

		if opt.Verbose {
			log.Infof("-------------------- [important parameters] --------------------")
			log.Infof("block size: %d", sBlock)
			log.Infof("block-max-kmers-threshold 1: %s", bytesize.ByteSize(kmerThreshold8))
			log.Infof("block-max-kmers-threshold 2: %s", bytesize.ByteSize(kmerThresholdS))
			log.Infof("-------------------- [important parameters] --------------------")
			log.Infof("buiding index ...")
		}

		if opt.Verbose {
			pbs = mpb.New(mpb.WithWidth(60), mpb.WithOutput(os.Stderr))
		}

		nIndexFiles := int((len(files) + sBlock - 1) / sBlock) // may be more if using -m and -M
		indexFiles := make([]string, 0, nIndexFiles)

		ch := make(chan string, nIndexFiles)
		done := make(chan int)
		go func() {
			for f := range ch {
				indexFiles = append(indexFiles, f)
			}
			done <- 1
		}()

		var fileSize float64
		chFileSize := make(chan float64, nIndexFiles)
		doneFileSize := make(chan int)
		go func() {
			for f := range chFileSize {
				fileSize += f
			}
			doneFileSize <- 1
		}()

		var prefix string

		var b int
		var wg0 sync.WaitGroup
		maxConc := opt.NumCPUs
		if dryRun {
			maxConc = 1 // just for loging in order
		}
		tokens0 := make(chan int, maxConc)
		tokensOpenFiles := make(chan int, maxOpenFiles)

		sBlock0 := sBlock

		batch := make([]*UnikFileInfo, 0, sBlock)
		blockSizes := make([]int, 0, 8)
		var flag8, flag bool
		var lastInfo *UnikFileInfo
		for i := 0; i <= nFiles; i++ {
			if i == nFiles { // process lastInfo
				// process files in batch or the last file
				if flag {
					batch = append(batch, lastInfo)
				}
			} else {
				info := fileInfos[i]

				if flag || flag8 {
					if flag {
						batch = append(batch, lastInfo)
						lastInfo = &info
					} else if info.Kmers > kmerThresholdS {
						// meet a very big file the first time
						flag = true      // mark
						lastInfo = &info // leave this file process in the next round
						// and we have to process files aleady in batch
					} else { // flag8 && !flag
						if lastInfo != nil { // right after found the first file > kmerThreshold8
							batch = append(batch, lastInfo)
							lastInfo = nil
						}

						batch = append(batch, &info)
						if len(batch) < sBlock { // not filled
							continue
						}
					}
				} else if info.Kmers > kmerThreshold8 {
					if info.Kmers > kmerThresholdS {
						// meet a very big file the first time
						flag = true      // mark
						lastInfo = &info // leave this file process in the next round
						// and we have to process files aleady in batch
					} else {
						// meet a big file > kmerThreshold8
						sBlock = 8
						flag8 = true     // mark
						lastInfo = &info // leave this file process in the next batch
						// and we have to process files aleady in batch
					}
				} else {
					batch = append(batch, &info)
					if len(batch) < sBlock { // not filled
						continue
					}
				}

			}

			if len(batch) == 0 {
				break
			}

			blockSizes = append(blockSizes, len(batch))

			b++

			prefix = fmt.Sprintf("[block #%03d]", b)

			var bar *mpb.Bar
			if opt.Verbose && !dryRun {
				bar = pbs.AddBar(int64((len(batch)+7)/8),
					mpb.PrependDecorators(
						decor.Name(prefix, decor.WC{W: len(prefix) + 1, C: decor.DidentRight}),
						decor.CountersNoUnit("%d / %d", decor.WCSyncWidth),
					),
					mpb.AppendDecorators(decor.Percentage(decor.WC{W: 5})),
					mpb.BarFillerClearOnComplete(),
				)
			}

			wg0.Add(1)
			tokens0 <- 1
			go func(files []*UnikFileInfo, b int, prefix string, bar *mpb.Bar) {
				var wg sync.WaitGroup
				tokens := make(chan int, opt.NumCPUs)

				var maxElements int64
				for _, info := range files {
					if maxElements < info.Kmers {
						maxElements = info.Kmers
					}
				}

				var nBatchFiles int
				nBatchFiles = int((len(files) + 7) / 8)

				sigsBlock := make([][]byte, 0, nBatchFiles)
				namesBlock := make([]string, 0, len(files))
				sizesBlock := make([]uint64, 0, len(files))

				chBatch8 := make(chan batch8, nBatchFiles)
				doneBatch8 := make(chan int)

				buf := make(map[int]batch8)

				go func() {
					var id int = 1
					for batch := range chBatch8 {
						if batch.id == id {
							sigsBlock = append(sigsBlock, batch.sigs)
							namesBlock = append(namesBlock, batch.names...)
							sizesBlock = append(sizesBlock, batch.sizes...)
							if opt.Verbose && !dryRun {
								bar.Increment()
							}
							id++
							continue
						}
						for {
							if _batch, ok := buf[id]; ok {
								sigsBlock = append(sigsBlock, _batch.sigs)
								namesBlock = append(namesBlock, _batch.names...)
								sizesBlock = append(sizesBlock, _batch.sizes...)
								if opt.Verbose && !dryRun {
									bar.Increment()
								}
								delete(buf, id)
								id++
							} else {
								break
							}
						}
						buf[batch.id] = batch
					}
					if len(buf) > 0 {
						ids := make([]int, 0, len(buf))
						for id := range buf {
							ids = append(ids, id)
						}
						sort.Ints(ids)
						for _, id := range ids {
							_batch := buf[id]

							sigsBlock = append(sigsBlock, _batch.sigs)
							namesBlock = append(namesBlock, _batch.names...)
							sizesBlock = append(sizesBlock, _batch.sizes...)
							if opt.Verbose && !dryRun {
								bar.Increment()
							}
						}
					}

					doneBatch8 <- 1
				}()

				numSigs := CalcSignatureSize(uint64(maxElements), numHashes, fpr)
				var eFileSize float64
				eFileSize = 24
				for _, info := range files {
					eFileSize += float64(len(info.Name) + 1 + 8)
				}
				eFileSize += float64(numSigs * uint64(nBatchFiles))

				if opt.Verbose && dryRun {
					log.Infof("%s #files: %d, max #k-mers: %d, #signatures: %d, file size: %s", prefix, len(files), maxElements, numSigs, bytesize.ByteSize(eFileSize))
				}

				// split into batches with 8 files
				var bb, jj int
				for ii := 0; ii < len(files); ii += 8 {
					if dryRun {
						continue
					}
					jj = ii + 8
					if jj > len(files) {
						jj = len(files)
					}
					wg.Add(1)
					tokens <- 1
					bb++

					var outFile string

					// 8 files
					go func(_files []*UnikFileInfo, bb int, maxElements int64, numSigs uint64, outFile string, id int) {
						defer func() {
							wg.Done()
							<-tokens
						}()

						names := make([]string, 0, 8)
						sizes := make([]uint64, 0, 8)
						for _, info := range _files {
							names = append(names, info.Name)
							sizes = append(sizes, uint64(info.Kmers))
						}

						sigs := make([]byte, numSigs)

						// every file in 8 files
						for _k, info := range _files {
							var infh *bufio.Reader
							var r *os.File
							var reader *unikmer.Reader
							var err error
							var code uint64
							var loc int

							tokensOpenFiles <- 1
							infh, r, _, err = inStream(info.Path)
							checkError(errors.Wrap(err, info.Path))

							reader, err = unikmer.NewReader(infh)
							checkError(errors.Wrap(err, info.Path))
							singleHash := numHashes == 1

							if reader.IsHashed() {
								if singleHash {
									for {
										code, _, err = reader.ReadCodeWithTaxid()
										if err != nil {
											if err == io.EOF {
												break
											}
											checkError(errors.Wrap(err, info.Path))
										}

										sigs[code%numSigs] |= 1 << (7 - _k)
									}
								} else {
									for {
										code, _, err = reader.ReadCodeWithTaxid()
										if err != nil {
											if err == io.EOF {
												break
											}
											checkError(errors.Wrap(err, info.Path))
										}

										for _, loc = range hashLocations(code, numHashes, numSigs) {
											sigs[loc] |= 1 << (7 - _k)
										}
									}
								}
							} else {
								if singleHash {
									for {
										code, _, err = reader.ReadCodeWithTaxid()
										if err != nil {
											if err == io.EOF {
												break
											}
											checkError(errors.Wrap(err, info.Path))
										}

										sigs[hash64(code)%numSigs] |= 1 << (7 - _k)
									}
								} else {
									for {
										code, _, err = reader.ReadCodeWithTaxid()
										if err != nil {
											if err == io.EOF {
												break
											}
											checkError(errors.Wrap(err, info.Path))
										}

										for _, loc = range hashLocations(hash64(code), numHashes, numSigs) {
											sigs[loc] |= 1 << (7 - _k)
										}
									}
								}
							}

							r.Close()

							<-tokensOpenFiles
						}

						chBatch8 <- batch8{
							id:    id,
							sigs:  sigs,
							names: names,
							sizes: sizes,
						}
						// if opt.Verbose {
						// 	log.Infof("%s batch #%03d/%d: %d signatures loaded\r", prefix, bb, nBatchFiles, len(sigs))
						// }

					}(files[ii:jj], bb, maxElements, numSigs, outFile, bb)
				}

				wg.Wait()
				close(chBatch8)
				<-doneBatch8

				blockFile := filepath.Join(outDir, fmt.Sprintf("_block%03d%s", b, extIndex))

				if !dryRun {

					outfh, gw, w, err := outStream(blockFile, false, opt.CompressionLevel)
					checkError(err)
					defer func() {
						outfh.Flush()
						if gw != nil {
							gw.Close()
						}
						w.Close()
					}()

					writer, err := index.NewWriter(outfh, k, canonical, uint8(numHashes), numSigs, namesBlock, sizesBlock)
					checkError(err)
					defer func() {
						checkError(writer.Flush())
					}()

					if nBatchFiles == 1 {
						checkError(writer.WriteBatch(sigsBlock[0], len(sigsBlock[0])))
					} else {
						row := make([]byte, nBatchFiles)
						for ii := 0; ii < int(numSigs); ii++ {
							for jj = 0; jj < nBatchFiles; jj++ {
								row[jj] = sigsBlock[jj][ii]
							}
							checkError(writer.Write(row))
						}
					}
					// if opt.Verbose {
					// 	log.Infof("%s signatures saved", prefix)
					// }
				}

				ch <- filepath.Base(blockFile)
				chFileSize <- eFileSize

				wg0.Done()
				<-tokens0
			}(batch, b, prefix, bar)

			batch = make([]*UnikFileInfo, 0, sBlock)
		}

		wg0.Wait()

		close(ch)
		close(chFileSize)
		<-done
		<-doneFileSize

		if opt.Verbose && !dryRun {
			pbs.Wait()
		}

		if !dryRun {
			sortutil.Strings(indexFiles)
			dbInfo := NewUnikIndexDBInfo(indexFiles)
			dbInfo.Alias = alias
			dbInfo.K = k
			dbInfo.Hashed = hashed
			dbInfo.Kmers = int(n)
			dbInfo.FPR = fpr
			dbInfo.BlockSize = sBlock0
			dbInfo.BlockSizes = blockSizes
			dbInfo.NumNames = len(names0)
			dbInfo.Names = names0
			dbInfo.Sizes = sizes0
			dbInfo.Paths = paths0
			dbInfo.NumHashes = numHashes
			dbInfo.Canonical = canonical
			dbInfo.Scaled = scaled
			dbInfo.Scale = scale
			dbInfo.Minimizer = minimizer
			dbInfo.MinimizerW = uint32(minimizerW)
			dbInfo.Syncmer = syncmer
			dbInfo.SyncmerS = uint32(syncmerS)
			checkError(dbInfo.WriteTo(filepath.Join(outDir, dbInfoFile)))

			// write name_mapping.tsv
			func() {
				outfh, gw, w, err := outStream(filepath.Join(outDir, dbNameMappingFile), false, opt.CompressionLevel)
				checkError(err)
				defer func() {
					outfh.Flush()
					if gw != nil {
						gw.Close()
					}
					w.Close()
				}()

				var name string
				for _, _name := range names0 {
					name, _ = filepathTrimExtension(_name)
					outfh.WriteString(fmt.Sprintf("%s\t%s\n", _name, name))
				}
			}()
		}

		// ------------------------------------------------------------------------------------

		if opt.Verbose {
			log.Infof("kmcp database with %d k-mers saved to %s", n, outDir)
			log.Infof("#index files: %d, total file size: %s", len(indexFiles), bytesize.ByteSize(fileSize))
		}
	},
}

func init() {
	RootCmd.AddCommand(indexCmd)

	indexCmd.Flags().StringP("in-dir", "I", "", `directory containing .unik files. if it's given and -O/--out-dir is empty, index files will be saved to -I/--in-dir`)
	indexCmd.Flags().StringP("file-regexp", "", ".unik$", `regular expression for matching files to index`)

	indexCmd.Flags().StringP("out-dir", "O", "", `output directory. if it's not given and -I/--in-dir is specified, index files will be saved to -I/--in-dir`)
	indexCmd.Flags().StringP("alias", "a", "", `database alias/name, default: basename of --out-dir/--in-dir. you can also manually edit in info file __db.yml`)

	indexCmd.Flags().Float64P("false-positive-rate", "f", 0.3, `false positive rate of single bloom filter`)
	indexCmd.Flags().IntP("num-hash", "n", 1, `number of hashes`)
	indexCmd.Flags().IntP("block-size", "b", 0, `block size, better be multiple of 64 for large number of input files. default: sqrt(#.files)`)
	indexCmd.Flags().StringP("block-max-kmers-t1", "m", "20M", `if kmers of single .unik file exceeds this threshold, change block size is changed to 8. unit supported: K, M, G`)
	indexCmd.Flags().StringP("block-max-kmers-t2", "M", "200M", `if kmers of single .unik file exceeds this threshold, an individual index is created for this file. unit supported: K, M, G`)

	indexCmd.Flags().BoolP("force", "", false, "overwrite tmp dir")
	indexCmd.Flags().StringP("name-regexp", "r", "", "regular expression for extracting name from .unik file name. if not given, base name are saved")
	indexCmd.Flags().IntP("max-open-files", "F", 256, "maximum number of opened files")
	indexCmd.Flags().BoolP("dry-run", "", false, "dry run, useful to adjust parameters")

	indexCmd.Flags().IntP("minimizer-w", "W", 0, `minimizer window size`)
	indexCmd.Flags().IntP("syncmer-s", "S", 0, `syncmer s`)

}

// batch8 contains data from 8 files, just for keeping order of all files of a block
type batch8 struct {
	id int

	sigs  []byte
	names []string
	sizes []uint64
}
