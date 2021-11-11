# Usage

KMCP is a command-line tool consisting of several subcommands.

```text

    Program: kmcp (K-mer-based Metagenomic Classification and Profiling)
    Version: v0.7.0
  Documents: https://shenwei356.github.io/kmcp
Source code: https://github.com/shenwei356/kmcp

KMCP is a tool for metagenomic classification and profiling.

KMCP can also be used for:
  1. Fast sequence search from large scales of genomic datasets
     as BIGSI and COBS do.
  2. Fast assembly/genome similarity estimation as Mash and sourmash do,
     by utilizing Minimizer, Scaled MinHash, or Closed Syncmers.

Usage:
  kmcp [command]

Available Commands:
  autocompletion Generate shell autocompletion script
  compute        Generate k-mers (sketch) from FASTA/Q sequences
  index          Construct database from k-mer files
  merge          Merge search results from multiple databases
  profile        Generate taxonomic profile from search results
  search         Search sequence against a database
  utils          Some utilities
  version        Print version information and check for update

Flags:
  -h, --help                 help for kmcp
  -i, --infile-list string   file of input files list (one file per line), if given, they are appended to files from cli arguments
      --log string           log file
  -q, --quiet                do not print any verbose information. you can write them to file with --log
  -j, --threads int          number of CPUs to use (default 16)

Use "kmcp [command] --help" for more information about a command.

```


## compute

```text
Generate k-mers (sketchs) from FASTA/Q sequences

Attentions:
  1. Input files can be given as list of FASTA/Q files via
     positional arguments or a directory containing sequence files
     via the flag -I/--in-dir. A regular expression for matching
     sequencing files is available by the flag -r/--file-regexp.
  2. Unwanted sequence like plasmid can be filtered out by
     the name via regular expressions (-B/--seq-name-filter).
  3. By default, kmcp computes k-mers (sketches) of every file,
     you can also use --by-seq to compute for every sequence,
     where sequence IDs in all input files better be distinct.
  4. It also supports splitting sequences into fragments, this
     could increase the specificity in profiling result in cost
     of slower searching speed.
  5. Multiple sizes of k-mers are supported.

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
     in a sequence file are concatenated with k-1 'N's before splitting.
  3. Both sequence IDs and fragments indices are saved for later use,
     in form of meta/description data in .unik files.

Meta data:
  1. Every outputted .unik file contains the sequence/reference ID,
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
  3. A summary file ("${outdir}/_info.txt") is generated for later use.
     Users need to check if the reference IDs (column "name") are what
     supposed to be.

Performance tips:
  1. Decrease value of -j/--threads for data in hard disk drives to
     reduce I/O pressure.

Usage:
  kmcp compute [flags]

Flags:
      --by-seq                    compute k-mers (sketches) for every sequence, instead of whole file
      --circular                  input sequence is circular
  -c, --compress                  output gzipped .unik files, it's slower and can saves little space
  -r, --file-regexp string        regular expression for matching files in -I/--in-dir to compute, case ignored (default "\\.(f[aq](st[aq])?|fna)(.gz)?$")
      --force                     overwrite output directory
  -h, --help                      help for compute
  -I, --in-dir string             directory containing FASTA/Q files. directory symlinks are followed
  -k, --kmer ints                 k-mer size(s) (default [21])
  -W, --minimizer-w int           minimizer window size
  -O, --out-dir string            output directory
  -N, --ref-name-regexp string    regular expression (must contains "(" and ")") for extracting reference name from file name (default "(?i)(.+)\\.(f[aq](st[aq])?|fna)(.gz)?$")
  -D, --scale int                 scale/down-sample factor (default 1)
  -B, --seq-name-filter strings   list of regular expressions for filtering out sequences by header/name, case ignored
  -m, --split-min-ref int         only splitting sequences >= M bp (default 1000)
  -n, --split-number int          fragment number, incompatible with -s/--split-size
  -l, --split-overlap int         fragment overlap for splitting sequences
  -s, --split-size int            fragment size for splitting sequences, incompatible with -n/--split-number
  -S, --syncmer-s int             closed syncmer length

```


## index

```text
Construct database from k-mer files

We build index for k-mers (sketches) with a modified compact bit-sliced
signature index (COBS). We totally rewrite the algorithms, data structure
and file format, and have improved the indexing and searching speed.

Attentions:
  1. All input .unik files should be generated by "kmcp compute".

Database size and searching accuracy:
  0. Use --dry-run to adjust parameters and check final number of 
     index files (#index-files) and the total file size.
  1. -f/--false-positive-rate: the default value 0.3 is enough for a
     query with tens of matched k-mers (see BIGSI/COBS paper).
     Small values could largely increase the size of database.
  2. -n/--num-hash: large values could reduce the database size,
     in cost of slower searching speed. Values <=4 is recommended.
  3. Value of block size -b/--block-size better be multiple of 64.
     The default value is:  (#unikFiles/#threads + 7) / 8 * 8
  4. Use flag -x/--block-sizeX-kmers-t, -8/--block-size8-kmers-t,
     and -1/--block-size1-kmers-t to separately create index for
     inputs with huge number of k-mers, for precise control of
     database size.

References:
  1. COBS: https://arxiv.org/abs/1905.09624

Taxonomy data:
  1. No taxonomy data are included in the database.
  2. Taxonomy information are only needed in "profile" command.
  
Performance tips:
  1. Number of blocks (.uniki files) better be smaller than or equal
     to number of CPU cores for faster searching speed. 
     We can set -j/--threads to control blocks number.
  2. #threads files are simultaneously opened, and max number
     of opened files is limited by the flag -F/--max-open-files.
     You may use a small value of -F/--max-open-files for 
     hard disk drive storage.
  3. When the database is used in a new computer with more CPU cores,
     'kmcp search' could automatically scale to utilize as many cores
     as possible.

Usage:
  kmcp index [flags]

Flags:
  -a, --alias string                 database alias/name, default: basename of --out-dir. you can also manually edit it in info file: ${outdir}/__db.yml
  -b, --block-size int               block size, better be multiple of 64 for large number of input files. default: min(#.files/#theads, 8)
  -1, --block-size1-kmers-t string   if k-mers of single .unik file exceeds this threshold, an individual index is created for this file. unit supported: K, M, G (default "200M")
  -8, --block-size8-kmers-t string   if k-mers of single .unik file exceeds this threshold, block size is changed to 8. unit supported: K, M, G (default "20M")
  -X, --block-sizeX int              if k-mers of single .unik file exceeds --block-sizeX-kmers-t, block size is changed to this value (default 256)
  -x, --block-sizeX-kmers-t string   if k-mers of single .unik file exceeds this threshold, block size is changed to --block-sizeX. unit supported: K, M, G (default "10M")
      --dry-run                      dry run, useful for adjusting parameters (recommended)
  -f, --false-positive-rate float    false positive rate of single bloom filter, range: (0, 1) (default 0.3)
      --file-regexp string           regular expression for matching files in -I/--in-dir to index, case ignored (default ".unik$")
      --force                        overwrite output directory
  -h, --help                         help for index
  -I, --in-dir string                directory containing .unik files. directory symlinks are followed
  -F, --max-open-files int           maximal number of opened files, please use a small value for hard disk drive storage (default 256)
  -n, --num-hash int                 number of hashes of bloom filters (default 1)
  -O, --out-dir string               output directory. default: ${indir}.kmcp-db

```

## search

```text
Search sequence against a database

Attentions:
  1. Input format should be (gzipped) FASTA or FASTQ from files or stdin.
     - Paired-end files should be given via -1/--read1 and -2/--read2.
        kmcp search -d db -1 read_1.fq.gz -2 read_2.fq.gz -o read.tsv.gz
     - Single-end can be given as positional arguments or -1/-2.
        kmcp search -d db file1.fq.gz file2.fq.gz -o result.tsv.gz
  2. A long query sequences may contain duplicated k-mers, which are
     not removed for short sequences by default. You may modify the
     value of -u/--kmer-dedup-threshold to remove duplicates.
  3. For long reads or contigs, you should split them in to short reads
     using "seqkit sliding", e.g.,
         seqkit sliding -s 100 -W 300

Shared flags between "search" and "profile":
  1. -t/--min-query-cov.
  2. -N/--name-map.

Special attentions:
  1. The values of tCov and jacc in results only apply to single size of k-mer.

Performance tips:
  1. Increase value of -j/--threads for acceleratation, but values larger
     than number of CPU cores won't bring extra speedup.
  2. Use --low-mem for database larger than RAM, but the searching would be
     very very slow for a large number of queries.

Usage:
  kmcp search [flags]

Flags:
  -d, --db-dir string              database directory created by "kmcp index"
  -D, --default-name-map           load ${db}/__name_mapping.tsv for mapping name first
  -S, --do-not-sort                do not sort matches of a query
  -h, --help                       help for search
  -n, --keep-top-scores int        keep matches with the top N score for a query, 0 for all
  -K, --keep-unmatched             keep unmatched query sequence information
  -u, --kmer-dedup-threshold int   remove duplicated kmers for a query with >= N k-mers (default 256)
      --low-mem                    do not load all index files into memory, the searching would be very very slow for a large number of queries
  -c, --min-kmers int              minimal number of matched k-mers (sketches) (default 30)
  -t, --min-query-cov float        minimal query coverage, i.e., proportion of matched k-mers and unique k-mers of a query (default 0.55)
  -m, --min-query-len int          minimal query length (default 70)
  -T, --min-target-cov float       minimal target coverage, i.e., proportion of matched k-mers and unique k-mers of a target
  -N, --name-map strings           tabular two-column file(s) mapping names to user-defined values
  -H, --no-header-row              do not print header row
  -o, --out-file string            out file, supports and recommends a ".gz" suffix ("-" for stdout) (default "-")
      --query-id string            custom query Id when using the whole file as a query
  -g, --query-whole-file           use the whole file as a query, e.g., for genome similarity estimation against k-mer sketch database
  -1, --read1 string               (gzipped) read1 file
  -2, --read2 string               (gzipped) read2 file
  -s, --sort-by string             sort hits by "qcov" (Containment Index), "tcov" or "jacc" (Jaccard Index) (default "qcov")
      --try-se                     if paired-end reads have no hits, re-search with read1, if still fails, try read2
  -G, --use-filename               use file name as query ID when using the whole file as a query

```

## merge

```text
Merge search results from multiple databases

Input:
  *. Searching results of the same reads on different databases.
  *. When only one input given, we just copy and write to the input file.
     This is friendly to workflows which assume multiple inputs are given.

Usage:
  kmcp merge [flags]

Flags:
  -n, --field-hits int       field of hits (default 5)
  -f, --field-queryIdx int   field of queryIdx (default 15)
  -h, --help                 help for merge
  -H, --no-header-row        do not print header row
  -o, --out-file string      out file, supports and recommends a ".gz" suffix ("-" for stdout) (default "-")
  -s, --sort-by string       sort hits by "qcov" (Containment Index), "tcov" or "jacc" (Jaccard Index) (default "qcov")

```

## profile


```text
Generate taxonomic profile from search results

Methods:
  1. Reference genomes can be splitted into fragments when computing
     k-mers (sketches), which could help to increase the specificity
     via a threshold, i.e., the minimal proportion of matched fragments
     (-p/--min-frags-prop). (***highly recommended***)
     Another flag -d/--max-frags-depth-stdev further reduces false positives.
  2. We require part of the uniquely matched reads of a reference
     having high similarity, i.e., with high confidence for decreasing
     the false positive rate.
  3. We also use the two-stage taxonomy assignment algorithm in MegaPath
     to reduce the false positive of ambiguous matches.
  4. Multi-aligned queries are proportionally assigned to references
     with a similar strategy in Metalign.
  5. Input files are parsed fours times, therefore STDIN is not supported.

Reference:
  1. MegaPath: https://doi.org/10.1186/s12864-020-06875-6
  2. Metalign: https://doi.org/10.1186/s13059-020-02159-0

Accuracy notes:
  *. Smaller -t/--min-qcov increase sensitivity in cost of higher false
     positive rate (-f/--max-fpr) of a query.
  *. We require part of the uniquely matched reads of a reference
     having high similarity, i.e., with high confidence for decreasing
     the false positive rate.
     E.g., -H >= 0.8 and -P >= 0.1 equals to 90th percentile >= 0.8
     *. -U/--min-hic-ureads,      minimal number, >= 1
     *. -H/--min-hic-ureads-qcov, minimal query coverage, >= -t/--min-qcov
     *. -P/--min-hic-ureads-prop, minimal proportion, higher values
        increase precision in cost of sensitivity.
  *. -R/--max-mismatch-err and -D/--min-dreads-prop is for determing
     the right reference for ambigous reads.
  *. --keep-perfect-match is not recommended, which decreases sensitivity. 
  *. -m/--keep-main-match is not recommended, which affects accuracy of
     abundance estimation.
  *. -n/--keep-top-qcovs  is not recommended, which affects accuracy of
     abundance estimation.

Taxonomy data:
  1. Mapping references IDs to TaxIds: -T/--taxid-map
  2. NCBI taxonomy dump files: -X/--taxdump

Performance notes:
  1. Searching results are parsed in parallel, and the number of
     lines proceeded by a thread can be set by the flag --chunk-size.
  2. However using a lot of threads does not always accelerate
     processing, 4 threads with chunk size of 500-5000 is fast enough.

Profiling output formats:
  1. KMCP      (-o/--out-prefix)
  2. CAMI      (-M/--metaphlan-report)
  3. MetaPhlAn (-C/--cami-report)

Taxonomic binning formats:
  1. CAMI      (-B/--binning-result)

Usage:
  kmcp profile [flags]

Flags:
  -B, --binning-result string         save extra binning result in CAMI report
  -C, --cami-report string            save extra CAMI-like report
      --chunk-size int                number of lines to process for each thread, and 4 threads is fast enough. Type "kmcp profile -h" for details (default 5000)
  -F, --filter-low-pct float          filter out predictions with the smallest relative abundances summing up N%. Range: [0,100)
  -h, --help                          help for profile
  -m, --keep-main-match               only keep main matches, abandon matches with sharply decreased qcov (> --max-qcov-gap)
      --keep-perfect-match            only keep the perfect matches (qcov == 1) if there are
  -n, --keep-top-qcovs int            keep matches with the top N qcovs for a query, 0 for all
      --level string                  level to estimate abundance at. available values: species, strain/assembly (default "species")
  -f, --max-fpr float                 maximal false positive rate of a read in search result (default 0.05)
  -d, --max-frags-depth-stdev float   maximal standard deviation of relative depths of all fragments (default 2)
  -R, --max-mismatch-err float        maximal error rate of a read being matched to a wrong reference, for determing the right reference for ambiguous reads. Range: (0, 1) (default 0.05)
      --max-qcov-gap float            max qcov gap between adjacent matches (default 0.2)
  -M, --metaphlan-report string       save extra metaphlan-like report
  -D, --min-dreads-prop float         minimal proportion of distinct reads, for determing the right reference for ambiguous reads. Range: (0, 1) (default 0.05)
  -p, --min-frags-prop float          minimal proportion of matched reference fragments with reads >= -r/--min-frags-reads (default 0.8)
  -r, --min-frags-reads int           minimal number of reads for a reference fragment (default 50)
  -U, --min-hic-ureads int            minimal number of high-confidence uniquely matched reads for a reference (default 1)
  -P, --min-hic-ureads-prop float     minimal proportion of high-confidence uniquely matched reads (default 0.1)
  -H, --min-hic-ureads-qcov float     minimal query coverage of high-confidence uniquely matched reads (default 0.75)
  -t, --min-query-cov float           minimal query coverage of a read in search result (default 0.55)
  -u, --min-uniq-reads int            minimal number of uniquely matched reads for a reference (default 10)
  -N, --name-map strings              tabular two-column file(s) mapping reference IDs to reference names
      --norm-abund string             method for normalize abundance of a reference by the mean/min/max abundance in all fragments, available values: mean, min, max (default "mean")
  -o, --out-prefix string             out file prefix ("-" for stdout) (default "-")
      --rank-prefix strings           prefixes of taxon name in certain ranks, used with --metaphlan-report  (default [k__,p__,c__,o__,f__,g__,s__,t__])
  -s, --sample-id string              sample ID in result file
  -S, --separator string              separator of TaxIds and taxonomy names (default ";")
      --show-rank strings             only show TaxIds and names of these ranks (default [superkingdom,phylum,class,order,family,genus,species,strain])
  -X, --taxdump string                directory of NCBI taxonomy dump files: names.dmp, nodes.dmp, optional with merged.dmp and delnodes.dmp
  -T, --taxid-map strings             tabular two-column file(s) mapping reference IDs to TaxIds

```

## utils

```text
Some utilities

Usage:
  kmcp utils [command]

Available Commands:
  filter        Filter search results and find species/assembly-specific queries
  index-info    Print information of index file
  merge-regions Merge species/assembly-specific regions
  unik-info     Print information of .unik file

```

## filter

```text
Filter search results and find species/assembly-specific queries

Taxonomy data:
  1. Mapping references IDs to TaxIds: -T/--taxid-map
  2. NCBI taxonomy dump files: -X/--taxdump

Performance notes:
  1. Searching results are parsed in parallel, and the number of
     lines proceeded by a thread can be set by the flag --chunk-size.
  2. However using a lot of threads does not always accelerate
     processing, 4 threads with chunk size of 500-5000 is fast enough.

Usage:
  kmcp utils filter [flags]

Flags:
      --chunk-size int        number of lines to process for each thread, and 4 threads is fast enough. Type "kmcp profile -h" for details (default 5000)
  -h, --help                  help for filter
      --level string          level to filter. available values: species, strain/assembly (default "species")
  -f, --max-fpr float         maximal false positive rate of a read in search result (default 0.05)
  -t, --min-query-cov float   minimal query coverage of a read in search result (default 0.55)
  -H, --no-header-row         do not print header row
  -o, --out-prefix string     out file prefix ("-" for stdout) (default "-")
  -X, --taxdump string        directory of NCBI taxonomy dump files: names.dmp, nodes.dmp, optional with merged.dmp and delnodes.dmp
  -T, --taxid-map strings     tabular two-column file(s) mapping reference IDs to TaxIds

```

## index-info

```text
Print information of index file

Usage:
  kmcp utils index-info [flags]

Flags:
  -a, --all                 all information
  -b, --basename            only output basenames of files
  -h, --help                help for index-info
  -o, --out-prefix string   out file prefix ("-" for stdout) (default "-")

```

## merge-regions

```text
Merge species/assembly-specific regions

Steps:
  # 1. Simulating reads and searching on one or more databases.
  seqkit sliding --step 10 --window 100 ref.fna.gz \
      | kmcp search -d db1.kmcp -o ref.fna.gz.kmcp@db1.tsv.gz
  seqkit sliding --step 10 --window 100 ref.fna.gz \
      | kmcp search -d db2.kmcp -o ref.fna.gz.kmcp@db2.tsv.gz
  
  # 2. Merging and filtering searching results
  kmcp merge ref.fna.gz.kmcp@*.tsv.gz \
      | kmcp utils filter -X taxdump -T taxid.map \
            -o ref.fna.gz.kmcp.uniq.tsv.gz
  
  # 3. Merging regions.
  # Here the value of --min-overlap should be k-1.
  kmcp utils merge-regions --min-overlap 20 ref.fna.gz.kmcp.uniq.tsv.gz \
      -o ref.fna.gz.kmcp.uniq.tsv.gz.bed

Output (BED6 format):
  1. chrom      - chromosome name
  2. chromStart - starting position (0-based)
  3. chromEnd   - ending position (0-based)
  4. name       - "species-specific" or "assembly-specific"
  5. score      - 0-1000, 1000 for "assembly-specific", others for ""species-specific"
  6. strand     - "."

Performance notes:
  1. Searching results are parsed in parallel, and the number of
     lines proceeded by a thread can be set by the flag --chunk-size.
  2. However using a lot of threads does not always accelerate
     processing, 4 threads with chunk size of 500-5000 is fast enough.

Usage:
  kmcp utils merge-regions [flags]

Flags:
      --chunk-size int         number of lines to process for each thread, and 4 threads is fast enough. Type "kmcp profile -h" for details (default 5000)
  -h, --help                   help for merge-regions
  -I, --ignore-type            merge species and assembly-specific regions
  -f, --max-fpr float          maximal false positive rate of a read in search result (default 0.05)
  -g, --max-gap int            maximal distance of starting positions of two adjacent regions, 0 for no limitation, 1 for no merging
  -l, --min-overlap int        minimal overlap of two adjacent regions, recommend K-1 (default 1)
  -t, --min-query-cov float    minimal query coverage of a read in search result (default 0.55)
  -a, --name-assembly string   name of assembly-specific regions (default "assembly-specific")
  -s, --name-species string    name of species-specific regions (default "species-specific")
  -o, --out-prefix string      out file prefix ("-" for stdout) (default "-")
  -r, --regexp string          regular expression for extract reference name and query locations (default "^(.+)_sliding:(\\d+)\\-(\\d+)$")
```

## unik-info

```text
Print information of .unik file

Tips:
  1. For lots of small files (especially on SDD), use big value of '-j' to
     parallelize counting.

Usage:
  kmcp utils unik-info [flags]

Flags:
  -a, --all                   all information, including number of k-mers
  -b, --basename              only output basename of files
  -h, --help                  help for unik-info
  -o, --out-file string       out file ("-" for stdout, suffix .gz for gzipped out) (default "-")
  -e, --skip-err              skip error, only show warning message
      --symbol-false string   smybol for false (default "✕")
      --symbol-true string    smybol for true (default "✓")
  -T, --tabular               output in machine-friendly tabular format

```


## autocomplete

```text
Generate shell autocompletion script

Supported shell: bash|zsh|fish|powershell

Bash:

    # generate completion shell
    kmcp autocompletion --shell bash

    # configure if never did.
    # install bash-completion if the "complete" command is not found.
    echo "for bcfile in ~/.bash_completion.d/* ; do source \$bcfile; done" >> ~/.bash_completion
    echo "source ~/.bash_completion" >> ~/.bashrc

Zsh:

    # generate completion shell
    kmcp autocompletion --shell zsh --file ~/.zfunc/_kmcp

    # configure if never did
    echo 'fpath=( ~/.zfunc "${fpath[@]}" )' >> ~/.zshrc
    echo "autoload -U compinit; compinit" >> ~/.zshrc

fish:

    kmcp autocompletion --shell fish --file ~/.config/fish/completions/kmcp.fish

Usage:
  kmcp autocompletion [flags]

Flags:
      --file string    autocompletion file (default "/home/shenwei/.bash_completion.d/kmcp.sh")
  -h, --help           help for autocompletion
      --shell string   autocompletion type (bash|zsh|fish|powershell) (default "bash")

```
