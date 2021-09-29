# Database

## Prebuilt databases

### A). Databases for metagenomic profiling

These databases are created following [steps below](#building-databases).
Users can also [build custom databases](#building-custom-databases), it's simple and fast.

|kingdoms                |source     |# species|# assembly|file                                    |size     |
|:-----------------------|:----------|:--------|:---------|:---------------------------------------|:--------|
|**Bacteria and Archaea**|GTDB r202  |43252    |47894     |[gtdb.kmcp.tar.gz]() (37.13 GB)         |55.12 GB |
|**Viruses**             |Refseq r207|7300     |11618     |[refseq-viruses.kmcp.tar.gz]() (1.09 GB)|4.14 GB  |
|**Fungi**               |Refseq r207|148      |390       |[refseq-fungi.kmcp.tar.gz]() (8.79GB)   |11.12 GB |


**Taxonomy data**: [taxdump.tar.gz]().

**Hardware requirements**:

- Prebuilt databases above were built for computers with >= 32 CPU cores
in consideration of better parallelization,
and computers should have at least 64 GB (> 55.12GB). 
- By default, `kmcp search` loads the whole database into main memory (RAM) for fast searching.
Optionally, the flag `--low-mem` can be set to avoid loading the whole database,
while it's much slower, >10X slower on SSD and should be much slower on HDD disks.
- To reduce memory requirements on computers without enough memory,
users can divide the reference genomes into several parts
and build smaller databases for all parts, so that the biggest
database can be loaded into RAM. After performing database searching,
search results on all small databases can be merged with `kmcp merge`
for downstream analysis.





### B). Databases for genome similarity estimation

Check [tutorial](/tutorial/searching).

|kingdoms                |source     |sketch         |parameters          |file                                           |size     |
|:-----------------------|:----------|:--------------|:-------------------|:----------------------------------------------|:--------|
|**Bacteria and Archaea**|GTDB r202  |Scaled MinHash |k=31, scale=1000    |[gtdb.minhash.kmcp.tar.gz]() (710 MB)          |1.52 GB  |
|**Bacteria and Archaea**|GTDB r202  |Closed Syncmers|k=31, s=15, scale=60|[gtdb.syncmer.kmcp.tar.gz]() (1.04 GB)         |2.28 GB  |
|**Viruses**             |Refseq r207|Scaled MinHash |K=31, scale=10      |[refseq-viruses.minhash.kmcp.tar.gz]() (181 MB)|407 MB   |
|**Viruses**             |Refseq r207|Closed Syncmers|k=31, s=21          |[refseq-viruses.syncmer.kmcp.tar.gz]() (142 MB)|312 MB   |

## Building databases

### GTDB

Tools:

- [brename](https://github.com/shenwei356/brename/releases) for batching renaming files.
- [rush](https://github.com/shenwei356/rush/releases) for executing jobs in parallel.
- [seqkit](https://github.com/shenwei356/seqkit/releases) for FASTA file processing.
- [kmcp](/download) for metagenomic profiling.

Files:

 - [gtdb_genomes_reps_r202.tar.gz](https://data.ace.uq.edu.au/public/gtdb/data/releases/release202/202.0/genomic_files_reps/gtdb_genomes_reps_r202.tar.gz)
 - [ar122_metadata_r202.tar.gz](https://data.ace.uq.edu.au/public/gtdb/data/releases/release202/202.0/ar122_metadata_r202.tar.gz)
 - [bac120_metadata_r202.tar.gz](https://data.ace.uq.edu.au/public/gtdb/data/releases/release202/202.0/bac120_metadata_r202.tar.gz)

Uncompressing and renaming:
 
    # uncompress
    mkdir -p gtdb
    tar -zxvf gtdb_genomes_reps_r202.tar.gz -O gtdb
    
    # rename
    brename -R -p '^(\w{3}_\d{9}\.\d+).+' -r '$1.fna.gz' gtdb
  
Mapping file:

    tar -zxvf ar122_metadata_r202.tar.gz  bac120_metadata_r202.tar.gz
    
    # assembly accesion -> full head
    find gtdb/ -name "*.fna.gz" \
        | rush -k 'echo -ne "{%@(.+).fna}\t$(seqkit head -n 1 {} | seqkit seq -n)\n" ' \
        > name.map
        
    # assembly accesion -> taxid
    (cat ar122_metadata_r202.tsv; sed 1d bac120_metadata_r202.tsv) \
        | csvtk cut -t -f accession,ncbi_taxid \
        | csvtk replace -t -p '^.._' \
        | csvtk grep -t -P <(cut -f 1 name.map) \
        | csvtk del-header \
        > taxid.map

    
    # stats (optional)

    # number of species/strains
    cat taxid.map \
        | taxonkit lineage -i 2 -r -n -L \
        | csvtk freq -Ht -f 4 -nr \
        | csvtk pretty -H -t \
        | tee taxid.map.stats
    species           43566
    strain            4108
    subspecies        111
    forma specialis   58
    no rank           26
    isolate           24
    serotype          1
    
    # number of unique species/strains
    cat taxid.map \
        | csvtk uniq -Ht -f 2 \
        | taxonkit lineage -i 2 -r -n -L \
        | csvtk freq -Ht -f 4 -nr \
        | csvtk pretty -H -t \
        | tee taxid.map.uniq.stats
    species           24739
    strain            4101
    subspecies        89
    forma specialis   58
    no rank           26
    isolate           24
    serotype          1
  
        
Building database:

    input=gtdb
    
    # compute k-mers
    #   sequence containing "plasmid" in name are ignored,
    #   reference genomes are splitted into 10 fragments with 100bp overlap
    #   k = 21
    kmcp compute -I $input -O gtdb-r202-k21-n10 -k 21 -n 10 -l 100 -B plasmid \
        --log gtdb-r202-k21-n10.log --force

    # build database
    #   number of index files: 32, for server with >= 32 CPU cores
    #   bloom filter parameter:
    #     number of hash function: 1
    #     false positive rate: 0.3
    kmcp index -j 32 -I gtdb-r202-k21-n10 -O gtdb.kmcp -n 1 -f 0.3 \
        --log gtdb.kmcp.log
    
    # cp taxid and name mapping file to database directory
    cp taxid.map name.map gtdb.kmcp/


### RefSeq

Tools

- [genome_updater](https://github.com/pirovc/genome_updater) for downloading genomes from NCBI.

Downloading viral and fungi sequences:

    # name=fungi
    name=viral
    
    # -k for dry-run
    # -i for fix
    time genome_updater.sh \
        -d "refseq"\
        -g $name \
        -c "all" \
        -l "all" \
        -f "genomic.fna.gz" \
        -o "refseq-$name" \
        -t 12 \
        -m -a -p

    # cd to 2021-07-30_21-54-19
        
    # taxdump
    mkdir -p taxdump
    tar -zxvf taxdump.tar.gz -C taxdump
    
    # assembly accesion -> taxid
    cut -f 1,6 assembly_summary.txt > taxid.map    
    # assembly accesion -> name
    cut -f 1,8 assembly_summary.txt > name.map
    
    
    # stats (optional)
    cat taxid.map  \
        | csvtk freq -Ht -f 2 -nr \
        | taxonkit lineage -r -n -L --data-dir taxdump/ \
        | taxonkit reformat -I 1 -f '{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}' --data-dir taxdump/ \
        | csvtk add-header -t -n 'taxid,count,name,rank,superkindom,phylum,class,order,family,genus,species' \
        > taxid.map.stats.tsv
    
    seqkit stats -T -j 12 --infile-list <(find files -name "*.fna.gz") > files.stats.tsv
    
    
    
    # create another directory linking to genome files
    
    input=files.renamed    
    
    mkdir -p $input
    # create soft links
    cd $input
    find ../files -name "*.fna.gz" \
        | rush 'ln -s {}'
    cd ..    
    # rename
    brename -R -p '^(\w{3}_\d{9}\.\d+).+' -r '$1.fna.gz' $input
        
Building database:
        
    # -----------------------------------------------------------------
    # for viral
    name=viral
    
    input=files.renamed
    
    kmcp compute -I $input -O refseq-$name-k21-n10 \
        -k 21 --seq-name-filter plasmid \
        --split-number 10 --split-overlap 100 \
        --log refseq-$name-k21-n10.log --force
    
    # viral genomes are small:
    #   using small false positive rate: 0.001
    #   using more hash functions: 3
    kmcp index -I refseq-$name-k21-n10/ -O refseq-viruses.kmcp \
        -j 32 -f 0.001 -n 3 \
        --log refseq-viruses.kmcp.log --force
    
    # cp taxid and name mapping file to database directory
    cp taxid.map name.map refseq-viruses.kmcp/

    # -----------------------------------------------------------------
    # for fungi
    name=fungi
    
    input=files.renamed
    
    kmcp compute -I $input -O refseq-$name-k21-n10 \
        -k 21 --seq-name-filter plasmid \
        --split-number 10 --split-overlap 100 \
        --log refseq-$name-k21-n10.log --force
      
    kmcp index -I refseq-$name-k21-n10/ -O refseq-fungi.kmcp \
        -j 32 -f 0.05 -n 2 \
        --log refseq-fungi.kmcp.log --force
    
    # cp taxid and name mapping file to database directory
    cp taxid.map name.map refseq-fungi.kmcp/

### HumGut


[HumGut](http://arken.nmbu.no/~larssn/humgut/) is a comprehensive Human Gut prokaryotic genomes collection filtered by metagenome data.

Dataset

- Genomes: [HumGut.tar.gz](http://arken.nmbu.no/~larssn/humgut/HumGut.tar.gz)
- Metadata: [HumGut.tsv](http://arken.nmbu.no/~larssn/humgut/HumGut.tsv)
- Taxdump files:
    - NCBI taxonomy: [ncbi_names.dump](http://arken.nmbu.no/~larssn/humgut/ncbi_names.dmp)
      and [ncbi_nodes.dump](http://arken.nmbu.no/~larssn/humgut/ncbi_names.dmp)
    - GTDB taxomomy: [gtdb_names.dump](http://arken.nmbu.no/~larssn/humgut/gtdb_names.dmp)
      and [gtdb_nodes.dump](http://arken.nmbu.no/~larssn/humgut/gtdb_names.dmp)
    
Uncompressing and renaming:
 
    # uncompress
    tar -zxvf HumGut.tar.gz
    
    # taxdump files
    mkdir taxdump
    mv ncbi_names.dmp taxdump/names.dmp
    mv ncbi_nodes.dmp taxdump/nodes.dmp
  
Mapping file:
    
    # assembly accesion -> taxid
    cat HumGut.tsv \
        | csvtk cut -t -f genome_file,ncbi_tax_id \
        | csvtk replace -f genome_file -t -p '\..+$' \
        | csvtk del-header \
        > taxid.map
        
    # assembly accesion -> name
    cat HumGut.tsv \
        | csvtk cut -t -f genome_file,ncbi_organism_name \
        | csvtk replace -f genome_file -t -p '\..+$' \
        | csvtk del-header \
        > name.map
    
Stats:

    # NCBI taxonomy

    # number of species/strains
    cat taxid.map \
        | taxonkit lineage --data-dir taxdump-ncbi -i 2 -r -n -L \
        | csvtk freq -Ht -f 4 -nr \
        | csvtk pretty -H -t \
        | tee taxid.map.stats
    species      23604
    strain       6816
    subspecies   161
    no rank      69
    serotype     38
    genus        3
    
    # number of unique species/strains
    cat taxid.map \
        | csvtk uniq -Ht -f 2 \
        | taxonkit lineage --data-dir taxdump-ncbi -i 2 -r -n -L \
        | csvtk freq -Ht -f 4 -nr \
        | csvtk pretty -H -t \
        | tee taxid.map.uniq.stats;
    species      1240
    strain       458
    subspecies   15
    genus        1
    no rank      1
    serotype     1
    
    # GTDB taxonomy
    
    cat HumGut.tsv \
        | csvtk cut -t -f genome_file,gtdbtk_tax_id \
        | csvtk replace -f genome_file -t -p '\..+$' \
        | csvtk del-header \
        > taxid-gtdb.map
        
    cat taxid-gtdb.map \
        | csvtk uniq -Ht -f 2 \
        | taxonkit lineage --data-dir taxdump-gtdb -i 2 -r -n -L \
        | csvtk freq -Ht -f 4 -nr \
        | csvtk pretty -H -t \
        | tee taxid-gtdb.map.uniq.stats
    species   2810
    genus     434
    family    52
    order     12
    class     2
        
Building database:
    
    # compute k-mers
    #   sequence containing "plasmid" in name are ignored,
    #   reference genomes are splitted into 10 fragments
    #   k = 21
    kmcp compute -I humgut/ -k 21 -n 10 -B plasmid -O humgut-k21-n10 --force

    # build database
    #   number of index files: 32, for server with >= 32 CPU cores
    #   bloom filter parameter:
    #     number of hash function: 1
    #     false positive rate: 0.3
    kmcp index -j 32 -I humgut-k21-n10 -O humgut.kmcp -n 1 -f 0.3


## Building custom databases

Files:

1. Genome files
    - (Gzip-compressed) FASTA/Q format.
    - One genome per file with the reference identifier in the file name.
2. TaxId mapping file (for metagenomic profiling)
    - Two-column (reference identifier and TaxId) tab-delimited.
3. [NCBI taxonomy dump files](ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz) (for metagenomic profiling)
    - `names.dmp`
    - `nodes.dmp`
    - `merged.dmp` (optional)
    - `delnodes.dmp` (optional)

Tools:

- [kmcp](/download) for metagenomic profiling.
- [rush](https://github.com/shenwei356/rush/releases) for executing jobs in parallel.
- [brename](https://github.com/shenwei356/brename/releases) for batching renaming files (optional)


**Memory notes**:

By default, `kmcp search` loads the whole database into main memory (RAM) for fast searching.
Optionally, the flag `--low-mem` can be set to avoid loading the whole database,
while it's much slower, >10X slower on SSD and should be much slower on HDD disks.

Another way is dividing the reference genomes into several parts
and building smaller databases for all parts, so that the biggest
database can be loaded into RAM. After performing database searching,
search results on all small databases can be merged with `kmcp merge`
for downstream analysis.


### Step 1. Computing k-mers

Plain or gzip-compressed input genome sequence files better be saved in one directory,
with multiple-level directories allowed. 
The file extension could be `fa`, `fasta`, `fna`, `fq`, `fastq`, `fa.gz`, `fasta.gz`,
 `fna.gz`, `fq.gz` or `fastq.gz`.

**Input files** can be given as a list of FASTA/Q files via
positional arguments or a directory (with multiple-level directories allowed) containing sequence files
via the flag `-I/--in-dir`. A regular expression for matching
sequencing files is available by the flag `-r/--file-regexp`.
The default pattern matches files with extension of 
`fa`, `fasta`, `fna`, `fq`, `fastq`, `fa.gz`, `fasta.gz`,
`fna.gz`, `fq.gz` or `fastq.gz`.

Unwanted sequence like plasmid sequences can be **filtered out** by
the name via regular expression(s) (`-B/--seq-name-filter`).

**How to compute k-mers**:
By default, `kmcp` computes k-mers (sketches) of every file,
you can also use `--by-seq` to compute for every sequence,
where sequence IDs in all input files better be distinct.
It also **supports splitting sequences into fragments, this
could increase the specificity in profiling result in cost
of slower searching speed**. 

**Splitting sequences**:

1. Sequences can be splitted into fragments by a fragment size 
    (`-s/--split-size`) or number of fragments (`-n/--split-number`)
    with overlap (-`l/--split-overlap`).
2. When splitting by number of fragments, **all sequences (except for
    these mathching any regular expression given by `-B/--seq-name-filter`)
    in a sequence file are concatenated before splitting**.
3. Both sequence/reference IDs and fragments indices are saved for later use,
    in form of meta/description data in `.unik` files, and will
    be reported in `kmcp search` results.

**Meta data**:

1. Every outputted `.unik` file contains the sequence/reference ID,
    fragment index, number of fragments, and genome size of reference.
2. When parsing whole sequence files or splitting by number of fragments,
    the **identifier of a reference** is the basename of the input file
    by default. It can also be extracted from the input file name via
    `-N/--ref-name-regexp`, e.g., `^(\w{3}_\d{9}\.\d+)` for **RefSeq assembly accessions**.

Multiple **sizes of k-mers** are supported, but a single k-mer size
is good enough. 

**Supported k-mer (sketches) types**:

1. K-mer:
    - ntHash of k-mer (`-k`)
2. K-mer sketchs (all using ntHash):
    - Scaled MinHash (`-k -D`)
    - Minimizer      (`-k -W`), optionally scaling/down-sampling (`-D`)
    - Closed Syncmer (`-k -S`), optionally scaling/down-sampling (`-D`)

**Output**:

1. All outputted `.unik` files are saved in `${outdir}`, with path

        ${outdir}/xxx/yyy/zzz/${infile}-id_${seqID}.unik

    where dirctory tree `/xxx/yyy/zzz/` is built for > 1000 output files.

2. For splitting sequence mode (`--split-size > 0` or `--split-number > 0`),
    output files are:
    
        ${outdir}//xxx/yyy/zzz/${infile}/{seqID}-frag_${fragIdx}.unik

**Performance tips**:

1. Decrease value of `-j/--threads` for data in hard disk drives (HDD) to
    reduce I/O pressure.

**Commands** ([usage](/usage/#compute)):

    # compute k-mers
    #   sequence containing "plasmid" in name are ignored,
    #   reference genomes are splitted into 10 fragments,
    #   k = 21
    kmcp compute --in-dir sequences-dir/ \
        --kmer 21 \
        --split-number 10 \
        --seq-name-filter plasmid \
        --out-dir refs-k21-n10

    # demo output
    16:01:56.269 [INFO] kmcp v0.6.0
    16:01:56.269 [INFO]   https://github.com/shenwei356/kmcp
    16:01:56.269 [INFO] 
    16:01:56.269 [INFO] checking input files ...
    16:01:56.269 [INFO]   9 input file(s) given
    16:01:56.269 [INFO] 
    16:01:56.269 [INFO] -------------------- [main parameters] --------------------
    16:01:56.269 [INFO] input and output:
    16:01:56.269 [INFO]   input directory: refs
    16:01:56.269 [INFO]     regular expression of input files: (?i)\.(f[aq](st[aq])?|fna)(.gz)?$
    16:01:56.269 [INFO]     *regular expression for extracting reference name from file name: 
    16:01:56.269 [INFO]     *regular expressions for filtering out sequences: [plasmid]
    16:01:56.269 [INFO]   output directory: refs-k21-n10
    16:01:56.269 [INFO] 
    16:01:56.269 [INFO] sequences splitting: true
    16:01:56.269 [INFO]   split parts: 10, overlap: 0 bp
    16:01:56.269 [INFO] 
    16:01:56.269 [INFO] k-mer (sketches) computing:
    16:01:56.269 [INFO]   k-mer size(s): 21
    16:01:56.269 [INFO]   circular genome: false
    16:01:56.269 [INFO]   saving exact number of k-mers: true
    16:01:56.269 [INFO] 
    16:01:56.269 [INFO] -------------------- [main parameters] --------------------
    16:01:56.269 [INFO] 
    16:01:56.269 [INFO] computing ...
    processed files:  9 / 9 [======================================] ETA: 0s. done
    16:01:56.706 [INFO] 
    16:01:56.706 [INFO] elapsed time: 437.515931ms
    16:01:56.706 [INFO]

A summary file (`_info.txt`) is generated for later use.
**Users needs check if the reference IDs (column `name`) are what supposed to be**.

|#path                                                              |name     |fragIdx|idxNum|genomeSize|kmers |
|:------------------------------------------------------------------|:--------|:------|:-----|:---------|:-----|
|refs-k21-n10/672/060/672/NC_010655.1.fasta.gz/NC_010655-frag_0.unik|NC_010655|0      |10    |2664102   |264229|
|refs-k21-n10/292/039/292/NC_000913.3.fasta.gz/NC_000913-frag_0.unik|NC_000913|0      |10    |4641652   |459259|
|refs-k21-n10/414/159/414/NC_011750.1.fasta.gz/NC_011750-frag_0.unik|NC_011750|0      |10    |5132068   |505905|
|refs-k21-n10/672/060/672/NC_010655.1.fasta.gz/NC_010655-frag_1.unik|NC_010655|1      |10    |2664102   |266219|
|refs-k21-n10/744/528/744/NC_018658.1.fasta.gz/NC_018658-frag_0.unik|NC_018658|0      |10    |5273097   |523691|

Meta data in the `.unik` file can be showed using [unikmer](https://github.com/shenwei356/unikmer/releases):

    unikmer info refs-k21-n10/072/380/072/NZ_CP028116.1.fasta.gz/NZ_CP028116-frag_0.unik -a


### Step 2. Building databases

KMCP builds index for k-mers (sketches) with a modified Compact Bit-sliced
Signature Index ([COBS](https://arxiv.org/abs/1905.09624)). 
We totally rewrite the algorithms, data structure and file format,
and have improved the indexing and searching speed
(check [benchmark](/benchmark)).

**Input**:

- The output directory generated by `kmcp compute`.

**Database size and searching accuracy**:

0. **Use `--dry-run` to adjust parameters and check final number of 
    index files (#index-files) and the total file size**.
1. `-f/--false-positive-rate`: the default value `0.3` is enough for a
    query with tens of matched k-mers (see BIGSI/COBS paper).
    Small values could largely increase the size of database.
2. `-n/--num-hash`: large values could reduce the database size,
    in cost of slower searching speed. Values <=4 is recommended.
3. Value of block size `-b/--block-size` better be multiple of 64.
    The default value is: 
    
        (#unikFiles/#threads + 7) / 8 * 8

4. Use flag `-x/--block-sizeX-kmers-t`, `-8/--block-size8-kmers-t`,
    and `-1/--block-size1-kmers-t` to separately create index for
    inputs with huge number of k-mers, for precise control of
    database size.

**Taxonomy data**:

1. No taxonomy data are included in the database.
2. Taxonomy information are only needed in `kmcp profile`.

**Performance tips**:

1. Number of blocks (.uniki files) better be smaller than or equal
   to number of CPU cores for faster searching speed. 
   **We can set `-j/--threads` to control blocks number**.
2. `#threads` files are simultaneously opened, and the max number
    of opened files is limited by the flag `-F/--max-open-files`.
    You may use a small value of `-F/--max-open-files` for 
    hard disk drive storage.
3. When the database is used in a new computer with more CPU cores,
   `kmcp search` could automatically scale to utilize as many cores
   as possible.

**Commands** ([usage](/usage/#compute)):

    # build database
    #   number of index files: 32, for server with >= 32 CPU cores
    #   bloom filter parameter:
    #     number of hash function: 1
    #     false positive rate: 0.3
    kmcp index -I refs-k21-n10 \
        --threads 32 \
        --num-hash 1 \
        --false-positive-rate 0.3 \
        --out-dir refs.kmcp 

    # demo output
    16:55:06.758 [INFO] kmcp v0.6.0
    16:55:06.758 [INFO]   https://github.com/shenwei356/kmcp
    16:55:06.758 [INFO] 
    16:55:06.758 [INFO] loading .unik file infos from file: refs-k21-n10/_info.txt
    16:55:06.767 [INFO]   90 cached file infos loaded
    16:55:06.767 [INFO] 
    16:55:06.767 [INFO] -------------------- [main parameters] --------------------
    16:55:06.767 [INFO]   number of hashes: 1
    16:55:06.767 [INFO]   false positive rate: 0.300000
    16:55:06.767 [INFO]   k-mer size(s): 21
    16:55:06.767 [INFO]   split seqequence size: 0, overlap: 0
    16:55:06.767 [INFO]   block-sizeX-kmers-t: 10.00 M
    16:55:06.767 [INFO]   block-sizeX        : 256.00
    16:55:06.767 [INFO]   block-size8-kmers-t: 20.00 M
    16:55:06.767 [INFO]   block-size1-kmers-t: 200.00 M
    16:55:06.767 [INFO] -------------------- [main parameters] --------------------
    16:55:06.767 [INFO] 
    16:55:06.767 [INFO] building index ...
    16:55:06.778 [WARN] ignore -X/--block-size (256) which is >= -b/--block-size (8)
    16:55:06.778 [INFO] 
    16:55:06.778 [INFO]   block size: 8
    16:55:06.778 [INFO]   number of index files: 32 (may be more)
    16:55:06.778 [INFO] 
    16:55:06.778 [block #001]   1 / 1  100 %
    16:55:06.778 [block #002]   1 / 1  100 %
    16:55:06.778 [block #003]   1 / 1  100 %
    16:55:06.778 [block #004]   1 / 1  100 %
    16:55:06.778 [block #005]   1 / 1  100 %
    16:55:06.778 [block #006]   1 / 1  100 %
    16:55:06.778 [block #007]   1 / 1  100 %
    16:55:06.778 [block #008]   1 / 1  100 %
    16:55:06.779 [block #009]   1 / 1  100 %
    16:55:06.779 [block #010]   1 / 1  100 %
    16:55:06.779 [block #011]   1 / 1  100 %
    16:55:06.779 [block #012]   1 / 1  100 %
    [saved index files] 12 / 12  ETA: 0s. done
    16:55:06.979 [INFO] 
    16:55:06.979 [INFO] kmcp database with 42713234 k-mers saved to refs.kmcp
    16:55:06.979 [INFO] total file size: 15.66 MB
    16:55:06.979 [INFO] total index files: 12
    16:55:06.979 [INFO] 
    16:55:06.979 [INFO] elapsed time: 221.79845ms
    16:55:06.979 [INFO]

Output:

    refs.kmcp/
    └── R001
        ├── _block001.uniki
        ├── _block002.uniki
        ├── _block003.uniki
        ├── _block004.uniki
        ├── _block005.uniki
        ├── _block006.uniki
        ├── _block007.uniki
        ├── _block008.uniki
        ├── _block009.uniki
        ├── _block010.uniki
        ├── _block011.uniki
        ├── _block012.uniki
        ├── __db.yml
        └── __name_mapping.tsv

`__db.yml` contains configuration of the database, and `_blockXXX.uniki` are index files.
`kmcp info` could show the basic information of index files:

    kmcp info  refs.kmcp/R001/_block001.uniki

|file                          |k  |canonical|num-hashes|num-sigs|num-names|
|:-----------------------------|:--|:--------|:---------|:-------|:--------|
|refs.kmcp/R001/_block001.uniki|21 |true     |1         |746392  |8        |

What's next? Check the [tutorials](/tutorial).
