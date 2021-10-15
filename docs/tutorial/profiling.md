# Metagenomic profiling

## Requirements

- Database
    - [Prebuilt databases](/database/#prebuilt-databases) are available.
    - Or [build custom databases](/database/#custom-database).
- Hardware.
    - CPU: ≥ 32 cores preferred.
    - RAM: ≥ 64 GB, depends on file size of the maximum database.

## Datasets

- Short reads, single or paired end, but **paired information is not used** in current version.
- Long reads (PacBio HIFI).

## Steps

### Step 1. Preprocessing reads

For example, removing adapters and trimming using [fastp](https://github.com/OpenGene/fastp):

    fastp -i in_1.fq.gz -I in_2.fq.gz  \
        -o out_1.fq.gz -O out_2.fq.gz \
        -l 75 -q 20 -W 4 -M 20 -3 20 --thread 32 \
        --html out.fastp.html

### Step 2. Removing host reads

Tools:

- [bowtie2](https://github.com/BenLangmead/bowtie2) is [recommended](https://doi.org/10.1099/mgen.0.000393) for removing host reads.
- [samtools](https://github.com/samtools/samtools) is also used for processing reads mapping file.
- [pigz](https://zlib.net/pigz/) is a parallel implementation of `gzip`, which is much faster than `gzip`.

Host reference genomes:

- Human: [CHM13](https://github.com/marbl/CHM13)
- Mouse: []()
- Rat: []()

Building the index:
    
    bowtie2-build --threads 32 chm13.draft_v1.1.fasta.gz chm13

Mapping and removing mapped reads:

    index=~/ws/db/bowtie2/chm13
    
    bowtie2 --threads 32 -x $index -1 in_1.fq.gz -2 in_2.fq.gz  \
        | samtools view -buS -f 4 - \
        | samtools fastq - \
        | gzip -c > sample.fq.gz

### Step 3. Searching

**Reads can be searched on different databases and then merged**.
Databases can be built with different parameters.

**Attentions**:

1. Input format should be (gzipped) FASTA or FASTQ from files or stdin.
2. A long query sequences may contain duplicated k-mers, which are
    not removed for short sequences by default. You may modify the
    value of `-u/--kmer-dedup-threshold` (default `256`) to remove duplicates.
3. For long reads or contigs, you should split them in to short reads
    using "seqkit sliding", e.g.,
        `seqkit sliding -s 100 -W 300`
4. The values of `tCov` and `jacc` in result only apply for single size of k-mer.

**`kmcp search` and `kmcp profile` share some flags**, therefore users
can use stricter criteria in `kmcp profile`.

1. `-t/--min-query-cov`, minimal query coverage, i.e., 
   proportion of matched k-mers and unique k-mers of a query (default `0.6`)
2. `-n/--keep-top-scores`, keep matches with the top N score for a query, 0 for all (default `5`),
   here it can reduce the output size, while     it does not effect the speed.
3. `-N/--name-map`, tabular two-column file(s) mapping names to user-defined values.

**Performance tips**:

1. Increase value of `-j/--threads` for acceleratation, but values larger
  than the number of CPU cores won't bring extra speedup.

**Commands**:

    # paired-ends: 
    #   file="sample_1.fq.gz sample_2.fq.gz"
    file=sample.fq.gz

    for db in refseq-viruses.kmcp gtdb.kmcp ; do
        dbname=$(basename $db)

        kmcp search \
            --threads            32 \
            --db-dir            $db \
            --min-kmers          30 \
            --min-query-len      70 \
            --min-query-cov    0.55 \
            --keep-top-scores     0 \
            $file                   \
            --out-file $file.kmcp@$dbname.tsv.gz \
            --log      $file.kmcp@$dbname.tsv.gz.log
    done

Merging searching results on multiple database:

    kmcp merge $file.kmcp@*.tsv.gz --out-file $file.kmcp.tsv.gz

### Step 4. Profiling

**Input**:

- TaxId mapping file(s).
- Taxdump files.
- KMCP search results.

**Methods**

1. Reference genomes can be splitted into fragments when computing
    k-mers (sketches), which could help to increase the specificity
    via a threshold, i.e., the minimal proportion of matched fragments
    (`-p`/--min-frags-prop). (***highly recommended***)
    Another flag `-d/--max-frags-cov-stdev` further reduces false positives.
2. We require part of the uniquely matched reads of a reference
    having high similarity, i.e., with high confidence for decreasing
    the false positive rate.
3. We also use the two-stage taxonomy assignment algorithm in [MegaPath](https://doi.org/10.1186/s12864-020-06875-6)
    to reduce the false positive of ambiguous matches.
4. Multi-aligned queries are proportionally assigned to references
    with a similar strategy in [Metalign](https://doi.org/10.1186/s13059-020-02159-0).
5. Input files are parsed 4 times, therefore STDIN is not supported.

<img src="/tutorial/profiling-steps.png" alt="" width="500"/>

**Accuracy notes**:

- Smaller `-t/--min-qcov` increase sensitivity in cost of higher false
    positive rate (`-f/--max-fpr`) of a query.
- And we require part of the uniquely matched reads of a reference
    having high similarity, i.e., with high confidence to decrease
    the false positive.
    E.g., `-H >= 0.8` and `-P >= 0.1` equals to `90th percentile >= 0.8`
    - `-U/--min-hic-ureads`,      minimal number, `>= 1`
    - `-H/--min-hic-ureads-qcov`, minimal query coverage, `>= -t/--min-qcov`
    - `-P/--min-hic-ureads-prop`, minimal proportion, higher values
    increase precision in cost of sensitivity.
- `-R/--max-mismatch-err` and `-D/--min-dreads-prop` is for determing
    the right reference for ambigous reads.
- `--keep-perfect-match` is not recommended, which decreases sensitivity.
- `-n/--keep-top-qcovs`  is not recommended, which affects accuracy of
     abundance estimation.

**Taxonomy data**:

1. Mapping references IDs to TaxIds: `-T/--taxid-map`
2. NCBI taxonomy dump files: `-X/--taxdump`

**Performance notes**:

1. Searching results are parsed in parallel, and the number of
    lines proceeded by a thread can be set by the flag `--chunk-size`.
2. However using a lot of threads does not always accelerate
    processing, 4 threads with chunk size of 500-5000 is fast enough.

Profiling output formats:

- KMCP      (`-o/--out-prefix`)
- CAMI      (`-M/--metaphlan-report`)
- MetaPhlAn (`-C/--cami-report`)

Taxonomic binning formats:
- CAMI      (`-B/--binning-result`)

**Commands**:

    # taxid mapping files, multiple files supported.
    taxid_map=refseq-viruses.kmcp/taxid.map,gtdb.kmcp/taxid.map

    # taxdump directory
    taxdump=taxdump/

    sfile=$file.kmcp.tsv.gz
    
    kmcp profile \
        --taxid-map      $taxid_map \
        --taxdump          $taxdump \
        --level             species \
        --min-query-cov        0.55 \
        --keep-top-qcovs          0 \
        --min-frags-reads        50 \
        --min-frags-prop        0.8 \
        --max-frags-cov-stdev     2 \
        --min-uniq-reads         10 \
        --min-hic-ureads          1 \
        --min-hic-ureads-qcov  0.75 \
        --min-hic-ureads-prop   0.1 \
        $sfile                      \
        --out-prefix       $sfile.kmcp.profile \
        --metaphlan-report $sfile.metaphlan.profile \
        --cami-report      $sfile.cami.profile \
        --sample-id        "0" \
        --binning-result   $sfile.binning

Default output:

|ref        |percentage|score |fragsProp|fragsRelCov                                      |fragsRelCovStd|reads |ureads|hicureads|refsize|refname|taxid |rank   |taxname                                  |taxpath                                                                                                                                              |taxpathsn                                        |
|:----------|:---------|:-----|:--------|:------------------------------------------------|:-------------|:-----|:-----|:--------|:------|:------|:-----|:------|:----------------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------|:------------------------------------------------|
|NC_000913.3|87.231001 |100.00|1.00     |1.04;0.99;0.99;1.00;0.99;0.99;0.99;0.97;1.04;0.99|0.02          |267738|191287|191287   |4641652|       |511145|no rank|Escherichia coli str. K-12 substr. MG1655|Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia coli;Escherichia coli K-12                   |2;1224;1236;91347;543;561;562;83333              |
|NC_002695.2|11.872107 |100.00|1.00     |0.97;0.99;0.86;1.06;1.02;0.93;1.04;1.05;1.08;1.01|0.06          |43166 |24210 |24210    |5498578|       |386585|strain |Escherichia coli O157:H7 str. Sakai      |Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia coli;Escherichia coli O157:H7 str. Sakai     |2;1224;1236;91347;543;561;562;386585             |
|NC_010655.1|0.896892  |100.00|1.00     |1.03;0.87;0.90;0.98;1.15;1.17;0.90;0.96;0.96;1.09|0.10          |1580  |1580  |1580     |2664102|       |349741|strain |Akkermansia muciniphila ATCC BAA-835     |Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;Akkermansia muciniphila;Akkermansia muciniphila ATCC BAA-835|2;74201;203494;48461;1647988;239934;239935;349741|

[CAMI format](https://github.com/CAMI-challenge/contest_information/blob/master/file_formats/CAMI_TP_specification.mkd):

    @SampleID:
    @Version:0.10.0
    @Ranks:superkingdom|phylum|class|order|family|genus|species|strain
    @TaxonomyID:ncbi-taxonomy
    @@TAXID	RANK	TAXPATH	TAXPATHSN	PERCENTAGE
    2	superkingdom	2	Bacteria	100.000000
    1224	phylum	2|1224	Bacteria|Proteobacteria	99.103108
    74201	phylum	2|74201	Bacteria|Verrucomicrobia	0.896892
    1236	class	2|1224|1236	Bacteria|Proteobacteria|Gammaproteobacteria	99.103108
    203494	class	2|74201|203494	Bacteria|Verrucomicrobia|Verrucomicrobiae	0.896892
    91347	order	2|1224|1236|91347	Bacteria|Proteobacteria|Gammaproteobacteria|Enterobacterales	99.103108
    48461	order	2|74201|203494|48461	Bacteria|Verrucomicrobia|Verrucomicrobiae|Verrucomicrobiales	0.896892
    543	family	2|1224|1236|91347|543	Bacteria|Proteobacteria|Gammaproteobacteria|Enterobacterales|Enterobacteriaceae	99.103108
    1647988	family	2|74201|203494|48461|1647988	Bacteria|Verrucomicrobia|Verrucomicrobiae|Verrucomicrobiales|Akkermansiaceae	0.896892
    561	genus	2|1224|1236|91347|543|561	Bacteria|Proteobacteria|Gammaproteobacteria|Enterobacterales|Enterobacteriaceae|Escherichia	99.103108
    239934	genus	2|74201|203494|48461|1647988|239934	Bacteria|Verrucomicrobia|Verrucomicrobiae|Verrucomicrobiales|Akkermansiaceae|Akkermansia	0.896892
    562	species	2|1224|1236|91347|543|561|562	Bacteria|Proteobacteria|Gammaproteobacteria|Enterobacterales|Enterobacteriaceae|Escherichia|Escherichia coli	99.103108
    239935	species	2|74201|203494|48461|1647988|239934|239935	Bacteria|Verrucomicrobia|Verrucomicrobiae|Verrucomicrobiales|Akkermansiaceae|Akkermansia|Akkermansia muciniphila	0.896892
    83333	strain	2|1224|1236|91347|543|561|562|83333	Bacteria|Proteobacteria|Gammaproteobacteria|Enterobacterales|Enterobacteriaceae|Escherichia|Escherichia coli|Escherichia coli K-12	87.230945
    386585	strain	2|1224|1236|91347|543|561|562|386585	Bacteria|Proteobacteria|Gammaproteobacteria|Enterobacterales|Enterobacteriaceae|Escherichia|Escherichia coli|Escherichia coli O157:H7 str. Sakai	11.872163
    349741	strain	2|74201|203494|48461|1647988|239934|239935|349741	Bacteria|Verrucomicrobia|Verrucomicrobiae|Verrucomicrobiales|Akkermansiaceae|Akkermansia|Akkermansia muciniphila|Akkermansia muciniphila ATCC BAA-835	0.896892


Metaphlan format:

    #SampleID	
    k__Bacteria	100.000000
    k__Bacteria|p__Proteobacteria	99.103108
    k__Bacteria|p__Verrucomicrobia	0.896892
    k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria	99.103108
    k__Bacteria|p__Verrucomicrobia|c__Verrucomicrobiae	0.896892
    k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales	99.103108
    k__Bacteria|p__Verrucomicrobia|c__Verrucomicrobiae|o__Verrucomicrobiales	0.896892
    k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae	99.103108
    k__Bacteria|p__Verrucomicrobia|c__Verrucomicrobiae|o__Verrucomicrobiales|f__Akkermansiaceae	0.896892
    k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Escherichia	99.103108
    k__Bacteria|p__Verrucomicrobia|c__Verrucomicrobiae|o__Verrucomicrobiales|f__Akkermansiaceae|g__Akkermansia	0.896892
    k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Escherichia|s__Escherichia coli	99.103108
    k__Bacteria|p__Verrucomicrobia|c__Verrucomicrobiae|o__Verrucomicrobiales|f__Akkermansiaceae|g__Akkermansia|s__Akkermansia muciniphila	0.896892
    k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Escherichia|s__Escherichia coli|t__Escherichia coli K-12	87.230945
    k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Escherichia|s__Escherichia coli|t__Escherichia coli O157:H7 str. Sakai	11.872163
    k__Bacteria|p__Verrucomicrobia|c__Verrucomicrobiae|o__Verrucomicrobiales|f__Akkermansiaceae|g__Akkermansia|s__Akkermansia muciniphila|t__Akkermansia muciniphila ATCC BAA-835	0.896892

[Binning result](https://github.com/CAMI-challenge/contest_information/blob/master/file_formats/CAMI_B_specification.mkd):

    # This is the bioboxes.org binning output format at
    # https://github.com/bioboxes/rfc/tree/master/data-format
    @Version:0.10.0
    @SampleID:
    @@SEQUENCEID	TAXID
    NC_000913.3_sliding:1244941-1245090	511145
    NC_002695.2_sliding:3465891-3466040	562
    NC_000913.3_sliding:3801041-3801190	511145
    NC_002695.2_sliding:3230881-3231030	562
    NC_000913.3_sliding:4080871-4081020	562
    NC_000913.3_sliding:3588091-3588240	511145
    NC_000913.3_sliding:2249621-2249770	562
    NC_000913.3_sliding:109271-109420	562
    NC_000913.3_sliding:2354841-2354990	511145
    NC_002695.2_sliding:4376441-4376590	386585
