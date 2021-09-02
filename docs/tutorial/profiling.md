# Metagenomic profiling

## Requirements

- Database
    - [Prebuilt databases](/database/#prebuilt-databases) are available.
    - Or [build custom database](/database/#custom-database).
- Hardware.
    - CPU: ≥ 32 cores preferred.
    - RAM: ≥ 64 GB, depends on file size of the maximum database.

## Datasets

- Short reads, single or paired end, but paired information is not used.
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

Reference genome:

- Human:
    - [CHM13](https://github.com/marbl/CHM13)
- Mouse:
    - []()
- Rat:
    - []()

Building an index:
    
    bowtie2-build --threads 32 chm13.draft_v1.1.fasta.gz chm13

Mapping and removing mapped reads:

    index=~/ws/db/bowtie2/chm13
    
    bowtie2 --threads 32 -x $index -1 in_1.fq.gz -2 in_2.fq.gz  \
        | samtools view -buS -f 4 - \
        | samtools fastq - \
        | gzip -c > sample.fq.gz

### Step 3. Searching

Reads can be searched on different databases and then merged.

    # paired-ends: 
    #   file="sample_1.fq.gz sample_2.fq.gz"
    file=sample.fq.gz

    for db in viruses.kmcp prokaryotes.kmcp ; do
        dbname=$(basename $db)

        kmcp search \
            --threads            32 \
            --db-dir            $db \
            --min-kmers          30 \
            --min-query-len      70 \
            --min-query-cov     0.6 \
            --keep-top-scores     5 \
            $file                   \
            --out-file $file.kmcp@$dbname.tsv.gz \
            --log      $file.kmcp@$dbname.tsv.gz.log
    done

Merging searching results on multiple database:

    kmcp merge $file.kmcp@*.tsv.gz --out-file $file.kmcp.tsv.gz

### Step 4. Profiling

Input:

- TaxId mapping file(s).
- Taxdump files.
- KMCP search results.

Commands:

    # taxid mapping files:
    taxid_map=taxid-prokaryotes.map,taxid-virus.map

    # taxdump directory
    taxdump=taxdump/

    sfile=$file.kmcp.tsv.gz
    
    kmcp profile \
        --taxid-map      $taxid_map \
        --taxdump          $taxdump \
        --level             species \
        --keep-top-qcovs          5 \
        --min-query-cov         0.6 \
        --min-reads              50 \
        --min-uniq-reads         10 \
        --min-hic-ureads          1 \
        --min-hic-ureads-qcov   0.8 \
        --min-hic-ureads-prop   0.1 \
        $sfile                      \
        --out-prefix       $sfile.kmcp.profile \
        --metaphlan-report $sfile.metaphlan.profile \
        --cami-report      $sfile.cami.profile \
        --sample-id        "0"
