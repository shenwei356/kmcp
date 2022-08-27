# Demo

## Dataset

### References

    # name.map
    ls refs/*.gz  | rush 'echo -ne "{%..}\t"; seqkit head -n 1 {} | seqkit seq -n | cut -d " " -f 2-' > name.map
    
    csvtk join -t -f id \
            <(seqkit stats -j 10 refs/*.gz -T -b \
                | csvtk mutate -t -n id -p "(.+)\.fa") \
            <(csvtk add-header -t -n id,name name.map) \
        | csvtk cut -t -f file,num_seqs,sum_len,name \
        | csvtk sort -t -k name \
        | csvtk csv2md -t
        
|file                 |num_seqs|sum_len|name                                                                              |
|:--------------------|:-------|:------|:---------------------------------------------------------------------------------|
|GCF_009759685.1.fa.gz|2       |3990388|Acinetobacter baumannii strain ATCC 19606 chromosome, complete genome             |
|GCF_000392875.1.fa.gz|3       |2881400|Enterococcus faecalis ATCC 19433 acAqW-supercont1.1, whole genome shotgun sequence|
|GCF_001544255.1.fa.gz|38      |2484851|Enterococcus faecium NBRC 100486, whole genome shotgun sequence                   |
|GCF_003697165.2.fa.gz|2       |5034834|Escherichia coli DSM 30083 = JCM 1649 = ATCC 11775 chromosome, complete genome    |
|GCF_000742135.1.fa.gz|5       |5545784|Klebsiella pneumoniae strain ATCC 13883 scaffold1, whole genome shotgun sequence  |
|GCF_000017205.1.fa.gz|1       |6588339|Pseudomonas aeruginosa PA7, complete genome                                       |
|GCF_000006945.2.fa.gz|2       |4951383|Salmonella enterica subsp. enterica serovar Typhimurium str. LT2, complete genome |
|GCF_002949675.1.fa.gz|2       |4578459|Shigella dysenteriae strain ATCC 13313 chromosome, complete genome                |
|GCF_002950215.1.fa.gz|3       |4938295|Shigella flexneri 2a strain ATCC 29903 chromosome, complete genome                |
|GCF_001027105.1.fa.gz|2       |2782562|Staphylococcus aureus subsp. aureus DSM 20231 chromosome, complete genome         |
|GCF_006742205.1.fa.gz|2       |2427041|Staphylococcus epidermidis NBRC 100911 DNA, complete genome                       |
|GCF_000148585.2.fa.gz|1       |1868883|Streptococcus mitis NCTC 12261 chromosome, complete genome                        |
|GCF_001096185.1.fa.gz|24      |2117177|Streptococcus pneumoniae strain SMRU824, whole genome shotgun sequence            |

### Taxonomy data

Please download and uncompress [taxdump.tar.gz](ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz),
and then copy `names.dmp`, `nodes.dmp`, `delnodes.dmp` and `merged.dmp` to directory `taxdump`.

Or create custom taxdump files with `taxonomy.tsv` using [taxonkit create-taxdump](https://bioinf.shenwei.me/taxonkit/usage/#create-taxdump):

    taxonkit create-taxdump -A 1 taxonomy.tsv -O taxdump-custom/

### Mock community 1

Generating mock dataset 1 using seqkit, no SNP, deletion, and insertion.

    # generating mock dataset
    (seqkit sliding -s 10 -W 150 refs/GCF_003697165.2.fa.gz | seqkit shuffle | seqkit sample -p 0.2   ; \
     seqkit sliding -s 10 -W 150 refs/GCF_002949675.1.fa.gz | seqkit shuffle | seqkit sample -p 0.2   ; \
     seqkit sliding -s 10 -W 150 refs/GCF_002950215.1.fa.gz | seqkit shuffle | seqkit sample -p 0.2   ; \
     seqkit sliding -s 10 -W 150 refs/GCF_000742135.1.fa.gz | seqkit shuffle | seqkit sample -p 0.2   ; \
     seqkit sliding -s 10 -W 150 refs/GCF_000006945.2.fa.gz | seqkit shuffle | seqkit sample -p 0.2   ; \
     seqkit sliding -s 10 -W 150 refs/GCF_000392875.1.fa.gz | seqkit shuffle | seqkit sample -p 0.02  ; \
     seqkit sliding -s 10 -W 150 refs/GCF_001544255.1.fa.gz | seqkit shuffle | seqkit sample -p 0.02  ; \
     seqkit sliding -s 10 -W 150 refs/GCF_001027105.1.fa.gz | seqkit shuffle | seqkit sample -p 0.01  ; \
     seqkit sliding -s 10 -W 150 refs/GCF_006742205.1.fa.gz | seqkit shuffle | seqkit sample -p 0.01  ; \
     seqkit sliding -s 10 -W 150 refs/GCF_000148585.2.fa.gz | seqkit shuffle | seqkit sample -p 0.002 ; \
     seqkit sliding -s 10 -W 150 refs/GCF_001096185.1.fa.gz | seqkit shuffle | seqkit sample -p 0.002 ; \
     seqkit sliding -s 10 -W 150 refs/GCF_000017205.1.fa.gz | seqkit shuffle | seqkit sample -p 0.001 ; \
     seqkit sliding -s 10 -W 150 refs/GCF_009759685.1.fa.gz | seqkit shuffle | seqkit sample -p 0.001 ) \
        | seqkit shuffle -o mock1.fastq.gz
        
### Mock community 2

Generating mock dataset 2 using https://github.com/iqbal-lab-org/simutator, 
every 2kb where a region of length 1.5kb gets 30 snps, 2 insertion and 4 deletions up to length 10bp.

Relative depths:

    $ cat depth.tsv 
    GCF_003697165.2.fa      1
    GCF_002949675.1.fa      1
    GCF_002950215.1.fa      1
    GCF_000742135.1.fa      1
    GCF_000006945.2.fa      1
    GCF_000392875.1.fa      0.1
    GCF_001544255.1.fa      0.1
    GCF_001027105.1.fa      0.05
    GCF_006742205.1.fa      0.05
    GCF_000148585.2.fa      0.01
    GCF_001096185.1.fa      0.01
    GCF_000017205.1.fa      0.005
    GCF_009759685.1.fa      0.005

Steps:
    
    # tools:
    #   - https://github.com/shenwei356/rush

    # unzip all references which are required by 'simutator'
    mkdir -p refs.plain
    ls refs/*.gz | rush 'seqkit seq {} -o refs.plain/{%.}'
    cd refs.plain
    
    # mutate
    ls *.fa | rush 'simutator mutate_fasta --complex 2000:1500:30:2:4:10 {} {}.m.fasta'
    rm *.vcf
    
    # simulate reads, since simutator only supports depth of integer,
    # we'll generate reads of depth 2 first and then perform downsampling.
    ls *m.fasta*.fa | rush 'simutator simulate_reads {} {@^(.+?\.fa)} --read_length 150 --read_depth 2 --fragment_length 350'

    # downsampling
    ls *.fa | grep -v m.fasta \
        | rush 'p=$(grep {%} ../depth.tsv | cut -f 2); \
            for f in {}.*.{1,2}.fq.gz; do seqkit sample -p $p $f -o $f.sample.fq.gz; done'
            
    # concatenate all sampled reads
    cd ..
    seqkit seq refs.plain/*.1.fq.gz.sample.fq.gz -o mock2_1.fastq.gz
    seqkit seq refs.plain/*.2.fq.gz.sample.fq.gz -o mock2_2.fastq.gz
    
## Metagenomic Profiling

### Building database

    # computing k-mers
    kmcp compute \
        --in-dir refs/ \
        --ref-name-regexp "^([\w\.\_]+\.\d+)" \
        --seq-name-filter "plasmid" \
        --kmer 21 \
        --split-number 10 \
        --split-overlap 150 \
        --out-dir refs-k21-n10 \
        --force

    # indexing k-mers
    kmcp index \
        --in-dir refs-k21-n10/\
        --num-hash 1 \
        --false-positive-rate 0.3 \
        --out-dir refs-k21-n10.kmcp \
        --force

### Searching and Profiling

Searching:
    
    # mock1
    kmcp search \
        --db-dir refs-k21-n10.kmcp/ \
        --min-query-cov 0.55 \
        mock1.fastq.gz \
        --out-file mock1.kmcp.gz \
        --log mock1.kmcp.gz.log
        
    # mock2
    # paired information are not used.
    kmcp search \
        --db-dir refs-k21-n10.kmcp/ \
        --min-query-cov 0.55 \
         mock2_1.fastq.gz mock2_2.fastq.gz \
        --out-file mock2.kmcp.gz \
        --log mock2.kmcp.gz.log
    

Profiling using mode 1 for low coverage data:

    for f in *.kmcp.gz; do
        kmcp profile \
            --taxid-map taxdump-custom/taxid.map \
            --taxdump taxdump-custom/ \
            $f \
            --mode 1 \
            --out-prefix $f.kmcp.profile \
            --metaphlan-report $f.metaphlan.profile \
            --cami-report $f.cami.profile \
            --binning-result $f.binning.gz
    done
    
Results (mock 1):

    grep "queries matched" mock1.kmcp.gz.log
    16:50:29.364 [INFO] 97.8931% (507195/518111) queries matched
    
    cat mock1.kmcp.gz.kmcp.profile \
        | csvtk cut -t -f ref,percentage,taxname \
        | csvtk csv2md -t
    
|ref            |percentage|taxname                   |
|:--------------|:---------|:-------------------------|
|GCF_000742135.1|19.123795 |Klebsiella pneumoniae     |
|GCF_003697165.2|18.887595 |Escherichia coli          |
|GCF_002949675.1|18.778146 |Shigella dysenteriae      |
|GCF_000006945.2|18.680933 |Salmonella enterica       |
|GCF_002950215.1|18.461495 |Shigella flexneri         |
|GCF_001544255.1|1.840951  |Enterococcus faecium      |
|GCF_000392875.1|1.829837  |Enterococcus faecalis     |
|GCF_001027105.1|0.922711  |Staphylococcus aureus     |
|GCF_006742205.1|0.918477  |Staphylococcus epidermidis|
|GCF_000148585.2|0.190200  |Streptococcus mitis       |
|GCF_001096185.1|0.187513  |Streptococcus pneumoniae  |
|GCF_000017205.1|0.090965  |Pseudomonas aeruginosa    |
|GCF_009759685.1|0.087383  |Acinetobacter baumannii   |

Results (mock 2):

    grep "queries matched" mock2.kmcp.gz.log
    16:46:43.561 [INFO] 89.3701% (307276/343824) queries matched

    cat mock2.kmcp.gz.kmcp.profile \
        | csvtk cut -t -f ref,percentage,taxname \
        | csvtk csv2md -t

|ref            |percentage|taxname                   |
|:--------------|:---------|:-------------------------|
|GCF_002949675.1|19.307308 |Shigella dysenteriae      |
|GCF_002950215.1|18.960647 |Shigella flexneri         |
|GCF_000742135.1|18.924976 |Klebsiella pneumoniae     |
|GCF_000006945.2|18.612709 |Salmonella enterica       |
|GCF_003697165.2|18.118441 |Escherichia coli          |
|GCF_001544255.1|1.883821  |Enterococcus faecium      |
|GCF_000392875.1|1.791524  |Enterococcus faecalis     |
|GCF_006742205.1|0.918205  |Staphylococcus epidermidis|
|GCF_001027105.1|0.911963  |Staphylococcus aureus     |
|GCF_000148585.2|0.192247  |Streptococcus mitis       |
|GCF_001096185.1|0.186028  |Streptococcus pneumoniae  |
|GCF_009759685.1|0.099627  |Acinetobacter baumannii   |
|GCF_000017205.1|0.092503  |Pseudomonas aeruginosa    |

