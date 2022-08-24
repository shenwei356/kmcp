# Demo of taxonomic profiling

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

Or create custom taxdump files with `taxonomy.tsv` using [taxonkit create-taxdump](https://bioinf.shenwei.me/taxonkit/usage/#create-taxdump) (v0.12.1 or later versions required):

    taxonkit create-taxdump -A 1 taxonomy.tsv -O taxdump-custom/

## Metagenomic Profiling

Building database:

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

Generating mock dataset. 

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
        | seqkit shuffle -o mock.fastq.gz

Searching

    # searching
    for f in *.fastq.gz; do
        kmcp search \
            --db-dir refs-k21-n10.kmcp/ \
            --min-query-cov 0.55 \
            $f \
            --out-file $f.kmcp.gz
    done

Profiling

    # profiling using mode 1 for low coverage data
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

    cat mock.fastq.gz.kmcp.gz.kmcp.profile \
        | csvtk cut -t -f ref,percentage,taxname \
        | csvtk csv2md -t
    
|ref            |percentage|taxname                   |
|:--------------|:---------|:-------------------------|
|GCF_003697165.2|23.974382 |Escherichia coli          |
|GCF_000742135.1|20.155870 |Klebsiella pneumoniae     |
|GCF_000006945.2|19.407535 |Salmonella enterica       |
|GCF_002950215.1|15.583846 |Shigella flexneri         |
|GCF_002949675.1|14.767599 |Shigella dysenteriae      |
|GCF_001544255.1|1.852659  |Enterococcus faecium      |
|GCF_000392875.1|1.844333  |Enterococcus faecalis     |
|GCF_001027105.1|0.929976  |Staphylococcus aureus     |
|GCF_006742205.1|0.923797  |Staphylococcus epidermidis|
|GCF_001096185.1|0.192491  |Streptococcus pneumoniae  |
|GCF_000148585.2|0.186310  |Streptococcus mitis       |
|GCF_000017205.1|0.091788  |Pseudomonas aeruginosa    |
|GCF_009759685.1|0.089414  |Acinetobacter baumannii   |
