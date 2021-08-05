## Demo dataset

    csvtk join -t -f id \
            <(seqkit stats -j 10 refs/*.fasta.gz -T -b \
                | csvtk mutate -t -n id -p "(.+)\.fasta") \
            <(csvtk add-header -t -n id,name name.map) \
        | csvtk cut -t -f file,format,type,num_seqs,sum_len,name \
        | csvtk csv2md -t

file                  |format|type|num_seqs|sum_len|name
:---------------------|:-----|:---|:-------|:------|:----------------------------------------
NC_000913.3.fasta.gz  |FASTA |DNA |1       |4641652|Escherichia coli str. K-12 substr. MG1655
NC_002695.2.fasta.gz  |FASTA |DNA |1       |5498578|Escherichia coli O157:H7 str. Sakai
NC_010655.1.fasta.gz  |FASTA |DNA |1       |2664102|Akkermansia muciniphila ATCC BAA-835
NC_011750.1.fasta.gz  |FASTA |DNA |1       |5132068|Escherichia coli IAI39
NC_012971.2.fasta.gz  |FASTA |DNA |1       |4558953|Escherichia coli BL21(DE3)
NC_013654.1.fasta.gz  |FASTA |DNA |1       |4717338|Escherichia coli SE15
NC_018658.1.fasta.gz  |FASTA |DNA |1       |5273097|Escherichia coli O104:H4 str. 2011C-3493
NZ_CP007592.1.fasta.gz|FASTA |DNA |1       |5104557|Escherichia coli O157:H16 strain Santai
NZ_CP028116.1.fasta.gz|FASTA |DNA |1       |5648177|Escherichia coli O26 str. RM8426

## Metagenomic Profiling

    # computing k-mers
    kmcp compute \
        --in-dir refs/ \
        --ref-name-regexp "^([\w\.\_]+\.\d+)" \
        --seq-name-filter "plasmid" \
        --kmer 31 --exact-number \
        --split-number 10 \
        --out-dir refs-k31-n10 \
        --force

    # indexing k-mers
    kmcp index \
        --in-dir refs-k31-n10/\
        --num-hash 3 \
        --false-positive-rate 0.01 \
        --out-dir refs-k31-n10.kmcp \
        --force

    # generating mock dataset
    (seqkit sliding -s 10 -W 150 refs/NC_000913.3.fasta.gz | seqkit sample -p 0.6 ; \
            seqkit sliding -s 10 -W 150 refs/NC_002695.2.fasta.gz | seqkit sample -p 0.06 ;  \
            seqkit sliding -s 10 -W 150 refs/NC_010655.1.fasta.gz | seqkit sample -p 0.006 ) \
        | seqkit shuffle -o mock.fastq.gz

    # searching
    for f in *.fastq.gz; do
        kmcp search \
            --db-dir refs-k31-n10.kmcp/ \
            --min-query-cov 0.8 \
            $f \
            --out-file $f.kmcp.gz
    done

    # profiling
    for f in *.kmcp.gz; do
        kmcp profile \
            --taxid-map taxid.map \
            --taxdump taxdump \
            $f \
            --level strain \
            --min-reads 20 \
            --min-uniq-reads 5 \
            --out-prefix $f.kmcp.profile \
            --metaphlan-report $f.metaphlan.profile \
            --cami-report $f.cami.profile
    done

    cat mock.fastq.gz.kmcp.gz.kmcp.profile \
        | csvtk cut -t -f ref,percentage,taxname \
        | csvtk csv2md -t
    
ref        |percentage|taxname
:----------|:---------|:----------------------------------------
NC_000913.3|87.230664 |Escherichia coli str. K-12 substr. MG1655
NC_002695.2|11.872444 |Escherichia coli O157:H7 str. Sakai
NC_010655.1|0.896892  |Akkermansia muciniphila ATCC BAA-835

## Sequence containment searching

    for f in *.fastq.gz; do
        kmcp search \
            --db-dir refs-k31-n10.kmcp/ \
            --min-query-cov 0.8 \
            $f \
            --out-file $f.kmcp.gz
    done

## Genome distance estimation

### Minimizer

    # computing k-mers
    kmcp compute \
        --in-dir refs/ \
        --ref-name-regexp "^([\w\.\_]+\.\d+)" \
        --seq-name-filter "plasmid" \
        --kmer 31 --exact-number \
        --minimizer-w 20 \
        --out-dir refs-k31-W20 \
        --force

    # indexing k-mers
    kmcp index \
        --in-dir refs-k31-W20/ \
        --num-hash 3 \
        --false-positive-rate 0.01 \
        --out-dir refs-k31-W20.kmcp \
        --force

    # searching
    kmcp search \
            --db-dir refs-k31-W20.kmcp \
            --query-whole-file \
            --min-query-cov 0.5 \
            --keep-top-scores 0 \
            --sort-by jacc \
            --name-map name.map \
            refs/NC_018658.1.fasta.gz \
        | csvtk cut -t -C '$' -f '#query,target,qCov,tCov,jacc' \
        | csvtk pretty -C '$' -t 

    #query        target                                      qCov     tCov     jacc
    -----------   -----------------------------------------   ------   ------   ------
    NC_018658.1   Escherichia coli O104:H4 str. 2011C-3493    1.0000   1.0000   1.0000
    NC_018658.1   Escherichia coli O26 str. RM8426            0.7372   0.7151   0.5698
    NC_018658.1   Escherichia coli str. K-12 substr. MG1655   0.5977   0.6728   0.4631
    NC_018658.1   Escherichia coli BL21(DE3)                  0.5884   0.6751   0.4585
    NC_018658.1   Escherichia coli O157:H16 strain Santai     0.5687   0.5816   0.4036
    NC_018658.1   Escherichia coli O157:H7 str. Sakai         0.5377   0.5268   0.3626

### Syncmer

    # computing k-mers
    kmcp compute \
        --in-dir refs/ \
        --ref-name-regexp "^([\w\.\_]+\.\d+)" \
        --seq-name-filter "plasmid" \
        --kmer 31 --exact-number \
        --syncmer-s 11 \
        --out-dir refs-k31-S11 \
        --force

    # indexing k-mers
    kmcp index \
        --in-dir refs-k31-S11/\
        --num-hash 3 \
        --false-positive-rate 0.01 \
        --out-dir refs-k31-S11.kmcp \
        --force

    # searching
    kmcp search \
            --db-dir refs-k31-S11.kmcp \
            --query-whole-file \
            --min-query-cov 0.5 \
            --keep-top-scores 0 \
            --sort-by jacc \
            --name-map name.map \
            refs/NC_018658.1.fasta.gz \
        | csvtk cut -t -C '$' -f '#query,target,qCov,tCov,jacc' \
        | csvtk pretty -C '$' -t 

    #query        target                                      qCov     tCov     jacc
    -----------   -----------------------------------------   ------   ------   ------
    NC_018658.1   Escherichia coli O104:H4 str. 2011C-3493    1.0000   1.0000   1.0000
    NC_018658.1   Escherichia coli O26 str. RM8426            0.7435   0.7212   0.5775
    NC_018658.1   Escherichia coli str. K-12 substr. MG1655   0.6061   0.6817   0.4724
    NC_018658.1   Escherichia coli BL21(DE3)                  0.5967   0.6841   0.4678
    NC_018658.1   Escherichia coli O157:H16 strain Santai     0.5779   0.5917   0.4132
    NC_018658.1   Escherichia coli O157:H7 str. Sakai         0.5487   0.5373   0.3726

### Scaled MinHash

    # computing k-mers
    kmcp compute \
        --in-dir refs/ \
        --ref-name-regexp "^([\w\.\_]+\.\d+)" \
        --seq-name-filter "plasmid" \
        --kmer 31 --exact-number \
        --scale 20 \
        --out-dir refs-k31-D20 \
        --force

    # indexing k-mers
    kmcp index \
        --in-dir refs-k31-D20/\
        --num-hash 3 \
        --false-positive-rate 0.01 \
        --out-dir refs-k31-D20.kmcp \
        --force

    # searching
    kmcp search \
            --db-dir refs-k31-D20.kmcp \
            --query-whole-file \
            --min-query-cov 0.5 \
            --keep-top-scores 0 \
            --sort-by jacc \
            --name-map name.map \
            refs/NC_018658.1.fasta.gz \
        | csvtk cut -t -C '$' -f '#query,target,qCov,tCov,jacc' \
        | csvtk pretty -C '$' -t 

    #query        target                                      qCov     tCov     jacc
    -----------   -----------------------------------------   ------   ------   ------
    NC_018658.1   Escherichia coli O104:H4 str. 2011C-3493    1.0000   1.0000   1.0000
    NC_018658.1   Escherichia coli O26 str. RM8426            0.7485   0.7257   0.5834
    NC_018658.1   Escherichia coli str. K-12 substr. MG1655   0.6117   0.6879   0.4788
    NC_018658.1   Escherichia coli BL21(DE3)                  0.6022   0.6903   0.4741
    NC_018658.1   Escherichia coli O157:H16 strain Santai     0.5853   0.5986   0.4203
    NC_018658.1   Escherichia coli O157:H7 str. Sakai         0.5573   0.5456   0.3806
