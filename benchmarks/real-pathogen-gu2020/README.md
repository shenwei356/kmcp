# Benchmarks on 87 metagenomic samples of infected body fluids

## Softwares

- KMCP [v0.8.0](https://github.com/shenwei356/kmcp/releases/tag/v0.8.0)
- Kraken [v2.1.2 (2021-05-10)](https://github.com/DerrickWood/kraken2/releases/tag/v2.1.2),
  Bracken [v2.6.2 (2021-03-22)](https://github.com/jenniferlu717/Bracken/releases/tag/v2.6.2)
- Centrifuge [v1.0.4 (2021-08-17)](https://github.com/DaehwanKimLab/centrifuge/releases/tag/v1.0.4)

## Databases and taxonomy version

- KMCP,  GTDB-RS202 (2021-04-27) + Genbank-viral (r246, 2021-12-06) + Refseq-fungi (r208, 2021-09-30), 2021-12-06
- Kraken, PlusPF (2021-05-17), 2021-05-17
- Centrifuge, built with the genomes same to KMCP.
- Kraken, built with the genomes same to KMCP.

**We create databases of GTDB and Refseq-fungi with a smaller false-positive rate `0.1` instead of `0.3`,
and use `2` hash functions instead of `1`.
The size of GTDB database increase fom 58 to 109GB, and that of Refseq-fungi from 4.2 to 7.9GB.
We use a small query coverage threshhold `0.4` instead of `0.55` during searching and profiling,
and use the re-built mode 0 (pathogen detection) in profiling.**

In this benchmark, we generate metagenomic profiles with the same NCBI Taxonomy version 2021-12-06,
including the gold-standard profiles.

## Datasets

> Gu, W., Deng, X., Lee, M. et al. Rapid pathogen detection by metagenomic next-generation 
> sequencing of infected body fluids.
> Nat Med 27, 115–124 (2021). https://doi.org/10.1038/s41591-020-1105-z

> A total of 182 body fluid samples from
> 160 patients, including 25 abscess, 21 joint, 32 pleural, 27 peritoneal,
> 35 cerebrospinal, 13 bronchoalveolar lavage (BAL) and 29 other
> body fluids (Table 1 and Supplementary Table 1), were collected as
> residual samples after routine clinical testing in the microbiology
> laboratory. Among these 182 samples, 170 were used to evaluate
> the accuracy of mNGS testing by Illumina sequencing (Fig. 1a and
> Supplementary Table 1). These accuracy samples included 127 posi-
> tive by culture (with pathogen(s) identified to genus or species
> level), nine culture negative but positive by 16S or 28S–ITS PCR
> and 34 negative controls from patients with alternative noninfec-
> tious diagnoses (for example, cancer, trauma) (Fig. 1b).

We download the 99 Illumina short-reads datasets from [PRJNA558701](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA558701),
in which 92 samples has detailed information in the supplementary table 1.
And 87 samples are verified by culture or 16S rRNA gene qPCR.

Run accession and corresponding sample names are listed in [acc2sample.tsv](https://github.com/shenwei356/kmcp/blob/main/benchmarks/pathogen-gu2020/acc2sample.tsv).

Download reads using SRAtoolkit:

    # download
    rush -j 12 'prefetch {1}' acc2sample.tsv
    
    # dump fastq
    ls SRR* | rush 'fasterq-dump -p {}'

    # compress fastq
    ls *.fastq | rush 'pigz {}'
    
    # brename
    brename -f ".fastq.gz$" -p '^(\w+)' -r '{kv}' -k acc2sample.tsv 
    
    mkdir reads
    mv *.fastq.gz reads/
    
    mkdir sra
    mv SRR* sra/
    
The files look like:

    ls reads/ | head -n 10
    P10.fastq.gz
    P42.fastq.gz
    P62.fastq.gz
    P78.fastq.gz
    P86.fastq.gz
    P87.fastq.gz
    P92.fastq.gz
    S10.fastq.gz
    S11.fastq.gz
    S12.fastq.gz

## KMCP

We search against GTDB, Genbank-viral, and Refseq-fungi respectively, and merge the results.

    # ------------------------------------------------------------------------
    # prepare folder and files
    
    # prepare folder and files.
    mkdir -p kmcp-se
    cd kmcp-se
    fd fastq.gz$ ../reads | rush 'ln -s {}'
    cd ..

   
    # ------------------------------------------------------------------------
    # search
      
    reads=kmcp-se
    j=4
    J=40
    
    # -------------------------------------
    # gtdb

    db=gtdb.kmcp/
    X=taxdump/
    T=gtdb.kmcp/taxid.map    
    dbname=gtdb
    
    fd fastq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -v db=$db -v dbname=$dbname -j $j -v j=$J -v 'p={:}' \
            'kmcp search -t 0.4 -d {db} {} \
                -o {p}.kmcp@{dbname}.tsv.gz \
                --log {p}.kmcp@{dbname}.tsv.gz.log -j {j}' \
            -c -C $reads@$dbname.rush
 
    # -------------------------------------
    # genbank-viral

    # genbank-viral
    db=genbank-viral.kmcp/
    X=taxdump/
    T=genbank-viral.kmcp/taxid.map    
    dbname=genbank-viral
    
    fd fastq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -v db=$db -v dbname=$dbname -j $j -v j=$J -v 'p={:}' \
            'kmcp search -t 0.4 -d {db} {} \
                -o {p}.kmcp@{dbname}.tsv.gz \
                --log {p}.kmcp@{dbname}.tsv.gz.log -j {j}' \
            -c -C $reads@$dbname.rush
    
    # -------------------------------------
    # refseq-fungi
    
    # refseq-fungi
    db=refseq-fungi.kmcp/
    X=taxdump/
    T=refseq-fungi.kmcp/taxid.map    
    dbname=refseq-fungi
    
    fd fastq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -v db=$db -v dbname=$dbname -j $j -v j=$J -v 'p={:}' \
            'kmcp search -t 0.4 -d {db} {} \
                -o {p}.kmcp@{dbname}.tsv.gz \
                --log {p}.kmcp@{dbname}.tsv.gz.log -j {j}' \
            -c -C $reads@$dbname.rush

    
    # ------------------------------------------------------------------------
    # merge results
    
    reads=kmcp-se
    j=16
    fd fastq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j $j -v 'p={:}' \
            'kmcp merge {p}.kmcp@*.tsv.gz -o {p}.kmcp.tsv.gz --log {p}.kmcp.tsv.gz.log'
    
    
    # ------------------------------------------------------------------------
    # [for merged search results] multiple profiling modes
     
    X=taxdump/
    cat genbank-viral.kmcp/taxid.map gtdb.kmcp/taxid.map refseq-fungi.kmcp/taxid.map > taxid.map
    cat genbank-viral.kmcp/name.map gtdb.kmcp/name.map refseq-fungi.kmcp/name.map > name.map
    T=taxid.map
    
    for m in $(seq 0 0); do
        fd kmcp.tsv.gz$ $reads/ \
            | csvtk sort -H -k 1:N \
            | rush -v X=$X -v T=$T -v m=$m \
                'kmcp profile -t 0.4 -m {m} -X {X} -T {T} {} -o {}.k-m{m}.profile -C {}.c-m{m}.profile -s {%:} --log {}.k-m{m}.profile.log' 
        
        profile=$reads.c-m$m.profile
        fd kmcp.tsv.gz.c-m$m.profile$ $reads/ \
            | csvtk sort -H -k 1:N \
            | rush -j 1 'cat {}' \
            > $profile
    done

    # ------------------------------------------------------------------------
    # [for merged search results] multiple profiling modes
    #  --no-amb-corr                   do not correct ambiguous reads (just for benchmark)
    
    for m in $(seq 0 0); do
        fd kmcp.tsv.gz$ $reads/ \
            | csvtk sort -H -k 1:N \
            | rush -v X=$X -v T=$T -v m=$m \
                'kmcp profile -t 0.4 -m {m} -X {X} -T {T} {} --no-amb-corr -o {}.k-m{m}.profile.no-amb-corr -s {%:}' 
    done


## Bracken

Preparing tocami.py which convert Bracken output to CAMI format

    # wget https://raw.githubusercontent.com/hzi-bifo/cami2_pipelines/master/bin/tocami.py
    chmod a+x tocami.py

    # install pacakge ete3
    pip install ete3
    
    # preparing taxdump.tar.gz for tocami.py
    tar -zcvf taxdump.tar.gz taxdump/*
    
    # creating database for ete3 (don't worry the error reports, just ignore):
    tocami.py -t taxdump.tar.gz -f bracken -s 1 -d . bracken-pe/Build_sample1.bracken

Steps

    # --------------------------------------------------
    # using kraken's PLUSPF database

    reads=bracken-se
    
    # prepare folder and files.
    mkdir -p $reads
    cd $reads
    fd fastq.gz$ ../reads | rush 'ln -s {}'
    cd ..

    reads=bracken-se    
    j=4
    J=40
    
    db=/home/shenwei/ws/db/kraken/pluspf    
    readlen=150
    threshold=2
    
    
    # --------------------------------------------------
    # using kraken database built with GTDB, Genbank-viral, Refseq-fungi
    
    reads = bracken-kmcp
    
    # prepare folder and files.
    mkdir -p $reads
    cd $reads
    fd fastq.gz$ ../reads | rush 'ln -s {}'
    cd ..

    reads=bracken-kmcp    
    j=4
    J=40
    
    db=/home/shenwei/ws/db/kraken/kmcp/kmcp
    readlen=150
    threshold=2
    
    
    fd fastq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j $j -v j=$J -v 'p={:}' -v db=$db -v db=$db -v readlen=$readlen -v threshold=$threshold \
            'memusg -t -s \
                "kraken2 --db {db} --threads {j} --memory-mapping --gzip-compressed  \
                    {} --report {p}.kreport > /dev/null; \
                for r in \"S\" \"G\" \"F\" \"O\" \"C\" \"P\"; do \
                    est_abundance.py -k {db}/database${readlen}mers.kmer_distrib -l \$r -t {threshold} \
                    -i {p}.kreport -o {p}.bracken.level-\$r ; \
                done; \
                cat {p}.bracken.level-* > {p}.bracken " \
                >{p}.log 2>&1 '

    # ------------------------------------------------------
    # convert to CAMI format
    fd .bracken$ $reads/ \
        | rush -j 12 'python3 ./tocami.py -d ./ -f bracken {} -s {%:} -o {}.profile'
    
    profile=$reads.profile
    fd bracken.profile$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j 1 'cat {}' \
        > $profile
        
    # sort by reads number
    fd level-S$ $reads/ \
        | rush 'csvtk sort -t -k new_est_reads:nr {} | csvtk pretty -t > {}.txt'

## Centrifuge

    # --------------------------------------------------
    # using centrifuge database built with GTDB, Genbank-viral, Refseq-fungi

    reads=centrifuge-pe
    
    # prepare folder and files.
    mkdir -p $reads
    cd $reads
    fd fastq.gz$ ../reads | rush 'ln -s {}'
    cd ..

    reads=centrifuge-pe    
    j=4
    J=40
    
    db=~/ws/db/centrifuge/kmcp  
    

    fd fastq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j $j -v j=$J -v 'p={:}' -v db=$db \
            'memusg -t -s \
                "centrifuge -t -p {j} --mm -x {db} \
                    -q -U {} \
                    --report-file {p}.cf-report.tsv -S {}.cf.tsv" \
                >{p}.log 2>&1 '

    # Kraken-style report
    fd fastq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j $j -v j=$J -v 'p={:}' -v db=$db \
            'centrifuge-kreport -x {db} {}.cf.tsv > {}.kreport'
                
    # ------------------------------------------------------
    # convert to CAMI format
    fd .kreport$ $reads/ \
        | rush -j 12 'python3 ./tocami.py -d ./ -f centrifuge {} -s {%:} -o {}.profile'
    
    profile=$reads.profile
    fd kreport.profile$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j 1 'cat {}' \
        > $profile

    # filter taxa with rank at or below species in reports and sort by numReads
    fd .cf-report.tsv$ $reads/ \
        | rush 'csvtk grep -t -f taxRank -v -p root,superkingdom,kingdom,phylum,class,family,order,genus {} \
            | csvtk filter2 -t -f "\$numReads >= 2" \
            | csvtk sort -t -k numReads:nr \
            | csvtk pretty -t > {}.txt'
