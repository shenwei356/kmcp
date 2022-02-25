# Benchmarks on 25 simulated prokaryotic communities from Sun et al

## Softwares

- KMCP [v0.8.0](https://github.com/shenwei356/kmcp/releases/tag/v0.8.0)
- mOTUs [3.0.1 (2021-07-28)](https://github.com/motu-tool/mOTUs/releases/tag/3.0.1)
- MetaPhlAn [3.0.13 (2021-07-27)](https://github.com/biobakery/MetaPhlAn/releases/tag/3.0.13)
- Kraken [v2.1.2 (2021-05-10)](https://github.com/DerrickWood/kraken2/releases/tag/v2.1.2),
  Bracken [v2.6.2 (2021-03-22)](https://github.com/jenniferlu717/Bracken/releases/tag/v2.6.2)
- Centrifuge [v1.0.4 (2021-08-17)](https://github.com/DaehwanKimLab/centrifuge/releases/tag/v1.0.4)
- DUDes [v0.08 (2017-11-08)](https://github.com/pirovc/dudes/releases/tag/dudes_v0.08)
- SLIMM [v0.3.4 (2018-09-04)](https://github.com/seqan/slimm/releases/tag/v0.3.4)

## Databases and taxonomy version

- KMCP,  GTDB-RS202 (2021-04-27) + Refseq-fungi (r208, 2021-09-30), 2021-12-06
- Centrifuge, built with the genomes same to KMCP.
- mOTUs, 3.0.1 (2021-06-28), 2019-01
- MetaPhlAn, mpa_v30_CHOCOPhlAn_201901 (?), 2019-01
- Kraken, PlusPF (2021-05-17), 2021-05-17
- Kraken, built with the genomes same to KMCP.
- DUDes, built with the genomes same to KMCP.
- SLIMM, built with the genomes same to KMCP.

In this benchmark, we generate metagenomic profiles with the same NCBI Taxonomy version 2021-12-06,
including the gold-standard profiles.

## Datasets

> We
> simulated metagenomic sequencing reads for 25 communities
> from distinct habitats (for example, gastrointestinal, oral, dermal,
> vaginal and building, five communities for each habitat; Methods).
> To avoid reference database bias of different metagenomic profil-
> ers, the genomes used to generate simulated communities were
> selected from the intersection among the reference databases of
> MetaPhlAn2, mOTUs2 and Kraken2.
> 
> Sun, Z., Huang, S., Zhang, M. et al. Challenges in benchmarking metagenomic profilers. Nat Methods > 18, 618–626 (2021). https://doi.org/10.1038/s41592-021-01141-3

Sun *et al.* simulated [25 metagenomic reads](https://figshare.com/projects/Pitfalls_and_Opportunities_in_Benchmarking_Metagenomic_Classifiers/79916) for benchmarking metagenomic profilers,
while the ground truth profiles format is not convenient for interpretation
with tools like [opal](https://github.com/CAMI-challenge/OPAL).
I've convert the ground truth profiles to CAMI format with Taxonomy version 2021-12-06,
available [here](https://github.com/shenwei356/sun2021-cami-profiles/releases/tag/v2021-12-06).

All FASTQ files were downloaded and saved in one directory `reads`.
        
## KMCP

We search against GTDB, Genbank-viral, and Refseq-fungi respectively, and merge the results.

    # ------------------------------------------------------------------------
    # prepare folder and files
    
    # prepare folder and files.
    mkdir -p kmcp-se
    cd kmcp-se
    fd fq.gz$ ../reads | rush 'ln -s {}'
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
    
    fd left.fq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -v db=$db -v dbname=$dbname -j $j -v j=$J -v 'p={:}' \
            'kmcp search -d {db} {p}.left.fq.gz {p}.right.fq.gz -o {p}.kmcp@{dbname}.tsv.gz \
                --log {p}.kmcp@{dbname}.tsv.gz.log -j {j}' \
            -c -C $reads@$dbname.rush
 
    # -------------------------------------
    # genbank-viral

    # genbank-viral
    db=genbank-viral.kmcp/
    X=taxdump/
    T=genbank-viral.kmcp/taxid.map    
    dbname=genbank-viral
    
    fd left.fq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -v db=$db -v dbname=$dbname -j $j -v j=$J -v 'p={:}' \
            'kmcp search -d {db} {p}.left.fq.gz {p}.right.fq.gz -o {p}.kmcp@{dbname}.tsv.gz \
                --log {p}.kmcp@{dbname}.tsv.gz.log -j {j}' \
            -c -C $reads@$dbname.rush
    
    # -------------------------------------
    # refseq-fungi
    
    # refseq-fungi
    db=refseq-fungi.kmcp/
    X=taxdump/
    T=refseq-fungi.kmcp/taxid.map    
    dbname=refseq-fungi
    
    fd left.fq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -v db=$db -v dbname=$dbname -j $j -v j=$J -v 'p={:}' \
            'kmcp search -d {db} {p}.left.fq.gz {p}.right.fq.gz -o {p}.kmcp@{dbname}.tsv.gz \
                --log {p}.kmcp@{dbname}.tsv.gz.log -j {j}' \
            -c -C $reads@$dbname.rush

    
    # ------------------------------------------------------------------------
    # merge results
    
    reads=kmcp-se
    j=16
    fd left.fq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j $j -v 'p={:}' \
            'kmcp merge {p}.kmcp@*.tsv.gz -o {p}.kmcp.tsv.gz --log {p}.kmcp.tsv.gz.log'
    
    
    # ------------------------------------------------------------------------
    # [for merged search results] multiple profiling modes
     
    X=taxdump/
    cat genbank-viral.kmcp/taxid.map gtdb.kmcp/taxid.map refseq-fungi.kmcp/taxid.map > taxid.map
    cat genbank-viral.kmcp/name.map gtdb.kmcp/name.map refseq-fungi.kmcp/name.map > name.map
    T=taxid.map
    
    for m in $(seq 1 5); do
        fd kmcp.tsv.gz$ $reads/ \
            | csvtk sort -H -k 1:N \
            | rush -v X=$X -v T=$T -v m=$m \
                'kmcp profile -m {m} -X {X} -T {T} {} -o {}.k-m{m}.profile -C {}.c-m{m}.profile -s {%:} \
                    --show-rank superkingdom,phylum,class,order,family,genus,species \
                    --log {}.k-m{m}.profile.log' 
        
        profile=$reads.c-m$m.profile
        fd kmcp.tsv.gz.c-m$m.profile$ $reads/ \
            | csvtk sort -H -k 1:N \
            | rush -j 1 'cat {}' \
            > $profile
    done


    
## mOTUs
    
    # prepare folder and files.
    reads=motus-pe
    
    mkdir -p $reads
    cd $reads
    fd fq.gz$ ../reads | rush 'ln -s {}'
    cd ..
    
    reads=motus-pe
    j=4
    J=40
    
    fd left.fq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j $j -v j=$J -v 'p={:}' \
            'memusg -t -s "motus profile -n {%:} -f {p}.left.fq.gz -r {p}.right.fq.gz -t {j} \
                -p -o {p}.motus.profile -C precision " >{p}.log 2>&1 '
     
    profile=$reads.profile
    fd motus.profile$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j 1 'cat {}' \
        > $profile
        

    # ------------------------------------------------------
    # change taxonomy version
    
    taxdump=taxdump/
    fd motus.profile$ $reads/ \
        | rush -v taxdump=$taxdump \
            'grep -E "^#|@" -v {} \
                | csvtk grep -Ht -f 2 -p species \
                | csvtk cut -Ht -f 1,5 \
                | taxonkit profile2cami --data-dir {taxdump} -s {%:} \
                > {.}.new.profile'
                
    newprofile=$reads.new.profile
    fd motus.new.profile$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j 1 'cat {}' \
        > $newprofile
    
## MetaPHlAn

    # prepare folder and files.
    mkdir -p mpa3-pe
    
    cd mpa3-pe
    fd fq.gz$ ../reads | rush 'ln -s {}'
    cd ..
    
    reads=mpa3-pe
    j=4
    J=40
    
    dbdir=~/ws/db/mpa3/
    db=mpa_v30_CHOCOPhlAn_201901
    
    /bin/rm -rf $reads/*.bz2
    fd left.fq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j $j -v j=$J -v 'p={:}' -v db=$db -v dbdir=$dbdir \
            'memusg -t -s "metaphlan --input_type fastq {p}.left.fq.gz,{p}.right.fq.gz -o {p}.mpa3.profile \
                -x {db} --bowtie2db {dbdir} \
                --bowtie2out {p}.bowtie2.bz2 --nproc {j} --CAMI_format_output" >{p}.log 2>&1 '

    profile=$reads.profile
    fd mpa3.profile$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j 1 'cat {}' \
        > $profile
        

    # ------------------------------------------------------
    # change taxonomy version
    
    taxdump=taxdump/
    fd mpa3.profile$ $reads/ \
        | rush -v taxdump=$taxdump \
            'grep -E "^#|@" -v {} \
                | csvtk grep -Ht -f 2 -p species \
                | csvtk cut -Ht -f 1,5 \
                | taxonkit profile2cami --data-dir {taxdump} -s {%:} \
                > {.}.new.profile'
                
    newprofile=$reads.new.profile
    fd mpa3.new.profile$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j 1 'cat {}' \
        > $newprofile
        
## MetaPHlAn_add_viruses

    # prepare folder and files.
    mkdir -p mpa3-pe-add-viruses
    
    cd mpa3-pe-add-viruses
    fd fq.gz$ ../reads | rush 'ln -s {}'
    cd ..
    
    reads=mpa3-pe-add-viruses
    j=4
    J=40
    
    dbdir=~/ws/db/mpa3/
    db=mpa_v30_CHOCOPhlAn_201901
    
    /bin/rm -rf $reads/*.bz2
    fd left.fq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j $j -v j=$J -v 'p={:}' -v db=$db -v dbdir=$dbdir \
            'memusg -t -s "metaphlan --add_viruses --input_type fastq {p}.left.fq.gz,{p}.right.fq.gz -o {p}.mpa3.profile \
                -x {db} --bowtie2db {dbdir} \
                --bowtie2out {p}.bowtie2.bz2 --nproc {j} --CAMI_format_output" >{p}.log 2>&1 '

    profile=$reads.profile
    fd mpa3.profile$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j 1 'cat {}' \
        > $profile
        

    # ------------------------------------------------------
    # change taxonomy version
    
    taxdump=taxdump/
    fd mpa3.profile$ $reads/ \
        | rush -v taxdump=$taxdump \
            'grep -E "^#|@" -v {} \
                | csvtk grep -Ht -f 2 -p species \
                | csvtk cut -Ht -f 1,5 \
                | taxonkit profile2cami --data-dir {taxdump} -s {%:} \
                > {.}.new.profile'
                
    newprofile=$reads.new.profile
    fd mpa3.new.profile$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j 1 'cat {}' \
        > $newprofile
    
    
## Bracken

Preparing tocami.py which convert Bracken output to CAMI format

    # wget https://raw.githubusercontent.com/hzi-bifo/cami2_pipelines/master/bin/tocami.py
    chmod a+x tocami.py

    # install pacakge ete3
    pip install ete3
    
    # preparing taxdump.tar.gz for tocami.py
    tar -zcvf taxdump.tar.gz taxdump/*
    
    # creating database for ete3 (don't worry the error reports, just ignore):
    ./tocami.py -t taxdump.tar.gz -f bracken -s 1 -d . bracken-pe/Build_sample1.bracken

Steps

    # --------------------------------------------------
    # using kraken's PLUSPF database

    reads=bracken-pe
    
    # prepare folder and files.
    mkdir -p $reads
    cd $reads
    fd fq.gz$ ../reads | rush 'ln -s {}'
    cd ..

    reads=bracken-pe    
    j=4
    J=40
    
    db=/home/shenwei/ws/db/kraken/pluspf    
    readlen=150
    threshold=10
    
    
    # --------------------------------------------------
    # using kraken database built with GTDB, Genbank-viral, Refseq-fungi
    
    reads=bracken-kmcp
    
    # prepare folder and files.
    mkdir -p $reads
    cd $reads
    fd fq.gz$ ../reads | rush 'ln -s {}'
    cd ..

    reads=bracken-kmcp    
    j=4
    J=40
    
    db=/home/shenwei/ws/db/kraken/kmcp/kmcp
    readlen=150
    threshold=10
    
    
    fd left.fq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j $j -v j=$J -v 'p={:}' -v db=$db -v readlen=$readlen -v threshold=$threshold \
            'memusg -t -s \
                "kraken2 --db {db} --threads {j} --memory-mapping --gzip-compressed --paired  \
                    {p}.left.fq.gz {p}.right.fq.gz --report {p}.kreport > /dev/null; \
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

## Centrifuge

    # --------------------------------------------------
    # using centrifuge database built with GTDB, Genbank-viral, Refseq-fungi

    reads=centrifuge-pe
    
    # prepare folder and files.
    mkdir -p $reads
    cd $reads
    fd fq.gz$ ../reads | rush 'ln -s {}'
    cd ..

    reads=centrifuge-pe    
    j=4
    J=40
    
    db=~/ws/db/centrifuge/kmcp  
    

    fd left.fq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j $j -v j=$J -v 'p={:}' -v db=$db \
            'memusg -t -s \
                "centrifuge -t -p {j} --mm -x {db} \
                    -q -1 {p}.left.fq.gz -2 {p}.right.fq.gz \
                    --report-file {p}.cf-report.tsv -S {}.cf.tsv" \
                >{p}.log 2>&1 '

    # Kraken-style report
    fd left.fq.gz$ $reads/ \
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

## DUDes

    # --------------------------------------------------
    # using dudes database built with GTDB, Genbank-viral, Refseq-fungi

    reads=dudes-pe
    
    # prepare folder and files.
    mkdir -p $reads
    cd $reads
    fd fq.gz$ ../reads | rush 'ln -s {}'
    cd ..

    reads=dudes-pe    
    j=4
    J=40
    n=100 # number of alignments for a read to keep. The paper use -k60
    
    db=~/ws/db/dudes/dudes-kmcp
    db2=~/ws/db/dudes/dudes-kmcp.npz
    
    # mapping with bowtie2
    fd left.fq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j $j -v j=$J -v 'p={:}' -v db=$db -v n=$n \
            'memusg -t -s \
                "bowtie2 --mm -p {j} -x {db} --no-unal --very-fast -k {n} -q \
                    -1 {p}.left.fq.gz -2 {p}.right.fq.gz -S {p}.sam" \
                >{p}.a.log 2>&1 '
    
    # profiling with dudes
    j=25
    J=1
    fd left.fq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j $j -v j=$J -v 'p={:}' -v db=$db2 -v n=$n \
            'memusg -t -s \
                "DUDes.py -t {j} -s {p}.sam -d {db} -o {p}" \
                >{p}.b.log 2>&1 '

    # rename sample id
    fd left.fq.gz$ $reads/ \
        | rush 'sed -i "s/{/}\///; s/.sam$//;" {:}.out'
                
    profile=$reads.profile
    fd .out$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j 1 'cat {}' \
        > $profile

## SIMM

    # --------------------------------------------------
    # using slimm database built with GTDB, Genbank-viral, Refseq-fungi

    reads=slimm-pe
    
    # prepare folder and files.
    mkdir -p $reads
    cd $reads
    fd fq.gz$ ../reads | rush 'ln -s {}'
    cd ..

    reads=slimm-pe    
    j=4
    J=40
    n=100 # number of alignments for a read to keep
    
    db=~/ws/db/slimm/dudes-kmcp    
    db2=~/ws/db/slimm/slimm-kmcp.sldb
    
    # mapping with bowtie2
    fd left.fq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j $j -v j=$J -v 'p={:}' -v db=$db -v n=$n \
            'memusg -t -s \
                "bowtie2 --mm -p {j} -x {db} --no-unal --very-fast -k {n} -q \
                    -1 {p}.left.fq.gz -2 {p}.right.fq.gz -S {p}.sam" \
                >{p}.a.log 2>&1 '
    
    # profiling with slimm
    j=25
    J=1
    fd left.fq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j $j -v j=$J -v 'p={:}' -v db=$db2 -v n=$n \
            'memusg -t -s \
                "slimm -w 1000 -o {p} {db} {p}.sam" \
                >{p}.b.log 2>&1 '

    # ------------------------------------------------------
    # convert profile table to cami format
    
    taxdump=taxdump/
    fd left.fq.gz$ $reads/ \
        | rush -v taxdump=$taxdump \
            'csvtk grep -t -f taxa_level -p species {:}_profile.tsv \
                | csvtk cut -t -f taxa_id,abundance \
                | csvtk grep -t -f taxa_id -r -p "\*$" -v \
                | sed 1d \
                | taxonkit profile2cami --data-dir {taxdump} -s {%:} \
                | taxonkit cami-filter \
                > {:}.profile'
                
    newprofile=$reads.profile
    fd profile$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j 1 'cat {}' \
        > $newprofile
    
