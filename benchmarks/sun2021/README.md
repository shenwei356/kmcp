# Benchmarks on 25 simulated communities from Sun et al

## Softwares

- KMCP ([v0.7.0](https://github.com/shenwei356/kmcp/releases/tag/v0.7.0))
- mOTUs [3.0.1 (Jul 28, 2021)](https://github.com/motu-tool/mOTUs/releases/tag/3.0.1)
- MetaPhlAn [3.0.13 (27 Jul, 2021)](https://github.com/biobakery/MetaPhlAn/releases/tag/3.0.13)
- Kraken [v2.1.2](https://github.com/DerrickWood/kraken2/releases/tag/v2.1.2),
  Bracken [v2.6.2](https://github.com/jenniferlu717/Bracken/releases/tag/v2.6.2)

## Databases and taxonomy version

- KMCP,  GTDB-RS202 (2021-04-27), 2021-10
- mOTUs, 3.0.1 (2021-06-28), 2019-01
- MetaPhlAn, mpa_v30_CHOCOPhlAn_201901 (?), 2019-01
- Kraken, PlusPF (2021-05-17), 2021-05-17

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
I've convert the ground truth profiles to CAMI format with Taxonomy version 2021-10-06,
available [here](https://github.com/shenwei356/sun2021-cami-profiles/releases/tag/v2021-12-06).

All FASTQ files were downloaded and saved in one directory `reads`.
        
## KMCP

        
    db=gtdb.kmcp/
    X=taxdump/
    T=gtdb.kmcp/taxid.map    
    dbname=gtdb
    
    
    # ------------------------------------------------------------------------
    # single-end mode    
    
    # prepare folder and files.
    mkdir -p kmcp-se
    cd kmcp-se
    fd fq.gz$ ../reads | rush 'ln -s {}'
    cd ..

    
    reads=kmcp-se
    j=4
    J=40
    
    fd left.fq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -v db=$db -v dbname=$dbname -j $j -v j=$J -v 'p={:}' \
            'kmcp search -d {db} {p}.left.fq.gz {p}.right.fq.gz -o {p}.kmcp@{dbname}.tsv.gz \
                --log {p}.kmcp@{dbname}.tsv.gz.log -j {j}' \
            -c -C $reads@$dbname.rush
            
            
    # ------------------------------------------------------------------------           
    # paired-end mode
    
    # prepare folder and files.
    mkdir -p kmcp-pe
    cd kmcp-pe
    fd fq.gz$ ../reads | rush 'ln -s {}'
    cd ..
    
    
    reads=kmcp-pe
    j=4
    J=40
    
    fd left.fq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -v db=$db -v dbname=$dbname -j $j -v j=$J -v 'p={:}' \
            'kmcp search -d {db} -1 {p}.left.fq.gz -2 {p}.right.fq.gz -o {p}.kmcp@{dbname}.tsv.gz \
                --log {p}.kmcp@{dbname}.tsv.gz.log -j {j}' \
            -c -C $reads@$dbname.rush
    
    
    # ------------------------------------------------------------------------
    # default profiling mode
    
    fd kmcp@$dbname.tsv.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -v db=$db -v dbname=$dbname -v X=$X -v T=$T \
            'kmcp profile -X {X} -T {T} {} -o {}.k.profile -C {}.c.profile -s {%:} --log {}.k.profile.log' 
    
    profile=$reads@$dbname.c.profile
    fd kmcp@$dbname.tsv.gz.c.profile$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j 1 'cat {}' \
        > $profile
        
    # ------------------------------------------------------------------------
    # multiple profiling modes
        
    for m in $(seq 1 5); do
        fd kmcp@$dbname.tsv.gz$ $reads/ \
            | csvtk sort -H -k 1:N \
            | rush -v db=$db -v dbname=$dbname -v X=$X -v T=$T -v m=$m \
                'kmcp profile -m {m} -X {X} -T {T} {} -o {}.k-m{m}.profile -C {}.c-m{m}.profile -s {%:} --log {}.k-m{m}.profile.log' 
        
        profile=$reads@$dbname.c-m$m.profile
        fd kmcp@$dbname.tsv.gz.c-m$m.profile$ $reads/ \
            | csvtk sort -H -k 1:N \
            | rush -j 1 'cat {}' \
            > $profile
    done
    
## mOTUs
    
    # prepare folder and files.
    mkdir -p motus-pe
    cd motus-pe
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

    # prepare folder and files.
    mkdir -p bracken-pe
    cd bracken-pe
    fd fq.gz$ ../reads | rush 'ln -s {}'
    cd ..
    
    reads=bracken-pe
    j=4
    J=40
    
    db=/home/shenwei/ws/db/kraken/
    readlen=150
    threshold=10
    
    fd left.fq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j $j -v j=$J -v 'p={:}' -v db=$db -v db=$db -v readlen=$readlen -v level=$level -v threshold=$threshold \
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
        | rush 'python3 ./tocami.py -d ./ -f bracken {} -s {%:} -o {}.profile'
    
    profile=$reads.profile
    fd bracken.profile$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j 1 'cat {}' \
        > $profile
