# Benchmarking analysis time and resource requirement

## Softwares

- KMCP [v0.9.0](https://github.com/shenwei356/kmcp/releases/tag/v0.9.0)
- mOTUs [3.0.1 (2021-07-28)](https://github.com/motu-tool/mOTUs/releases/tag/3.0.1)
- MetaPhlAn [3.0.13 (2021-07-27)](https://github.com/biobakery/MetaPhlAn/releases/tag/3.0.13)
- Kraken [v2.1.2 (2021-05-10)](https://github.com/DerrickWood/kraken2/releases/tag/v2.1.2),
  Bracken [v2.6.2 (2021-03-22)](https://github.com/jenniferlu717/Bracken/releases/tag/v2.6.2)
- Centrifuge [v1.0.4 (2021-08-17)](https://github.com/DaehwanKimLab/centrifuge/releases/tag/v1.0.4)
- DUDes [v0.08 (2017-11-08)](https://github.com/pirovc/dudes/releases/tag/dudes_v0.08)
- SLIMM [v0.3.4 (2018-09-04)](https://github.com/seqan/slimm/releases/tag/v0.3.4)
- Ganon [1.1.2](https://github.com/pirovc/ganon/releases/tag/1.1.2)

## Databases and taxonomy version

- KMCP,  GTDB-RS202 (2021-04-27) + Genbank-viral (r246, 2021-12-06) + Refseq-fungi (r208, 2021-09-30), 2021-12-06
- Centrifuge, built with the genomes same to KMCP.
- mOTUs, 3.0.1 (2021-06-28), 2019-01
- MetaPhlAn, mpa_v30_CHOCOPhlAn_201901 (?), 2019-01
- Kraken, PlusPF (2021-05-17), 2021-05-17
- Kraken, built with the genomes same to KMCP.
- DUDes, built with the genomes same to KMCP.
- SLIMM, built with the genomes same to KMCP.
- Ganon, built with the genomes same to KMCP.

## Database of KMCP (local machine)

GTDB

    input=gtdb
    
    memusg -t -s "kmcp compute -I $input -O gtdb-r202-k21-n10 -k 21 -n 10 -l 150 -B plasmid \
        --log gtdb-r202-k21-n10.log0 -j 40 --force" > gtdb-r202-k21-n10.log 2>&1

    memusg -t -s "kmcp index -j 40 -I gtdb-r202-k21-n10 -O gtdb.kmcp -n 1 -f 0.3 \
        --log gtdb.kmcp.log0 --force" > gtdb.kmcp.log 2>&1
    
    cp taxid.map name.map gtdb.kmcp/

Genbank Viral

    name=viral    
    input=files.renamed.slim
    
    memusg -t -s "kmcp compute -I $input -O genbank-$name-k21-n10 \
        -k 21 --seq-name-filter plasmid \
        --split-number 10 --split-overlap 150 \
        --log genbank-$name-k21-n10.log0 -j 40 --force " > genbank-$name-k21-n10.log 2>&1
    memusg -t -s "kmcp index -I genbank-$name-k21-n10/ -O genbank-viral.kmcp \
        -j 40 -f 0.05 -n 1 -x 100K -8 1M \
        --log genbank-viral.kmcp.log0 --force" > genbank-viral.kmcp.log 2>&1
    
    cp taxid.map name.map genbank-viral.kmcp/
    
Refseq fungi

    name=fungi    
    input=files.renamed
    
    memusg -t -s "kmcp compute -I $input -O refseq-$name-k21-n10 \
        -k 21 --seq-name-filter plasmid \
        --split-number 10 --split-overlap 150 \
        --log refseq-$name-k21-n10.log0 -j 40 --force" > refseq-$name-k21-n10.log 2>&1
      
    memusg -t -s "kmcp index -I refseq-$name-k21-n10/ -O refseq-fungi.kmcp \
        -j 40 -f 0.3 -n 1 \
        --log refseq-fungi.kmcp.log0 --force" > refseq-fungi.kmcp.log 2>&1
    
    cp taxid.map name.map refseq-fungi.kmcp/
    
## Database of KMCP (HPC)
    
https://bioinf.shenwei.me/kmcp/database/
    
## Datasets

Eight simulated short-read mouse gut metagenome datasets from the CAMI challenge, 
with 5Gb 150-bp paired-end reads per sample, are used to assess analysis time and resource requirement.

    # linking files
    mkdir reads0; cd reads0
    fd fq.gz$ ~/ws/data/cami2/mouse-gut/short-reads/19122017_mousegut_scaffolds/ \
        | csvtk sort -Ht -k 1:N \
        | head -n 8 \
        | rush 'ln -s {}'
    cd ..
    
    # split merged paired-end reads.
    mkdir reads
    cd reads
    fd .fq.gz$ ../reads0 | rush -j 8 'seqkit split2 -p 2 {} -O .'
    brename -p .part_00 -r _
    brename -p '_1\.fq.gz' -r .left.fq.gz
    brename -p '_2\.fq.gz' -r .right.fq.gz
    cd ..

    
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
            'memusg -t -s "kmcp search -d {db} {p}.left.fq.gz {p}.right.fq.gz -o {p}.kmcp@{dbname}.tsv.gz \
                --log {p}.kmcp@{dbname}.tsv.gz.log0 -j {j}" > {p}.kmcp@{dbname}.tsv.gz.log 2>&1' \
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
            'memusg -t -s "kmcp search -d {db} {p}.left.fq.gz {p}.right.fq.gz -o {p}.kmcp@{dbname}.tsv.gz \
                --log {p}.kmcp@{dbname}.tsv.gz.log0 -j {j}" > {p}.kmcp@{dbname}.tsv.gz.log 2>&1' \
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
            'memusg -t -s "kmcp search -d {db} {p}.left.fq.gz {p}.right.fq.gz -o {p}.kmcp@{dbname}.tsv.gz \
                --log {p}.kmcp@{dbname}.tsv.gz.log0 -j {j}" > {p}.kmcp@{dbname}.tsv.gz.log 2>&1' \
            -c -C $reads@$dbname.rush

    
    # ------------------------------------------------------------------------
    # merge results
    
    reads=kmcp-se
    j=16
    fd left.fq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j $j -v 'p={:}' \
            'memusg -t -s "kmcp merge {p}.kmcp@*.tsv.gz -o {p}.kmcp.tsv.gz \
                --log {p}.kmcp.tsv.gz.log0" > {p}.kmcp.tsv.gz.log 2>&1'
    
    
    # ------------------------------------------------------------------------
    # [for merged search results] multiple profiling modes
     
    X=taxdump/
    cat genbank-viral.kmcp/taxid.map gtdb.kmcp/taxid.map refseq-fungi.kmcp/taxid.map > taxid.map
    cat genbank-viral.kmcp/name.map gtdb.kmcp/name.map refseq-fungi.kmcp/name.map > name.map
    T=taxid.map
    
    fd kmcp.tsv.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -v X=$X -v T=$T \
            'memusg -t -s "kmcp profile -X {X} -T {T} {} -o {}.k.profile -C {}.c.profile -s {%:} \
                --log {}.k.profile.log0" > {}.k.profile.log 2>&1 ' 

          
    # ------------------------------------------------------------------------
    # time and peak memory
    
    # retrive time and memory from the log files
    stats(){
        app=$1
        dir=$2
        
        fd .log$ $dir \
            | rush -k -v n=$app \
                'echo -ne "{n}\t"{%:}; \
                echo -en "\t"$(tail -n 3 {} | grep "elapsed time" | sed -r "s/.+: //"); \
                echo -en "\t"$(tail -n 3 {} | grep "peak rss" | sed -r "s/.+: //")"\n"; ' \
            | csvtk add-header -Ht -n app,sample,time0,memory \
            | csvtk mutate -t -f 3 -n h --na -p '([\d\.]+)h' \
            | csvtk replace -t -f h -p '^$' -r 0 \
            | csvtk mutate -t -f 3 -n m --na -p '([\d\.]+)m' \
            | csvtk replace -t -f m -p '^$' -r 0 \
            | csvtk mutate -t -f 3 -n s -p '([\d\.]+)s' \
            | csvtk mutate2 -t -n time -e '$h*3600 + $m*60 + $s' -w 0 \
            | csvtk mutate -t -f 4 -n gb --na -p '([\d\.]+) GB' \
            | csvtk replace -t -f gb -p '^$' -r 0 \
            | csvtk mutate -t -f 4 -n mb --na -p '([\d\.]+) MB' \
            | csvtk replace -t -f mb -p '^$' -r 0 \
            | csvtk mutate2 -t -n mem -e '$gb + $mb/1024' -w 1 \
            | csvtk cut -t -f app,sample,time,mem \
            | csvtk summary -t -g app,sample -f time:sum -f mem:max \
            | csvtk rename -t -f time:sum -n time \
            | csvtk rename -t -f mem:max -n mem
    }
    
    stats kmcp kmcp-se > kmcp.stats
    
## KMCP in clusters

    # ------------------------------------------------------------------------
    # prepare folder and files
    
    # prepare folder and files.
    mkdir -p kmcp-se-cluster
    cd kmcp-se-cluster
    fd fq.gz$ ../reads | rush 'ln -s {}'
    cd ..

    # ------------------------------------
    # searching
    
    
    j=32
    reads=kmcp-se-cluster
    
    # -----------------
    # gtdb
    
    dbprefix=~/ws/db/kmcp/gtdb.n16-
    
    for file in $reads/*.left.fq.gz; do
        prefix=$(echo $file | sed 's/.left.fq.gz//')
        read1=$file
        read2=$(echo $file | sed 's/left.fq.gz/right.fq.gz/')
        
        ls -d $dbprefix*.kmcp \
            | easy_sbatch \
                -c $j -J $(basename $prefix) \
                "memusg -t -s 'kmcp search   \
                            --load-whole-db  \
                            --threads   $j   \
                            --db-dir    {}   \
                            $read1 $read2    \
                            --out-file  $prefix.kmcp@\$(basename {}).tsv.gz \
                            --log       $prefix.kmcp@\$(basename {}).tsv.gz.log0 \
                            --quiet ' \
                > $prefix.kmcp@\$(basename {}).tsv.gz.log 2>&1"
    done
    
    # -----------------
    # viral
    
    dbprefix=~/ws/db/kmcp/genbank-viral.n4-
    
    for file in $reads/*.left.fq.gz; do
        prefix=$(echo $file | sed 's/.left.fq.gz//')
        read1=$file
        read2=$(echo $file | sed 's/left.fq.gz/right.fq.gz/')
        
        ls -d $dbprefix*.kmcp \
            | easy_sbatch \
                -c $j -J $(basename $prefix) \
                "memusg -t -s 'kmcp search   \
                            --load-whole-db  \
                            --threads   $j   \
                            --db-dir    {}   \
                            $read1 $read2    \
                            --out-file  $prefix.kmcp@\$(basename {}).tsv.gz \
                            --log       $prefix.kmcp@\$(basename {}).tsv.gz.log0 \
                            --quiet ' \
                > $prefix.kmcp@\$(basename {}).tsv.gz.log 2>&1"
    done
    
    
    # -----------------
    # fungi
    
    dbprefix=~/ws/db/kmcp/refseq-fungi

    for file in $reads/*.left.fq.gz; do
        prefix=$(echo $file | sed 's/.left.fq.gz//')
        read1=$file
        read2=$(echo $file | sed 's/left.fq.gz/right.fq.gz/')
        
        ls -d $dbprefix*.kmcp \
            | easy_sbatch \
                -c $j -J $(basename $prefix) \
                "memusg -t -s 'kmcp search   \
                            --load-whole-db  \
                            --threads   $j   \
                            --db-dir    {}   \
                            $read1 $read2    \
                            --out-file  $prefix.kmcp@\$(basename {}).tsv.gz \
                            --log       $prefix.kmcp@\$(basename {}).tsv.gz.log0 \
                            --quiet ' \
                > $prefix.kmcp@\$(basename {}).tsv.gz.log 2>&1"
    done
    

    # ----------------------------------
    # wait all job being done
    
    
    
    # ----------------------------------
    # merge result and profiling
    
    # merge results
    # there's no need to submit to slurm, which could make it slower, cause the bottleneck is file IO
    for file in $reads/*.left.fq.gz; do
        prefix=$(echo $file | sed 's/.left.fq.gz//')        
        
        echo $prefix; date
        memusg -t -s "kmcp merge $prefix.kmcp@*.tsv.gz --out-file $prefix.kmcp.tsv.gz \
            --quiet --log $prefix.kmcp.tsv.gz.merge.log0" > $prefix.kmcp.tsv.gz.merge.log 2>&1
    done
    
    # profiling
    X=taxdump/
    T=taxid.map
    
    fd kmcp.tsv.gz$ $reads/ \
        | grep @ -v \
        | rush -v X=$X -v T=$T \
            'memusg -t -s "kmcp profile -X {X} -T {T} {} -o {}.k.profile -C {}.c.profile -s {%:} \
                --quiet --log {}.k.profile.log0" > {}.k.profile.log 2>&1 ' 

    # ----------------------------------
    # wait all job being done

    
    
    # ----------------------------------
    # time and peak memory
    
    app=kmcp_c
    
    # 1) search
    fd .log$ $reads \
        | grep @ \
        | grep merge -v \
        | grep profile -v \
        | rush -k -v n=$app \
            'echo -ne "{n}\t"{%:}; \
            echo -en "\t"$(tail -n 3 {} | grep "elapsed time" | sed -r "s/.+: //"); \
            echo -en "\t"$(tail -n 3 {} | grep "peak rss" | sed -r "s/.+: //")"\n"; ' \
        | csvtk add-header -Ht -n app,sample,time0,memory \
        | csvtk mutate -t -f 3 -n h --na -p '([\d\.]+)h' \
        | csvtk replace -t -f h -p '^$' -r 0 \
        | csvtk mutate -t -f 3 -n m --na -p '([\d\.]+)m' \
        | csvtk replace -t -f m -p '^$' -r 0 \
        | csvtk mutate -t -f 3 -n s -p '([\d\.]+)s' \
        | csvtk mutate2 -t -n time -e '$h*3600 + $m*60 + $s' -w 0 \
        | csvtk mutate -t -f 4 -n gb --na -p '([\d\.]+) GB' \
        | csvtk replace -t -f gb -p '^$' -r 0 \
        | csvtk mutate -t -f 4 -n mb --na -p '([\d\.]+) MB' \
        | csvtk replace -t -f mb -p '^$' -r 0 \
        | csvtk mutate2 -t -n mem -e '$gb + $mb/1024' -w 1 \
        | csvtk cut -t -f app,sample,time,mem \
        | csvtk summary -t -g app,sample -f time:max -f mem:max \
        | csvtk rename -t -f time:max -n time \
        | csvtk rename -t -f mem:max -n mem \
        > $reads/$app.search.stats
    
    # 2) merge
    fd .log$ $reads \
        | grep merge  \
        | rush -k -v n=$app \
            'echo -ne "{n}\t"{%:}; \
            echo -en "\t"$(tail -n 3 {} | grep "elapsed time" | sed -r "s/.+: //"); \
            echo -en "\t"$(tail -n 3 {} | grep "peak rss" | sed -r "s/.+: //")"\n"; ' \
        | csvtk add-header -Ht -n app,sample,time0,memory \
        | csvtk mutate -t -f 3 -n h --na -p '([\d\.]+)h' \
        | csvtk replace -t -f h -p '^$' -r 0 \
        | csvtk mutate -t -f 3 -n m --na -p '([\d\.]+)m' \
        | csvtk replace -t -f m -p '^$' -r 0 \
        | csvtk mutate -t -f 3 -n s -p '([\d\.]+)s' \
        | csvtk mutate2 -t -n time -e '$h*3600 + $m*60 + $s' -w 0 \
        | csvtk mutate -t -f 4 -n gb --na -p '([\d\.]+) GB' \
        | csvtk replace -t -f gb -p '^$' -r 0 \
        | csvtk mutate -t -f 4 -n mb --na -p '([\d\.]+) MB' \
        | csvtk replace -t -f mb -p '^$' -r 0 \
        | csvtk mutate2 -t -n mem -e '$gb + $mb/1024' -w 1 \
        | csvtk cut -t -f app,sample,time,mem \
        | csvtk summary -t -g app,sample -f time:max -f mem:max \
        | csvtk rename -t -f time:max -n time \
        | csvtk rename -t -f mem:max -n mem \
        > $reads/$app.merge.stats
    
    # 3) profile
    fd .log$ $reads \
        | grep profile  \
        | rush -k -v n=$app \
            'echo -ne "{n}\t"{%:}; \
            echo -en "\t"$(tail -n 3 {} | grep "elapsed time" | sed -r "s/.+: //"); \
            echo -en "\t"$(tail -n 3 {} | grep "peak rss" | sed -r "s/.+: //")"\n"; ' \
        | csvtk add-header -Ht -n app,sample,time0,memory \
        | csvtk mutate -t -f 3 -n h --na -p '([\d\.]+)h' \
        | csvtk replace -t -f h -p '^$' -r 0 \
        | csvtk mutate -t -f 3 -n m --na -p '([\d\.]+)m' \
        | csvtk replace -t -f m -p '^$' -r 0 \
        | csvtk mutate -t -f 3 -n s -p '([\d\.]+)s' \
        | csvtk mutate2 -t -n time -e '$h*3600 + $m*60 + $s' -w 0 \
        | csvtk mutate -t -f 4 -n gb --na -p '([\d\.]+) GB' \
        | csvtk replace -t -f gb -p '^$' -r 0 \
        | csvtk mutate -t -f 4 -n mb --na -p '([\d\.]+) MB' \
        | csvtk replace -t -f mb -p '^$' -r 0 \
        | csvtk mutate2 -t -n mem -e '$gb + $mb/1024' -w 1 \
        | csvtk cut -t -f app,sample,time,mem \
        | csvtk summary -t -g app,sample -f time:max -f mem:max \
        | csvtk rename -t -f time:max -n time \
        | csvtk rename -t -f mem:max -n mem \
        > $reads/$app.profile.stats
    
    # sum up
    csvtk concat $reads/$app.*.stats \
        | csvtk summary -t -g app,sample -f time:sum -f mem:max \
        | csvtk rename -t -f time:sum -n time \
        | csvtk rename -t -f mem:max -n mem \
        > $app.stats
    
    
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
   
    stats motus motus-pe > motus.stats
                
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
    
    stats metaphlan mpa3-pe > metaphlan.stats
   
## Bracken

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
    
    stats bracken-std bracken-pe > bracken-std.stats
    
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

    stats bracken-kmcp bracken-kmcp > bracken-kmcp.stats

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
                    --report-file {p}.cf-report.tsv -S {}.cf.tsv; \
                centrifuge-kreport -x {db} {}.cf.tsv > {}.kreport " \
                >{p}.log 2>&1 '

    stats centrifuge centrifuge-pe > centrifuge.stats

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

    stats dudes dudes-pe > dudes.stats
    
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

    stats slimm slimm-pe > slimm.stats

## ganon

    # --------------------------------------------------
    # using ganon database built with GTDB, Genbank-viral, Refseq-fungi

    reads=ganon-pe
    
    # prepare folder and files.
    mkdir -p $reads
    cd $reads
    fd fq.gz$ ../reads | rush 'ln -s {}'
    cd ..

    reads=ganon-pe    
    j=4
    J=40
    
    db=~/ws/db/ganon/ganon-kmcp    
    
    fd left.fq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j $j -v j=$J -v 'p={:}' -v db=$db \
            'memusg -t -s \
                "ganon classify -d {db} -t {j} \
                    -p {p}.left.fq.gz {p}.right.fq.gz -o {p}; \
                ganon table -l percentage --header taxid -r species -i {p}.tre -o {p}.tsv " \
                >{p}.log 2>&1 '

    stats ganon ganon-pe > ganon.stats
    
