# Taxonomy profiling accuracies with different parameters

## Databases and reads

Reference genomes (Bacteria and Archaea):

1. Microbial genomes were extracted from CAMI2 RefSeq snapshot (`2019-01-08`) using
corresponding taxonomy information .
2. **For every species, at most 5 assemblies (sorted by assembly accession) were kept**.

Only reads of the sample 0-7 were used for benchmark.

Building databases:

    genomes=refseq-cami2-slim
    prefix=refseq-cami2
    
    for chunks in 1 5 10 20; do
        j=40    
        k=21
        kmcp compute -j $j -I $genomes/ -O $prefix-k$k-n$chunks -k $k -n $chunks -l 150 -B plasmid \
            --log $prefix-k$k-n$chunks.log --force
            
        n=1
        f=0.3
        kmcp index -I $prefix-k$k-n$chunks/ -O $prefix-k$k-n$chunks.db -j $j -n $n -f $f \
            --log $prefix-k$k-n$chunks.db.log --force
        
        # remove k-mer files
        /bin/rm -rf $prefix-k$k-n$chunks/
    done


## The number of genome chunks

Searching and profiling

    # Searching --------------------

    for chunks in 1 5 10 20; do
        db=refseq-cami2-k21-n$chunks.db
        dbname=refseq-cami2-k21-n$chunks
        
        reads=single
        j=4
        J=40
        
        fd fq.gz$ $reads/ \
            | csvtk sort -H -k 1:N \
            | rush -v db=$db -v dbname=$dbname -j $j -v j=$J \
                'memusg -H -t -s "kmcp search -j {j} -w -d {db} {} -o {}.kmcp@{dbname}.tsv.gz \
                    --log {}.kmcp@{dbname}.tsv.gz.log0" > {}.kmcp@{dbname}.tsv.gz.log 2>&1' \
                -c -C $reads@$dbname.rush
    done
    
    # Profiling --------------------

    for chunks in 1 5 10 20; do
        db=refseq-cami2-k21-n$chunks.db
        dbname=refseq-cami2-k21-n$chunks
        
        reads=single
        
        X=taxdump
        T=taxid.map
    
        fd kmcp@$dbname.tsv.gz$ $reads/ \
            | csvtk sort -H -k 1:N \
            | rush -v db=$db -v dbname=$dbname -v X=$X -v T=$T \
                'memusg -H -t -s "kmcp profile -X {X} -T {T} {} -o {}.k.profile -C {}.c.profile -s {@sample_(\d+)} \
                    --log {}.k.profile.log0" > {}.k.profile.log 2>&1' 
        
        profile=$reads@$dbname.c.profile
        fd kmcp@$dbname.tsv.gz.c.profile$ $reads/ \
            | csvtk sort -H -k 1:N \
            | rush -j 1 'cat {}' \
            > $profile
    done

Metrics of profiling accuracy

    # OPAL
    opal.py -g cami2_mouse_gut_gs.profile \
        single@refseq-cami2-k21-n1.c.profile \
        single@refseq-cami2-k21-n5.c.profile \
        single@refseq-cami2-k21-n10.c.profile \
        single@refseq-cami2-k21-n20.c.profile \
        -l chunks-1,chunks-5,chunks-10,chunks-20 \
        -o opal

    cat opal/results.tsv \
        | csvtk grep -t -f tool -p 'Gold standard' -v \
        | csvtk grep -t -f rank -p na -p species -p genus \
        | csvtk grep -t -f metric -p 'Completeness' -p 'Purity' -p 'F1 score' \
            -p 'L1 norm error' -p 'Weighted UniFrac error' \
        | tee accuracy.tsv \
        | csvtk summary -t -g tool,metric,rank -f value:mean -w 6 \
        | csvtk rename -t -f value:mean -n value \
        | csvtk sort -t -k rank -k metric -k tool:N \
        | csvtk grep -t -f rank -p genus -v \
        | csvtk csv2md -t

|tool     |rank   |metric                |value   |
|:--------|:------|:---------------------|:-------|
|chunks-1 |na     |Weighted UniFrac error|0.334769|
|chunks-5 |na     |Weighted UniFrac error|0.336242|
|chunks-10|na     |Weighted UniFrac error|0.337509|
|chunks-20|na     |Weighted UniFrac error|0.339137|
|chunks-1 |species|Completeness          |0.725203|
|chunks-5 |species|Completeness          |0.713924|
|chunks-10|species|Completeness          |0.698402|
|chunks-20|species|Completeness          |0.671732|
|chunks-1 |species|F1 score              |0.629712|
|chunks-5 |species|F1 score              |0.745917|
|chunks-10|species|F1 score              |0.771154|
|chunks-20|species|F1 score              |0.773635|
|chunks-1 |species|L1 norm error         |0.520477|
|chunks-5 |species|L1 norm error         |0.525122|
|chunks-10|species|L1 norm error         |0.524143|
|chunks-20|species|L1 norm error         |0.527089|
|chunks-1 |species|Purity                |0.570523|
|chunks-5 |species|Purity                |0.795017|
|chunks-10|species|Purity                |0.869483|
|chunks-20|species|Purity                |0.920199|

Computing time (searching + profiling) and peak memory occupation.

    # retrive time (Minute) and peak memory (MB) from the log files
    time_mem(){
        dir=$1
        pattern=$2
        name=$3
        
        fd .log$ $dir \
            | grep $pattern \
            | rush -k -v n=$name \
                'echo -ne "{n}\t"{%:}; \
                echo -en "\t"$(tail -n 3 {} | grep "elapsed time" | sed -r "s/.+: //"); \
                echo -en "\t"$(tail -n 3 {} | grep "peak rss" | sed -r "s/.+: //")"\n"; ' \
            | csvtk add-header -Ht -n name,sample,time0,mem0\
            | csvtk mutate2 -t -n time -e '$time0/60' -w 1 \
            | csvtk mutate2 -t -n mem -e '$mem0/1024' -w 1 \
            | csvtk cut -t -f name,sample,time,mem \
            | csvtk summary -t -g name,sample -f time:sum -f mem:max \
            | csvtk rename -t -f time:sum -n time \
            | csvtk rename -t -f mem:max -n mem
    }
    
    for chunks in 1 5 10 20; do
        time_mem single/ "n$chunks\." chunks-$chunks > stats_chunks-$chunks.tsv
    done
    
    csvtk concat -t stats_chunks-*.tsv \
        | csvtk sort -t -k name:N \
        | tee stats_chunks.tsv \
        | csvtk summary -t -g name -f time:mean -f time:stdev -f mem:mean -f mem:stdev \
        | csvtk sort -t -k name:N \
        | csvtk csv2md -t
        
|name     |time:mean|time:stdev|mem:mean|mem:stdev|
|:--------|:--------|:---------|:-------|:--------|
|chunks-1 |19.91    |4.68      |26528.66|1622.59  |
|chunks-5 |42.25    |5.67      |22907.95|2940.77  |
|chunks-10|71.33    |7.09      |19954.61|3045.18  |
|chunks-20|127.91   |11.71     |17577.61|1849.73  |

## Paired-end and single-end reads

Single-end reads

    # re-use the search and profiling results
    mkdir -p single-n10; cd single-n10
    fd n10 ../single | rush 'ln -s {}'
    cd ..
    
    
    # Profiling  --------------------
    
    chunks=10
    db=refseq-cami2-k21-n$chunks.db
    dbname=refseq-cami2-k21-n$chunks
    
    reads=single-n10
    
    profile=$reads@$dbname.c.profile
    fd kmcp@$dbname.tsv.gz.c.profile$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j 1 'cat {}' \
        > $profile

Paired-end reads

    # Searching --------------------
    
    chunks=10
    db=refseq-cami2-k21-n$chunks.db
    dbname=refseq-cami2-k21-n$chunks
    
    reads=paired
    j=4
    J=40
    
    fd _1.fq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -v db=$db -v dbname=$dbname -j $j -v j=$J -v 'p={@(.+)_1.fq.gz}' \
            'memusg -H -t -s "kmcp search -j {j} -w -d {db} -1 {} -2 {p}_2.fq.gz -o {p}.kmcp@{dbname}.tsv.gz \
                --log {p}.kmcp@{dbname}.tsv.gz.log0" > {p}.kmcp@{dbname}.tsv.gz.log 2>&1' \
            -c -C $reads@$dbname.rush
    
    
    # Profiling --------------------
    
    reads=paired
    
    X=taxdump
    T=taxid.map

    fd kmcp@$dbname.tsv.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -v db=$db -v dbname=$dbname -v X=$X -v T=$T \
            'memusg -H -t -s "kmcp profile -X {X} -T {T} {} -o {}.k.profile -C {}.c.profile -s {@sample_(\d+)} \
                --log {}.k.profile.log0" > {}.k.profile.log 2>&1' 
    
    profile=$reads@$dbname.c.profile
    fd kmcp@$dbname.tsv.gz.c.profile$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j 1 'cat {}' \
        > $profile

Paired-end reads (with flag `--try-se`: if paired-end reads have no hits, re-search with read1, if still
fails, try read2)

    # Searching --------------------
    
    chunks=10
    db=refseq-cami2-k21-n$chunks.db
    dbname=refseq-cami2-k21-n$chunks
    
    reads=paired-se
    j=4
    J=40
    
    fd _1.fq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -v db=$db -v dbname=$dbname -j $j -v j=$J -v 'p={@(.+)_1.fq.gz}' \
            'memusg -H -t -s "kmcp search -j {j} -w -d {db} -1 {} -2 {p}_2.fq.gz -o {p}.kmcp@{dbname}.tsv.gz \
                --try-se \
                --log {p}.kmcp@{dbname}.tsv.gz.log0" > {p}.kmcp@{dbname}.tsv.gz.log 2>&1' \
            -c -C $reads@$dbname.rush
    
    
    # Profiling --------------------
    
    reads=paired-se
    
    X=taxdump
    T=taxid.map

    fd kmcp@$dbname.tsv.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -v db=$db -v dbname=$dbname -v X=$X -v T=$T \
            'memusg -H -t -s "kmcp profile -X {X} -T {T} {} -o {}.k.profile -C {}.c.profile -s {@sample_(\d+)} \
                --log {}.k.profile.log0" > {}.k.profile.log 2>&1' 
    
    profile=$reads@$dbname.c.profile
    fd kmcp@$dbname.tsv.gz.c.profile$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j 1 'cat {}' \
        > $profile

        
Metrics of profiling accuracy

    # OPAL
    opal.py -g cami2_mouse_gut_gs.profile \
        single-n10@refseq-cami2-k21-n10.c.profile \
        paired@refseq-cami2-k21-n10.c.profile \
        paired-se@refseq-cami2-k21-n10.c.profile \
        -l "SE,PE,PE/SE" \
        -o opal

    cat opal/results.tsv \
        | csvtk grep -t -f tool -p 'Gold standard' -v \
        | csvtk grep -t -f rank -p na -p species -p genus \
        | csvtk grep -t -f metric -p 'Completeness' -p 'Purity' -p 'F1 score' \
            -p 'L1 norm error' -p 'Weighted UniFrac error' \
        | tee accuracy.tsv \
        | csvtk summary -t -g tool,metric,rank -f value:mean -w 6 \
        | csvtk rename -t -f value:mean -n value \
        | csvtk sort -t -k rank -k metric -k tool:N \
        | csvtk grep -t -f rank -p genus -v \
        | csvtk csv2md -t


## Profiling modes

    # re-use the search results
    
    mkdir -p single-modes; cd single-modes
    fd n10.tsv.gz$ ../single | rush 'ln -s {}'
    cd ..
    
    
    # Profiling  --------------------
    
    chunks=10
    db=refseq-cami2-k21-n$chunks.db
    dbname=refseq-cami2-k21-n$chunks    
    
    reads=single-modes    
    X=taxdump
    T=taxid.map
    
    for mode in 0 1 2 3 4 5; do
        fd kmcp@$dbname.tsv.gz$ $reads/ \
            | csvtk sort -H -k 1:N \
            | rush -v db=$db -v dbname=$dbname -v X=$X -v T=$T -v mode=$mode \
                'memusg -H -t -s "kmcp profile -m {mode} -X {X} -T {T} {} -o {}.k-m{mode}.profile -C {}.c-m{mode}.profile -s {@sample_(\d+)} \
                    --log {}.k-m{mode}.profile.log0" > {}.k-m{mode}.profile.log 2>&1' 
        
        profile=$reads@$dbname.c-m$mode.profile
        fd kmcp@$dbname.tsv.gz.c-m$mode.profile$ $reads/ \
            | csvtk sort -H -k 1:N \
            | rush -j 1 'cat {}' \
            > $profile
            
        # do not correct ambiguous reads
        
        fd kmcp@$dbname.tsv.gz$ $reads/ \
            | csvtk sort -H -k 1:N \
            | rush -v db=$db -v dbname=$dbname -v X=$X -v T=$T -v mode=$mode \
                'memusg -H -t -s "kmcp profile -m {mode} -X {X} -T {T} {} -o {}.k-m{mode}-nc.profile -C {}.c-m{mode}-nc.profile -s {@sample_(\d+)} \
                    --no-amb-corr \
                    --log {}.k-m{mode}-nc.profile.log0" > {}.k-m{mode}-nc.profile.log 2>&1'
        
        profile=$reads@$dbname.c-m$mode-nc.profile
        fd kmcp@$dbname.tsv.gz.c-m$mode-nc.profile$ $reads/ \
            | csvtk sort -H -k 1:N \
            | rush -j 1 'cat {}' \
            > $profile
    done
    
Metrics of profiling accuracy

    # OPAL
    opal.py -g cami2_mouse_gut_gs.profile \
        single-modes@refseq-cami2-k21-n10.c-m0.profile \
        single-modes@refseq-cami2-k21-n10.c-m1.profile \
        single-modes@refseq-cami2-k21-n10.c-m2.profile \
        single-modes@refseq-cami2-k21-n10.c-m3.profile \
        single-modes@refseq-cami2-k21-n10.c-m4.profile \
        single-modes@refseq-cami2-k21-n10.c-m5.profile \
        single-modes@refseq-cami2-k21-n10.c-m0-nc.profile \
        single-modes@refseq-cami2-k21-n10.c-m1-nc.profile \
        single-modes@refseq-cami2-k21-n10.c-m2-nc.profile \
        single-modes@refseq-cami2-k21-n10.c-m3-nc.profile \
        single-modes@refseq-cami2-k21-n10.c-m4-nc.profile \
        single-modes@refseq-cami2-k21-n10.c-m5-nc.profile \
        -l m0,m1,m2,m3,m4,m5,m0-nc,m1-nc,m2-nc,m3-nc,m4-nc,m5-nc \
        -o opal

    cat opal/results.tsv \
        | csvtk grep -t -f tool -p 'Gold standard' -v \
        | csvtk grep -t -f rank -p na -p species -p genus \
        | csvtk grep -t -f metric -p 'Completeness' -p 'Purity' -p 'F1 score' \
            -p 'L1 norm error' -p 'Weighted UniFrac error' \
        | tee accuracy.tsv \
        | csvtk summary -t -g tool,metric,rank -f value:mean -w 6 \
        | csvtk rename -t -f value:mean -n value \
        | csvtk sort -t -k rank -k metric -k tool:N \
        | csvtk csv2md -t
