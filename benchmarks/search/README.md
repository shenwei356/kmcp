## Introduction

Softwares

- [cobs](https://github.com/bingmann/cobs) ([1915fc0](https://github.com/bingmann/cobs/commit/1915fc061bbe47946116b4a051ed7b4e3f3eca15))
- [sourmash](https://github.com/dib-lab/sourmash) (v4.2.1)
- kmcp (v0.6.0)
- [memusg](https://github.com/shenwei356/memusg), version: [91a19ab](https://github.com/shenwei356/memusg/commit/91a19abaf041c3046b91ef3a35ed28aade1e05fc)

## Reference sequences

GTDB r202

- [gtdb_genomes_reps_r202.tar.gz](https://data.ace.uq.edu.au/public/gtdb/data/releases/release202/202.0/genomic_files_reps/gtdb_genomes_reps_r202.tar.gz)
- file size: 46.26 GB
- files: 47,894
- bases: 151.94 Gb


Uncompressing and renaming

    # uncompress
    mkdir -p gtdb202
    tar -zxvf gtdb_genomes_reps_r202.tar.gz -O gtdb202

    # rename
    brename -R -p '^(\w{3}_\d{9}\.\d+).+' -r '$1.fa.gz' gtdb202    

## Whole sequence

Building database

    # indexing ---------------------------------------------------------------------------------
    
    seqs=gtdb202 # directory of gtdb sequences
    db=gtdb
    k=31
    threads=40
    
    dbCOBS=gtdb-cobs-k$k.cobs_compact
    dbKMCPtmp=gtdb-kmcp-k$k
    dbKMCP=gtdb-kmcp-k$k.db
    
    
    # --------------- cobs ---------------
     
    # time: 51m:01s, Peak mem: 76.14G
    /bin/rm -rf $dbCOBS $seqs/*.cobs_cache 
    memusg -t -s "cobs compact-construct -T $threads -k $k -f 0.3 --num-hashes 1 -p 1024 --file-type fasta $seqs $dbCOBS --clobber" \
        2>$dbCOBS.time
    
    # 86.96 G
    du -sh $dbCOBS > $dbCOBS.size
    
    
    # --------------- kmcp ---------------
    
    # 24m:52s, 49.31G
    /bin/rm -rf $dbKMCPtmp $dbKMCP
    memusg -t -s "kmcp compute -e -j $threads -k $k -I $seqs -O $dbKMCPtmp --force --quiet \
        && kmcp index -j $threads -f 0.3 -n 1 -b 1024 -I $dbKMCPtmp -O $dbKMCP --force --quiet " \
        2>$dbKMCP.time

    # 54.63 G
    du -sh $dbKMCP > $dbKMCP.size
    
-              |cobs      |kmcp
:--------------|:---------|:-------
database size  | 86.96 G  | 55.15 G
building time  | 51min01s | 18min25s
    
Searching (bacterial genome)
    
    # searching  ---------------------------------------------------------------------------------

    # cobs query needs plain text
    for f in refs/*.fasta.gz; do
        gzip -d $f
    done
    
    
    # whole genome ----------------
    t=0.8
    
    # emptying the buffers cache
    # su -c "free && sync && echo 3 > /proc/sys/vm/drop_caches && free"
    
    for f in refs/*.fasta; do
        echo $f
        # cobs
        memusg -t -H -s "cobs query  -T $threads -i $dbCOBS -f $f -t $t > $f.cobs@$db.txt" 2>$f.cobs@$db.txt.time
        
        # kmcp
        memusg -t -H -s "kmcp search -g -j $threads -d $dbKMCP    $f -t $t --quiet > $f.kmcp@$db.txt" 2>$f.kmcp@$db.txt.time
    done

    # time
    find refs/ -name "*.fasta" \
        | rush -k 'echo -ne {%}; \
            echo -en "\t"$(grep "elapsed time" {}.cobs@gtdb.txt.time | sed -r "s/.+: //"); \
            echo -en "\t"$(grep "peak rss" {}.cobs@gtdb.txt.time | sed -r "s/.+: //"); \
            echo -en "\t"$(grep "elapsed time" {}.kmcp@gtdb.txt.time | sed -r "s/.+: //"); \
            echo -e "\t"$(grep "peak rss" {}.kmcp@gtdb.txt.time | sed -r "s/.+: //"); \
            ' \
        | sed 's/ /\t/g' \
        | csvtk add-header -t -n "query,cobs_t,cobs_m0,kmcp_t,kmcp_m0" \
        | csvtk mutate2 -t -n cobs_m -e '$cobs_m0 / 1048576' \
        | csvtk mutate2 -t -n kmcp_m -e '$kmcp_m0 / 1048576' \
        | csvtk cut -t -f query,cobs_t,cobs_m,kmcp_t,kmcp_m \
        | csvtk rename -t -f 1-5 -n "query,cobs:time(s),cobs:memory(GB),kmcp:time(s),kmcp:memory(GB)" \
        | csvtk csv2md -t
        
    
Query time and peak memory

query              |cobs:time(s)|cobs:memory(GB)|kmcp:time(s)|kmcp:memory(GB)
:------------------|:-----------|:--------------|:-----------|:--------------
NC_000913.3.fasta  |585.024     |58.24          |165.475     |54.91
NC_002695.2.fasta  |262.735     |59.00          |17.583      |55.07
NC_010655.1.fasta  |197.728     |40.02          |12.379      |54.52
NC_011750.1.fasta  |161.446     |57.29          |13.800      |55.07
NC_012971.2.fasta  |60.677      |60.10          |13.695      |55.07
NC_013654.1.fasta  |108.834     |60.74          |13.542      |55.10
NC_018658.1.fasta  |74.624      |66.96          |13.304      |55.13
NZ_CP007592.1.fasta|72.536      |63.80          |13.985      |55.04
NZ_CP028116.1.fasta|62.162      |72.72          |13.505      |55.05



Searching (short reads)

    # short sequence ----------------
    t=0.8    
    for f in refs/*.fasta; do
        cat $f | seqkit sliding -s 20 -W 150  > $f.short
    done
    
    for f in refs/*.fasta.short; do
        echo $f
        # cobs
        memusg -t -H -s "cobs query  -T $threads -i $dbCOBS -f $f -t $t > $f.cobs@$db.txt" 2>$f.cobs@$db.txt.time
        
        # kmcp
        memusg -t -H -s "kmcp search -g -j $threads -d $dbKMCP    $f -t $t --quiet > $f.kmcp@$db.txt" 2>$f.kmcp@$db.txt.time
    done

    find refs/ -name "*.fasta.short" \
        | rush -k 'echo -ne {%}; \
            echo -en "\t"$(grep "elapsed time" {}.cobs@gtdb.txt.time | sed -r "s/.+: //"); \
            echo -en "\t"$(grep "peak rss" {}.cobs@gtdb.txt.time | sed -r "s/.+: //"); \
            echo -en "\t"$(grep "elapsed time" {}.kmcp@gtdb.txt.time | sed -r "s/.+: //"); \
            echo -e "\t"$(grep "peak rss" {}.kmcp@gtdb.txt.time | sed -r "s/.+: //"); \
            ' \
        | sed 's/ /\t/g' \
        | csvtk add-header -t -n "query,cobs_t,cobs_m0,kmcp_t,kmcp_m0" \
        | csvtk mutate2 -t -n cobs_m -e '$cobs_m0 / 1048576' \
        | csvtk mutate2 -t -n kmcp_m -e '$kmcp_m0 / 1048576' \
        | csvtk cut -t -f query,cobs_t,cobs_m,kmcp_t,kmcp_m \
        | csvtk rename -t -f 1-5 -n "query,cobs:time(s),cobs:memory(GB),kmcp:time(s),kmcp:memory(GB)" \
        | csvtk csv2md -t
    
query                    |cobs:time(s)|cobs:memory(GB)|kmcp:time(s)|kmcp:memory(GB)
:------------------------|:-----------|:--------------|:-----------|:--------------
NC_000913.3.fasta.short  |147.911     |61.02          |17.991      |55.93
NC_002695.2.fasta.short  |173.254     |66.30          |19.631      |56.13
NC_010655.1.fasta.short  |81.341      |43.49          |15.903      |55.33
NC_011750.1.fasta.short  |159.786     |63.92          |19.090      |55.70
NC_012971.2.fasta.short  |141.801     |60.35          |19.053      |55.70
NC_013654.1.fasta.short  |150.562     |61.73          |18.837      |55.92
NC_018658.1.fasta.short  |165.385     |65.26          |19.471      |55.74
NZ_CP007592.1.fasta.short|159.737     |64.34          |18.471      |55.73
NZ_CP028116.1.fasta.short|176.580     |66.83          |21.340      |56.22


## Scaled minhash


Building database

    # indexing ---------------------------------------------------------------------------------

    seqs=gtdb202
    db=gtdb
    k=31
    threads=16
    scale=1000
    
    dbSOURMASHtmp=gtdb-sourmash-k$k-D$scale
    dbSOURMASH=gtdb-sourmash-k$k-D$scale/_db.sbt.json
    dbKMCPtmp=gtdb-kmcp-k$k-D$scale
    dbKMCP=gtdb-kmcp-k$k-D$scale.db
    
    
    # --------------- sourmash ---------------
    
    # 28m12s
    # 30G
    mkdir -p $dbSOURMASHtmp
    indexSourmash() {
        fd .fa.gz$ $seqs \
            | rush -j $threads -v d=$dbSOURMASHtmp -v s=$scale -v k=$k \
                'sourmash -q sketch dna -p k={k},scaled={s} {} -o {d}/{%}.sig'     
        sourmash -q index $dbSOURMASH --from-file <(ls $dbSOURMASHtmp/*.sig)
    }
    
    { time indexSourmash ; } 2> $dbSOURMASH.time
    
    # --------------- kmcp ---------------
    
    # 
    # 
    memusg -t -s "kmcp compute -e -j $threads -k $k -I $seqs -O $dbKMCPtmp -D $scale --force --quiet \
        && kmcp index -j $threads -f 0.001 -n 3 -I $dbKMCPtmp -O $dbKMCP --force --quiet " \
        2>$dbKMCP.time
    
-              |sourmash      |kmcp
:--------------|:---------|:-------
database size  |       G  |  G
buiding time   |          | 

    # searching  ---------------------------------------------------------------------------------
    
    t=0.8
    for f in refs/*.fasta; do
        memusg -t -s "sourmash -q sketch dna -p k=$k,scaled=$scale $f -o $f.sig; \
            sourmash -q search --containment $f.sig $dbSOURMASH  --threshold $t > $f.sourmash@$db.txt" 2>$f.sourmash@$db.txt.time
        
        # 0.558s, 1.32G
        # single-thread 1.2s
        memusg -t -s "kmcp search -g -j $threads -d $dbKMCP    $f -t $t --quiet > $f.kmcp.scaled@$db.txt" 2>$f.kmcp.scaled@$db.txt.time
    done
