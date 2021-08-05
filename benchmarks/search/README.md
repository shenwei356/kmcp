## Softwares

Softwares

- [cobs](https://github.com/bingmann/cobs) ([1915fc0](https://github.com/bingmann/cobs/commit/1915fc061bbe47946116b4a051ed7b4e3f3eca15))
- [sourmash](https://github.com/dib-lab/sourmash) (v4.2.1)
- kmcp (v0.6.0)

Utilities

- [brename](https://github.com/shenwei356/brename/releases) for batching renaming files.
- [rush](https://github.com/shenwei356/rush/releases) for executing jobs in parallel.
- [memusg](https://github.com/shenwei356/memusg) ([91a19ab](https://github.com/shenwei356/memusg/commit/91a19abaf041c3046b91ef3a35ed28aade1e05fc)) for monitoring peak memory usage and executing time.

## Datasets

GTDB representative genomes

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

## KMCP vs COBS

### Building database
    
    seqs=gtdb202 # directory of gtdb sequences
    db=gtdb
    k=31
    threads=40
    
    dbCOBS=gtdb-cobs-k$k.cobs_compact
    dbKMCPtmp=gtdb-kmcp-k$k
    dbKMCP=gtdb-kmcp-k$k.db
        
    # --------------- cobs ---------------
     
    /bin/rm -rf $dbCOBS
    find $seqs -name "*.cobs_cache" | rush "/bin/rm {}"

    memusg -t -s "cobs compact-construct -T $threads -k $k -f 0.3 --num-hashes 1 -p 1024 --file-type fasta $seqs $dbCOBS --keep-temporary --clobber" \
        2>$dbCOBS.time
    
    du -sh $dbCOBS > $dbCOBS.size
        
    # --------------- kmcp ---------------
    
    /bin/rm -rf $dbKMCPtmp $dbKMCP
    memusg -t -s "kmcp compute -e -j $threads -k $k -I $seqs -O $dbKMCPtmp --force --quiet \
        && kmcp index -j $threads -f 0.3 -n 1 -b 1024 -I $dbKMCPtmp -O $dbKMCP --force --quiet " \
        2>$dbKMCP.time

    du -sh $dbKMCP > $dbKMCP.size
    
Database size and building time.

-              |cobs      |kmcp
:--------------|:---------|:-------
database size  | 86.96GB  | 55.15GB
building time  | 29m:55s  | 24min52s
temporary files| 160.76GB | 1.19TB
    
### Searching with bacterial genomes

    # cobs query needs plain text
    for f in refs/*.fasta.gz; do
        gzip -d $f
    done
    
    # whole genome ----------------

    t=0.5
    
    # emptying the buffers cache
    # su -c "free && sync && echo 3 > /proc/sys/vm/drop_caches && free"
    
    for f in refs/*.fasta; do
        echo $f
        # cobs
        memusg -t -H -s "cobs query --load-complete -T $threads -i $dbCOBS -f $f -t $t > $f.cobs@$db.txt" 2>$f.cobs@$db.txt.time
        
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
        
Query time and peak memory (cold start, loading all data in memory)

query              |cobs:time(s)|cobs:memory(GB)|kmcp:time(s)|kmcp:memory(GB)
:------------------|:-----------|:--------------|:-----------|:--------------
NC_000913.3.fasta  |66.781      |113.56         |158.206     |54.93
NC_002695.2.fasta  |53.470      |118.48         |16.568      |55.10
NC_010655.1.fasta  |49.788      |102.24         |12.286      |54.63
NC_011750.1.fasta  |54.853      |116.37         |12.893      |55.11
NC_012971.2.fasta  |53.294      |113.09         |12.680      |55.02
NC_013654.1.fasta  |50.987      |114.00         |13.160      |55.16
NC_018658.1.fasta  |55.782      |117.18         |12.983      |55.15
NZ_CP007592.1.fasta|51.755      |116.22         |13.262      |55.14
NZ_CP028116.1.fasta|51.402      |119.34         |13.116      |55.03


### Searching with short reads

    
    for f in refs/*.fasta; do
        cat $f | seqkit sliding -s 20 -W 150  > $f.short
    done

    seqkit stats -j 10 -T refs/*.short | csvtk csv2md -t
    
file                          |format|type|num_seqs|sum_len |min_len|avg_len|max_len
:-----------------------------|:-----|:---|:-------|:-------|:------|:------|:------
refs/NC_000913.3.fasta.short  |FASTA |DNA |232076  |34811400|150    |150.0  |150
refs/NC_002695.2.fasta.short  |FASTA |DNA |274922  |41238300|150    |150.0  |150
refs/NC_010655.1.fasta.short  |FASTA |DNA |133198  |19979700|150    |150.0  |150
refs/NC_011750.1.fasta.short  |FASTA |DNA |256596  |38489400|150    |150.0  |150
refs/NC_012971.2.fasta.short  |FASTA |DNA |227941  |34191150|150    |150.0  |150
refs/NC_013654.1.fasta.short  |FASTA |DNA |235860  |35379000|150    |150.0  |150
refs/NC_018658.1.fasta.short  |FASTA |DNA |263648  |39547200|150    |150.0  |150
refs/NZ_CP007592.1.fasta.short|FASTA |DNA |255221  |38283150|150    |150.0  |150
refs/NZ_CP028116.1.fasta.short|FASTA |DNA |282402  |42360300|150    |150.0  |150

    t=0.8
    
    for f in refs/*.fasta.short; do
        echo $f
        # cobs
        memusg -t -H -s "cobs query --load-complete -T $threads -i $dbCOBS -f $f -t $t > $f.cobs@$db.txt" 2>$f.cobs@$db.txt.time
        
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
    
Query time and peak memory (hot start)

query                    |cobs:time(s)|cobs:memory(GB)|kmcp:time(s)|kmcp:memory(GB)
:------------------------|:-----------|:--------------|:-----------|:--------------
NC_000913.3.fasta.short  |256.995     |97.61          |15.551      |55.96
NC_002695.2.fasta.short  |220.510     |99.57          |17.097      |56.80
NC_010655.1.fasta.short  |136.154     |93.07          |13.100      |55.28
NC_011750.1.fasta.short  |209.332     |98.73          |16.675      |55.71
NC_012971.2.fasta.short  |189.145     |97.41          |15.916      |55.89
NC_013654.1.fasta.short  |196.493     |97.78          |15.675      |55.93
NC_018658.1.fasta.short  |215.541     |99.05          |16.475      |55.88
NZ_CP007592.1.fasta.short|205.543     |98.66          |16.624      |55.70
NZ_CP028116.1.fasta.short|227.358     |99.92          |17.197      |56.27


## KMCP vs sourmash


### Building database (MinHash sketches)

    seqs=gtdb202
    db=gtdb
    k=31
    threads=8
    scale=1000
    
    dbSOURMASHtmp=gtdb-sourmash-k$k-D$scale
    dbSOURMASH=gtdb-sourmash-k$k-D$scale/_db.sbt.json
    dbKMCPtmp=gtdb-kmcp-k$k-D$scale
    dbKMCP=gtdb-kmcp-k$k-D$scale.db
    
    
    # --------------- sourmash ---------------
    
    # 89m59.612s
    # 5.19G
    mkdir -p $dbSOURMASHtmp
    indexSourmash() {
        find $seqs -name "*.fa.gz" \
            | rush -j $threads -v d=$dbSOURMASHtmp -v s=$scale -v k=$k \
                'sourmash -q sketch dna -p k={k},scaled={s} {} -o {d}/{%}.sig'     
        sourmash -q index $dbSOURMASH --from-file <(find $dbSOURMASHtmp -name "*.sig")
    }
    
    { time indexSourmash ; } 2> $dbSOURMASH.time
    

    # --------------- kmcp ---------------
    
    memusg -t -s "kmcp compute -e -j $threads -k $k -I $seqs -O $dbKMCPtmp -D $scale --force --quiet \
        && kmcp index -j $threads -f 0.001 -n 3 -I $dbKMCPtmp -O $dbKMCP --force --quiet " \
        2>$dbKMCP.time


-              |sourmash  |kmcp
:--------------|:---------|:-------
database size  |  5.19GB  | 1.52GB
buiding time   |  89m59s  | 7min02s
temporary files| -        | 3.41GB


### Searching with bacterial genomes
    
    # emptying the buffers cache
    # su -c "free && sync && echo 3 > /proc/sys/vm/drop_caches && free"
    
    t=0.5
    for f in refs/*.fasta; do
        echo $f
        # sourmash
        memusg -t -H -s "sourmash -q sketch dna -p k=$k,scaled=$scale $f -o $f.sig; \
            sourmash -q search $f.sig $dbSOURMASH  --threshold $t > $f.sourmash@$db.txt" 2>$f.sourmash@$db.txt.time
        
        # kmcp, single-thread
        memusg -t -H -s "kmcp search -g -j 8 -d $dbKMCP    $f -t $t --quiet > $f.kmcp.scaled@$db.txt" 2>$f.kmcp.scaled@$db.txt.time
    done

    # time
    find refs/ -name "*.fasta" \
        | rush -k 'echo -ne {%}; \
            echo -en "\t"$(grep "elapsed time" {}.sourmash@gtdb.txt.time | sed -r "s/.+: //"); \
            echo -en "\t"$(grep "peak rss" {}.sourmash@gtdb.txt.time | sed -r "s/.+: //"); \
            echo -en "\t"$(grep "elapsed time" {}.kmcp.scaled@gtdb.txt.time | sed -r "s/.+: //"); \
            echo -e "\t"$(grep "peak rss" {}.kmcp.scaled@gtdb.txt.time | sed -r "s/.+: //"); \
            ' \
        | sed 's/ /\t/g' \
        | csvtk add-header -t -n "query,sourmash_t,sourmash_m0,kmcp_t,kmcp_m0" \
        | csvtk mutate2 -t -n sourmash_m -e '$sourmash_m0 / 1024' \
        | csvtk mutate2 -t -n kmcp_m -e '$kmcp_m0 / 1024' \
        | csvtk cut -t -f query,sourmash_t,sourmash_m,kmcp_t,kmcp_m \
        | csvtk rename -t -f 1-5 -n "query,sourmash:time(s),sourmash:memory(MB),kmcp:time(s),kmcp:memory(MB)" \
        | csvtk csv2md -t
        
    
Query time and peak memory (kmcp utilizes single thread, cold start)

query              |sourmash:time(s)|sourmash:memory(MB)|kmcp:time(s)|kmcp:memory(MB)
:------------------|:---------------|:------------------|:-----------|:--------------
NC_000913.3.fasta  |12.758          |230.31             |19.993      |697.38
NC_002695.2.fasta  |6.584           |242.73             |4.828       |738.39
NC_010655.1.fasta  |3.830           |225.41             |1.467       |363.21
NC_011750.1.fasta  |4.677           |232.41             |1.468       |676.07
NC_012971.2.fasta  |3.977           |221.39             |1.068       |546.64
NC_013654.1.fasta  |4.028           |221.39             |1.493       |673.60
NC_018658.1.fasta  |4.235           |221.39             |1.066       |729.62
NZ_CP007592.1.fasta|4.233           |221.40             |1.274       |609.60
NZ_CP028116.1.fasta|4.232           |221.40             |1.068       |737.34


Query time and peak memory (kmcp utilizes 8 threads, cold start)

query              |sourmash:time(s)|sourmash:memory(MB)|kmcp:time(s)|kmcp:memory(MB)
:------------------|:---------------|:------------------|:-----------|:--------------
NC_000913.3.fasta  |12.474          |242.72             |4.819       |689.82
NC_002695.2.fasta  |6.932           |229.43             |2.298       |746.47
NC_010655.1.fasta  |3.797           |221.36             |0.852       |476.36
NC_011750.1.fasta  |4.568           |222.50             |0.633       |300.12
NC_012971.2.fasta  |4.005           |225.62             |0.605       |70.46
NC_013654.1.fasta  |4.196           |225.18             |0.617       |155.55
NC_018658.1.fasta  |4.366           |232.67             |0.616       |159.48
NZ_CP007592.1.fasta|4.651           |244.54             |0.639       |278.61
NZ_CP028116.1.fasta|4.439           |221.40             |0.632       |128.37
