## Softwares

Softwares

- [cobs](https://github.com/bingmann/cobs) ([1915fc0](https://github.com/bingmann/cobs/commit/1915fc061bbe47946116b4a051ed7b4e3f3eca15))
- [sourmash](https://github.com/dib-lab/sourmash) (v4.2.1)
- kmcp ([v0.6.0](https://github.com/shenwei356/kmcp/releases/tag/v0.6.0))

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
    memusg -t -s "kmcp compute -j $threads -k $k -I $seqs -O $dbKMCPtmp --force --quiet \
        && kmcp index -j $threads -f 0.3 -n 1 -b 1024 -I $dbKMCPtmp -O $dbKMCP --force --quiet " \
        2>$dbKMCP.time

    du -sh $dbKMCP > $dbKMCP.size
    
Database size and building time.

|               |cobs      |kmcp     |
|:--------------|:---------|:--------|
|database size  | 86.96GB  | 55.15GB |
|building time  | 29m:55s  | 24min52s|
|temporary files| 160.76GB | 1.19TB  |
    
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
NC_000913.3.fasta  |68.113      |113.56         |161.158     |54.74
NC_002695.2.fasta  |51.283      |118.47         |15.429      |54.86
NC_010655.1.fasta  |44.998      |102.24         |11.647      |54.51
NC_011750.1.fasta  |52.408      |116.37         |11.928      |54.87
NC_012971.2.fasta  |51.739      |113.09         |12.038      |54.84
NC_013654.1.fasta  |50.885      |114.00         |12.297      |54.87
NC_018658.1.fasta  |55.553      |117.19         |12.288      |54.84
NZ_CP007592.1.fasta|52.696      |116.22         |12.641      |54.85
NZ_CP028116.1.fasta|52.767      |119.33         |12.210      |54.87


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
NC_000913.3.fasta.short  |184.737     |97.60          |15.813      |55.70
NC_002695.2.fasta.short  |212.143     |99.57          |17.209      |56.08
NC_010655.1.fasta.short  |132.232     |93.07          |13.830      |55.33
NC_011750.1.fasta.short  |199.384     |98.73          |16.925      |55.69
NC_012971.2.fasta.short  |176.219     |97.42          |16.189      |55.65
NC_013654.1.fasta.short  |180.648     |97.78          |16.169      |55.67
NC_018658.1.fasta.short  |198.074     |99.06          |17.354      |55.64
NZ_CP007592.1.fasta.short|194.866     |98.67          |16.223      |55.71
NZ_CP028116.1.fasta.short|207.710     |99.91          |17.642      |56.12


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
    
    memusg -t -s "kmcp compute -j $threads -k $k -I $seqs -O $dbKMCPtmp -D $scale --force --quiet \
        && kmcp index -j $threads -f 0.001 -n 3 -I $dbKMCPtmp -O $dbKMCP --force --quiet " \
        2>$dbKMCP.time


|               |sourmash  |kmcp    |
:---------------|:---------|:-------|
|database size  |  5.19GB  | 1.52GB |
|buiding time   |  89m59s  | 7min02s|
|temporary files| -        | 3.41GB |


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
        memusg -t -H -s "kmcp search -g -j 1 -d $dbKMCP    $f -t $t --quiet > $f.kmcp.scaled@$db.txt" 2>$f.kmcp.scaled@$db.txt.time
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
NC_000913.3.fasta  |12.491          |187.77             |20.866      |693.22
NC_002695.2.fasta  |6.499           |185.69             |5.054       |753.05
NC_010655.1.fasta  |3.611           |205.79             |1.275       |477.62
NC_011750.1.fasta  |4.416           |185.69             |1.461       |565.40
NC_012971.2.fasta  |3.812           |183.69             |1.236       |516.36
NC_013654.1.fasta  |3.956           |185.68             |1.271       |590.29
NC_018658.1.fasta  |4.157           |183.70             |1.235       |721.58
NZ_CP007592.1.fasta|3.959           |183.70             |1.278       |664.78
NZ_CP028116.1.fasta|4.333           |187.68             |1.250       |416.22



Query time and peak memory (kmcp utilizes 8 threads, cold start)

query              |sourmash:time(s)|sourmash:memory(MB)|kmcp:time(s)|kmcp:memory(MB)
:------------------|:---------------|:------------------|:-----------|:--------------
NC_000913.3.fasta  |12.515          |185.86             |4.891       |692.79
NC_002695.2.fasta  |6.559           |181.71             |2.495       |745.04
NC_010655.1.fasta  |3.806           |208.96             |0.792       |467.58
NC_011750.1.fasta  |4.029           |195.29             |0.639       |575.79
NC_012971.2.fasta  |3.939           |195.64             |0.641       |299.26
NC_013654.1.fasta  |3.969           |183.69             |0.630       |186.05
NC_018658.1.fasta  |4.196           |185.69             |0.636       |152.86
NZ_CP007592.1.fasta|4.027           |187.68             |0.658       |59.88
NZ_CP028116.1.fasta|4.194           |187.68             |0.644       |58.36
