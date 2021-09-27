# Benchmarks on sequence and genome searching

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

    memusg -t -s "cobs compact-construct -T $threads -k $k -f 0.3 --num-hashes 1 -p 1024 \
        --file-type fasta $seqs $dbCOBS --keep-temporary --clobber" \
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
        memusg -t -H -s "cobs query --load-complete -T $threads -i $dbCOBS -f $f -t $t > $f.cobs@$db.txt" \
            2>$f.cobs@$db.txt.time
        
        # kmcp
        memusg -t -H -s "kmcp search -g -j $threads -d $dbKMCP    $f -t $t --quiet > $f.kmcp@$db.txt" \
            2>$f.kmcp@$db.txt.time
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

|query              |cobs:time(s)|cobs:memory(GB)|kmcp:time(s)|kmcp:memory(GB)|
|:------------------|:-----------|:--------------|:-----------|:--------------|
|NC_000913.3.fasta  |67.733      |113.56         |170.532     |54.74          |
|NC_002695.2.fasta  |53.863      |118.48         |16.069      |54.91          |
|NC_010655.1.fasta  |49.738      |102.23         |12.001      |54.40          |
|NC_011750.1.fasta  |52.687      |116.38         |12.355      |54.87          |
|NC_012971.2.fasta  |51.707      |113.09         |12.268      |54.83          |
|NC_013654.1.fasta  |54.652      |114.00         |12.536      |54.84          |
|NC_018658.1.fasta  |52.455      |117.18         |12.632      |54.87          |
|NZ_CP007592.1.fasta|50.348      |116.22         |12.691      |54.88          |
|NZ_CP028116.1.fasta|54.619      |119.33         |12.999      |54.90          |


### Searching with short reads

    
    for f in refs/*.fasta; do
        cat $f | seqkit sliding -s 20 -W 150  > $f.short
    done

    seqkit stats -j 10 -T refs/*.short | csvtk csv2md -t
    
|file                          |format|type|num_seqs|sum_len |min_len|avg_len|max_len|
|:-----------------------------|:-----|:---|:-------|:-------|:------|:------|:------|
|refs/NC_000913.3.fasta.short  |FASTA |DNA |232076  |34811400|150    |150.0  |150    |
|refs/NC_002695.2.fasta.short  |FASTA |DNA |274922  |41238300|150    |150.0  |150    |
|refs/NC_010655.1.fasta.short  |FASTA |DNA |133198  |19979700|150    |150.0  |150    |
|refs/NC_011750.1.fasta.short  |FASTA |DNA |256596  |38489400|150    |150.0  |150    |
|refs/NC_012971.2.fasta.short  |FASTA |DNA |227941  |34191150|150    |150.0  |150    |
|refs/NC_013654.1.fasta.short  |FASTA |DNA |235860  |35379000|150    |150.0  |150    |
|refs/NC_018658.1.fasta.short  |FASTA |DNA |263648  |39547200|150    |150.0  |150    |
|refs/NZ_CP007592.1.fasta.short|FASTA |DNA |255221  |38283150|150    |150.0  |150    |
|refs/NZ_CP028116.1.fasta.short|FASTA |DNA |282402  |42360300|150    |150.0  |150    |

    t=0.8
    
    for f in refs/*.fasta.short; do
        echo $f
        # cobs
        memusg -t -H -s "cobs query --load-complete -T $threads -i $dbCOBS -f $f -t $t > $f.cobs@$db.txt" \
            2>$f.cobs@$db.txt.time
        
        # kmcp
        memusg -t -H -s "kmcp search -g -j $threads -d $dbKMCP    $f -t $t --quiet > $f.kmcp@$db.txt" \
            2>$f.kmcp@$db.txt.time
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
    
Query time and peak memory (hot start, loading all data in memory)

|query                    |cobs:time(s)|cobs:memory(GB)|kmcp:time(s)|kmcp:memory(GB)|
|:------------------------|:-----------|:--------------|:-----------|:--------------|
|NC_000913.3.fasta.short  |181.941     |97.61          |16.781      |55.79          |
|NC_002695.2.fasta.short  |210.005     |99.57          |16.960      |56.18          |
|NC_010655.1.fasta.short  |121.183     |93.08          |14.384      |55.19          |
|NC_011750.1.fasta.short  |202.856     |98.73          |16.355      |55.79          |
|NC_012971.2.fasta.short  |184.675     |97.42          |16.041      |55.65          |
|NC_013654.1.fasta.short  |183.921     |97.77          |16.694      |55.79          |
|NC_018658.1.fasta.short  |203.736     |99.05          |17.266      |55.79          |
|NZ_CP007592.1.fasta.short|193.875     |98.67          |17.195      |55.79          |
|NZ_CP028116.1.fasta.short|213.337     |99.92          |17.253      |56.18          |


## KMCP vs Mash and Sourmash


### Building database (MinHash sketches)

    seqs=gtdb202
    db=gtdb
    k=31
    threads=8
    scale=1000
    scaleMash=3400 // average genome size: 3.4Mb.
    
    dbMASH=gtdb-mash-k$k-S$scaleMash.msh
    dbSOURMASHtmp=gtdb-sourmash-k$k-D$scale
    dbSOURMASH=gtdb-sourmash-k$k-D$scale/_db.sbt.json
    dbKMCPtmp=gtdb-kmcp-k$k-D$scale
    dbKMCP=gtdb-kmcp-k$k-D$scale.db
    
    
    # --------------- mash ---------------
    
    find $seqs -name *.fa.gz > $seqs.list
    memusg -t -s "mash sketch -p $threads -s $scaleMash -k $k -o $dbMASH -l $seqs.list " \
        2>$dbMASH.time
    
    
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


|               |mash   |sourmash  |kmcp    |
|:--------------|:------|:---------|:-------|
|database size  |743MB  | 5.19GB   | 1.52GB |
|buiding time   |11m39s | 89m59s   | 7min02s|
|temporary files|-      |-         | 3.41GB |


### Searching with bacterial genomes
    
    # emptying the buffers cache
    # su -c "free && sync && echo 3 > /proc/sys/vm/drop_caches && free"
    
    
    threads=8
    
    t=0.4
    for f in refs/*.fasta; do
        echo $f
        # mash
        memusg -t -H -s "mash dist -p $threads -s $scaleMash -v 0.01 -d 0.05 $dbMASH $f > $f.mash@$db.txt " \
            2>$f.mash@$db.txt.time
        
        # sourmash
        memusg -t -H -s "sourmash -q sketch dna -p k=$k,scaled=$scale $f -o $f.sig; \
            sourmash -q search $f.sig $dbSOURMASH  --threshold $t > $f.sourmash@$db.txt" \
            2>$f.sourmash@$db.txt.time
        
        # kmcp, single-thread
        memusg -t -H -s "kmcp search -g -j $threads -d $dbKMCP    $f -t $t --quiet > $f.kmcp.scaled@$db.txt" \
            2>$f.kmcp.scaled@$db.txt.time
    done

    # time
    find refs/ -name "*.fasta" \
        | rush -k 'echo -ne {%}; \
            echo -en "\t"$(grep "elapsed time" {}.mash@gtdb.txt.time | sed -r "s/.+: //"); \
            echo -en "\t"$(grep "peak rss" {}.mash@gtdb.txt.time | sed -r "s/.+: //"); \
            echo -en "\t"$(grep "elapsed time" {}.sourmash@gtdb.txt.time | sed -r "s/.+: //"); \
            echo -en "\t"$(grep "peak rss" {}.sourmash@gtdb.txt.time | sed -r "s/.+: //"); \
            echo -en "\t"$(grep "elapsed time" {}.kmcp.scaled@gtdb.txt.time | sed -r "s/.+: //"); \
            echo -e "\t"$(grep "peak rss" {}.kmcp.scaled@gtdb.txt.time | sed -r "s/.+: //"); \
            ' \
        | sed 's/ /\t/g' \
        | csvtk add-header -t -n "query,mash_t,mash_m0,sourmash_t,sourmash_m0,kmcp_t,kmcp_m0" \
        | csvtk mutate2 -t -n mash_m -e '$mash_m0 / 1024' \
        | csvtk mutate2 -t -n sourmash_m -e '$sourmash_m0 / 1024' \
        | csvtk mutate2 -t -n kmcp_m -e '$kmcp_m0 / 1024' \
        | csvtk cut -t -f query,mash_t,mash_m,sourmash_t,sourmash_m,kmcp_t,kmcp_m \
        | csvtk rename -t -f 1-7 -n "query,mash:time(s),mash:mem(MB),sourmash:time(s),sourmash:mem(MB),kmcp:time(s),kmcp:mem(MB)" \
        | csvtk csv2md -t
        
    
Query time and peak memory (mash and kmcp utilize single thread, cold start)

|query              |mash:time(s)|mash:mem(MB)|sourmash:time(s)|sourmash:mem(MB)|kmcp:time(s)|kmcp:mem(MB)|
|:------------------|:-----------|:-----------|:---------------|:---------------|:-----------|:-----------|
|NC_000913.3.fasta  |3.382       |1492.70     |12.960          |241.41          |20.112      |680.85      |
|NC_002695.2.fasta  |3.615       |1495.46     |7.182           |228.50          |4.862       |746.45      |
|NC_010655.1.fasta  |3.150       |1488.69     |3.912           |232.46          |1.686       |454.02      |
|NC_011750.1.fasta  |3.369       |1494.16     |4.664           |228.50          |1.278       |655.17      |
|NC_012971.2.fasta  |3.597       |1492.43     |4.216           |230.56          |0.861       |615.64      |
|NC_013654.1.fasta  |3.914       |1493.46     |4.254           |234.48          |1.454       |642.32      |
|NC_018658.1.fasta  |3.413       |1493.92     |4.307           |228.51          |1.051       |639.74      |
|NZ_CP007592.1.fasta|3.485       |1494.42     |5.099           |234.48          |1.272       |584.30      |
|NZ_CP028116.1.fasta|3.783       |1496.68     |5.099           |228.51          |1.297       |614.72      |



Query time and peak memory (mash and kmcp utilize 8 threads, cold start)

|query              |mash:time(s)|mash:mem(MB)|sourmash:time(s)|sourmash:mem(MB)|kmcp:time(s)|kmcp:mem(MB)|
|:------------------|:-----------|:-----------|:---------------|:---------------|:-----------|:-----------|
|NC_000913.3.fasta  |2.292       |1495.78     |12.994          |228.59          |4.868       |682.43      |
|NC_002695.2.fasta  |2.517       |1494.37     |7.227           |238.47          |2.219       |713.38      |
|NC_010655.1.fasta  |2.292       |1488.93     |3.906           |224.49          |0.852       |460.16      |
|NC_011750.1.fasta  |2.343       |1496.29     |4.671           |230.50          |0.563       |705.36      |
|NC_012971.2.fasta  |2.347       |1495.78     |4.220           |230.49          |0.646       |227.00      |
|NC_013654.1.fasta  |2.544       |1495.78     |4.234           |230.49          |0.644       |360.12      |
|NC_018658.1.fasta  |2.340       |1496.29     |4.468           |228.51          |0.620       |41.04       |
|NZ_CP007592.1.fasta|2.341       |1496.30     |4.843           |234.49          |0.632       |81.90       |
|NZ_CP028116.1.fasta|2.547       |1496.84     |4.866           |232.50          |0.647       |76.75       |
