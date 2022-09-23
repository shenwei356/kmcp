# Benchmarks on sequence and genome searching

<img src="bench.searching.jpg" alt="" width="600"/>

## Softwares

Softwares

- [COBS](https://github.com/bingmann/cobs) ([1915fc0](https://github.com/bingmann/cobs/commit/1915fc061bbe47946116b4a051ed7b4e3f3eca15))
- [Sourmash](https://github.com/dib-lab/sourmash) (v4.5.0)
- [Mash](https://github.com/marbl/Mash) (v2.3)
- KMCP ([v0.9.0](https://github.com/shenwei356/kmcp/releases/tag/v0.9.0))

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
|building time  | 29m55s   | 21min04s|
|temporary files| 160.76GB | 935.11G |
    
### Searching with bacterial genomes

Generating the queries

    # cobs query needs plain text
    for f in refs/*.fasta.gz; do
        gzip -d $f
    done
    
    seqkit stats -j 10 -T refs/*.fasta | csvtk csv2md -t
    
|file                    |format|type|num_seqs|sum_len|min_len|avg_len  |max_len|
|:-----------------------|:-----|:---|:-------|:------|:------|:--------|:------|
|refs/NC_000913.3.fasta  |FASTA |DNA |1       |4641652|4641652|4641652.0|4641652|
|refs/NC_002695.2.fasta  |FASTA |DNA |1       |5498578|5498578|5498578.0|5498578|
|refs/NC_011750.1.fasta  |FASTA |DNA |1       |5132068|5132068|5132068.0|5132068|
|refs/NC_012971.2.fasta  |FASTA |DNA |1       |4558953|4558953|4558953.0|4558953|
|refs/NC_013654.1.fasta  |FASTA |DNA |1       |4717338|4717338|4717338.0|4717338|
|refs/NC_018658.1.fasta  |FASTA |DNA |1       |5273097|5273097|5273097.0|5273097|
|refs/NZ_CP007592.1.fasta|FASTA |DNA |1       |5104557|5104557|5104557.0|5104557|
|refs/NZ_CP028116.1.fasta|FASTA |DNA |1       |5648177|5648177|5648177.0|5648177|

Searching
    
    # whole genome ----------------

    t=0.5
    
    # emptying the buffers cache
    # su -c "free && sync && echo 3 > /proc/sys/vm/drop_caches && free"
    
    for f in refs/*.fasta; do
        echo $f
        # cobs
        memusg -t -H -s "cobs query -T $threads -i $dbCOBS -f $f -t $t > $f.cobs@$db.txt" \
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
        > bench.kmcp-cobs.long.tsv
        
    csvtk csv2md -t bench.kmcp-cobs.long.tsv
        
Query time and peak memory (cold start)

|query              |cobs:time(s)|cobs:memory(GB)|kmcp:time(s)|kmcp:memory(GB)|
|:------------------|:-----------|:--------------|:-----------|:--------------|
|NC_000913.3.fasta  |594.866     |57.70          |180.354     |54.72          |
|NC_002695.2.fasta  |255.000     |59.00          |16.792      |54.81          |
|NC_011750.1.fasta  |188.998     |57.14          |13.428      |54.81          |
|NC_012971.2.fasta  |60.119      |61.79          |12.871      |54.78          |
|NC_013654.1.fasta  |122.684     |57.75          |13.094      |54.79          |
|NC_018658.1.fasta  |82.371      |66.80          |13.171      |54.81          |
|NZ_CP007592.1.fasta|71.601      |68.53          |13.263      |54.81          |
|NZ_CP028116.1.fasta|62.426      |72.64          |13.507      |54.81          |

Query time and peak memory (hot start)

|query              |cobs:time(s)|cobs:memory(GB)|kmcp:time(s)|kmcp:memory(GB)|
|:------------------|:-----------|:--------------|:-----------|:--------------|
|NC_000913.3.fasta  |23.289      |77.02          |12.682      |54.76          |
|NC_002695.2.fasta  |27.296      |85.22          |13.526      |54.81          |
|NC_011750.1.fasta  |25.922      |81.61          |13.439      |54.81          |
|NC_012971.2.fasta  |24.484      |73.90          |13.056      |54.75          |
|NC_013654.1.fasta  |24.075      |76.57          |13.492      |54.78          |
|NC_018658.1.fasta  |26.822      |83.44          |13.709      |54.84          |
|NZ_CP007592.1.fasta|26.832      |81.84          |13.687      |54.86          |
|NZ_CP028116.1.fasta|26.892      |84.15          |12.951      |54.83          |

### Searching with short reads

Generating the queries
    
    for f in refs/*.fasta; do
        cat $f | seqkit sliding -s 4 -W 150  > $f.short
    done

    seqkit stats -j 10 -T refs/*.short | csvtk csv2md -t
    
|file                          |format|type|num_seqs|sum_len  |min_len|avg_len|max_len|
|:-----------------------------|:-----|:---|:-------|:--------|:------|:------|:------|
|refs/NC_000913.3.fasta.short  |FASTA |DNA |1160376 |174056400|150    |150.0  |150    |
|refs/NC_002695.2.fasta.short  |FASTA |DNA |1374608 |206191200|150    |150.0  |150    |
|refs/NC_011750.1.fasta.short  |FASTA |DNA |1282980 |192447000|150    |150.0  |150    |
|refs/NC_012971.2.fasta.short  |FASTA |DNA |1139701 |170955150|150    |150.0  |150    |
|refs/NC_013654.1.fasta.short  |FASTA |DNA |1179298 |176894700|150    |150.0  |150    |
|refs/NC_018658.1.fasta.short  |FASTA |DNA |1318237 |197735550|150    |150.0  |150    |
|refs/NZ_CP007592.1.fasta.short|FASTA |DNA |1276102 |191415300|150    |150.0  |150    |
|refs/NZ_CP028116.1.fasta.short|FASTA |DNA |1412007 |211801050|150    |150.0  |150    |

Query time and peak memory (hot start)

    t=0.8
    
    for f in refs/*.fasta.short; do
        echo $f
        # cobs
        # --load-complete 
        memusg -t -H -s "cobs query -T $threads -i $dbCOBS -f $f -t $t > $f.cobs@$db.txt" \
            2>$f.cobs@$db.txt.time
        
        # kmcp
        # --load-whole-db
        memusg -t -H -s "kmcp search -j $threads -d $dbKMCP    $f -t $t --quiet > $f.kmcp@$db.txt" \
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
        > bench.kmcp-cobs.short.tsv
        
    csvtk csv2md -t bench.kmcp-cobs.short.tsv

|query                    |cobs:time(s)|cobs:memory(GB)|kmcp:time(s)|kmcp:memory(GB)|
|:------------------------|:-----------|:--------------|:-----------|:--------------|
|NC_000913.3.fasta.short  |678.906     |103.54         |61.289      |54.87          |
|NC_002695.2.fasta.short  |815.897     |116.64         |62.255      |54.89          |
|NC_011750.1.fasta.short  |754.642     |110.94         |66.638      |54.89          |
|NC_012971.2.fasta.short  |658.405     |102.15         |53.403      |54.89          |
|NC_013654.1.fasta.short  |705.408     |104.92         |62.649      |54.86          |
|NC_018658.1.fasta.short  |763.529     |113.57         |60.111      |54.87          |
|NZ_CP007592.1.fasta.short|760.646     |111.12         |58.620      |54.90          |
|NZ_CP028116.1.fasta.short|857.118     |118.57         |72.758      |54.90          |



Query time and peak memory (hot start, loading all data in memory)

    t=0.8
    
    for f in refs/*.fasta.short; do
        echo $f
        # cobs
        # --load-complete 
        memusg -t -H -s "cobs query --load-complete  -T $threads -i $dbCOBS -f $f -t $t > $f.cobs@$db.txt2" \
            2>$f.cobs@$db.txt2.time
        
        # kmcp
        # --load-whole-db
        memusg -t -H -s "kmcp search --load-whole-db -j $threads -d $dbKMCP    $f -t $t --quiet > $f.kmcp@$db.txt2" \
            2>$f.kmcp@$db.txt2.time
    done

    find refs/ -name "*.fasta.short" \
        | rush -k 'echo -ne {%}; \
            echo -en "\t"$(grep "elapsed time" {}.cobs@gtdb.txt2.time | sed -r "s/.+: //"); \
            echo -en "\t"$(grep "peak rss" {}.cobs@gtdb.txt2.time | sed -r "s/.+: //"); \
            echo -en "\t"$(grep "elapsed time" {}.kmcp@gtdb.txt2.time | sed -r "s/.+: //"); \
            echo -e "\t"$(grep "peak rss" {}.kmcp@gtdb.txt2.time | sed -r "s/.+: //"); \
            ' \
        | sed 's/ /\t/g' \
        | csvtk add-header -t -n "query,cobs_t,cobs_m0,kmcp_t,kmcp_m0" \
        | csvtk mutate2 -t -n cobs_m -e '$cobs_m0 / 1048576' \
        | csvtk mutate2 -t -n kmcp_m -e '$kmcp_m0 / 1048576' \
        | csvtk cut -t -f query,cobs_t,cobs_m,kmcp_t,kmcp_m \
        | csvtk rename -t -f 1-5 -n "query,cobs:time(s),cobs:memory(GB),kmcp:time(s),kmcp:memory(GB)" \
        > bench.kmcp-cobs-loadall.short.tsv
        
    csvtk csv2md -t bench.kmcp-cobs-loadall.short.tsv
    
For a database of large number of reference genomes, loading complete database does not help increase the speed.

|query                    |cobs:time(s)|cobs:memory(GB)|kmcp:time(s)|kmcp:memory(GB)|
|:------------------------|:-----------|:--------------|:-----------|:--------------|
|NC_000913.3.fasta.short  |727.551     |140.11         |65.974      |70.08          |
|NC_002695.2.fasta.short  |855.974     |149.92         |67.339      |69.61          |
|NC_011750.1.fasta.short  |752.955     |145.73         |65.929      |69.10          |
|NC_012971.2.fasta.short  |673.242     |139.16         |60.343      |70.77          |
|NC_013654.1.fasta.short  |612.932     |140.97         |61.690      |63.97          |
|NC_018658.1.fasta.short  |771.647     |147.33         |69.084      |71.66          |
|NZ_CP007592.1.fasta.short|734.883     |145.41         |63.677      |65.70          |
|NZ_CP028116.1.fasta.short|887.062     |151.63         |73.977      |68.63          |

## KMCP vs Mash and Sourmash


### Building database (MinHash sketches)

    seqs=gtdb202
    db=gtdb
    k=31
    threads=8
    scale=1000
    scaleMash=3400 # average genome size: 3.4Mb.
    
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
|database size  |1.22GB | 5.19GB   | 1.52GB |
|buiding time   |13m30s | 40m39s   | 7min59s|
|temporary files|-      |-         | 1.85GB |


### Searching with bacterial genomes
    
    # emptying the buffers cache
    # su -c "free && sync && echo 3 > /proc/sys/vm/drop_caches && free"
    
    threads=8
    # threads=1
    
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
        
        # kmcp
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
        > bench.kmcp-mash-sourmash.thread$threads.tsv
        
    csvtk csv2md -t bench.kmcp-mash-sourmash.thread$threads.tsv
        
    
Query time and peak memory (mash and kmcp utilize single thread, cold start)

|query              |mash:time(s)|mash:mem(MB)|sourmash:time(s)|sourmash:mem(MB)|kmcp:time(s)|kmcp:mem(MB)|
|:------------------|:-----------|:-----------|:---------------|:---------------|:-----------|:-----------|
|NC_000913.3.fasta  |6.578       |2511.22     |11.669          |126.69          |20.990      |681.25      |
|NC_002695.2.fasta  |6.691       |2518.22     |6.271           |126.57          |4.869       |737.11      |
|NC_011750.1.fasta  |6.923       |2517.48     |3.979           |126.57          |1.896       |707.25      |
|NC_012971.2.fasta  |6.130       |2515.52     |3.344           |126.56          |1.280       |668.68      |
|NC_013654.1.fasta  |5.892       |2515.63     |3.383           |126.56          |1.250       |539.38      |
|NC_018658.1.fasta  |6.126       |2518.09     |3.966           |126.58          |1.472       |640.68      |
|NZ_CP007592.1.fasta|6.661       |2517.76     |3.535           |126.57          |1.462       |690.38      |
|NZ_CP028116.1.fasta|6.507       |2516.80     |3.760           |126.57          |1.248       |657.16      |

Query time and peak memory (mash and kmcp utilize single thread, hot start)

|query              |mash:time(s)|mash:mem(MB)|sourmash:time(s)|sourmash:mem(MB)|kmcp:time(s)|kmcp:mem(MB)|
|:------------------|:-----------|:-----------|:---------------|:---------------|:-----------|:-----------|
|NC_000913.3.fasta  |6.280       |2516.58     |3.373           |126.57          |0.998       |425.38      |
|NC_002695.2.fasta  |6.224       |2517.02     |3.539           |126.58          |1.236       |494.00      |
|NC_011750.1.fasta  |6.869       |2519.32     |3.539           |126.57          |1.045       |484.11      |
|NC_012971.2.fasta  |6.881       |2516.20     |3.335           |126.56          |1.045       |592.29      |
|NC_013654.1.fasta  |6.534       |2516.38     |3.370           |126.57          |1.068       |636.11      |
|NC_018658.1.fasta  |6.931       |2517.14     |3.590           |126.57          |1.039       |477.69      |
|NZ_CP007592.1.fasta|5.685       |2518.06     |3.549           |126.57          |1.252       |550.79      |
|NZ_CP028116.1.fasta|6.259       |2518.54     |3.545           |126.58          |1.227       |481.67      |

Query time and peak memory (mash and kmcp utilize 8 threads, cold start)

|query              |mash:time(s)|mash:mem(MB)|sourmash:time(s)|sourmash:mem(MB)|kmcp:time(s)|kmcp:mem(MB)|
|:------------------|:-----------|:-----------|:---------------|:---------------|:-----------|:-----------|
|NC_000913.3.fasta  |4.258       |2512.18     |11.631          |126.69          |5.038       |683.47      |
|NC_002695.2.fasta  |4.214       |2519.62     |6.306           |126.58          |2.255       |723.79      |
|NC_011750.1.fasta  |4.154       |2519.36     |3.983           |126.57          |0.831       |664.24      |
|NC_012971.2.fasta  |4.157       |2512.34     |3.340           |126.56          |0.621       |82.98       |
|NC_013654.1.fasta  |4.166       |2511.69     |3.351           |126.56          |0.569       |680.69      |
|NC_018658.1.fasta  |4.250       |2519.62     |4.014           |126.57          |0.619       |139.04      |
|NZ_CP007592.1.fasta|4.127       |2519.41     |3.538           |126.57          |0.585       |41.66       |
|NZ_CP028116.1.fasta|4.197       |2519.88     |3.757           |126.57          |0.622       |213.89      |

Query time and peak memory (mash and kmcp utilize 8 threads, hot start)

|query              |mash:time(s)|mash:mem(MB)|sourmash:time(s)|sourmash:mem(MB)|kmcp:time(s)|kmcp:mem(MB)|
|:------------------|:-----------|:-----------|:---------------|:---------------|:-----------|:-----------|
|NC_000913.3.fasta  |4.126       |2518.59     |3.333           |126.56          |0.619       |385.54      |
|NC_002695.2.fasta  |4.193       |2519.67     |3.553           |126.58          |0.620       |230.27      |
|NC_011750.1.fasta  |4.220       |2519.45     |3.537           |126.57          |0.595       |349.68      |
|NC_012971.2.fasta  |4.110       |2518.58     |3.375           |126.56          |0.606       |301.81      |
|NC_013654.1.fasta  |4.415       |2518.90     |3.532           |126.56          |0.621       |185.03      |
|NC_018658.1.fasta  |4.193       |2519.42     |3.550           |126.57          |0.532       |639.11      |
|NZ_CP007592.1.fasta|4.192       |2519.41     |3.372           |126.57          |0.599       |412.88      |
|NZ_CP028116.1.fasta|4.144       |2519.93     |3.583           |126.57          |0.572       |390.53      |
