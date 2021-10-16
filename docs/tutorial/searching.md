# Sequence and genome searching

## Sequence search

KMCP can be used for fast sequence search from large scales of genomic dataset
as [BIGSI](https://github.com/Phelimb/BIGSI) and [COBS](https://github.com/bingmann/cobs) do.

KMCP reimplemented and modified the Compact Bit-Sliced Signature index (COBS) algorithm,
bringing a small database size and much faster searching speed
 (check the [benchmark](/benchmark/searching)).


### Step 1. Building databases

The database building process is similar to that of metagenomic profiling,
with one difference:

1. Reference genomes are not splitted into fragments which slow down searching speed.

Taken GTDB for example:

    # mask low-complexity region
    mkdir -p gtdb.masked
    find gtdb/ -name "*.fna.gz" \
        | rush 'dustmasker -in <(zcat {}) -outfmt fasta \
            | sed -e "/^>/!s/[a-z]/n/g" \
            | gzip -c > gtdb.masked/{%}'

    # compute k-mers
    #   sequence containing "plasmid" in name are ignored,
    #   k = 21
    kmcp compute -I gtdb.masked/ -k 21 -B plasmid -O gtdb-r202-k21 --force

    # build database
    #   number of index files: 32, for server with >= 32 CPU cores
    #   bloom filter parameter:
    #     number of hash function: 1
    #     false positive rate: 0.3
    kmcp index -j 32 -I gtdb-r202-k21 -O gtdb-r202-k21-n 1 -f 0.3
    
    # cp name mapping file to database directory
    cp name.map gtdb-r202-k21.kmcp/

The size of database is 56GB. 

By default, `kmcp search` loads the whole database into main memory (RAM) for fast searching.
Optionally, the flag `--low-mem` can be set to avoid loading the whole database,
while it's much slower, >10X slower on SSD and should be much slower on HDD disks.

### Step 2. Searching

`kmcp search` supports FASTA/Q format from STDIN or files ([usage](/usage/#search)).

    kmcp search -d gtdb-r202-k21.kmcp/ test.fq.gz -o result.tsv.gz

    13:36:14.522 [INFO] kmcp v0.7.0
    13:36:14.522 [INFO]   https://github.com/shenwei356/kmcp
    13:36:14.522 [INFO] 
    13:36:14.522 [INFO] checking input files ...
    13:36:14.522 [INFO]   1 input file(s) given
    13:36:14.522 [INFO] loading database with mmap enabled ...
    13:36:14.523 [INFO] number of extra workers for every index file: 3
    13:36:14.571 [INFO] database loaded: gtdb-r202-k21.kmcp/
    13:36:14.571 [INFO] 
    13:36:14.571 [INFO] -------------------- [main parameters] --------------------
    13:36:14.571 [INFO]   minimum    query length: 70
    13:36:14.571 [INFO]   minimum  matched k-mers: 30
    13:36:14.571 [INFO]   minimum  query coverage: 0.550000
    13:36:14.571 [INFO]   minimum target coverage: 0.000000
    13:36:14.571 [INFO]   minimum target coverage: 0.000000
    13:36:14.571 [INFO] -------------------- [main parameters] --------------------
    13:36:14.571 [INFO] 
    13:36:14.571 [INFO] searching ...
    13:36:14.575 [INFO] reading sequence file: test.fq.gz
    processed queries: 999424, speed: 2.47 million queries per minute
    13:36:38.843 [INFO] 
    13:36:38.843 [INFO] processed queries: 1000000, speed: 2.47 million queries per minute
    13:36:38.843 [INFO] 68.2500% (682500/1000000) queries matched
    13:36:38.843 [INFO] done searching
    13:36:42.852 [INFO] 
    13:36:42.853 [INFO] elapsed time: 28.331134572s
    13:36:42.853 [INFO]

The result is tab-delimited.

|#query|qLen|qKmers|FPR       |hits|target         |fragIdx|frags|tLen   |kSize|mKmers|qCov  |tCov  |jacc  |queryIdx|
|:-----|:---|:-----|:---------|:---|:--------------|:------|:----|:------|:----|:-----|:-----|:-----|:-----|:-------|
|S0R0/1|150 |130   |3.0168e-03|2   |GCF_002872255.1|0      |1    |2582291|21   |88    |0.6769|0.0000|0.0000|0       |
|S0R0/1|150 |130   |3.0168e-03|2   |GCF_001434585.1|0      |1    |2219511|21   |74    |0.5692|0.0000|0.0000|0       |
|S0R0/2|150 |130   |3.0168e-03|5   |GCF_002872255.1|0      |1    |2582291|21   |81    |0.6231|0.0000|0.0000|1       |
|S0R0/2|150 |130   |3.0168e-03|5   |GCF_001434585.1|0      |1    |2219511|21   |80    |0.6154|0.0000|0.0000|1       |
|S0R0/2|150 |130   |3.0168e-03|5   |GCF_001437405.1|0      |1    |2260906|21   |77    |0.5923|0.0000|0.0000|1       |

Reference IDs (column `target`) can be optionally mapped to their names, let's print the main columns only:

    csvtk head -t -C $ -n 5 result.tsv.gz \
        | csvtk rename -t -C $ -f 1 -n query \
        | csvtk cut -t -f query,FPR,qCov,target \
        | csvtk csv2md -t 

|query |FPR       |qCov  |target         |
|:-----|:---------|:-----|:--------------|
|S0R0/1|3.0168e-03|0.6769|GCF_002872255.1|
|S0R0/1|3.0168e-03|0.5692|GCF_001434585.1|
|S0R0/2|3.0168e-03|0.6231|GCF_002872255.1|
|S0R0/2|3.0168e-03|0.6154|GCF_001434585.1|
|S0R0/2|3.0168e-03|0.5923|GCF_001437405.1|

## Genome similarity estimation

KMCP can be used for **fast similarity estimation of newly assembled genome against known reference genomes**.
BLAST is an option while it focuses on local similarity.
Besides, the database is big and the searching speed is relative slow.

Genome sketching is a method of utilizing small and approximate summaries of
genomic data for fast searching and comparison.
[Mash](https://github.com/marbl/Mash) and [Sourmash](https://github.com/sourmash-bio/sourmash)
provide fast genome distance estimation using MinHash (Mash) or Scaled MinHash (Sourmash).
Here KMCP utilizes multiple sketches 
([Minimizer](https://academic.oup.com/bioinformatics/article/20/18/3363/202143), 
[Scaled MinHash](https://f1000research.com/articles/8-1006) and
[Syncmers](https://peerj.com/articles/10805/)) for genome similarity estimation.

[Prebuilt databases](/kmcp/database) are available and users can also build custom databases following steps below.

### Step 1. Building databases

The database building process is similar to that of metagenomic profiling,
with few differences:

1. K-mer sketches instead of all k-mers are computed.
2. Reference genomes are not splitted into fragments.
3. Smaller false positive rate (`-f`) and more than one hash functions (`-n`) are used to improve accuracy.
4. Using bigger k-mer size: 31.

Supported k-mer (sketches) types:

1. K-mer:
    - ntHash of k-mer (`-k`)
2. K-mer sketchs (all using ntHash):
    - Scaled MinHash (`-k -D`)
    - Minimizer      (`-k -W`), optionally scaling/down-sampling (`-D`)
    - Closed Syncmer (`-k -S`), optionally scaling/down-sampling (`-D`)

Taken GTDB for example:
    
    # compute Scaled MinHash with scale 1000
    #   sequence containing "plasmid" in name are ignored,
    #   k = 31
    kmcp compute -I gtdb/ -k 31 -D 1000 -B plasmid -O gtdb-r202-minhash

    # build database
    #   number of index files: 8, for server with >= 8 CPU cores
    #   bloom filter parameter:
    #     number of hash function: 3
    #     false positive rate: 0.001
    kmcp index -j 8 -I gtdb-r202-minhash -O gtdb.minhash.kmcp -n 3 -f 0.001
    
    # cp name mapping file to database directory
    cp taxid.map name.map gtdb.minhash.kmcp/

Attention:

- For small genomes like viruses, sketching parameters should be adjusted. 
For examples, setting a smaller scale like 10.
    
### Step 2. Searching

The searching process is simple and [very fast](https://bioinf.shenwei.me/kmcp/benchmark) (<1 second).

    kmcp search --query-whole-file -d gtdb.minhash.kmcp/ contigs.fasta -o result.tsv

The output is in tab-delimited format, full search result:

|#query                             |qLen   |qKmers|FPR         |hits|target         |fragIdx|frags|tLen   |kSize|mKmers|qCov  |tCov  |jacc  |queryIdx|
|:----------------------------------|:------|:-----|:-----------|:---|:--------------|:------|:----|:------|:----|:-----|:-----|:-----|:-----|:-------|
|NODE_1_length_106749_cov_161.610981|4563649|9052  |0.000000e+00|1   |GCF_900460465.1|0      |1    |4777134|31   |8942  |0.9878|0.9933|0.9813|0       |

Reference IDs can be optionally mapped to their names, let's print the main columns only:

    kmcp search --query-whole-file -d gtdb.minhash.kmcp/ -N gtdb.minhash.kmcp/name.map \
        contigs.fasta --quiet \
        | csvtk rename -t -C $ -f 1 -n query \
        | csvtk cut -t -f query,jacc,target \
        > result.tsv
    
|query                              |jacc  |target                                                                          |
|:----------------------------------|:-----|:-------------------------------------------------------------------------------|
|NODE_1_length_106749_cov_161.610981|0.9813|NZ_UAVH01000012.1 Yersinia pestis strain NCTC5923, whole genome shotgun sequence|

Using closed syncmer:

    kmcp search --query-whole-file -d gtdb.syncmer.kmcp/ -N gtdb.syncmer.kmcp/name.map \
        contigs.fasta --quiet \
        | csvtk rename -t -C $ -f 1 -n query \
        | csvtk cut -t -f query,jacc,target
        
|query                              |jacc  |target                                                                          |
|:----------------------------------|:-----|:-------------------------------------------------------------------------------|
|NODE_1_length_106749_cov_161.610981|0.9383|NZ_UAVH01000012.1 Yersinia pestis strain NCTC5923, whole genome shotgun sequence|
