# Sequence and genome searching


## Sequence search

Similar tools

- [BIGSI](https://github.com/Phelimb/BIGSI) and [COBS](https://github.com/bingmann/cobs)
- [Mantis](https://github.com/splatlab/mantis)


## Genome similarity estimation

KMCP can be used for **fast similarity estimation of newly assembled genome against known reference genomes**.
BLAST is an option while it focus on local similarity,
while the database is big and searching speed is relative slow.

Genome sketching is a method of utilizing small and approximate summaries of
genomic data for faster searching and comparison.
[Mash](https://github.com/marbl/Mash) and [Sourmash](https://github.com/sourmash-bio/sourmash)
provide fast genome distance estimation using MinHash (Mash) or Scaled MinHash (Sourmash).
Here KMCP utilizes multiple sketches (Minimizer, Scaled MinHash and Syncmers) for
genome similarity estimation.

[Prebuilt databases] are available and users can also build custom databases.

### Step 1. Building database of genome sketches.

The database building process is similar to that of metagenomic profiling,
with few differences:

1. K-mer sketches instead of all k-mers are computed.
2. Reference genomes are not splitted into fragements.
3. Smaller false positive rate (`-f`) and more than one hash functions (`-n`) are used to increase accuracy.
4. Using bigger k-mer size: 31.

Taken GTDB for example:

    # mask low-complexity region
    mkdir -p gtdb.masked
    find gtdb/ -name "*.fna.gz" \
        | rush 'dustmasker -in <(zcat {}) -outfmt fasta \
            | sed -e "/^>/!s/[a-z]/n/g" \
            | gzip -c > gtdb.masked/{%}'
    
    # compute Scaled MinHash with scale 1000
    #   sequence containing "plasmid" in name are ignored,
    #   k = 31
    kmcp compute -I gtdb.masked/ -k 31 -D 1000 -B plasmid -O gtdb-r202-minhash

    # build database
    #   number of index files: 8, for server with >= 8 CPU cores
    #   bloom filter parameter:
    #     number of hash function: 3
    #     false positive rate: 0.001
    kmcp index -j 8 -I gtdb-r202-minhash -O gtdb.minhash.kmcp -n 3 -f 0.001
    
    # cp name mapping file to database directory
    cp name.map gtdb.minhash.kmcp/

Attention: for small genomes like viruses, sketching parameters should be adjusted. 
For examples, setting smaller scale 10.
    
### Step 2. Searching

The searching process is simple and veryfast (<1 second).

    kmcp search --query-whole-file -d gtdb.minhash.kmcp/ contigs.fasta -o result.tsv

The output is in tab-delimited format, full search result:

|#query                             |qLen   |qKmers|FPR         |hits|target         |fragIdx|frags|tLen   |kSize|mKmers|qCov  |tCov  |jacc  |queryIdx|
|:----------------------------------|:------|:-----|:-----------|:---|:--------------|:------|:----|:------|:----|:-----|:-----|:-----|:-----|:-------|
|NODE_1_length_106749_cov_161.610981|4563649|9052  |0.000000e+00|1   |GCF_900460465.1|0      |1    |4777134|31   |8942  |0.9878|0.9933|0.9813|0       |

Reference IDs can be optional mapped to their names, let's print the main columns only:

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
