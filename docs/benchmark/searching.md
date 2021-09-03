# Searching benchmarks

[Software, datasets and commands details](https://github.com/shenwei356/kmcp/tree/main/benchmarks/searching).

[GTDB r202 representative genomes](https://data.ace.uq.edu.au/public/gtdb/data/releases/release202/202.0/genomic_files_reps/gtdb_genomes_reps_r202.tar.gz) are used for tests:

- file size: 46.26 GB
- files: 47,894
- bases: 151.94 Gb

## KMCP vs COBS

All k-mers are indexed and searched.

Database size and building time:

|               |cobs      |kmcp     |
|:--------------|:---------|:--------|
|database size  | 86.96GB  | 55.15GB |
|building time  | 29m:55s  | 24min52s|
|temporary files| 160.76GB | 1.19TB  |

Searching with bacterial genomes (cold start, loading all data in memory):

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


Searching with short reads (~200K reads) (hot start, loading all data in memory)):

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

Only Scaled MinHash (scale=1000) are indexed and searched.

Database size and building time:

|               |sourmash  |kmcp    |
:---------------|:---------|:-------|
|database size  |  5.19GB  | 1.52GB |
|buiding time   |  89m59s  | 7min02s|
|temporary files| -        | 3.41GB |


Searching with bacterial genomes (kmcp utilizes single thread, cold start)

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


Searching with bacterial genomes (kmcp utilizes 8 threads, cold start)

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

