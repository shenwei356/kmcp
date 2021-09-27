# Searching benchmarks

[Software, datasets and commands details](https://github.com/shenwei356/kmcp/tree/main/benchmarks/searching).

Softwares

- [cobs](https://github.com/bingmann/cobs) ([1915fc0](https://github.com/bingmann/cobs/commit/1915fc061bbe47946116b4a051ed7b4e3f3eca15))
- [sourmash](https://github.com/dib-lab/sourmash) (v4.2.1)
- kmcp ([v0.6.0](https://github.com/shenwei356/kmcp/releases/tag/v0.6.0))

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



Searching with short reads (~200K reads) (hot start, loading all data in memory)):

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

Only Scaled MinHash (scale=1000) are indexed and searched.

Database size and building time:

|               |mash   |sourmash  |kmcp    |
|:--------------|:------|:---------|:-------|
|database size  |743MB  | 5.19GB   | 1.52GB |
|buiding time   |11m39s | 89m59s   | 7min02s|
|temporary files|-      |-         | 3.41GB |


Searching with bacterial genomes (mash and kmcp utilize single thread, cold start)

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


Searching with bacterial genomes (mash and kmcp utilize 8 threads, cold start)

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
