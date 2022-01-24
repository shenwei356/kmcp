# Searching benchmarks

[Software, datasets and commands details](https://github.com/shenwei356/kmcp/tree/main/benchmarks/searching).

Softwares

- [COBS](https://github.com/bingmann/cobs) ([1915fc0](https://github.com/bingmann/cobs/commit/1915fc061bbe47946116b4a051ed7b4e3f3eca15))
- [Sourmash](https://github.com/dib-lab/sourmash) (v4.2.2)
- KMCP ([v0.7.0](https://github.com/shenwei356/kmcp/releases/tag/v0.7.0))

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

Searching with bacterial genomes or short reads (~1M reads).

## KMCP vs Mash and Sourmash

Only FracMinHash (Scaled MinHash) (scale=1000 for Sourmash and KMCP) or MinHash (scale=3400 for Mash) are indexed and searched.

Database size and building time:

|               |mash   |sourmash  |kmcp    |
|:--------------|:------|:---------|:-------|
|database size  |743MB  | 5.19GB   | 1.52GB |
|building time  |11m39s | 89m59s   | 7min02s|
|temporary files|-      |-         | 3.41GB |


Searching with bacterial genomes.

## Result

<img src="/benchmark/bench.searching.jpg" alt="" width="600"/>
