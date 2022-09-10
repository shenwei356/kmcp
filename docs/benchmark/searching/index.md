# Searching benchmarks

[Software, datasets and commands details](https://github.com/shenwei356/kmcp/tree/main/benchmarks/searching).

Softwares

- [COBS](https://github.com/bingmann/cobs) ([1915fc0](https://github.com/bingmann/cobs/commit/1915fc061bbe47946116b4a051ed7b4e3f3eca15))
- [Sourmash](https://github.com/dib-lab/sourmash) (v4.5.0)
- KMCP ([v0.9.0](https://github.com/shenwei356/kmcp/releases/tag/v0.9.0))

[GTDB r202 representative genomes](https://data.ace.uq.edu.au/public/gtdb/data/releases/release202/202.0/genomic_files_reps/gtdb_genomes_reps_r202.tar.gz) are used for tests:

- file size: 46.26 GB
- files: 47,894
- bases: 151.94 Gb

## KMCP vs COBS

All k-mers are indexed and searched.

Database size and building time:

|               |cobs      |kmcp     |
|:--------------|:---------|:--------|
|database size  | 86.96GB  | **55.15GB** |
|building time  | 29m55s   | **21min04s**|
|temporary files| **160.76GB** | 935.11G |

Searching with bacterial genomes or short reads (~1M reads).

## KMCP vs Mash and Sourmash

Only FracMinHash (Scaled MinHash) (scale=1000 for Sourmash and KMCP) or MinHash (scale=3400 for Mash) are indexed and searched.

Database size and building time:

|               |mash   |sourmash  |kmcp    |
|:--------------|:------|:---------|:-------|
|database size  |**1.22GB** | 5.19GB   | 1.52GB |
|buiding time   |13m30s | 40m39s   | **7min59s**|
|temporary files|-      |-         | 1.85GB |


Searching with bacterial genomes.

## Result

<img src="bench.searching.jpg" alt="" width="600"/>
