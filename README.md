# KMCP - K-mer-based Metagenomic Classification and Profiling

<img src="kmcp.png" alt="" width="800"/>

## What can we do?

### 1. Accurate metagenomic profiling and binning

xxx

(check the [benchmark](https://bioinf.shenwei.me/kmcp/benchmark/profiling))


### 2. Fast sequence search from large scales of genomic datasets

KMCP can be used for fast sequence search from large scales of genomic dataset
as [BIGSI](https://github.com/Phelimb/BIGSI) and [COBS](https://github.com/bingmann/cobs) do.
We reimplemented and modified the Compact Bit-Sliced Signature index (COBS) algorithm,
bringing a small database size and much faster searching speed
 (check the [tutorial](https://bioinf.shenwei.me/kmcp/tutorial/searching) and [benchmark](https://bioinf.shenwei.me/kmcp/benchmark/searching)).
 
### 3. Fast genome similarity estimation

KMCP can be used for fast similarity estimation of newly assembled genomes against known reference genomes.

Genome sketching is a method of utilizing small and approximate summaries of
genomic data for fast searching and comparison.
[Mash](https://github.com/marbl/Mash) and [Sourmash](https://github.com/sourmash-bio/sourmash)
provide fast genome distance estimation using MinHash (Mash) or Scaled MinHash (Sourmash).
Here KMCP utilizes multiple sketches 
([Minimizer](https://academic.oup.com/bioinformatics/article/20/18/3363/202143), 
[Scaled MinHash](https://f1000research.com/articles/8-1006) and
[Closed Syncmers](https://peerj.com/articles/10805/)) for genome similarity estimation
 (check the [tutorial](https://bioinf.shenwei.me/kmcp/tutorial/searching) and [benchmark](https://bioinf.shenwei.me/kmcp/benchmark/searching)).


## Tutorials and resources:

- [Installation](https://bioinf.shenwei.me/kmcp/download)
- [Databases](https://bioinf.shenwei.me/kmcp/database)
- Tutorials
    - [Taxonomic profiling](https://bioinf.shenwei.me/kmcp/tutorial/profiling)
    - [Sequence and genome searching](https://bioinf.shenwei.me/kmcp/tutorial/searching)
- [Usage](https://bioinf.shenwei.me/kmcp/usage)
- [Benchmarks](https://bioinf.shenwei.me/kmcp/benchmark)
- [FAQs](https://bioinf.shenwei.me/kmcp/faq)

## Installation

Download [executable binaries](https://github.com/shenwei356/kmcp/releases),
or install using conda:

    conda install -c bioconda kmcp

## Quick Start

    # compute k-mers
    kmcp compute -k 21 --split-number 10 --split-overlap 100 \
        --in-dir genomes/ --out-dir genomes-k21-n10

    # index k-mers
    kmcp index --in-dir genomes-k21-n10/ --out-dir genomes.kmcp
    
    # delete temporary files
    # rm -rf genomes-k21-n10/
    
    # search    
    kmcp search --db-dir genomes.kmcp/ test.fa.gz --out-file search.tsv.gz

    # profile and binning
    kmcp profile search.tsv.gz \
        --taxid-map        taxid.map \
        --taxdump          taxdump/ \
        --out-prefix       search.tsv.gz.k.profile \
        --metaphlan-report search.tsv.gz.m.profile \
        --cami-report      search.tsv.gz.c.profile \
        --binning-result   search.tsv.gz.binning.gz

## Support

Please [open an issue](https://github.com/shenwei356/kmcp/issues) to report bugs,
propose new functions or ask for help.

## License

[MIT License](https://github.com/shenwei356/kmcp/blob/master/LICENSE)

## Acknowledgements

- [Zhi-Luo Deng](https://dawnmy.github.io/CV/) (Helmholtz Centre for Infection Research, Germany)
  gave many valuable advice on metagenomic profiling and benchmarking.
- [Robert Clausecker](https://github.com/clausecker/) (Zuse Institute Berlin, Germany)
  wrote the high-performance vectorized positional popcount package 
  ([pospop](https://github.com/clausecker/pospop)) 
  [during my development of KMCP](https://stackoverflow.com/questions/63248047/),
  which greatly accelerated the bit-matrix searching.
