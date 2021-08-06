# KMCP

A K-mer-based Metagenomics Classification and Profiling tool.

Tutorials and resources:

- [Installation](https://bioinf.shenwei.me/kmcp/download)
- [Database]((https://bioinf.shenwei.me/kmcp/database))
- Tutorials
    - [Taxonomic Profiling](https://bioinf.shenwei.me/kmcp/profiling)
    - [Sequence and genome searching](https://bioinf.shenwei.me/kmcp/searching)
- [Usage](https://bioinf.shenwei.me/kmcp/usage)
- [Benchmarks](https://bioinf.shenwei.me/kmcp/benchmark)
- [FAQs](https://bioinf.shenwei.me/kmcp/faq)

## Installation

Download [executable binaries](https://github.com/shenwei356/kmcp/releases),
or install using conda:

    conda install -c conda-forge -c bioconda kmcp

## Quick Start

    # compute k-mers
    kmcp compute -I genomes/ -O genomes-k21 -k 21 --force

    # index k-mers
    kmcp index -I genomes-k21/ -O genomes-k21.kmcp --force
    
    # search    
    kmcp search -d genomes-k21.kmcp/ test.fa.gz -o search.tsv.gz

    # profile
    kmcp profile search.tsv.gz -T taxid.map -X taxdump \
        --out-prefix search.tsv.k.profile \
        --metaphlan-report search.tsv.m.profile \
        --cami-report search.tsv.c.profile

## Support

Please [open an issue](https://github.com/shenwei356/kmcp/issues) to report bugs,
propose new functions or ask for help.

## License

[MIT License](https://github.com/shenwei356/kmcp/blob/master/LICENSE)

## Credits

- [Zhi-Luo Deng](https://dawnmy.github.io/CV/) (Helmholtz Centre for Infection Research, Germany)
  gave many valuable advice on metagenomic profiling tool development and benchmarking.
- [Robert Clausecker](https://github.com/clausecker/) (Zuse Institute Berlin, Germany)
  wrote the high-performance vectorized positional popcount package 
  ([pospop](https://github.com/clausecker/pospop)) 
  [during my development KMCP](https://stackoverflow.com/questions/63248047/),
  which greatly accelerated bit-matrix searching.
