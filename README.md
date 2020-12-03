# kmcp

Kmer-based Metagenomics Classification and Profiling

## Quick Start

    # compute k-mer sketch
    time kmcp compute -i <(find genomes/ -name "*fa.gz") -O kmcp-db -W 16 --force --verbose

    # index sketch
    time kmcp index -I kmcp-db/ -W 16 --verbose
    
    # search    
    time kmcp search -m -t 1 -d kmcp-db/ t_t2.fq.gz --verbose > result.tsv
