# kmcp

Kmer-based Metagenomics Classification and Profiling

## Quick Start

    # compute k-mer sketch
    time kmcp compute -i <(find genomes/ -name "*fa.gz") -O kmcp-db -W 16 --force --verbose

    # index sketch
    time kmcp index -I kmcp-db/ -a db-alias --verbose
    
    # search    
    time kmcp search -d kmcp-db/ -d kmcp-db2/ -t 1 t_t2.fq.gz -o result.tsv --verbose 
