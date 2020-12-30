# kmcp

Kmer-based Metagenomics Classification and Profiling

## Quick Start

    # compute minimizer sketch
    time kmcp compute -I genomes -O genomes.m16 -k 31 -W 16 --force --verbose

    # index sketch
    time kmcp index -I genomes.m16 -O genomes.m16.db -a minimizer-16 --verbose
    
    # search    
    time kmcp search -d genomes.m16.db -t 0.9 test.fa.gz -o result.tsv --verbose 
