# kmcp

K-mer-based Metagenomic Classification and Profiling

## Quick Start

    # compute k-mers
    kmcp compute -I genomes/ -O genomes-k21 -k 21 --force

    # index sketch
    kmcp index -I genomes-k21/ -O genomes-k21.kmcp --force
    
    # search    
    kmcp search -d genomes-k21.kmcp/ -t 0.7 test.fa.gz -o search.tsv.gz

    # profile
    kmcp profile search.tsv.gz -T taxid.map -X taxdump \
        -o search.tsv.k.profile \
        --metaphlan-report search.tsv.m.profile \
        --cami-report search.tsv.c.profile

