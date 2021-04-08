# kmcp

K-mer-based Metagenomic Classification and Profiling

## Quick Start

    # compute minimizer sketch
    time kmcp compute -I genomes/ -O genomes.m16 -k 31 -W 16 --force

    # index sketch
    time kmcp index -I genomes.m16/ -O genomes.m16.db --force
    
    # search    
    time kmcp search -d genomes.m16.db/ -t 0.7 test.fa.gz -o search.tsv

    # profile
    time kmcp profile search.tsv -T taxid_mapping.tsv -X taxdata-dir \
        -o search.tsv.k.profile \
        --metaphlan-report search.tsv.m.profile \
        --cami-report search.tsv.c.profile
