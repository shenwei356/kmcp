## GTDB

Dataset: [gtdb_genomes_reps_r95.tar.gz](https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/genomic_files_reps/gtdb_genomes_reps_r95.tar.gz)

    # uncompress
    mkdir -p gtdb
    tar -zxvf gtdb_genomes_reps_r95.tar.gz -O gtdb
    
    # rename
    brename -R -p _genomic.fna -r .fa gtdb
    
Name mapping file:

    # sequence accesion -> full head
    find gtdb -name *.fa.gz \
        | rush -k 'echo -ne "{%@(.+).fa}\t$(seqkit head -n 1 {} | seqkit seq -n)\n" ' \
        > gtdb.ass2name.tsv
    
    # assembly accession -> full head
    find gtdb -name *.fa.gz \
        | rush -k 'seqkit seq -n {}'\
        | csvtk mutate -Ht -p '^(\S+)' \
        | csvtk cut -Ht -f 2,1 > gtdb.acc2name.tsv
    
