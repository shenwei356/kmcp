## Generating datasets 

Generating reads

    len=100
    
    file=GCF_003697165.2.fasta.gz
    
    # splitting genomes
    seqkit sliding -s 10 -W $len -w 0 $file -o $file.fa.gz
    
Searching

    kmcp search -d gtdb.kmcp/ $file.fa.gz -o $file.fa.gz.kmcp@gtdb.tsv.gz
    
    kmcp search -d gtdb.kmcp/ $file.fa.gz -o $file.fa.gz.kmcp@gtdb.tsv.gz
