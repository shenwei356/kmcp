## For one genome

Generating reads

    len=100
    
    file=GCF_003697165.2.fasta.gz
    
    # splitting genomes
    seqkit sliding -s 10 -W $len -w 0 $file -o $file.fa.gz
    
Searching

    kmcp search -d gtdb.kmcp/ $file.fa.gz -o $file.fa.gz.kmcp@gtdb.tsv.gz    
    kmcp search -d refseq-fungi.kmcp/ $file.fa.gz -o $file.fa.gz.kmcp@refseq-fungi.tsv.gz    
    kmcp search -d refseq-viruses.kmcp/ $file.fa.gz -o $file.fa.gz.kmcp@refseq-viruses.tsv.gz
    
    kmcp merge $file.fa.gz.kmcp@* -o $file.fa.gz.kmcp.tsv.gz
    
Filtering

    kmcp filter -T gtdb.kmcp/taxid.map -X ~/.taxonkit $file.fa.gz.kmcp.tsv.gz -o $file.fa.gz.kmcp.uniq.tsv.gz

## For all genomes


    
