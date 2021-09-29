## For one genome


    len=100
    
    file=GCF_003697165.2.fasta.gz
    
    file=GCF_002402265.1.fna.gz # EB virus
    
    file=GCF_000861825.2.fna.gz # HBV
        
    # search
    seqkit sliding -s 10 -W $len -w 0 $file \
        | kmcp search -d gtdb.kmcp/ -o $file.fa.gz.kmcp@gtdb.tsv.gz
        
    seqkit sliding -s 10 -W $len -w 0 $file \
        | kmcp search -d refseq-fungi.kmcp/ -o $file.fa.gz.kmcp@refseq-fungi.tsv.gz
        
    seqkit sliding -s 10 -W $len -w 0 $file \
        | kmcp search -d refseq-viruses.kmcp/ -o $file.fa.gz.kmcp@refseq-viruses.tsv.gz
    
    # merge
    kmcp merge $file.fa.gz.kmcp@*.tsv.gz -o $file.fa.gz.kmcp.tsv.gz
    
    cat gtdb.kmcp/taxid.map refseq-fungi.kmcp/taxid.map refseq-viruses.kmcp/taxid.map \
        > taxid.map        
    cat gtdb.kmcp/name.map refseq-fungi.kmcp/name.map refseq-viruses.kmcp/name.map \
        > name.map
    
    # filter
    kmcp filter -T taxid.map -X ~/.taxonkit \
        $file.fa.gz.kmcp.tsv.gz -o $file.fa.gz.kmcp.uniq.tsv.gz
        
    # export regions
    kmcp utils merge-regions $file.fa.gz.kmcp.uniq.tsv.gz \
        -o $file.fa.gz.kmcp.uniq.tsv.gz.bed


## For all genomes


    
