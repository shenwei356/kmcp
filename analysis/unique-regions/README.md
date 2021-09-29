## For one genome


    step=10
    len=100
    
    file=GCF_003697165.2.fasta.gz # E.coli
    
    file=GCF_002402265.1.fna.gz # EB virus
    
    file=GCF_000861825.2.fna.gz # HBV
        
    # simulating reads and searching
    
    seqkit grep -n -i -v -p plasmid  $file \
        | seqkit sliding --step $step --window $len -w 0 \
        | kmcp search -d gtdb.kmcp/ -o $file.fa.gz.kmcp@gtdb.tsv.gz
        
    seqkit grep -n -i -v -p plasmid  $file \
        | seqkit sliding --step $step --window $len -w 0 \
        | kmcp search -d refseq-fungi.kmcp/ -o $file.fa.gz.kmcp@refseq-fungi.tsv.gz
        
    seqkit grep -n -i -v -p plasmid  $file \
        | seqkit sliding --step $step --window $len -w 0 \
        | kmcp search -d refseq-viruses.kmcp/ -o $file.fa.gz.kmcp@refseq-viruses.tsv.gz
    
    # merging searching results
    kmcp merge $file.fa.gz.kmcp@*.tsv.gz -o $file.fa.gz.kmcp.tsv.gz
    
    cat gtdb.kmcp/taxid.map refseq-fungi.kmcp/taxid.map refseq-viruses.kmcp/taxid.map \
        > taxid.map        
    cat gtdb.kmcp/name.map refseq-fungi.kmcp/name.map refseq-viruses.kmcp/name.map \
        > name.map
    
    # filtering searching results
    kmcp filter -T taxid.map -X ~/.taxonkit \
        $file.fa.gz.kmcp.tsv.gz -o $file.fa.gz.kmcp.uniq.tsv.gz
        
    # merge regions
    kmcp utils merge-regions -g $step $file.fa.gz.kmcp.uniq.tsv.gz \
        -o $file.fa.gz.kmcp.uniq.tsv.gz.bed


## For all genomes


    
