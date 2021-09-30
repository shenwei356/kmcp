## For one genome


    step=10
    len=100
    
    file=GCF_003697165.2.fna.gz # E.coli
    
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
    
    # filtering searching results
    kmcp filter -T taxid.map -X ~/.taxonkit \
        $file.fa.gz.kmcp.tsv.gz -o $file.fa.gz.kmcp.uniq.tsv.gz
        
    # merge regions
    kmcp utils merge-regions -g $step $file.fa.gz.kmcp.uniq.tsv.gz \
        -o $file.fa.gz.kmcp.uniq.tsv.gz.bed


## For all genomes

    # --------------------------------------------------
    # preprare files

    mkdir -p uniq
    
    mkdir -p uniq/gtdb
    cd uniq/gtdb
    find ../../gtdb/gtdb -name "*.fna.gz" | rush 'ln -s {}'
    cd ../../
    
    mkdir -p uniq/viruses
    cd uniq/viruses
    find ../../refseq-viral/2021-07-30_21-51-31/files.renamed/ -name "*.fna.gz" | rush 'ln -s {}'
    cd ../../
    
    mkdir -p uniq/fungi
    cd uniq/fungi
    find ../../refseq-fungi/2021-07-30_21-54-19/files.renamed/ -name "*.fna.gz" | rush 'ln -s {}'
    cd ../../
    
    # --------------------------------------------------
    # search and filter
    
    step=10
    len=100
    dbs="gtdb.kmcp refseq-fungi.kmcp refseq-viruses.kmcp"
    taxidmap=taxid.map
    taxdump=taxdump/
    
    j=4
    J=40
    find uniq/gtdb -name "*.fna.gz" \
        | rush -j $j -v j=$J -v step=$step -v len=$len -v "dbs=$dbs" \
            -v taxidmap=$taxidmap -v taxdump=$taxdump \
            'for db in {dbs}; do \
                seqkit grep -n -i -v -p plasmid {} \
                    | seqkit sliding --step {step} --window {len} -w 0 \
                    | kmcp search -j {j} -d $db -o {}.kmcp@$db.tsv.gz; \
            done; \
            kmcp merge {}.kmcp@*.tsv.gz \
                | kmcp filter -T {taxidmap} -X {taxdump} -o {}.kmcp.uniq.tsv.gz; \
            /bin/rm {}.kmcp@*.tsv.gz; ' -c -C search.rush
    
    
