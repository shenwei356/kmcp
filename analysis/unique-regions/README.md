## Create database of Scaled MinHash sketch

Using scale `5`.

GTDB

    k=21
    scale=5
    
    in=gtdb
    out=gtdb.smin$scale.kmcp
    tmp=gtdb-r202-k21-n10-d$scale
    
    kmcp compute -I $in -O $tmp -k $k -n 10 -l 100 -B plasmid -D $scale \
        --log $tmp.log --force        
    kmcp index -j 40 -I $tmp -O $out -n 1 -f 0.3 \
        --log $out.log
    
    cp taxid.map name.map $out/

Refseq-fungi

    k=21
    scale=5
    
    in=files.renamed
    out=refseq-fungi.smin$scale.kmcp
    tmp=refseq-fungi-k21-n10-d$scale
    
    kmcp compute -I $in -O $tmp -k $k -n 10 -l 100 -B plasmid -D $scale \
        --log $tmp.log --force        
    kmcp index -j 40 -I $tmp -O $out -n 1 -f 0.3 \
        --log $out.log
    cp taxid.map name.map $out/
    
Refseq-viral

    k=21
    scale=5
    
    in=files.renamed
    out=refseq-viral.smin$scale.kmcp
    tmp=refseq-viral-k21-n5-d$scale
    
    kmcp compute -I $in -O $tmp -k $k -n 5 -l 100 -B plasmid -D $scale \
        --log $tmp.log --force        
    kmcp index -j 40 -I $tmp -O $out -n 3 -f 0.001 \
        --log $out.log
    cp taxid.map name.map $out/



## For one genome

    step=200
    len=500
    
    
    file=GCF_003697165.2.fna.gz # E.coli
    
    file=GCF_002402265.1.fna.gz # EB virus
    
    file=GCF_000861825.2.fna.gz # HBV
        
    # simulating reads and searching
    
    seqkit grep -n -i -v -r -p plasmid  $file \
        | seqkit sliding --greedy --step $step --window $len -w 0 \
        | kmcp search -d gtdb.kmcp/ -o $file.kmcp@gtdb.tsv.gz
        
    seqkit grep -n -i -v -r -p plasmid  $file \
        | seqkit sliding --greedy --step $step --window $len -w 0 \
        | kmcp search -d refseq-fungi.kmcp/ -o $file.kmcp@refseq-fungi.tsv.gz
        
    seqkit grep -n -i -v -r -p plasmid  $file \
        | seqkit sliding --greedy --step $step --window $len -w 0 \
        | kmcp search -d refseq-viruses.kmcp/ -o $file.kmcp@refseq-viruses.tsv.gz
    
    # merging searching results
    kmcp merge $file.kmcp@*.tsv.gz -o $file.kmcp.tsv.gz
    
    # filtering searching results
    kmcp filter -T taxid.map -X ~/.taxonkit \
        $file.kmcp.tsv.gz -o $file.kmcp.uniq.tsv.gz
        
    # merge regions
    kmcp utils merge-regions -I -l 20 $file.kmcp.uniq.tsv.gz \
        -o $file.kmcp.uniq.tsv.gz.bed


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
    
    step=200
    len=500
    dbs="gtdb.kmcp refseq-fungi.kmcp refseq-viruses.kmcp"
    taxidmap=taxid.map
    taxdump=taxdump/
    
    j=4
    J=40
    find uniq/gtdb -name "*.fna.gz" \
        | rush -j $j -v j=$J -v step=$step -v len=$len -v "dbs=$dbs" \
            -v taxidmap=$taxidmap -v taxdump=$taxdump \
            'for db in {dbs}; do \
                seqkit grep -n -i -v -r -p plasmid {} \
                    | seqkit sliding --greedy --step {step} --window {len} -w 0 \
                    | kmcp search -j {j} -d $db -o {}.kmcp@$db.tsv.gz; \
            done; \
            kmcp merge {}.kmcp@*.tsv.gz \
                | kmcp filter -T {taxidmap} -X {taxdump} -o {}.kmcp.uniq.tsv.gz; \
            /bin/rm {}.kmcp@*.tsv.gz; ' -c -C search.rush
    
    
