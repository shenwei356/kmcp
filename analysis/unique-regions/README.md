## Create database of Scaled MinHash sketch

Using scale `10`.

GTDB

    k=21
    scale=10
    
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
    scale=10
    
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
    scale=10
    
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
    
    seqkit grep -n -i -v -r -p plasmid  $file -w 0 \
        | seqkit sliding --greedy --step $step --window $len -w 0 \
        | kmcp search -j 40 -d gtdb.smin10.kmcp/ -o $file.kmcp@gtdb.tsv.gz
        
    seqkit grep -n -i -v -r -p plasmid  $file -w 0 \
        | seqkit sliding --greedy --step $step --window $len -w 0 \
        | kmcp search -j 40 -d refseq-fungi.smin10.kmcp/ -o $file.kmcp@refseq-fungi.tsv.gz
        
    seqkit grep -n -i -v -r -p plasmid  $file -w 0 \
        | seqkit sliding --greedy --step $step --window $len -w 0 \
        | kmcp search -j 40 -d refseq-viral.smin10.kmcp/ -o $file.kmcp@refseq-viral.tsv.gz
    
    # merging searching results
    kmcp merge $file.kmcp@*.tsv.gz -o $file.kmcp.tsv.gz
    
    # filtering searching results
    kmcp utils filter -T taxid.map -X ~/.taxonkit \
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
    
    mkdir -p uniq/viral
    cd uniq/viral
    find ../../refseq-viral/files.renamed/ -name "*.fna.gz" | rush 'ln -s {}'
    cd ../../
    
    mkdir -p uniq/fungi
    cd uniq/fungi
    find ../../refseq-fungi/files.renamed/ -name "*.fna.gz" | rush 'ln -s {}'
    cd ../../
    
    
    
    
    input=uniq/viral
    
    # --------------------------------------------------
    # search and filter (2days for gtdb, 3hours for fungi)
    
    step=200
    len=500
    dbs="gtdb.smin10.kmcp refseq-fungi.smin10.kmcp refseq-viral.smin10.kmcp"
    taxidmap=taxid.map
    taxdump=taxdump/
            
    j=4
    J=40
    time find $input -name "*fna.gz" \
        | rush -j $j -v j=$J -v step=$step -v len=$len -v "dbs=$dbs" \
            -v taxidmap=$taxidmap -v taxdump=$taxdump \
            'for db in {dbs}; do \
                seqkit grep -n -i -v -r -p plasmid {} -w 0 \
                    | seqkit sliding --greedy --step {step} --window {len} -w 0 \
                    | kmcp search -q -j {j} -d $db -o {}.kmcp@$db.tsv.gz; \
            done; \
            kmcp merge -q {}.kmcp@*.tsv.gz \
                | kmcp utils filter -q -T {taxidmap} -X {taxdump} -o {}.kmcp.uniq.tsv.gz; \
            /bin/rm {}.kmcp@*.tsv.gz; ' -c -C search.rush
    
    # merge regions
    time find $input -name "*.fna.gz" \
        | rush 'kmcp utils merge-regions -q -I -l 20 {}.kmcp.uniq.tsv.gz -o {}.kmcp.uniq.tsv.gz.bed' \
            -c -C merge.rush
    
    # -----------------------
    # length summary
    time find $input -name "*.fna.gz" \
        | rush -k 'glen=$(seqkit stats -T {} | csvtk cut -t -f sum_len | csvtk del-header); \
                len=$(awk "{print \$3-\$2}" {}.kmcp.uniq.tsv.gz.bed | csvtk summary -Ht -n 0  -f 1:sum); \
                echo -ne "{%:},$len,$glen\n" ' \
        | csvtk add-header -n file,uniqs,genome -o $input.len.csv.gz
    
    # --------------------------------------------------
    # extract subsequences

    output=uniqs/gtdb
    mkdir -p $output
    
    find $input/ -name "*.fna.gz" \
        | rush -v out=$output \
            'seqkit subseq --quiet --bed {}.kmcp.uniq.tsv.gz.bed {} -o {out}/{%}' \
            -c -C subseq.rush
