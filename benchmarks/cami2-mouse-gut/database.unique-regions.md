## Database

    k=21
    scale=5
    
    in=refseq-cami2-slim 
    out=refseq-cami2.smin$scale.kmcp
    tmp=refseq-cami2-k21-n10-d$scale
    
    kmcp compute -I $in -O $tmp -k $k -n 10 -l 100 -B plasmid -D $scale \
        --log $tmp.log --force        
    kmcp index -j 40 -I $tmp -O $out -n 1 -f 0.3 \
        --log $out.log
    
    cp taxid.map name.map $out/

## Steps

    # --------------------------------------------------
    # preprare files
    input=uniq
    
    mkdir -p $input
    cd $input
    
    find ../refseq-cami2-slim -name "*.fna.gz" | rush 'ln -s {}'
    cd ../
    

    # --------------------------------------------------
    # search and filter
    
    step=200
    len=500
    db=refseq-cami2.smin5.kmcp
    taxidmap=taxid.map
    taxdump=taxdump/
    
    j=4
    J=40
    find $input/ -name "*.fna.gz" \
        | rush -j $j -v j=$J -v step=$step -v len=$len -v "db=$db" \
            -v taxidmap=$taxidmap -v taxdump=$taxdump \
            'seqkit grep -n -i -v -r -p plasmid {} \
                | seqkit sliding --greedy --step {step} --window {len} -w 0 \
                | kmcp search -j {j} -d {db} \
                | kmcp filter -T {taxidmap} -X {taxdump} -o {}.kmcp.uniq.tsv.gz;' \
          -c -C $input.search.rush
    
    # merge regions
    find $input/ -name "*.fna.gz" \
        | rush 'kmcp utils merge-regions -q -I -l 20 {}.kmcp.uniq.tsv.gz -o {}.kmcp.uniq.tsv.gz.bed'

    
    # -----------------------
    # length summary
    time find $input -name "*.fna.gz" \
        | rush -k 'glen=$(seqkit stats -T {} | csvtk cut -t -f sum_len | csvtk del-header); \
                len=$(awk "{print \$3-\$2}" {}.kmcp.uniq.tsv.gz.bed | csvtk summary -Ht -n 0  -f 1:sum); \
                echo -ne "{%..},$len,$glen\n" ' \
        | csvtk add-header -n file,uniqs,genome -o $input.len.csv.gz

    # --------------------------------------------------
    # extract subsequences
    
    input=uniq
    output=refseq-cami2-slim.uniq
    mkdir -p $output
    
    find $input/ -name "*.fna.gz" \
        | rush -v out=$output \
            'seqkit subseq --quiet --bed {}.kmcp.uniq.tsv.gz.bed {} -o {out}/{%}'

## check some genome without unique regions found

    input=uniq2
    
    step=10
    len=150
    db=refseq-cami2-k31-n10.db
    taxidmap=taxid.map
    taxdump=taxdump/
    
    
    j=4
    J=40
    find $input/ -name "*.fna.gz" \
        | rush -j $j -v j=$J -v step=$step -v len=$len -v "db=$db" \
            -v taxidmap=$taxidmap -v taxdump=$taxdump \
            'seqkit grep -n -i -v -r -p plasmid {} \
                | seqkit sliding --greedy --step {step} --window {len} -w 0 \
                | kmcp search -j {j} -d {db} -o {}.kmcp.tsv.gz; \
             kmcp filter -T {taxidmap} -X {taxdump} {}.kmcp.tsv.gz -o {}.kmcp.uniq.tsv.gz;' \
          -c -C $input.search.rush
