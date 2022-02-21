
## KMCP

    # ------------------------------------------------------------------------
    # prepare folder and files
    
    mkdir -p kmcp-se
    cd kmcp-se
    fd fq.gz$ ../reads | rush 'ln -s {}'
    cd ..

    # ------------------------------------------------------------------------
    # search
    
    reads=kmcp-se
    j=4
    J=40
    
    # -------------------------------------
    # prokaryotic    
    
    db=../refseq-cami2-k21-n10.db  
    dbname=refseq-cami2
    
    fd reads.fq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -v db=$db -v dbname=$dbname -j $j -v j=$J -v 'p={@(.+)_reads}' \
            'kmcp search -d {db} {} -o {p}.kmcp@{dbname}.tsv.gz \
                --log {p}.kmcp@{dbname}.tsv.gz.log -j {j}' \
            -c -C $reads@$dbname.rush
 
    # -------------------------------------
    # viral
    
    db=../refseq-cami2-viral-k21-n5.db  
    dbname=refseq-cami2-viral
    
    fd reads.fq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -v db=$db -v dbname=$dbname -j $j -v j=$J -v 'p={@(.+)_reads}' \
            'kmcp search -d {db} {} -o {p}.kmcp@{dbname}.tsv.gz \
                --log {p}.kmcp@{dbname}.tsv.gz.log -j {j}' \
            -c -C $reads@$dbname.rush
            
 
    # ------------------------------------------------------------------------
    # merge results
    
    reads=kmcp-se
    j=10
    fd reads.fq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j $j -v 'p={@(.+)_reads}' \
            'kmcp merge {p}.kmcp@*.tsv.gz -o {p}.kmcp.tsv.gz --log {p}.kmcp.tsv.gz.log'
    
    
    # ------------------------------------------------------------------------
    # [for merged search results] multiple profiling modes
    
    
    X=../taxdump/
    cat ../taxid.map ../taxid-viral.map > ../_taxid.map
    T=../_taxid.map
        
    for m in $(seq 1 5); do
        fd kmcp.tsv.gz$ $reads/ \
            | csvtk sort -H -k 1:N \
            | rush -v X=$X -v T=$T -v m=$m \
                'kmcp profile -m {m} -X {X} -T {T} {} -o {}.k-m{m}.profile -C {}.c-m{m}.profile -s {%:} --log {}.k-m{m}.profile.log' 
        
        profile=$reads.c-m$m.profile
        fd kmcp.tsv.gz.c-m$m.profile$ $reads/ \
            | csvtk sort -H -k 1:N \
            | rush -j 1 'cat {}' \
            > $profile
    done
