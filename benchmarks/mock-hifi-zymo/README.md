## Benchmarks on the ZymoBIOMICS Gut Microbiome Standard D6331


## Datasets

[The ZymoBIOMICS Gut Microbiome Standard D6331](https://www.zymoresearch.com/collections/zymobiomics-microbial-community-standards/products/zymobiomics-gut-microbiome-standard )
is comprised of 21 different strains to mimic the human gut microbiome.

The abundance table comes from the `Genome Copy` in [ds1712_zymobiomics_gut_microbiome_standard_data_sheet.pdf](https://files.zymoresearch.com/datasheets/ds1712_zymobiomics_gut_microbiome_standard_data_sheet.pdf)

Three PacBio HiFi Reads are from [PRJNA680590](https://www.ebi.ac.uk/ena/browser/view/PRJNA680590):

- SRR13128012: Zymo D6331 PacBio Ultra-Low Input Library
- SRR13128013: Zymo D6331 PacBio Low Input Library
- SRR13128014: Zymo D6331 PacBio Standard Input Library

Download and rename. The files look like:

    $ ls reads/
    ZymoBIOMICS-D6331-Low.fastq.gz
    ZymoBIOMICS-D6331-Standard.fastq.gz
    ZymoBIOMICS-D6331-Ultra-Low.fastq.gz

    # stats
    $ seqkit stats reads/*.fastq.gz
    file                                   format   type   num_seqs    sum_len          min_len   avg_len   max_len   Q1      Q2      Q3       sum_gap   N50      Q20(%)   Q30(%)
    ZymoBIOMICS-D6331-Low.fastq.gz         FASTQ    DNA    2,770,833   25,783,223,248   47        9,305.2   45,807    5,859   8,348   11,908   0         10,753   99.01    97.61
    ZymoBIOMICS-D6331-Standard.fastq.gz    FASTQ    DNA    1,978,852   17,993,566,711   45        9,092.9   39,601    5,592   8,077   11,716   0         10,628   99.03    97.66
    ZymoBIOMICS-D6331-Ultra-Low.fastq.gz   FASTQ    DNA    2,480,208   21,331,978,947   199       8,600.9   32,014    6,492   8,306   10,332   0         9,273    99       97.63


## KMCP (k-mer database)


We search against GTDB, Genbank-viral, and Refseq-fungi respectively, and merge the results.

    # ------------------------------------------------------------------------
    # prepare folder and files.
    
    mkdir -p kmcp-all-kmer
    cd kmcp-all-kmer
    fd fastq.gz$ ../reads | rush 'ln -s {}'
    cd ..
    
    
    # ------------------------------------------------------------------------
    # search
    
    reads=kmcp-all-kmer
    j=4
    J=40
    s=1000   # sliding step   
    w=212   # sliding window, 212=64*3+21-1, for better searching speed.
  

    # -------------------------------------
    # gtdb
    
    db=gtdb.kmcp/
    X=taxdump/
    T=gtdb.kmcp/taxid.map    
    dbname=gtdb
    
    fd fastq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -v db=$db -v dbname=$dbname -v s=$s -v w=$w -j $j -v j=$J -v 'p={:}' \
            'seqkit sliding -s {s} -W {w} -g {} \
                | kmcp search -d {db} \
                    -o {p}.kmcp@{dbname}.tsv.gz \
                    --log {p}.kmcp@{dbname}.tsv.gz.log -j {j}' \
            -c -C $reads@$dbname.rush
     
    
    # -------------------------------------
    # refseq-fungi
    
    db=refseq-fungi.kmcp/
    X=taxdump/
    T=refseq-fungi.kmcp/taxid.map    
    dbname=refseq-fungi
    
    fd fastq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -v db=$db -v dbname=$dbname -v s=$s -v w=$w -j $j -v j=$J -v 'p={:}' \
            'seqkit sliding -s {s} -W {w} -g {} \
                | kmcp search -d {db} \
                    -o {p}.kmcp@{dbname}.tsv.gz \
                    --log {p}.kmcp@{dbname}.tsv.gz.log -j {j}' \
            -c -C $reads@$dbname.rush
    
        
    # -------------------------------------
    # genbank-viral
    
    db=genbank-viral.kmcp/
    X=taxdump/
    T=genbank-viral.kmcp/taxid.map    
    dbname=genbank-viral    
    
    fd fastq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -v db=$db -v dbname=$dbname -v s=$s -v w=$w -j $j -v j=$J -v 'p={:}' \
            'seqkit sliding -s {s} -W {w} -g {} \
                | kmcp search -d {db} \
                    -o {p}.kmcp@{dbname}.tsv.gz \
                    --log {p}.kmcp@{dbname}.tsv.gz.log -j {j}' \
            -c -C $reads@$dbname.rush
    
    
    
    # ------------------------------------------------------------------------
    # merge results
    
    reads=kmcp-all-kmer
    j=16
    fd fastq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j $j -v 'p={:}' \
            'kmcp merge {p}.kmcp@*.tsv.gz -o {p}.kmcp.tsv.gz --log {p}.kmcp.tsv.gz.log'
    
    
    # ------------------------------------------------------------------------
    # [for merged search results] multiple profiling modes        
    
    X=taxdump/
    # cat genbank-viral.kmcp/taxid.map gtdb.kmcp/taxid.map refseq-fungi.kmcp/taxid.map > taxid.map
    # cat genbank-viral.kmcp/name.map gtdb.kmcp/name.map refseq-fungi.kmcp/name.map > name.map
    T=taxid.map
    
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

## KMCP (k-mer sketches database)

    # ------------------------------------------------------------------------
    # prepare folder and files.
    
    mkdir -p kmcp-syncmers
    cd kmcp-syncmers
    fd fastq.gz$ ../reads | rush 'ln -s {}'
    cd ..
    
    
    # ------------------------------------------------------------------------
    # search
    
    reads=kmcp-syncmers
    j=4
    J=40
    s=500   # sliding step   
    w=500   # sliding window
  

    # -------------------------------------
    # gtdb
    
    db=gtdb.sync16.kmcp/
    X=taxdump/
    T=gtdb.sync16.kmcp/taxid.map    
    dbname=gtdb
    
    fd fastq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -v db=$db -v dbname=$dbname -v s=$s -v w=$w -j $j -v j=$J -v 'p={:}' \
            'seqkit sliding -s {s} -W {w} -g {} \
                | kmcp search -d {db} \
                    -o {p}.kmcp@{dbname}.tsv.gz \
                    --log {p}.kmcp@{dbname}.tsv.gz.log -j {j}' \
            -c -C $reads@$dbname.rush
     
    
    # -------------------------------------
    # refseq-fungi
    
    db=refseq-fungi.sync16.kmcp/
    X=taxdump/
    T=refseq-fungi.sync16.kmcp/taxid.map    
    dbname=refseq-fungi
    
    fd fastq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -v db=$db -v dbname=$dbname -v s=$s -v w=$w -j $j -v j=$J -v 'p={:}' \
            'seqkit sliding -s {s} -W {w} -g {} \
                | kmcp search -d {db} \
                    -o {p}.kmcp@{dbname}.tsv.gz \
                    --log {p}.kmcp@{dbname}.tsv.gz.log -j {j}' \
            -c -C $reads@$dbname.rush
    
        
    # -------------------------------------
    # genbank-viral
    
    db=genbank-viral.sync16.kmcp/
    X=taxdump/
    T=genbank-viral.sync16.kmcp/taxid.map    
    dbname=genbank-viral    
    
    fd fastq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -v db=$db -v dbname=$dbname -v s=$s -v w=$w -j $j -v j=$J -v 'p={:}' \
            'seqkit sliding -s {s} -W {w} -g {} \
                | kmcp search -d {db} \
                    -o {p}.kmcp@{dbname}.tsv.gz \
                    --log {p}.kmcp@{dbname}.tsv.gz.log -j {j}' \
            -c -C $reads@$dbname.rush
    
    # ------------------------------------------------------------------------
    # merge results
    
    reads=kmcp-syncmers
    j=16
    fd fastq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j $j -v 'p={:}' \
            'kmcp merge {p}.kmcp@*.tsv.gz -o {p}.kmcp.tsv.gz --log {p}.kmcp.tsv.gz.log'
            

    # ------------------------------------------------------------------------
    # [for merged search results] multiple profiling modes        
    
    X=taxdump/
    # cat genbank-viral.kmcp/taxid.map gtdb.kmcp/taxid.map refseq-fungi.kmcp/taxid.map > taxid.map
    # cat genbank-viral.kmcp/name.map gtdb.kmcp/name.map refseq-fungi.kmcp/name.map > name.map
    T=taxid.map
    
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

