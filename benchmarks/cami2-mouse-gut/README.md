# Benchmarks on CAMI2 toy mouse gut short reads dataset


## Softwares

- kmcp ([v0.7.0](https://github.com/shenwei356/kmcp/releases/tag/v0.7.0)

## Databases

Prebuilt databases

- DB for bacteria: [refseq-cami2-k21-n10.db.tar.gz]()
- DB for viruses: [refseq-cami2-viral-k21-n5.db.tar.gz]()
- TaxId mapping: [taxid.map](), and [taxid-viral.map]()
- [taxdump.tar.gz]()

**Attention**: the CAMI2 RefSeq snapshot did not include viruses,
and CAMI2 toy mouse gut dataset did not contain viral reads either.

Reference genomes (Bacteria and Archaea):

1. Microbial genomes were extracted from CAMI2 RefSeq snapshot (`2019-01-08`) using
corresponding taxonomy information .
2. For every species, at most 5 assemblies (sorted by assembly accession) were kept.

Reference genomes (Viruses):

1. Current RefSeq virual `assembly_summary.txt` was downloaded, and filter
records before `2019-01-08`.
2. All viral reference genomes were downloaded for database building.

[Databases building steps](./database.md).

## Datasets

https://data.cami-challenge.org/participate 

2nd CAMI Toy Mouse Gut Dataset

    # short reads
    java -jar camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMISIM_MOUSEGUT . \
        -p "mousegut_scaffolds.+sample_\d+\/.+.fq.gz"

    # rename and re-organizing files
    # ...

    $ tree 19122017_mousegut_scaffolds
    19122017_mousegut_scaffolds
    ├── sample_0.fq.gz
    ├── sample_10.fq.gz
        
## Searching and Profiling
    
    db=refseq-cami2-k21-n10.db
    dbname=refseq-cami2-k21-n10
    
    
    # ------------------------------------------------------------------------
    # single-end mode

    # ln -s 19122017_mousegut_scaffolds single

    reads=single
    j=4
    J=40
    fd fq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -v db=$db -v dbname=$dbname -j $j -v j=$J \
            'kmcp search -d {db} {} -o {}.kmcp@{dbname}.tsv.gz --log {}.kmcp@{dbname}.tsv.gz.log -j {j}' \
            -c -C $reads@$dbname.rush
            
    # ------------------------------------------------------------------------           
    # paired-end mode

    # # split merged paired-end reads.
    # mkdir paired
    # cd paired
    # fd .fq.gz$ ../$reads | rush -j 12 'seqkit split2 -p 2 {} -O .'
    # brename -p .part_00 -r _
    # cd ..

    reads=paired
    j=4
    J=40
    fd _1.fq.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -v db=$db -v dbname=$dbname -j $j -v j=$J -v 'p={@(.+)_1.fq.gz}' \
            'kmcp search -d {db} --try-se -1 {} -2 {p}_2.fq.gz -o {p}.kmcp@{dbname}.tsv.gz \
            --log {p}.kmcp@{dbname}.tsv.gz.log -j {j}' \
            -c -C $reads@$dbname.rush

    # ------------------------------------------------------------------------
            
    X=taxdump
    T=taxid.map
    fd kmcp@$dbname.tsv.gz$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -v db=$db -v dbname=$dbname -v X=$X -v T=$T \
            'kmcp profile -X {X} -T {T} {} -o {}.k.profile -C {}.c.profile -s {@sample_(\d+)} --log {}.k.profile.log' 
    
    profile=$reads@$dbname.c.profile
    fd kmcp@$dbname.tsv.gz.c.profile$ $reads/ \
        | csvtk sort -H -k 1:N \
        | rush -j 1 'cat {}' \
        > $profile
        
    # ------------------------------------------------------------------------
    # profiling modes
    
    X=taxdump
    T=taxid.map
    
    for m in $(seq 1 5); do
        fd kmcp@$dbname.tsv.gz$ $reads/ \
            | csvtk sort -H -k 1:N \
            | rush -v db=$db -v dbname=$dbname -v X=$X -v T=$T -v m=$m \
                'kmcp profile -m {m} -X {X} -T {T} {} -o {}.k-m{m}.profile -C {}.c-m{m}.profile -s {@sample_(\d+)} --log {}.k-m{m}.profile.log' 
    
        profile=$reads@$dbname.c-m$m.profile
        fd kmcp@$dbname.tsv.gz.c-m$m.profile$ $reads/ \
            | csvtk sort -H -k 1:N \
            | rush -j 1 'cat {}' \
            > $profile
    done
    
## Test with prebuilt database

Searching

    # search on multiple database
    file=sample_0.fq.gz
    sampleId=0
    for db in refseq-cami2-k21-n10.db refseq-cami2-viral-k21-n5.db; do
        kmcp search -j 40 -d $db $file -o $file.kmcp@$db.tsv.gz
    done
    
    # merge search result
    kmcp merge $file.kmcp@*.tsv.gz -o $file.kmcp.tsv.gz
    
Profiling

    # profile
    kmcp profile -T taxid.map -T taxid-virus.map -X taxdump/ \
        $file.kmcp.tsv.gz \
        -o $file.kmcp.tsv.gz.k.profile \
        -C $file.kmcp.tsv.gz.cami.profile -s $sampleId        
    
