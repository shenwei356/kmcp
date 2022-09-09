# Benchmarks on datasets from CAMI2 challenge

We created KMCP databaes with the RefSeq and Taxonomy snapshot provided by CAMI2 (2019-01-08).

Gold standard and profiles of other tools were  downloaded from: https://zenodo.org/record/5006866.
Sequence datasets were download from: https://data.cami-challenge.org/participate.

## Softwares

- kmcp [v0.9.0](https://github.com/shenwei356/kmcp/releases/tag/v0.9.0)

## Databases

[Prebuilt databases and the reference genomes](https://1drv.ms/u/s!Ag89cZ8NYcqtjVVADr8r--fnKFt-?e=ivNZNK):

- DB for bacteria: [refseq-cami2-k21-n10.db.tar.gz](https://1drv.ms/u/s!Ag89cZ8NYcqtjV62KmQmOojxwBRr?e=lp5a9F), [md5](https://1drv.ms/t/s!Ag89cZ8NYcqtjWISqJGcxQD39FCv?e=CQ0E8d)
- DB for viruses: [refseq-cami2-viral-k21-n5.db.tar.gz](https://1drv.ms/u/s!Ag89cZ8NYcqtjVyYFIHY01PtDMcx?e=AO7xkY), [md5](https://1drv.ms/t/s!Ag89cZ8NYcqtjWDTIXL4eMpZNVA0?e=1YXKkk)
- TaxId mapping: [taxid.map](https://1drv.ms/u/s!Ag89cZ8NYcqtjVvZBPDumqTp0LLX?e=2OGhTe), [md5](https://1drv.ms/t/s!Ag89cZ8NYcqtjWXOc2bP9cmE2H9C?e=yyZnaB);
  [taxid-viral.map](https://1drv.ms/u/s!Ag89cZ8NYcqtjVclRm9-rd-K2MA3?e=6aOgqm), [md5](https://1drv.ms/t/s!Ag89cZ8NYcqtjWZR9Zfs7m33k_lV?e=ALSUe0)
- [taxdump.tar.gz](https://1drv.ms/u/s!Ag89cZ8NYcqtjVjXKnxzq-sUb8Cw?e=7AwTnG), [md5](https://1drv.ms/t/s!Ag89cZ8NYcqtjWQf06gXeM0rLHJ9?e=dxSW9g)

**Attention**: the CAMI2 RefSeq snapshot did not include viruses,
and CAMI2 toy mouse gut dataset did not contain viral reads either.

Reference genomes (Bacteria and Archaea):

1. Microbial genomes were extracted from CAMI2 RefSeq snapshot (`2019-01-08`) using
corresponding taxonomy information .
2. **For every species, at most 5 assemblies (sorted by assembly accession) were kept**.

Reference genomes (Viruses):

1. Current RefSeq virual `assembly_summary.txt` was downloaded, and filter
records before `2019-01-08`.
2. All viral reference genomes were downloaded for database building.

[Databases building steps](./database.md).

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
