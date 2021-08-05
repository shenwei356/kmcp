# Database

Prebuilt databases are available.

kingdom             |source     |#species|#assembly|file|file size
:-------------------|:----------|:-------|:--------|:---|:--------
Bacteria and Archaea|GTDB r202  |43252   |47894    |    |55.12 GB
Viruses             |Refseq r207|7300    |11618    |    |4.14 GB
Fungi               |Refseq r207|148     |390      |    |11.12 GB

Taxonomy data

- [taxdump.tar.gz]()


These databases are created with steps below.

## GTDB

Tools

- [brename](https://github.com/shenwei356/brename/releases) for batching renaming files.
- [rush](https://github.com/shenwei356/rush/releases) for executing jobs in parallel.
- [dustmasker](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) for masking low-complexity regions
- [kmcp]()

Files

 - [gtdb_genomes_reps_r202.tar.gz](https://data.ace.uq.edu.au/public/gtdb/data/releases/release202/202.0/genomic_files_reps/gtdb_genomes_reps_r202.tar.gz)
 - [ar122_metadata_r202.tar.gz](https://data.ace.uq.edu.au/public/gtdb/data/releases/release202/202.0/ar122_metadata_r202.tar.gz)
 - [bac120_metadata_r202.tar.gz](https://data.ace.uq.edu.au/public/gtdb/data/releases/release202/202.0/bac120_metadata_r202.tar.gz)

Uncompressing and renaming
 
    # uncompress
    mkdir -p gtdb
    tar -zxvf gtdb_genomes_reps_r202.tar.gz -O gtdb
    
    # rename
    brename -R -p '^(\w{3}_\d{9}\.\d+).+' -r '$1.fna.gz' gtdb    
  
Mapping file

    tar -zxvf ar122_metadata_r202.tar.gz  bac120_metadata_r202.tar.gz
    
    (cat ar122_metadata_r202.tsv; sed 1d bac120_metadata_r202.tsv) \
        | csvtk cut -t -f accession,ncbi_taxid \
        | csvtk replace -t -p '^.._' \
        | csvtk del-header \
        > taxid.map
        
    # sequence accesion -> full head
    find gtdb/ -name *.fna.gz \
        | rush -k 'echo -ne "{%@(.+).fna}\t$(seqkit head -n 1 {} | seqkit seq -n)\n" ' \
        > name.map
    
    # number of species
    cat taxid.map \
        | csvtk grep -Ht -P <(cut -f 1 name.map) \
        | taxonkit filter -i 2 -E species \
        | wc -l
        
Building database

    # mask low-complexity region
    mkdir -p gtdb.masked
    find gtdb/ -name *.fna.gz \
        | rush 'dustmasker -in <(zcat {}) -outfmt fasta \
            | sed -e "/^>/!s/[a-z]/n/g" \
            | gzip -c > gtdb.masked/{%}'
    
    # 
    kmcp compute -I gtdb.masked/ -e -k 21 -n 10 -B plasmid -O gtdb-r202-k21-n10 --force

    
    kmcp index -j 32 -I gtdb-r202-k21-n10 -O gtdb-r202-k21-n10.db -n 1 -f 0.3 --dry-run


## RefSeq

Downloading viral ad fungi sequences

    # name=fungi
    name=viral
    
    # https://github.com/pirovc/genome_updater
    # -k for dry-run
    # -i for fix
    time genome_updater.sh \
        -d "refseq"\
        -g $name \
        -c "all" \
        -l "all" \
        -f "genomic.fna.gz" \
        -o "refseq-$name" \
        -t 12 \
        -m -a -p -k

    # cd to 2021-07-30_21-54-19
        
    # taxdump
    mkdir -p taxdump
    tar -zxvf taxdump.tar.gz -C taxdump
    
    # ass2taxid
    cut -f 1,6 assembly_summary.txt > taxid.map    
    # ass2name
    cut -f 1,8 assembly_summary.txt > name.map
        
    # optional
    # stats
    cat taxid.map  \
        | csvtk freq -Ht -f 2 -nr \
        | taxonkit lineage -r -n -L --data-dir taxdump/ \
        | taxonkit reformat -I 1 -f '{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}' --data-dir taxdump/ \
        | csvtk add-header -t -n 'taxid,count,name,rank,superkindom,phylum,class,order,family,genus,species' \
        > taxid.map.stats.tsv
    
    seqkit stats -T -j 12 --infile-list <(find files -name *.fna.gz) > files.stats.tsv
        
Building database

    # mask
    mkdir -p files.masked
    fd fna.gz files \
        | rush 'dustmasker -in <(zcat {}) -outfmt fasta \
            | sed -e "/^>/!s/[a-z]/n/g" \
            | gzip -c > files.masked/{%}'
            
    # rename
    brename -R -p '^(\w{3}_\d{9}\.\d+).+' -r '$1.fna.gz' files.masked   
    
    
    
    # -----------------------------------------------------------------
    # for viral
    name=viral
    
    kmcp compute -I files.masked/ -O refseq-$name-k21-n10 \
        -k 21 --exact-number --seq-name-filter plasmid \
        --split-number 10 --split-overlap 100 --force
      
    kmcp index -I refseq-$name-k21-n10/ -O refseq-$name-k21-n10.db \
        -j 32 -f 0.001 -n 3 --force
    
    # -----------------------------------------------------------------
    # for fungi
    name=fungi
    
    kmcp compute -I files.masked/ -O refseq-$name-k21-n10 \
        -k 21 --exact-number --seq-name-filter plasmid \
        --split-number 10 --split-overlap 100 --force
      
    kmcp index -I refseq-$name-k21-n10/ -O refseq-$name-k21-n10.db \
        -j 32 -f 0.05 -n 2 --force

## HumGut

http://arken.nmbu.no/~larssn/humgut/


## Custom database

