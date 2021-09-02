# Database

KMCP is a reference based taxonomic profiling tool.

## Prebuilt Databases

### A). Databases for metagenomic profiling

Prebuilt databases are available, you can also [build custom databases](#custom-database).

|kingdoms                |source     |# species|# assembly|file                          |file size|
|:-----------------------|:----------|:--------|:---------|:-----------------------------|:--------|
|**Bacteria and Archaea**|GTDB r202  |43252    |47894     |[gtdb.kmcp.tar.gz]()          |55.12 GB |
|**Viruses**             |Refseq r207|7300     |11618     |[refseq-viruses.kmcp.tar.gz]()|4.14 GB  |
|**Fungi**               |Refseq r207|148      |390       |[refseq-fungi.kmcp.tar.gz]()  |11.12 GB |

Taxonomy data

- [taxdump.tar.gz]()

These databases are created with steps below.

### B). Databases for genome similarity estimation

Check [tutorial](https://bioinf.shenwei.me/kmcp/tutorial/searching).

|kingdoms                |source     |sketch         |parameters          |file                                  |file size|
|:-----------------------|:----------|:--------------|:-------------------|:-------------------------------------|:--------|
|**Bacteria and Archaea**|GTDB r202  |Scaled MinHash |k=31, scale=1000    |[gtdb.minhash.kmcp.tar.gz]()          |850 MB   |
|**Bacteria and Archaea**|GTDB r202  |Closed Syncmers|k=31, s=15, scale=60|[gtdb.syncmer.kmcp.tar.gz]()          |1.04 GB  |
|**Viruses**             |Refseq r207|Scaled MinHash |K=31, scale=10      |[refseq-viruses.minhash.kmcp.tar.gz]()|         |
|**Viruses**             |Refseq r207|Closed Syncmers|k=31, s=21          |[refseq-viruses.syncmer.kmcp.tar.gz]()|         |

## Building Databases

### GTDB

Tools

- [brename](https://github.com/shenwei356/brename/releases) for batching renaming files.
- [rush](https://github.com/shenwei356/rush/releases) for executing jobs in parallel.
- [dustmasker](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) for masking low-complexity regions.
- [seqkit](https://github.com/shenwei356/seqkit/releases) for FASTA file processing.
- [kmcp](/download) for metagenomic profiling.

Files

 - [gtdb_genomes_reps_r202.tar.gz](https://data.ace.uq.edu.au/public/gtdb/data/releases/release202/202.0/genomic_files_reps/gtdb_genomes_reps_r202.tar.gz)
 - [ar122_metadata_r202.tar.gz](https://data.ace.uq.edu.au/public/gtdb/data/releases/release202/202.0/ar122_metadata_r202.tar.gz)
 - [bac120_metadata_r202.tar.gz](https://data.ace.uq.edu.au/public/gtdb/data/releases/release202/202.0/bac120_metadata_r202.tar.gz)

Uncompressing and renaming:
 
    # uncompress
    mkdir -p gtdb
    tar -zxvf gtdb_genomes_reps_r202.tar.gz -O gtdb
    
    # rename
    brename -R -p '^(\w{3}_\d{9}\.\d+).+' -r '$1.fna.gz' gtdb    
  
Mapping file:

    tar -zxvf ar122_metadata_r202.tar.gz  bac120_metadata_r202.tar.gz
    
    # assembly accesion -> full head
    find gtdb/ -name "*.fna.gz" \
        | rush -k 'echo -ne "{%@(.+).fna}\t$(seqkit head -n 1 {} | seqkit seq -n)\n" ' \
        > name.map
        
    # assembly accesion -> taxid
    (cat ar122_metadata_r202.tsv; sed 1d bac120_metadata_r202.tsv) \
        | csvtk cut -t -f accession,ncbi_taxid \
        | csvtk replace -t -p '^.._' \
        | csvtk grep -Ht -P <(cut -f 1 name.map) \
        | csvtk del-header \
        > taxid.map

    
    # stats (optional)

    # number of species/strains
    cat taxid.map \
        | taxonkit lineage -i 2 -r -n -L \
        | csvtk freq -Ht -f 4 -nr \
        | csvtk pretty -H -t \
        | tee taxid.map.stats
    species           43566
    strain            4107
    subspecies        111
    forma specialis   58
    no rank           26
    isolate           24
    serotype          1
    
    # number of unique species/strains
    cat taxid.map \
        | csvtk uniq -Ht -f 2 \
        | taxonkit lineage -i 2 -r -n -L \
        | csvtk freq -Ht -f 4 -nr \
        | csvtk pretty -H -t \
        | tee taxid.map.uniq.stats
    species           24739
    strain            4100
    subspecies        89
    forma specialis   58
    no rank           26
    isolate           24
    serotype          1
  
        
Building database:

    # mask low-complexity region
    mkdir -p gtdb.masked
    find gtdb/ -name "*.fna.gz" \
        | rush 'dustmasker -in <(zcat {}) -outfmt fasta \
            | sed -e "/^>/!s/[a-z]/n/g" \
            | gzip -c > gtdb.masked/{%}'
    
    # compute k-mers
    #   sequence containing "plasmid" in name are ignored,
    #   reference genomes are splitted into 10 fragments
    #   k = 21
    kmcp compute -I gtdb.masked/ -k 21 -n 10 -B plasmid -O gtdb-r202-k21-n10 --force

    # build database
    #   number of index files: 32, for server with >= 32 CPU cores
    #   bloom filter parameter:
    #     number of hash function: 1
    #     false positive rate: 0.3
    kmcp index -j 32 -I gtdb-r202-k21-n10 -O gtdb.kmcp -n 1 -f 0.3
    
    # cp taxid and name mapping file to database directory
    cp taxid.map name.map gtdb.kmcp/


### RefSeq

Tools

- [genome_updater](https://github.com/pirovc/genome_updater) for downloading genomes from NCBI.

Downloading viral and fungi sequences:

    # name=fungi
    name=viral
    
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
        -m -a -p

    # cd to 2021-07-30_21-54-19
        
    # taxdump
    mkdir -p taxdump
    tar -zxvf taxdump.tar.gz -C taxdump
    
    # assembly accesion -> taxid
    cut -f 1,6 assembly_summary.txt > taxid.map    
    # assembly accesion -> name
    cut -f 1,8 assembly_summary.txt > name.map
    
    
    # stats (optional)
    cat taxid.map  \
        | csvtk freq -Ht -f 2 -nr \
        | taxonkit lineage -r -n -L --data-dir taxdump/ \
        | taxonkit reformat -I 1 -f '{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}' --data-dir taxdump/ \
        | csvtk add-header -t -n 'taxid,count,name,rank,superkindom,phylum,class,order,family,genus,species' \
        > taxid.map.stats.tsv
    
    seqkit stats -T -j 12 --infile-list <(find files -name "*.fna.gz") > files.stats.tsv
        
Building database:

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
        -k 21 --seq-name-filter plasmid \
        --split-number 10 --split-overlap 100 --force
    
    # viral genomes are small:
    #   using small false positive rate: 0.001
    #   using more hash functions: 3
    kmcp index -I refseq-$name-k21-n10/ -O refseq-$name-k21-n10.db \
        -j 32 -f 0.001 -n 3 --force
    
    mv refseq-$name-k21-n10.db refseq-viruses.kmcp
    
    # cp taxid and name mapping file to database directory
    cp taxid.map name.map gtdb.kmcp/

    # -----------------------------------------------------------------
    # for fungi
    name=fungi
    
    kmcp compute -I files.masked/ -O refseq-$name-k21-n10 \
        -k 21 --seq-name-filter plasmid \
        --split-number 10 --split-overlap 100 --force
      
    kmcp index -I refseq-$name-k21-n10/ -O refseq-$name-k21-n10.db \
        -j 32 -f 0.05 -n 2 --force

    mv refseq-$name-k21-n10.db refseq-fungi.kmcp
    
    # cp taxid and name mapping file to database directory
    cp taxid.map name.map gtdb.kmcp/

### HumGut


[HumGut](http://arken.nmbu.no/~larssn/humgut/) is a comprehensive Human Gut prokaryotic genomes collection filtered by metagenome data.

Dataset

- [HumGut.tar.gz](http://arken.nmbu.no/~larssn/humgut/HumGut.tar.gz)
- [HumGut.tsv](http://arken.nmbu.no/~larssn/humgut/HumGut.tsv)
- Taxdump files:
    - NCBI taxonomy: [ncbi_names.dump](http://arken.nmbu.no/~larssn/humgut/ncbi_names.dmp)
      and [ncbi_nodes.dump](http://arken.nmbu.no/~larssn/humgut/ncbi_names.dmp)
    - GTDB taxomomy: [gtdb_names.dump](http://arken.nmbu.no/~larssn/humgut/gtdb_names.dmp)
      and [gtdb_nodes.dump](http://arken.nmbu.no/~larssn/humgut/gtdb_names.dmp)
    
Uncompressing and renaming:
 
    # uncompress
    tar -zxvf HumGut.tar.gz
    
    # taxdump files
    mkdir taxdump
    mv ncbi_names.dmp taxdump/names.dmp
    mv ncbi_nodes.dmp taxdump/nodes.dmp
  
Mapping file:
    
    # assembly accesion -> taxid
    cat HumGut.tsv \
        | csvtk cut -t -f genome_file,ncbi_tax_id \
        | csvtk replace -f genome_file -t -p '\..+$' \
        | csvtk del-header \
        > taxid.map
        
    # assembly accesion -> name
    cat HumGut.tsv \
        | csvtk cut -t -f genome_file,ncbi_organism_name \
        | csvtk replace -f genome_file -t -p '\..+$' \
        | csvtk del-header \
        > name.map
    
Stats:

    # NCBI taxonomy

    # number of species/strains
    cat taxid.map \
        | taxonkit lineage --data-dir taxdump-ncbi -i 2 -r -n -L \
        | csvtk freq -Ht -f 4 -nr \
        | csvtk pretty -H -t \
        | tee taxid.map.stats
    species      23604
    strain       6816
    subspecies   161
    no rank      69
    serotype     38
    genus        3
    
    # number of unique species/strains
    cat taxid.map \
        | csvtk uniq -Ht -f 2 \
        | taxonkit lineage --data-dir taxdump-ncbi -i 2 -r -n -L \
        | csvtk freq -Ht -f 4 -nr \
        | csvtk pretty -H -t \
        | tee taxid.map.uniq.stats;
    species      1240
    strain       458
    subspecies   15
    genus        1
    no rank      1
    serotype     1
    
    # GTDB taxonomy
    
    cat HumGut.tsv \
        | csvtk cut -t -f genome_file,gtdbtk_tax_id \
        | csvtk replace -f genome_file -t -p '\..+$' \
        | csvtk del-header \
        > taxid-gtdb.map
        
    cat taxid-gtdb.map \
        | csvtk uniq -Ht -f 2 \
        | taxonkit lineage --data-dir taxdump-gtdb -i 2 -r -n -L \
        | csvtk freq -Ht -f 4 -nr \
        | csvtk pretty -H -t \
        | tee taxid-gtdb.map.uniq.stats
    species   2810
    genus     434
    family    52
    order     12
    class     2
        
Building database:

    # mask low-complexity region
    mkdir -p humgut.masked
    find fna/ -name "*.fna.gz" \
        | rush 'dustmasker -in <(zcat {}) -outfmt fasta \
            | sed -e "/^>/!s/[a-z]/n/g" \
            | gzip -c > humgut.masked/{%}'
    
    # compute k-mers
    #   sequence containing "plasmid" in name are ignored,
    #   reference genomes are splitted into 10 fragments
    #   k = 21
    kmcp compute -I humgut.masked/ -k 21 -n 10 -B plasmid -O humgut-k21-n10 --force

    # build database
    #   number of index files: 32, for server with >= 32 CPU cores
    #   bloom filter parameter:
    #     number of hash function: 1
    #     false positive rate: 0.3
    kmcp index -j 32 -I humgut-k21-n10 -O humgut.kmcp -n 1 -f 0.3


### Custom database

Files:

1. Genome files
    - (Gzip-compressed) FASTA/Q format
    - One genome per file with the reference identifier in the file name.
2. TaxId mapping file (for metagenomic profiling)
    - Two-column (reference identifier and TaxId) tab-delimited.
3. [NCBI taxonomy dump files](ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz) (for metagenomic profiling)
    - `names.dmp`
    - `nodes.dmp`
    - `merged.dmp` (optional)
    - `delnodes.dmp` (optional)

Steps

    # directory containing genome files
    genomes=genomes

    # mask low-complexity region
    mkdir -p masked
    find $genomes -name "*" \
        | rush 'dustmasker -in <(zcat {}) -outfmt fasta \
            | sed -e "/^>/!s/[a-z]/n/g" \
            | gzip -c > masked/{%}'
    
    # compute k-mers
    #   sequence containing "plasmid" in name are ignored,
    #   reference genomes are splitted into 10 fragments
    #   k = 21
    kmcp compute --in-dir masked/ \
        --kmer 21 \
        --split-number 10 \
        --seq-name-filter plasmid \
        --out-dir custom-k21-n10 \
        --force

    # build database
    #   number of index files: 32, for server with >= 32 CPU cores
    #   bloom filter parameter:
    #     number of hash function: 1
    #     false positive rate: 0.3
    kmcp index -I custom-k21-n10 \
        --threads 32 \
        --num-hash 1 \
        --false-positive-rate 0.3 \
        --out-dir custom.kmcp 
