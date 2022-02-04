# Database

## Prebuilt databases

All prebuilt database are available at:

- [OneDrive](https://1drv.ms/u/s!Ag89cZ8NYcqtjHwpe0ND3SUEhyrp?e=QDRbEC) for global users.
- [CowTransfer](https://shenwei356.cowtransfer.com/s/c7220dd5901c42) for Chinese users, [a command-line tool](https://github.com/Mikubill/cowtransfer-uploader) is recommended.

**Hardware requirements**

- Prebuilt databases were built for computers with >= 32 CPU cores
in consideration of better parallelization,
and computers should have at least 64 GB. 
- By default, `kmcp search` loads the whole database into main memory (via [mmap](https://en.wikipedia.org/wiki/Mmap)) for fast searching.
Optionally, the flag `--low-mem` can be set to avoid loading the whole database,
while it's much slower, >10X slower on SSD and should be much slower on HDD disks.
- ***To reduce memory requirements on computers without enough memory,
users can divide the reference genomes into several parts
and build smaller databases for all parts, so that the biggest
database can be loaded into RAM***. 
This can also accelerate searching on a **computer cluster**, where every node searches a part of the database.
After performing database searching,
search results on all small databases can be merged with `kmcp merge`
for downstream analysis.



### A). Databases for metagenomic profiling

These databases are created following [steps below](#building-databases).
Users can also [build custom databases](#building-custom-databases), it's simple and fast.


|DB                      |source     |#species|#assemblies|parameters      |archive file                                                                                                                                                      |size     |
|:-----------------------|:----------|:-------|:----------|:---------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------|:--------|
|**Bacteria and Archaea**|GTDB r202  |43252   |47894      |k=21, frags=10  |[gtdb.kmcp.tar.gz](https://1drv.ms/u/s!Ag89cZ8NYcqtjRewpV1B37CO1Ghe?e=g0cwiI) (50.16 GB, [md5](https://1drv.ms/t/s!Ag89cZ8NYcqtjQmv5Nn4bt3hUSpN?e=H0PxRa))        |58.02 GB |
|**Bacteria and Archaea**|HumGut     |23604   |30691      |k=21, frags=10  |[humgut.kmcp.tar.gz](https://1drv.ms/u/s!Ag89cZ8NYcqtjUxZymOTLu1qJyDI?e=ZPWhDt) (18.70 GB, [md5](https://1drv.ms/t/s!Ag89cZ8NYcqtjUVZu1Y-Vtussvdc?e=wHlWdm))      |21.52 GB |
|**Fungi**               |Refseq r208|161     |403        |k=21, frags=10  |[refseq-fungi.kmcp.tar.gz](https://1drv.ms/t/s!Ag89cZ8NYcqtjROm3VsX6PVrxHe5?e=PO1N78) (3.67 GB, [md5](https://1drv.ms/t/s!Ag89cZ8NYcqtjQM4tAbFU2bFS07e?e=0CwT1E)) |4.18 GB  |
|**Viruses**             |GenBank 246|19584   |27936      |k=21, frags=5   |[genbank-viral.kmcp.tar.gz](https://1drv.ms/u/s!Ag89cZ8NYcqtjkI5lQI-ygIiDe-B?e=Viaev8) (1.15 GB, [md5](https://1drv.ms/t/s!Ag89cZ8NYcqtjkEqcjvzQczfIAMr?e=LBsj4X))|3.75 GB  |
|**Viruses**             |Refseq r208|7189    |11618      |k=21, frags=5   |[refseq-viral.kmcp.tar.gz](https://1drv.ms/u/s!Ag89cZ8NYcqtjQj5zEDzlN9kCYzT?e=bZNzAk) (967 MB, [md5](https://1drv.ms/t/s!Ag89cZ8NYcqtjQBrtR3Ol5GsJ6e3?e=wAgY1e))  |2.00 GB  |
|**Human**               |CHM13      |1       |1          |k=21, frags=1024|[human-chm13.kmcp.tar.gz](https://1drv.ms/u/s!Ag89cZ8NYcqtjVQgKPCZ7jciZqEp?e=jAO76U) (818 MB, [md5](https://1drv.ms/t/s!Ag89cZ8NYcqtjU1nGeOJaFf70y_K?e=bzJPcE))   |946 MB   |


**Taxonomy data**:

- Taxonomy dump file: [taxdump.tar.gz](https://1drv.ms/u/s!Ag89cZ8NYcqtjhCcJiTgJ7-dZPg3?e=vUeMmJ) (2021-12-06, [md5](https://1drv.ms/t/s!Ag89cZ8NYcqtjg4lwBMJt6ryNrvS?e=lAKZjU))

**Taxonomy data for HumGut**:

- Taxonomy dump file: [taxdump-humgut.tar.gz](https://1drv.ms/u/s!Ag89cZ8NYcqtjUeiSOIBztt87yfi?e=LqKhiC) ([md5](https://1drv.ms/t/s!Ag89cZ8NYcqtjUhiSa_tbjKHtRUl?e=k07pzf))
- Taxid mapping file: [taxid-humgut.map](https://1drv.ms/u/s!Ag89cZ8NYcqtjUusPDpqb2qfNKtj?e=PY9dxA) ([md5](https://1drv.ms/t/s!Ag89cZ8NYcqtjUqDQMxCC_PwJC6K?e=VowaM2))
- Name mapping file: [name-humgut.map](https://1drv.ms/u/s!Ag89cZ8NYcqtjUnwt8woH-HjN8TI?e=hM4iec) ([md5](https://1drv.ms/t/s!Ag89cZ8NYcqtjUZegCvp42p52yFv?e=6iqman))


### B). Databases for genome similarity estimation

Check the [tutorial](/tutorial/searching).

FracMinHash (Scaled MinHash):

|kingdoms                |source     |parameters      |file                                                                                                                                                                     |size     |
|:-----------------------|:----------|:---------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:--------|
|**Bacteria and Archaea**|GTDB r202  |k=31, scale=1000|[gtdb.minhash.kmcp.tar.gz](https://1drv.ms/u/s!Ag89cZ8NYcqtjQpRVpcW1zZJHWTR?e=A9Oh1q) (710 MB, [md5](https://1drv.ms/t/s!Ag89cZ8NYcqtjQuoumn3pK_jFMDq?e=6wUa4r))         |1.52 GB  |
|**Fungi**               |Refseq r208|k=31, scale=1000|[refseq-fungi.minhash.kmcp.tar.gz](https://1drv.ms/u/s!Ag89cZ8NYcqtjQTxGnMeL7n66_iY?e=kDwJi8) (49 MB, [md5](https://1drv.ms/t/s!Ag89cZ8NYcqtjQV_O4GTRoOUhc7W?e=w0p3l1))  |98 MB    |
|**Viruses**             |Genbank 246|k=31, scale=10  |[genbank-viral.minhash.kmcp.tar.gz](https://1drv.ms/u/s!Ag89cZ8NYcqtjhNMoB2v8OcXPwlB?e=dOCBL3) (580 MB, [md5](https://1drv.ms/t/s!Ag89cZ8NYcqtjhH5MRvKtbb-OW4v?e=dY6Yh9))|1.19 GB  |
|**Viruses**             |Refseq r208|k=31, scale=10  |[refseq-viral.minhash.kmcp.tar.gz](https://1drv.ms/u/s!Ag89cZ8NYcqtjH-2z4mrt4H3pcP8?e=irv43B) (205 MB, [md5](https://1drv.ms/t/s!Ag89cZ8NYcqtjQGBiyrNm8LOJcVZ?e=PzlBLl)) |555 MB   |

Closed Syncmers:

|kingdoms                |source     |parameters          |file                                                                                                                                                                     |size     |
|:-----------------------|:----------|:-------------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:--------|
|**Bacteria and Archaea**|GTDB r202  |k=31, s=15, scale=60|[gtdb.syncmer.kmcp.tar.gz](https://1drv.ms/t/s!Ag89cZ8NYcqtjQuoumn3pK_jFMDq?e=6wUa4r) (1.03 GB, [md5](https://1drv.ms/t/s!Ag89cZ8NYcqtjQ06tZBHsxwPvTN5?e=XDKAJk))        |2.28 GB  |
|**Fungi**               |Refseq r208|k=31, s=15, scale=60|[refseq-fungi.syncmer.kmcp.tar.gz](https://1drv.ms/u/s!Ag89cZ8NYcqtjQZqTZrnLyiEJ0IX?e=aORH7a) (73 MB, [md5](https://1drv.ms/t/s!Ag89cZ8NYcqtjQdqsq1rGf0YuCdz?e=g0gOtx))  |145 MB   |
|**Viruses**             |Genbank 246|k=31, s=10          |[genbank-viral.syncmer.kmcp.tar.gz](https://1drv.ms/u/s!Ag89cZ8NYcqtjhQiHSTYQNFUgqPx?e=xlJm35) (473 MB, [md5](https://1drv.ms/t/s!Ag89cZ8NYcqtjhL1YSXJGeBmI_c7?e=lWVkON))|1.06 GB  |
|**Viruses**             |Refseq r208|k=31, s=21          |[refseq-viral.syncmer.kmcp.tar.gz](https://1drv.ms/t/s!Ag89cZ8NYcqtjQGBiyrNm8LOJcVZ?e=PzlBLl) (162 MB, [md5](https://1drv.ms/t/s!Ag89cZ8NYcqtjQKyoa8lcPTO5yt_?e=e7Ki65)) |441 MB   |

### C). Databases of plasmid

|source     |# assembly|type          |parameters    |file                                                                                                                                                                       |size   |
|:----------|:---------|:-------------|:-------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:------|
|Refseq r208|37318     |All k-mers    |k=21          |[refseq-plasmid.kmcp.tar.gz](https://1drv.ms/u/s!Ag89cZ8NYcqtjUQpIuKC5ju2TIdD?e=yKdRs2) (5.29 GB, [md5](https://1drv.ms/t/s!Ag89cZ8NYcqtjT-TO94gsfpUpe8_?e=1RAvr1))        |7.80 GB|
|Refseq r208|37318     |FracMinHash   |K=31, scale=10|[refseq-plasmid.minhash.kmcp.tar.gz](https://1drv.ms/u/s!Ag89cZ8NYcqtjUF10X19vzRpOD86?e=Gfnau8) (1.01 GB, [md5](https://1drv.ms/t/s!Ag89cZ8NYcqtjUB8HzWSO-hu8iI3?e=KbVXKp))|2.00 GB|
|Refseq r208|37318     |Closed Syncmer|K=31, s=21    |[refseq-plasmid.syncmer.kmcp.tar.gz](https://1drv.ms/u/s!Ag89cZ8NYcqtjUPv9o2SA4kiSxwU?e=UvkSDV) (806 MB, [md5](https://1drv.ms/t/s!Ag89cZ8NYcqtjUKt3A_p7Q7BE7qT?e=dApTX8)) |1.54 GB|


## Building databases

### GTDB

Tools:

- [brename](https://github.com/shenwei356/brename/releases) for batching renaming files.
- [rush](https://github.com/shenwei356/rush/releases) for executing jobs in parallel.
- [seqkit](https://github.com/shenwei356/seqkit/releases) for FASTA file processing.
- [kmcp](/download) for metagenomic profiling.

Files:

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
        | csvtk grep -t -P <(cut -f 1 name.map) \
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
    strain            4108
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
    strain            4101
    subspecies        89
    forma specialis   58
    no rank           26
    isolate           24
    serotype          1
  
        
Building database (all k-mers, for profiling on short-reads):
    
    input=gtdb
    
    # compute k-mers
    #   sequence containing "plasmid" in name are ignored,
    #   reference genomes are splitted into 10 fragments with 100bp overlap
    #   k = 21
    kmcp compute -I $input -O gtdb-r202-k21-n10 -k 21 -n 10 -l 100 -B plasmid \
        --log gtdb-r202-k21-n10.log -j 32 --force

    # build database
    #   number of index files: 32, for server with >= 32 CPU cores
    #   bloom filter parameter:
    #     number of hash function: 1
    #     false positive rate: 0.3
    kmcp index -j 32 -I gtdb-r202-k21-n10 -O gtdb.kmcp -n 1 -f 0.3 \
        --log gtdb.kmcp.log
    
    # cp taxid and name mapping file to database directory
    cp taxid.map name.map gtdb.kmcp/
    
Building database (k-mer sketches, for profiling on long-reads):
    
    # -------------------------------------------------------------------------------------
    # Closed Syncmers with s=16
    
    input=gtdb
    
    # compute k-mers
    #   sequence containing "plasmid" in name are ignored,
    #   reference genomes are splitted into 10 fragments with 100bp overlap
    #   k = 21
    #   s = 16 # Closed Syncmers
    kmcp compute -I $input -O gtdb-r202-k21-n10-S16 -k 21 -S 16 -n 10 -l 100 -B plasmid \
        --log gtdb-r202-k21-n10-S16.log -j 32 --force

    # build database
    #   number of index files: 32, for server with >= 32 CPU cores
    #   bloom filter parameter:
    #     number of hash function: 1
    #     false positive rate: 0.2
    kmcp index -j 32 -I gtdb-r202-k21-n10-S16 -O gtdb.sync16.kmcp -n 1 -f 0.2 \
        --log gtdb.sync16.kmcp.log
    
    # cp taxid and name mapping file to database directory
    cp taxid.map name.map gtdb.sync16.kmcp/
    
    
    # -------------------------------------------------------------------------------------
    # FracMinhash/Scaled MinHash with d=5
    
    input=gtdb
       
    # compute k-mers
    #   sequence containing "plasmid" in name are ignored,
    #   reference genomes are splitted into 10 fragments with 100bp overlap
    #   k = 21
    #   D = 5 # FracMinhash
    kmcp compute -I $input -O gtdb-r202-k21-n10-D5 -k 21 -D 5 -n 10 -l 100 -B plasmid \
        --log gtdb-r202-k21-n10-D5.log -j 32 --force

    # build database
    #   number of index files: 32, for server with >= 32 CPU cores
    #   bloom filter parameter:
    #     number of hash function: 1
    #     false positive rate: 0.2
    kmcp index -j 32 -I gtdb-r202-k21-n10-D5 -O gtdb.minh5.kmcp -n 1 -f 0.2 \
        --log gtdb.minh5.kmcp.log
    
    # cp taxid and name mapping file to database directory
    cp taxid.map name.map gtdb.minh5.kmcp/
    
    
    # -------------------------------------------------------------------------------------
    # Minimizer with W=5
    
    input=gtdb
       
    # compute k-mers
    #   sequence containing "plasmid" in name are ignored,
    #   reference genomes are splitted into 10 fragments with 100bp overlap
    #   k = 21
    #   W = 5 # Minimizer
    kmcp compute -I $input -O gtdb-r202-k21-n10-W5 -k 21 -W 5 -n 10 -l 100 -B plasmid \
        --log gtdb-r202-k21-n10-W5.log -j 32 --force

    # build database
    #   number of index files: 32, for server with >= 32 CPU cores
    #   bloom filter parameter:
    #     number of hash function: 1
    #     false positive rate: 0.2
    kmcp index -j 32 -I gtdb-r202-k21-n10-W5 -O gtdb.mini5.kmcp -n 1 -f 0.2 \
        --log gtdb.mini5.kmcp.log
    
    # cp taxid and name mapping file to database directory
    cp taxid.map name.map gtdb.mini5.kmcp/
    
Building small databases (all k-mers, for profiling with a computer cluster):
    
    input=gtdb
    
    find $input -name "*.fna.gz" > $input.files.txt
    
    
    # number of databases
    n=16
        
    # split into $n chunks
    split -n l/$n $chunksize -d  $input.files.txt $input.n$n-
    
    # create database for every chunks
    for f in $input.n$n-*; do 
        echo $f
        
        # compute k-mers
        #   sequence containing "plasmid" in name are ignored,
        #   reference genomes are splitted into 10 fragments with 100bp overlap
        #   k = 21
        kmcp compute -i $f -O $f-k21-n10 -k 21 -n 10 -l 100 -B plasmid \
            --log $f-k21-n10.log -j 24 --force

        # build database
        #   number of index files: 24, for server with >= 24 CPU cores
        #   bloom filter parameter:
        #     number of hash function: 1
        #     false positive rate: 0.3
        kmcp index -j 24 -I $f-k21-n10 -O $f.kmcp -n 1 -f 0.3 \
            --log $f.kmcp.log
        
        # cp taxid and name mapping file to database directory
        cp taxid.map name.map $f.kmcp/
    done

### RefSeq viral or fungi

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

    # cd to 2021-09-30_19-35-19
        
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
    
    
    
    # ------------------------------------------------------ 
    # create another directory linking to genome files
    
    input=files.renamed    
    
    mkdir -p $input
    # create soft links
    cd $input
    find ../files -name "*.fna.gz" \
        | rush 'ln -s {}'
    cd ..    
    # rename
    brename -R -p '^(\w{3}_\d{9}\.\d+).+' -r '$1.fna.gz' $input
        
Building database (all k-mers, for profiling on short-reads):
        
    # -----------------------------------------------------------------
    # for viral, only splitting into 5 fragments
    name=viral
    
    input=files.renamed
    
    # -------------
    # all kmers
    
    kmcp compute -I $input -O refseq-$name-k21-n5 \
        -k 21 --seq-name-filter plasmid \
        --split-number 5 --split-overlap 100 \
        --log refseq-$name-k21-n5.log -j 32 --force
    
    # viral genomes are small:
    #   using small false positive rate: 0.001
    #   using more hash functions: 3
    kmcp index -I refseq-$name-k21-n5/ -O refseq-viral.kmcp \
        -j 32 -f 0.001 -n 3 -x 100K \
        --log refseq-viral.kmcp.log --force
    
    # cp taxid and name mapping file to database directory
    cp taxid.map name.map refseq-viral.kmcp/

    
    # -----------------------------------------------------------------
    # for fungi
    name=fungi
    
    input=files.renamed
    
    
    # -------------
    # all kmers
    
    kmcp compute -I $input -O refseq-$name-k21-n10 \
        -k 21 --seq-name-filter plasmid \
        --split-number 10 --split-overlap 100 \
        --log refseq-$name-k21-n10.log -j 32 --force
      
    kmcp index -I refseq-$name-k21-n10/ -O refseq-fungi.kmcp \
        -j 32 -f 0.3 -n 1 \
        --log refseq-fungi.kmcp.log --force
    
    # cp taxid and name mapping file to database directory
    cp taxid.map name.map refseq-fungi.kmcp/
    
   
    
Building database (k-mer sketches, for profiling on long-reads):


    # -----------------------------------------------------------------
    # for viral, only splitting into 5 fragments
    name=viral
    
    input=files.renamed

     
    # ---------------------------------------------
    # here we compute Closed Syncmers with s=16
    
    kmcp compute -I $input -O refseq-$name-k21-n5-S16 \
        -k 21 -S 16 --seq-name-filter plasmid \
        --split-number 5 --split-overlap 100 \
        --log refseq-$name-k21-n5-S16.log -j 32 --force
    
    # viral genomes are small:
    #   using small false positive rate: 0.001
    #   using more hash functions: 3
    kmcp index -I refseq-$name-k21-n5-S16/ -O refseq-viral.sync16.kmcp \
        -j 32 -f 0.001 -n 3 -x 100K \
        --log refseq-viral.sync16.kmcp.log --force
    
    # cp taxid and name mapping file to database directory
    cp taxid.map name.map refseq-viral.sync16.kmcp/
    
    
    # ---------------------------------------------
    # here we compute FracMinHash with D=5
    
    kmcp compute -I $input -O refseq-$name-k21-n5-D5 \
        -k 21 -D 5 --seq-name-filter plasmid \
        --split-number 5 --split-overlap 100 \
        --log refseq-$name-k21-n5-D5.log -j 32 --force
    
    # viral genomes are small:
    #   using small false positive rate: 0.001
    #   using more hash functions: 3
    kmcp index -I refseq-$name-k21-n5-D5/ -O refseq-viral.minh5.kmcp \
        -j 32 -f 0.001 -n 3 -x 100K \
        --log refseq-viral.minh5.kmcp.log --force
    
    # cp taxid and name mapping file to database directory
    cp taxid.map name.map refseq-viral.minh5.kmcp/
    
        
    # -----------------------------------------------------------------
    # for fungi
    name=fungi
    
    input=files.renamed
     
    # ---------------------------------------------
    # here we compute Closed Syncmers with s=16
    
    kmcp compute -I $input -O refseq-$name-k21-n10-S16 \
        -k 21 -S 16 --seq-name-filter plasmid \
        --split-number 10 --split-overlap 100 \
        --log refseq-$name-k21-n10-S16.log -j 32 --force
      
    kmcp index -I refseq-$name-k21-n10-S16/ -O refseq-fungi.sync16.kmcp \
        -j 32 -f 0.2 -n 1 \
        --log refseq-fungi.sync16.kmcp.log --force
    
    # cp taxid and name mapping file to database directory
    cp taxid.map name.map refseq-fungi.sync16.kmcp/
    
  
    # ---------------------------------------------
    # here we compute FracMinHash with D=5
    
    kmcp compute -I $input -O refseq-$name-k21-n10-D5 \
        -k 21 -D 5 --seq-name-filter plasmid \
        --split-number 10 --split-overlap 100 \
        --log refseq-$name-k21-n10-D5.log -j 32 --force
      
    kmcp index -I refseq-$name-k21-n10-D5/ -O refseq-fungi.minh5.kmcp \
        -j 32 -f 0.2 -n 1 \
        --log refseq-fungi.minh5.kmcp.log --force
    
    # cp taxid and name mapping file to database directory
    cp taxid.map name.map refseq-fungi.minh5.kmcp/
    
    
    # ---------------------------------------------
    # here we compute Minimizer with W=5
    
    kmcp compute -I $input -O refseq-$name-k21-n10-W5 \
        -k 21 -W 5 --seq-name-filter plasmid \
        --split-number 10 --split-overlap 100 \
        --log refseq-$name-k21-n10-W5.log -j 32 --force
      
    kmcp index -I refseq-$name-k21-n10-W5/ -O refseq-fungi.mini5.kmcp \
        -j 32 -f 0.2 -n 1 \
        --log refseq-fungi.mini5.kmcp.log --force
    
    # cp taxid and name mapping file to database directory
    cp taxid.map name.map refseq-fungi.mini5.kmcp/
    
### Genbank viral


Tools

- [genome_updater](https://github.com/pirovc/genome_updater) for downloading genomes from NCBI.

Downloading viral sequences:

    name=viral

    # -k for dry-run
    # -i for fix
    time genome_updater.sh \
        -d "genbank"\
        -g $name \
        -c "all" \
        -l "all" \
        -f "genomic.fna.gz" \
        -o "genbank-$name" \
        -t 12 \
        -m -a -p

    # cd genbank-viral/2021-12-06_15-27-37/
        
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
    
    
    
    # ------------------------------------------------------ 
    # create another directory linking to genome files
    
    input=files.renamed    
    
    mkdir -p $input
    # create soft links
    cd $input
    find ../files -name "*.fna.gz" \
        | rush 'ln -s {}'
    cd ..    
    # rename
    brename -R -p '^(\w{3}_\d{9}\.\d+).+' -r '$1.fna.gz' $input
        
        
keep at most 5 genomes for a taxid:


    # -----------------------------------------------------------------
    # keep at most 5 genomes for a taxid
    
    # backup
    cat taxid.map > taxid.all.map
    cat name.map > name.all.map
    
    # keep at most 5 ids for a taxid
    cat taxid.all.map \
        | csvtk sort -Ht -k 2:n -k 1:n \
        | csvtk uniq -Ht -f 2 -n 5 \
        > taxid.map
    
    cat name.all.map \
        | csvtk grep -Ht -P <(cut -f 1 taxid.map) \
        > name.map
        
    # organize files
    input=files.renamed
    output=files.renamed.slim
    mkdir -p $output 
    cd $output
    find ../$input -name "*.fna.gz" \
        | csvtk mutate -Ht -p '/([^\/]+).fna.gz' \
        | csvtk grep -Ht -f 2 -P <(cut -f 1 ../taxid.map) \
        | rush 'ln -s {1}'
    cd ..
        
    # check number of genomes
    ls $output | wc -l
    wc -l taxid.map

Building database (all k-mers, for profiling on short-reads):

    # -----------------------------------------------------------------
    # for viral, only splitting into 5 fragments
    name=viral
    
    input=files.renamed.slim
    
    
    # ----------------
    # all kmers
    
    kmcp compute -I $input -O genbank-$name-k21-n5 \
        -k 21 --seq-name-filter plasmid \
        --split-number 5 --split-overlap 100 \
        --log genbank-$name-k21-n5.log -j 32 --force
    
    # viral genomes are small:
    #   using small false positive rate: 0.05
    #   still using one hash function: 1
    kmcp index -I genbank-$name-k21-n5/ -O genbank-viral.kmcp \
        -j 32 -f 0.05 -n 1 -x 100K -8 1M \
        --log genbank-viral.kmcp.log --force
    
    # cp taxid and name mapping file to database directory
    cp taxid.map name.map genbank-viral.kmcp/


Building database (k-mer sketches, for profiling on long-reads):

    name=viral
    
    input=files.renamed.slim
    
    
    # ----------------------------------------------
    # here we compute Closed Syncmers with s=16
    
    kmcp compute -I $input -O genbank-$name-k21-n5-S16 \
        -k 21 -S 16 --seq-name-filter plasmid \
        --split-number 5 --split-overlap 100 \
        --log genbank-$name-k21-n5-S16.log -j 32 --force
    
    # viral genomes are small:
    #   using small false positive rate: 0.001
    #   using more hash functions: 3
    kmcp index -I genbank-$name-k21-n5-S16/ -O genbank-viral.sync16.kmcp \
        -j 32 -f 0.001 -n 3 -x 50K -8 1M \
        --log genbank-viral.sync16.kmcp.log --force
    
    # cp taxid and name mapping file to database directory
    cp taxid.map name.map genbank-viral.sync16.kmcp/
    
    
    # ----------------------------------------------
    # here we compute FracMinHash with D=5
    
    kmcp compute -I $input -O genbank-$name-k21-n5-D5 \
        -k 21 -D 5 --seq-name-filter plasmid \
        --split-number 5 --split-overlap 100 \
        --log genbank-$name-k21-n5-D5.log -j 32 --force
    
    # viral genomes are small:
    #   using small false positive rate: 0.001
    #   using more hash functions: 3
    kmcp index -I genbank-$name-k21-n5-D5/ -O genbank-viral.minh5.kmcp \
        -j 32 -f 0.001 -n 3 -x 50K -8 1M \
        --log genbank-viral.minh5.kmcp.log --force
    
    # cp taxid and name mapping file to database directory
    cp taxid.map name.map genbank-viral.minh5.kmcp/
    
    
    # ----------------------------------------------
    # here we compute Minimizer with W=5
    
    kmcp compute -I $input -O genbank-$name-k21-n5-W5 \
        -k 21 -W 5 --seq-name-filter plasmid \
        --split-number 5 --split-overlap 100 \
        --log genbank-$name-k21-n5-W5.log -j 32 --force
    
    # viral genomes are small:
    #   using small false positive rate: 0.001
    #   using more hash functions: 3
    kmcp index -I genbank-$name-k21-n5-W5/ -O genbank-viral.mini5.kmcp \
        -j 32 -f 0.001 -n 3 -x 50K -8 1M \
        --log genbank-viral.mini5.kmcp.log --force
    
    # cp taxid and name mapping file to database directory
    cp taxid.map name.map genbank-viral.mini5.kmcp/
    
Building small databases (all k-mers, for profiling with a computer cluster):
    
    name=genbank-viral
    
    input=files.renamed.slim
    
    find $input -name "*.fna.gz" > $input.files.txt
    
    
    # number of databases
    n=4
        
    # split into $n chunks
    split -n l/$n $chunksize -d  $input.files.txt $name.n$n-
    
    # create database for every chunks
    for f in $name.n$n-*; do 
        echo $f
  
        kmcp compute -i $f -O $f-k21-n5 \
            -k 21 --seq-name-filter plasmid \
            --split-number 5 --split-overlap 100 \
            --log $f-k21-n5.log -j 32 --force
        
        # viral genomes are small:
        #   using small false positive rate: 0.001
        #   using more hash functions: 3
        kmcp index -I $f-k21-n5/ -O $f.kmcp \
            -j 32 -f 0.05 -n 1 -x 100K -8 1M \
            --log $f.kmcp.log --force
        
        # cp taxid and name mapping file to database directory
        cp taxid.map name.map $f.kmcp/
    done
    

    
### Human genome


Downloading human genome file from [CHM13](https://github.com/marbl/CHM13):

    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.3_CHM13_T2T_v1.1/GCA_009914755.3_CHM13_T2T_v1.1_genomic.fna.gz
    
Building database (all k-mers, < 6min):

    input=GCA_009914755.3_CHM13_T2T_v1.1_genomic.fna.gz
    
    # splitting human genome into 1024 fragments.
    # The regular expression '^(\w{3}_\d{9}\.\d+).+' is for extracting 'GCA_009914755.3' from the file name.
    kmcp compute $input -O human-chm13-k21-n1024 \
        --ref-name-regexp '^(\w{3}_\d{9}\.\d+).+' \
        -k 21 \
        --split-number 1024 --split-overlap 100 \
        --log human-chm13-k21-n1024.log -j 32 --force
    
    #   using small false positive rate: 0.3
    #   using more hash functions: 1
    kmcp index -I human-chm13-k21-n1024/ -O human-chm13.kmcp \
        -j 8 -f 0.3 -n 1 \
        --log human-chm13.kmcp.log --force
    
    # taxid.map
    echo -ne "GCA_009914755.3\t9606\n" > taxid.map
    
    # name.mapp
    echo -ne "GCA_009914755.3\tHomo sapiens isolate CHM13\n" > name.map
    
    # cp name mapping file to database directory
    cp taxid.map name.map human-chm13.kmcp/
    
### Refseq plasmid

Downloading plasmid sequences:

    wget https://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/
    cat index.html \
        | perl -ne 'next unless /"(plasmid.*genomic.fna.gz)"/; print "$1\n"' > files.txt
        
    baseurl=https://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/
    outdir=archives; mkdir -p $outdir
    cat files.txt \
        | rush -j 12 -v b=$baseurl -v out=$outdir \
            'wget -c {b}/{} -o /dev/null -O {out}/{}' -c -C download.rush

            
    # name mapping
    seqkit seq -n archives/*.fna.gz \
        | sed -E 's/\s+/\t/' \
        > name.map
        
    # length stats
    seqkit stats -b -a -j 12 -T archives/*.fna.gz > archives.stats.tsv
    
    
    # split to individual files
    for f in archives/*.fna.gz; do \
        seqkit split2 -s 1 --quiet $f -O refseq-plasmid; \
    done
    # rename files with sequence ID
    find refseq-plasmid -name "*.fna.gz" \
        | rush 'mv {} {/}/$(seqkit seq -ni {}).fna.gz'
        
Building database (all k-mers):

    name=plasmid
    
    input=refseq-plasmid
    
    kmcp compute -I $input -O refseq-$name-k21-n5 \
        -k 21 --circular \
        --split-number 5 --split-overlap 100 \
        --log refseq-$name-k21-n5.log -j 32 --force
    
    #   using small false positive rate: 0.01
    #   using more hash functions: 3
    kmcp index -I refseq-$name-k21-n5/ -O refseq-$name.kmcp \
        -j 32 -f 0.01 -n 3 -x 200K -X 1024 \
        --log refseq-$name.kmcp.log --force
    
    # cp name mapping file to database directory
    cp name.map refseq-$name.kmcp/

Building database (FracMinHash/Scaled MinHash):

    name=plasmid
    
    input=refseq-plasmid
    
    kmcp compute -I $input -O refseq-$name-k31-D10 \
        -k 31 --circular --scale 10 \
        --log refseq-$name-k31-D10.log -j 32 --force
    
    #   using small false positive rate: 0.01
    #   using more hash functions: 3
    kmcp index -I refseq-$name-k31-D10/ -O refseq-$name.minhash.kmcp \
        -j 8 -f 0.01 -n 3 -x 100K -X 1024 \
        --log refseq-$name.minhash.kmcp.log --force
    
    # cp name mapping file to database directory
    cp name.map refseq-$name.minhash.kmcp/

Building database (Closed Syncmer):

    name=plasmid
    
    input=refseq-plasmid
    
    kmcp compute -I $input -O refseq-$name-k31-S21 \
        -k 31 --circular --syncmer-s 21 \
        --log refseq-$name-k31-S21.log -j 32 --force
    
    #   using small false positive rate: 0.01
    #   using more hash functions: 3
    kmcp index -I refseq-$name-k31-S21/ -O refseq-$name.syncmer.kmcp \
        -j 8 -f 0.01 -n 3 -x 100K -X 1024 \
        --log refseq-$name.syncmer.kmcp.log --force
    
    # cp name mapping file to database directory
    cp name.map refseq-$name.syncmer.kmcp/

### HumGut


[HumGut](http://arken.nmbu.no/~larssn/humgut/) is a comprehensive Human Gut prokaryotic genomes collection filtered by metagenome data.

Dataset

- Genomes: [HumGut.tar.gz](http://arken.nmbu.no/~larssn/humgut/HumGut.tar.gz)
- Metadata: [HumGut.tsv](http://arken.nmbu.no/~larssn/humgut/HumGut.tsv)
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

    # -------------------------------------------------------------
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
    
    # -------------------------------------------------------------
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
    
    # compute k-mers
    #   sequence containing "plasmid" in name are ignored,
    #   reference genomes are splitted into 10 fragments
    #   k = 21
    kmcp compute -I fna/ -k 21 -n 10 -B plasmid -O humgut-k21-n10 -j 32 --force

    # build database
    #   number of index files: 32, for server with >= 32 CPU cores
    #   bloom filter parameter:
    #     number of hash function: 1
    #     false positive rate: 0.3
    kmcp index -j 32 -I humgut-k21-n10 -O humgut.kmcp -n 1 -f 0.3
    
    cp taxid.map taxid-gtdb.map humgut.kmcp/

### proGenomes2

[proGenomes](https://progenomes.embl.de/) v2.1 provides 84,096 consistently annotated bacterial
and archaeal genomes from over 12000 species

Dataset:

- Representative genomes: [contigs.representatives.fasta.gz](https://progenomes.embl.de/data/repGenomes/freeze12.contigs.representatives.fasta.gz)
- NCBI taxonomy: [proGenomes2.1_specI_lineageNCBI.tab](https://progenomes.embl.de/data/proGenomes2.1_specI_lineageNCBI.tab)
    
Organize sequences:
    
    # download
    wget https://progenomes.embl.de/data/repGenomes/freeze12.contigs.representatives.fasta.gz
    
    # unzip for splitting by genome ID
    time seqkit seq freeze12.contigs.representatives.fasta.gz -o freeze12.contigs.representatives.fasta
    
    # split by genome ID
    seqkit split --by-id --two-pass --id-regexp '(^\d+\.\w+)\.' freeze12.contigs.representatives.fasta --out-dir genomes
    
    # batch rename
    brename -p '^.+id_' genomes
    
    # compress genomes for saving space
    find genomes/ -name "*.fasta" | rush 'gzip {}'

Mapping file:

    # download
    wget https://progenomes.embl.de/data/proGenomes2.1_specI_lineageNCBI.tab
    
    ls genomes/ | sed 's/.fasta.gz//' > id.txt
    
    # id -> taxid
    cut -f 1 proGenomes2.1_specI_lineageNCBI.tab \
        | awk -F . '{print $0"\t"$1}' \
        | csvtk grep -Ht -P id.txt \
        > taxid.map
    
    cat taxid.map \
        | csvtk uniq -Ht -f 2 \
        | taxonkit lineage --data-dir taxdump -i 2 -r -n -L \
        | csvtk freq -Ht -f 4 -nr \
        | csvtk pretty -H -t \
        | tee taxid.map.uniq.stats
    species           8088
    strain            3262
    subspecies        37
    isolate           25
    no rank           16
    forma specialis   14
    biotype           1
    
    
    # use the taxonomy verion of the refseq-virus: 2021-10-01
    # id -> name
    cat taxid.map \
        | taxonkit lineage --data-dir taxdump -i 2 -n -L \
        | cut -f 1,3 \
        > name.map

Building database:
    
    # compute k-mers
    #   sequence containing "plasmid" in name are ignored,
    #   reference genomes are splitted into 10 fragments
    #   k = 21
    kmcp compute -I genomes/ -k 21 -n 10 -B plasmid -O progenomes-k21-n10 --force

    # build database
    #   number of index files: 32, for server with >= 32 CPU cores
    #   bloom filter parameter:
    #     number of hash function: 1
    #     false positive rate: 0.3
    kmcp index -j 32 -I progenomes-k21-n10 -O progenomes.kmcp -n 1 -f 0.3
    
    cp taxid.map name.map progenomes.kmcp/
        
## Building custom databases

Files:

1. Genome files
    - (Gzip-compressed) FASTA/Q format.
    - One genome per file with the reference identifier in the file name.
2. TaxId mapping file (for metagenomic profiling)
    - Two-column (reference identifier and TaxId) tab-delimited.
3. [NCBI taxonomy dump files](ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz) (for metagenomic profiling)
    - `names.dmp`
    - `nodes.dmp`
    - `merged.dmp` (optional)
    - `delnodes.dmp` (optional)

Tools:

- [kmcp](/download) for metagenomic profiling.
- [rush](https://github.com/shenwei356/rush/releases) for executing jobs in parallel.
- [brename](https://github.com/shenwei356/brename/releases) for batching renaming files (optional)


**Memory notes**:

By default, `kmcp search` loads the whole database into main memory (RAM) for fast searching.
Optionally, the flag `--low-mem` can be set to avoid loading the whole database,
while it's much slower, >10X slower on SSD and should be much slower on HDD disks.

Another way is dividing the reference genomes into several parts
and building smaller databases for all parts, so that the biggest
database can be loaded into RAM. After performing database searching,
search results on all small databases can be merged with `kmcp merge`
for downstream analysis.

Buiding small databases can also accelerate searching on a *computer cluster*,
where every node searches a part of the database.


### Step 1. Computing k-mers

Plain or gzip-compressed input genome sequence files better be saved in one directory,
with multiple-level directories allowed. 
The file extension could be `fa`, `fasta`, `fna`, `fq`, `fastq`, `fa.gz`, `fasta.gz`,
 `fna.gz`, `fq.gz` or `fastq.gz`.

**Input files** can be given as a list of FASTA/Q files via
positional arguments or a directory (with multiple-level directories allowed) containing sequence files
via the flag `-I/--in-dir`. A regular expression for matching
sequencing files is available by the flag `-r/--file-regexp`.
The default pattern matches files with extension of 
`fa`, `fasta`, `fna`, `fq`, `fastq`, `fa.gz`, `fasta.gz`,
`fna.gz`, `fq.gz` or `fastq.gz`.

Unwanted sequence like plasmid sequences can be **filtered out** by
the name via regular expression(s) (`-B/--seq-name-filter`).

**How to compute k-mers**:
By default, `kmcp` computes k-mers (sketches) of every file,
you can also use `--by-seq` to compute for every sequence,
where sequence IDs in all input files better be distinct.
It also **supports splitting sequences into fragments, this
could increase the specificity in profiling result in cost
of slower searching speed**. 

**Splitting sequences**:

1. Sequences can be splitted into fragments by a fragment size 
    (`-s/--split-size`) or number of fragments (`-n/--split-number`)
    with overlap (`-l/--split-overlap`).
    ***In this mode, the sequences of each genome should be saved in an
    individual file***.
2. When splitting by number of fragments, **all sequences (except for
    these mathching any regular expression given by `-B/--seq-name-filter`)
    in a sequence file are concatenated with k-1 Ns before splitting**.
3. Both sequence/reference IDs and fragments indices are saved for later use,
    in form of meta/description data in `.unik` files, and will
    be reported in `kmcp search` results.

**Meta data**:

1. Every outputted `.unik` file contains the sequence/reference ID,
    fragment index, number of fragments, and genome size of reference.
2. When parsing whole sequence files or splitting by number of fragments,
    the **identifier of a reference** is the basename of the input file
    by default. It can also be extracted from the input file name via
    `-N/--ref-name-regexp`, e.g., `^(\w{3}_\d{9}\.\d+)` for **RefSeq assembly accessions**.

Multiple **sizes of k-mers** are supported, but a single k-mer size
is good enough. 

**Supported k-mer (sketches) types**:

1. K-mer:
    - ntHash of k-mer (`-k`)
2. K-mer sketchs (all using ntHash):
    - Scaled MinHash (`-k -D`), reviously named Scaled MinHash
    - Minimizer      (`-k -W`), optionally scaling/down-sampling (`-D`)
    - Closed Syncmer (`-k -S`), optionally scaling/down-sampling (`-D`)

**Output**:

1. All outputted `.unik` files are saved in `${outdir}`, with path

        ${outdir}/xxx/yyy/zzz/${infile}-id_${seqID}.unik

    where dirctory tree `/xxx/yyy/zzz/` is built for > 1000 output files.

2. For splitting sequence mode (`--split-size > 0` or `--split-number > 0`),
    output files are:
    
        ${outdir}//xxx/yyy/zzz/${infile}/{seqID}-frag_${fragIdx}.unik
        
3. A summary file (`${outdir}/_info.txt`) is generated for later use.
     ***Users need to check if the reference IDs (column `name`) are what
     supposed to be***.

**Performance tips**:

1. Decrease value of `-j/--threads` for data in hard disk drives (HDD) to
    reduce I/O pressure.

**Commands** ([usage](/usage/#compute)):

    # compute k-mers
    #   sequence containing "plasmid" in name are ignored,
    #   reference genomes are splitted into 10 fragments,
    #   k = 21
    kmcp compute --in-dir refs/ \
        --kmer 21 \
        --split-number 10 \
        --seq-name-filter plasmid \
        --ref-name-regexp '(.+).fasta.gz' \
        --out-dir refs-k21-n10

    # demo output
    22:33:10.685 [INFO] kmcp v0.7.0
    22:33:10.685 [INFO]   https://github.com/shenwei356/kmcp
    22:33:10.685 [INFO] 
    22:33:10.685 [INFO] checking input files ...
    22:33:10.685 [INFO]   9 input file(s) given
    22:33:10.685 [INFO] 
    22:33:10.685 [INFO] -------------------- [main parameters] --------------------
    22:33:10.685 [INFO] input and output:
    22:33:10.685 [INFO]   input directory: refs/
    22:33:10.685 [INFO]     regular expression of input files: (?i)\.(f[aq](st[aq])?|fna)(.gz)?$
    22:33:10.685 [INFO]     *regular expression for extracting reference name from file name: 
    22:33:10.685 [INFO]     *regular expressions for filtering out sequences: [plasmid]
    22:33:10.685 [INFO]   output directory: refs-k21-n10
    22:33:10.685 [INFO] 
    22:33:10.685 [INFO] sequences splitting: true
    22:33:10.685 [INFO]   split parts: 10, overlap: 0 bp
    22:33:10.685 [INFO] 
    22:33:10.685 [INFO] k-mer (sketches) computing:
    22:33:10.685 [INFO]   k-mer size(s): 21
    22:33:10.685 [INFO]   circular genome: false
    22:33:10.685 [INFO]   saving exact number of k-mers: true
    22:33:10.685 [INFO] 
    22:33:10.685 [INFO] -------------------- [main parameters] --------------------
    22:33:10.685 [INFO] 
    22:33:10.685 [INFO] computing ...
    processed files:  9 / 9 [======================================] ETA: 0s. done
    22:33:11.121 [INFO] 
    22:33:11.121 [INFO] elapsed time: 436.367564ms
    22:33:11.121 [INFO]

A summary file (`_info.txt`) is generated for later use.
***Users need to check if the reference IDs (column `name`) are what supposed to be***.

|#path                                                                |name       |fragIdx|idxNum|genomeSize|kmers |
|:--------------------------------------------------------------------|:----------|:------|:-----|:---------|:-----|
|refs-k21-n10/672/060/672/NC_010655.1.fasta.gz/NC_010655.1-frag_0.unik|NC_010655.1|0      |10    |2664102   |264247|
|refs-k21-n10/672/060/672/NC_010655.1.fasta.gz/NC_010655.1-frag_1.unik|NC_010655.1|1      |10    |2664102   |266237|
|refs-k21-n10/595/162/595/NC_012971.2.fasta.gz/NC_012971.2-frag_0.unik|NC_012971.2|0      |10    |4558953   |450494|
|refs-k21-n10/292/039/292/NC_000913.3.fasta.gz/NC_000913.3-frag_0.unik|NC_000913.3|0      |10    |4641652   |459277|
|refs-k21-n10/934/859/934/NC_013654.1.fasta.gz/NC_013654.1-frag_0.unik|NC_013654.1|0      |10    |4717338   |470575|

Meta data in the `.unik` file can be showed using `kmcp utils unik-info`:

    kmcp utils unik-info refs-k21-n10/072/380/072/NZ_CP028116.1.fasta.gz/NZ_CP028116.1-frag_0.unik -a


### Step 2. Building databases

KMCP builds index for k-mers (sketches) with a modified Compact Bit-sliced
Signature Index ([COBS](https://arxiv.org/abs/1905.09624)). 
We totally rewrite the algorithms, data structure and file format,
and have improved the indexing and searching speed
(check [benchmark](/benchmark/searching)).

**Input**:

- The output directory generated by `kmcp compute`.

**Database size and searching accuracy**:

0. **Use `--dry-run` to adjust parameters and check final number of 
    index files (#index-files) and the total file size**.
1. `-f/--false-positive-rate`: the default value `0.3` is enough for a
    query with tens of matched k-mers (see BIGSI/COBS paper).
    Small values could largely increase the size of database.
2. `-n/--num-hash`: large values could reduce the database size,
    in cost of slower searching speed. Values <=4 is recommended.
3. Value of block size `-b/--block-size` better be multiple of 64.
    The default value is: 
    
        (#unikFiles/#threads + 7) / 8 * 8

4. Use flag `-x/--block-sizeX-kmers-t`, `-8/--block-size8-kmers-t`,
    and `-1/--block-size1-kmers-t` to separately create index for
    inputs with huge number of k-mers, for precise control of
    database size.

**Taxonomy data**:

1. No taxonomy data are included in the database.
2. Taxonomy information are only needed in `kmcp profile`.

**Performance tips**:

1. Number of blocks (`.uniki` files) better be smaller than or equal
   to number of CPU cores for faster searching speed. 
   **We can set `-j/--threads` to control the blocks number**.
   When more threads (>= 1.3 * #blocks) are given, extra workers are
   automatically created.
2. `#threads` files are simultaneously opened, and the max number
    of opened files is limited by the flag `-F/--max-open-files`.
    You may use a small value of `-F/--max-open-files` for 
    hard disk drive storage.
3. When the database is used in a new computer with more CPU cores,
   `kmcp search` could automatically scale to utilize as many cores
   as possible.

**Commands** ([usage](/usage/#compute)):

    # build database
    #   number of index files: 32, for server with >= 32 CPU cores
    #   bloom filter parameter:
    #     number of hash function: 1
    #     false positive rate: 0.3
    kmcp index -I refs-k21-n10 \
        --threads 32 \
        --num-hash 1 \
        --false-positive-rate 0.3 \
        --out-dir refs.kmcp 

    # demo output
    22:44:17.710 [INFO] kmcp v0.7.0
    22:44:17.710 [INFO]   https://github.com/shenwei356/kmcp
    22:44:17.710 [INFO] 
    22:44:17.710 [INFO] loading .unik file infos from file: refs-k21-n10/_info.txt
    22:44:17.716 [INFO]   90 cached file infos loaded
    22:44:17.716 [INFO] 
    22:44:17.716 [INFO] -------------------- [main parameters] --------------------
    22:44:17.716 [INFO]   number of hashes: 1
    22:44:17.716 [INFO]   false positive rate: 0.300000
    22:44:17.717 [INFO]   k-mer size(s): 21
    22:44:17.717 [INFO]   split seqequence size: 0, overlap: 0
    22:44:17.717 [INFO]   block-sizeX-kmers-t: 10.00 M
    22:44:17.717 [INFO]   block-sizeX        : 256.00
    22:44:17.717 [INFO]   block-size8-kmers-t: 20.00 M
    22:44:17.717 [INFO]   block-size1-kmers-t: 200.00 M
    22:44:17.717 [INFO] -------------------- [main parameters] --------------------
    22:44:17.717 [INFO] 
    22:44:17.717 [INFO] building index ...
    22:44:17.726 [WARN] ignore -X/--block-size (256) which is >= -b/--block-size (8)
    22:44:17.726 [INFO] 
    22:44:17.726 [INFO]   block size: 8
    22:44:17.726 [INFO]   number of index files: 32 (may be more)
    22:44:17.726 [INFO] 
    22:44:17.726 [block #001]   1 / 1  100 %
    22:44:17.726 [block #002]   1 / 1  100 %
    22:44:17.726 [block #003]   1 / 1  100 %
    22:44:17.726 [block #004]   1 / 1  100 %
    22:44:17.726 [block #005]   1 / 1  100 %
    22:44:17.726 [block #006]   1 / 1  100 %
    22:44:17.726 [block #007]   1 / 1  100 %
    22:44:17.726 [block #008]   1 / 1  100 %
    22:44:17.726 [block #009]   1 / 1  100 %
    22:44:17.727 [block #010]   1 / 1  100 %
    22:44:17.727 [block #011]   1 / 1  100 %
    22:44:17.727 [block #012]   1 / 1  100 %
    [saved index files] 12 / 12  ETA: 0s. done
    22:44:17.933 [INFO] 
    22:44:17.933 [INFO] kmcp database with 42713316 k-mers saved to refs.kmcp
    22:44:17.933 [INFO] total file size: 15.66 MB
    22:44:17.933 [INFO] total index files: 12
    22:44:17.933 [INFO] 
    22:44:17.933 [INFO] elapsed time: 223.524128ms
    22:44:17.933 [INFO]

Output:

    refs.kmcp/
    └── R001
        ├── _block001.uniki
        ├── _block002.uniki
        ├── _block003.uniki
        ├── _block004.uniki
        ├── _block005.uniki
        ├── _block006.uniki
        ├── _block007.uniki
        ├── _block008.uniki
        ├── _block009.uniki
        ├── _block010.uniki
        ├── _block011.uniki
        ├── _block012.uniki
        ├── __db.yml
        └── __name_mapping.tsv

`__db.yml` contains configuration of the database, and `_blockXXX.uniki` are index files.
`kmcp utils index-info` could show the basic information of index files:

    kmcp utils index-info  refs.kmcp/R001/_block001.uniki

|file                          |k  |canonical|num-hashes|num-sigs|num-names|
|:-----------------------------|:--|:--------|:---------|:-------|:--------|
|refs.kmcp/R001/_block001.uniki|21 |true     |1         |746442  |8        |

What's next? Check the [tutorials](/tutorial).
