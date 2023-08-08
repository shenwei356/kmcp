# Database

## Prebuilt databases

All prebuilt databases and the used reference genomes are available at:

- [OneDrive](https://1drv.ms/u/s!Ag89cZ8NYcqtjHwpe0ND3SUEhyrp?e=QDRbEC) for *global users*.
- [CowTransfer](https://shenwei356.cowtransfer.com/s/c7220dd5901c42) for *Chinese users and global users*.<br>
  **Please click the "kmcp" icon on the left to browse directories of the current version, and choose each individual file to download.**.<br>

<p style="color:Tomato;">Please check file integrity with `md5sum` after downloading the files:</p>

    md5sum -c *.kmcp.tar.gz.md5.txt
  

### A). Databases for metagenomic profiling

These databases are created following [steps below](#building-databases).
Users can also [build custom databases](#building-custom-databases), it's simple and fast.

**Current version**: v2023.05 ([OneDrive](https://1drv.ms/f/s!Ag89cZ8NYcqtlgAtJo--uKPNVT4t?e=KEDFrc), [CowTransfer](https://shenwei356.cowtransfer.com/s/1c5dfb28218748))

|kingdoms            |source      |#NCBI-species| #assemblies|db-parameters                     |size     |
|:-------------------|:-----------|:------------|:-----------|:---------------------------------|:--------|
|Bacteria and Archaea|GTDB r214   |34395+       |85205       |k=21, chunks=10; fpr=0.3, hashes=1|50+49 GB* |
|Fungi               |Refseq r217 |491          |496         |k=21, chunks=10; fpr=0.3, hashes=1|5.5 GB   |
|Viruses             |GenBank r255|26680        |33479       |k=21, chunks=5; fpr=0.05, hashes=1|5.9 GB   |

**Notes**: *GTDB representative genomes were split into 2 parts to build relatively small databases
which can be filled into workstations with small RAM (around 64GB).
Users need to search reads on all these small databases and merge the results before running `kmcp profile`.

**Taxonomy information**: Either NCBI taxonomy or a combination of GTDB and NCBI are available:

- NCBI: taxdmp_2023-05-01
- GTDB+NCBI: GTDB r214 + NCBI taxdmp_2023-05-01

**Hardware requirements**

- Prebuilt databases were built for computers with >= 32 CPU cores
in consideration of better parallelization,
and computers should have at least 64 GB.
- By default, `kmcp search` loads the whole database into main memory
(via [mmap](https://en.wikipedia.org/wiki/Mmap) by default) for fast searching.
Optionally, the flag `--low-mem` can be set to avoid loading the whole database,
while it's much slower, >10X slower on SSD and should be much slower on HDD disks.
- **To reduce memory requirements on computers without enough memory,
users can split the reference genomes into partitions
and build a smaller database for each partition, so that the biggest
database can be loaded into RAM**.
This can also **accelerate searching on a computer cluster, where every node searches against a small database.
After performing database searching,
search results from all small databases can be merged with `kmcp merge`
for downstream analysis**.


### B). Databases for genome similarity estimation

Check the [tutorial](/kmcp/tutorial/searching).

FracMinHash (Scaled MinHash), links: [OneDrive](https://1drv.ms/f/s!Ag89cZ8NYcqtlgEtskAjVVosItE7?e=CdkwhG), [CowTransfer](https://shenwei356.cowtransfer.com/s/21f68041334440)

|kingdoms            |source      |#NCBI-species| #assemblies|parameters      |size      |
|:-------------------|:-----------|:------------|:-----------|:---------------|:---------|
|Bacteria and Archaea|GTDB r214   |34395+       |85205       |k=31, scale=1000|2.7 GB    |
|Fungi               |Refseq r217 |491          |496         |k=31, scale=1000|131 MB    |
|**Viruses**         |GenBank r255|26680        |33479       |k=31, scale=10  |2.0 GB    |


## Building databases

### GTDB

Tools:

- [brename](https://github.com/shenwei356/brename/releases) for batching renaming files.
- [rush](https://github.com/shenwei356/rush/releases) for executing jobs in parallel.
- [seqkit](https://github.com/shenwei356/seqkit/releases) for FASTA file processing.
- [csvtk](https://github.com/shenwei356/csvtk/releases) for tsv/csv data manipulations.
- [taxonkit](https://github.com/shenwei356/taxonkit/releases) for NCBI taxonomy data manipulations.
- [kmcp](/kmcp/download) for metagenomic profiling.

Files (https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/)

- representative genomes: gtdb_genomes_reps.tar.gz
- metagdata: ar53_metadata.tar.gz
- metagdata: bac120_metadata.tar.gz
- taxdump file: https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_archive/

Uncompressing and renaming:
 
    # uncompress
    mkdir -p gtdb
    tar -zxvf gtdb_genomes_reps.tar.gz -C gtdb
    
    # rename
    brename -R -p '^(\w{3}_\d{9}\.\d+).+' -r '$1.fna.gz' gtdb
  
Mapping file:

    # assembly accession -> full head
    find gtdb/ -name "*.fna.gz" \
        | rush -k 'echo -ne "{%@(.+).fna}\t$(seqkit sort --quiet -lr {} | head -n 1 | seqkit seq -n)\n" ' \
        > name.map
        
    tar -zxvf ar53_metadata.tar.gz
    tar -zxvf bac120_metadata.tar.gz

    # assembly accession -> taxid
    cat *_metadata*.tsv \
        | csvtk cut -t -f accession,ncbi_taxid \
        | csvtk replace -t -p '^.._' \
        | csvtk grep -t -P <(cut -f 1 name.map) \
        | csvtk del-header \
        > taxid.map

    
    # stats (optional)
    # taxdump is a directory containing NCBI taxdump files, including names.dmp, nodes.dmp, delnodes.dmp and merged.dmp
    cat taxid.map  \
        | csvtk freq -Ht -f 2 -nr \
        | taxonkit lineage -r -n -L --data-dir taxdump/ \
        | taxonkit reformat -I 1 -f '{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}' --data-dir taxdump/ \
        | csvtk add-header -t -n 'taxid,count,name,rank,superkindom,phylum,class,order,family,genus,species' \
        > gtdb.taxid.map.stats.tsv
        
    
    # number of unique species/strains
    cat taxid.map \
        | csvtk uniq -Ht -f 2 \
        | taxonkit lineage --data-dir taxdump/ -i 2 -r -n -L \
        | csvtk freq -Ht -f 4 -nr \
        | csvtk pretty -H -t \
        | tee taxid.map.uniq.stats
    species           31165
    strain            4073
    subspecies        104
    forma specialis   59
    no rank           31
    isolate           26

Masking prophage regions and removing plasmid sequences with [genomad](https://github.com/apcamargo/genomad) (optional):

    conda activate genomad

    # https://github.com/shenwei356/cluster_files
    cluster_files -p '^(.+).fna.gz$' gtdb -o gtdb.genomad

    # indentify prophages and plasmids sequences
    # > 10 days
    db=~/ws/db/genomad/genomad_db
    find gtdb.genomad/ -name "*.fna.gz" \
        | rush -j 6 -v db=$db \
            'genomad end-to-end --relaxed --cleanup --threads 25 {} {/}/genomad {db}' \
            -c -C genomad.rush --eta

    # mask proviruses and remove plasmid sequences
    find gtdb.genomad/ -name "*.fna.gz" \
        | grep -v masked \
        | rush -j 10 -v 'p={/}/genomad/{/%}_summary/{/%}' \
            'csvtk grep -t -f coordinates -p NA {p}_virus_summary.tsv \
                | sed 1d | cut -f 1 \
                > {p}.v.id; \
             csvtk grep -t -f coordinates -p NA {p}_virus_summary.tsv -v \
                | sed 1d | cut -f 1 \
                | perl -ne "/^(.+)\|.+_(\d+)_(\d+)$/; \$s=\$2-1; print qq/\$1\t\$s\t\$3\n/;" \
                > {p}.v.bed; \
            cut -f 1 {p}_plasmid_summary.tsv | sed 1d > {p}.p.id; \
            bedtools maskfasta -fi <(seqkit grep --quiet -v -f <(cat {p}.p.id {p}.v.id) {}) \
                                -bed {p}.v.bed -fo {}.masked.fna; \
            gzip -f {}.masked.fna' \
            -c -C genomad.b.rush --eta

    # replace masked sequence with k-Ns, for more accurate genome size.
    find gtdb.genomad/ -name "*.fna.gz.masked.fna.gz" \
        | rush -j 20 'seqkit replace -s -i -p N+ -r NNNNNNNNNNNNNNNNNNNNN {} -o {...}.masked2.fna.gz' \
            -c -C genomad.c.rush --eta

    # stats
    time find gtdb.genomad/ -name "*.fna.gz" \
        | seqkit stats -j 20 --infile-list - -T -b -o gtdb.genomad.stats.tsv

    # create symbol links in another directory, and rename
    mkdir -p gtdb_masked
    cluster_files  -p '^(.+).fna.gz.masked2.fna.gz$' gtdb.genomad -o gtdb_masked
    brename -R -p .masked.fna.gz$ gtdb_masked -d
        
Building database (all k-mers, for profiling on short-reads):
    
    input=gtdb
    
    # compute k-mers
    #   sequence containing "plasmid" in name are ignored,
    #   reference genomes are split into 10 chunks with 100bp overlap
    #   k = 21
    kmcp compute -I $input -O gtdb-k21-n10 -k 21 -n 10 -l 150 -B plasmid \
        --log gtdb-k21-n10.log -j 32 --force

    # build database
    #   number of index files: 32, for server with >= 32 CPU cores
    #   bloom filter parameter:
    #     number of hash function: 1
    #     false positive rate: 0.3
    kmcp index -j 32 -I gtdb-k21-n10 -O gtdb.kmcp -n 1 -f 0.3 \
        --log gtdb.kmcp.log --force
    
    # cp taxid and name mapping file to database directory
    cp taxid.map name.map gtdb.kmcp/

    # clean up
    rm -rf gtdb-k21-n10
    
Building small databases (all k-mers, for profiling with a computer cluster or computer with limited RAM):
    
    input=gtdb
    
    find $input -name "*.fna.gz" > $input.files.txt
    
    # sort files by genome size, so we can split them into chunks with similar genome sizes
    cp $input.files.txt $input.files0.txt
    cat $input.files0.txt \
        | rush --eta -k 'echo -e {}"\t"$(seqkit stats -T {} | sed 1d | cut -f 5)' \
        > $input.files0.size.txt
    cat $input.files0.size.txt \
        | csvtk sort -Ht -k 2:nr \
        | csvtk cut -t -f 1 \
        > $input.files.txt
    
    # number of databases
    n=2
        
    # split into $n chunks using round robin distribution
    split -n r/$n -d  $input.files.txt $input.n$n-
    
    # create database for every chunks
    for f in $input.n$n-*; do 
        echo $f
        
        # compute k-mers
        #   sequence containing "plasmid" in name are ignored,
        #   reference genomes are split into 10 chunks with 100bp overlap
        #   k = 21
        kmcp compute -i $f -O $f-k21-n10 -k 21 -n 10 -l 150 -B plasmid \
            --log $f-k21-n10.log --force

        # build database
        #   number of index files: 32, for server with >= 32 CPU cores
        #   bloom filter parameter:
        #     number of hash function: 1
        #     false positive rate: 0.3
        kmcp index -j 32 -I $f-k21-n10 -O $f.kmcp -n 1 -f 0.3 \
            --log $f.kmcp.log --force
        
        # cp taxid and name mapping file to database directory
        cp taxid.map name.map $f.kmcp/
    done

    # rename small databases
    # for example:
    #   [INFO] checking: [ ok ] 'gtdb.n2-00.kmcp' -> 'gtdb.part_1.kmcp'
    #   [INFO] checking: [ ok ] 'gtdb.n2-01.kmcp' -> 'gtdb.part_2.kmcp'
    #   [INFO] 2 path(s) to be renamed

    brename -D --only-dir -p "n$n-\d+" -r "part_{nr}" -d

    
    # clean up
    for f in $input.n$n-*; do
        rm -rf $f-k21-n10
    done

### RefSeq viral or fungi

Tools

- [genome_updater](https://github.com/pirovc/genome_updater) (0.5.1) for downloading genomes from NCBI.

Downloading viral and fungi sequences:

    name=fungi
    # name=viral
    
    # -k for dry-run
    # -i for fix
    time genome_updater.sh \
        -d "refseq"\
        -g $name \
        -c "" \
        -l "" \
        -f "genomic.fna.gz" \
        -o "refseq-$name" \
        -t 12 \
        -m -a -p

    # cd to 2021-09-30_19-35-19
        
    # taxdump
    mkdir -p taxdump
    tar -zxvf taxdump.tar.gz -C taxdump
    
    # assembly accession -> taxid
    cut -f 1,6 assembly_summary.txt > taxid.map    
    # assembly accession -> name
    cut -f 1,8 assembly_summary.txt > name.map
    
    
    # stats (optional)
    cat taxid.map  \
        | csvtk freq -Ht -f 2 -nr \
        | taxonkit lineage -r -n -L --data-dir taxdump/ \
        | taxonkit reformat -I 1 -f '{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}' --data-dir taxdump/ \
        | csvtk add-header -t -n 'taxid,count,name,rank,superkindom,phylum,class,order,family,genus,species' \
        > refseq-$name.taxid.map.stats.tsv
    
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
    # for viral, only splitting into 10 chunks
    name=viral
    
    input=files.renamed
    
    # -------------
    # all kmers
    
    kmcp compute -I $input -O refseq-$name-k21-n10 \
        -k 21 --seq-name-filter plasmid \
        --split-number 10 --split-overlap 150 \
        --log refseq-$name-k21-n10.log -j 32 --force
    
    # viral genomes are small:
    #   using small false positive rate: 0.001
    #   using more hash functions: 3
    kmcp index -I refseq-$name-k21-n10/ -O refseq-viral.kmcp \
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
        --split-number 10 --split-overlap 150 \
        --log refseq-$name-k21-n10.log --force
      
    kmcp index -I refseq-$name-k21-n10/ -O refseq-fungi.kmcp \
        -j 32 -f 0.3 -n 1 \
        --log refseq-fungi.kmcp.log --force
    
    # cp taxid and name mapping file to database directory
    cp taxid.map name.map refseq-fungi.kmcp/
    

### Genbank viral


Tools

- [genome_updater](https://github.com/pirovc/genome_updater) for downloading genomes from NCBI.

Downloading viral sequences:

    name=viral

    # -k for dry-run
    # -i for fix
    # keep at most 5 genomes for a taxid:
    #    genome_updater v0.5.1: -A 5 -c "" -l ""
    #    genome_updater v0.2.5: -j taxids:5 -c "all" -l "all"
    time genome_updater.sh \
        -d "genbank"\
        -A 5 \
        -g $name \
        -c "" \
        -l "" \
        -f "genomic.fna.gz" \
        -o "genbank-$name" \
        -t 12 \
        -m -a -p

    # cd genbank-viral/2021-12-06_15-27-37/
        
    # taxdump
    mkdir -p taxdump
    tar -zxvf taxdump.tar.gz -C taxdump
    
    # assembly accession -> taxid
    cut -f 1,6 assembly_summary.txt > taxid.map    
    # assembly accession -> name
    cut -f 1,8 assembly_summary.txt > name.map
    
    
    # stats (optional)
    cat taxid.map  \
        | csvtk freq -Ht -f 2 -nr \
        | taxonkit lineage -r -n -L --data-dir taxdump/ \
        | taxonkit reformat -I 1 -f '{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}' --data-dir taxdump/ \
        | csvtk add-header -t -n 'taxid,count,name,rank,superkindom,phylum,class,order,family,genus,species' \
        > genbank-viral.taxid.map.stats.tsv
    
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
        
        
Keep at most 5 genomes for a taxid (optional)


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
    name=viral
    
    # input=files.renamed.slim
    input=files.renamed
    
    
    # ----------------
    # all kmers
    
    kmcp compute -I $input -O genbank-$name-k21-n10 \
        -k 21 --seq-name-filter plasmid \
        --split-number 10 --split-overlap 150 \
        --log genbank-$name-k21-n10.log --force
    
    # viral genomes are small:
    #   using small false positive rate: 0.05
    #   still using one hash function: 1
    kmcp index -I genbank-$name-k21-n10/ -O genbank-viral.kmcp \
        -j 32 -f 0.05 -n 1 -x 100K -8 1M \
        --log genbank-viral.kmcp.log --force
    
    # cp taxid and name mapping file to database directory
    cp taxid.map name.map genbank-viral.kmcp/

Building small databases (all k-mers, for profiling with a computer cluster or computer with limited RAM):
    
    name=genbank-viral
    
    input=files.renamed.slim
    
    find $input -name "*.fna.gz" > $input.files.txt
    
    # sort files by genome size, so we can split them into chunks with similar genome sizes
    cp $input.files.txt $input.files0.txt
    cat $input.files0.txt \
        | rush -k 'echo -e {}"\t"$(seqkit stats -T {} | sed 1d | cut -f 5)' \
        > $input.files0.size.txt
    cat $input.files0.size.txt \
        | csvtk sort -Ht -k 2:nr \
        | csvtk cut -t -f 1 \
        > $input.files.txt
        
    # number of databases
    n=4
        
    # split into $n chunks using round robin distribution
    split -n r/$n -d  $input.files.txt $name.n$n-
    
    # create database for every chunks
    for f in $name.n$n-*; do 
        echo $f
  
        kmcp compute -i $f -O $f-k21-n10 \
            -k 21 --seq-name-filter plasmid \
            --split-number 10 --split-overlap 150 \
            --log $f-k21-n10.log -j 24 --force
        
        # viral genomes are small:
        #   using small false positive rate: 0.001
        #   using more hash functions: 3
        kmcp index -I $f-k21-n10/ -O $f.kmcp \
            -j 24 -f 0.05 -n 1 -x 100K -8 1M \
            --log $f.kmcp.log --force
        
        # cp taxid and name mapping file to database directory
        cp taxid.map name.map $f.kmcp/
    done


### Human genome


Downloading human genome file from [CHM13](https://github.com/marbl/CHM13):

    # v1.1: wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.3_CHM13_T2T_v1.1/GCA_009914755.3_CHM13_T2T_v1.1_genomic.fna.gz
    
    # v2.0
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.4_T2T-CHM13v2.0/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz
    
Building database (all k-mers, < 10min):
    
    # v1.1: input=GCA_009914755.3_CHM13_T2T_v1.1_genomic.fna.gz
    
    # v2.0
    input=GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz
    
    # splitting human genome into 3072 chunks, each with ~1Mb
    # The regular expression '^(\w{3}_\d{9}\.\d+).+' is for extracting 'GCA_009914755.4' from the file name.
    kmcp compute $input -O human-chm13-k21-n3072 \
        --ref-name-regexp '^(\w{3}_\d{9}\.\d+).+' \
        -k 21 \
        --split-number 3072 --split-overlap 150 \
        --log human-chm13-k21-n3072.log -j 32 --force
    
    # using a small false positive rate: 0.05
    # using one hash function: 1
    kmcp index -I human-chm13-k21-n3072/ -O human-chm13.kmcp \
        -j 8 -f 0.05 -n 1 \
        --log human-chm13.kmcp.log --force
    
    # taxid.map
    echo -ne "GCA_009914755.4\t9606\n" > taxid.map
    
    # name.map
    echo -ne "GCA_009914755.4\tHomo sapiens isolate CHM13\n" > name.map
    
    # cp name mapping file to database directory
    cp taxid.map name.map human-chm13.kmcp/

    # database size
    $ du -sh human-chm13.kmcp
    6.7G    human-chm13.kmcp
    
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
    
    kmcp compute -I $input -O refseq-$name-k21-n10 \
        -k 21 --circular \
        --split-number 10 --split-overlap 150 \
        --log refseq-$name-k21-n10.log -j 32 --force
    
    #   using small false positive rate: 0.01
    #   using more hash functions: 3
    kmcp index -I refseq-$name-k21-n10/ -O refseq-$name.kmcp \
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


## Building databases (prokaryotic genome collections)

### HumGut (30,691 clusters)

[HumGut](http://arken.nmbu.no/~larssn/humgut/) is a comprehensive Human Gut prokaryotic genomes collection filtered by metagenome data.

>  In this work, we aimed to create a collection of the most prevalent healthy human gut prokaryotic genomes,
>  to be used as a reference database, including both MAGs from the human gut and ordinary RefSeq genomes.

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
    
    # assembly accession -> taxid
    cat HumGut.tsv \
        | csvtk cut -t -f genome_file,ncbi_tax_id \
        | csvtk replace -f genome_file -t -p '\.fna\.gz$' \
        | csvtk del-header \
        > taxid.map
        
    # assembly accession -> name
    cat HumGut.tsv \
        | csvtk cut -t -f genome_file,ncbi_organism_name \
        | csvtk replace -f genome_file -t -p '\.fna\.gz$' \
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
        | csvtk replace -f genome_file -t -p '\.fna\.gz$' \
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
    #   reference genomes are split into 10 chunks
    #   k = 21
    kmcp compute -I fna/ -k 21 -n 10 -B plasmid -O humgut-k21-n10 -j 32 --force

    # build database
    #   number of index files: 32, for server with >= 32 CPU cores
    #   bloom filter parameter:
    #     number of hash function: 1
    #     false positive rate: 0.3
    kmcp index -j 32 -I humgut-k21-n10 -O humgut.kmcp -n 1 -f 0.3
    
    cp taxid.map taxid-gtdb.map humgut.kmcp/

### proGenomes

[proGenomes](https://progenomes.embl.de/) v2.1 provides 84,096 consistently annotated bacterial
and archaeal genomes from over 12000 species.

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
    #   reference genomes are split into 10 chunks
    #   k = 21
    kmcp compute -I genomes/ -k 21 -n 10 -B plasmid -O progenomes-k21-n10 --force

    # build database
    #   number of index files: 32, for server with >= 32 CPU cores
    #   bloom filter parameter:
    #     number of hash function: 1
    #     false positive rate: 0.3
    kmcp index -j 32 -I progenomes-k21-n10 -O progenomes.kmcp -n 1 -f 0.3
    
    cp taxid.map name.map progenomes.kmcp/
        


## Building databases (viral genome collections)

### MGV

> Bacteriophages have important roles in the ecology of the human gut microbiome but are under-represented in reference data-
> bases. To address this problem, we assembled the Metagenomic Gut Virus catalogue that comprises 189,680 viral genomes
> from 11,810 publicly available human stool metagenomes. Over 75% of genomes represent double-stranded DNA phages that
> infect members of the Bacteroidia and Clostridia classes. Based on sequence clustering we identified 54,118 candidate viral spe-
> cies, 92% of which were not found in existing databases. The Metagenomic Gut Virus catalogue improves detection of viruses
> in stool metagenomes and accounts for nearly 40% of CRISPR spacers found in human gut Bacteria and Archaea. We also pro-
> duced a catalogue of 459,375 viral protein clusters to explore the functional potential of the gut virome. This revealed tens of
> thousands of diversity-generating retroelements, which use error-prone reverse transcription to mutate target genes and may
> be involved in the molecular arms race between phages and their bacterial hosts.
> 
> https://doi.org/10.1038/s41564-021-00928-6
> https://portal.nersc.gov/MGV/

Basic information (optional)

    $ seqkit  stats mgv_contigs.fna 
    file             format  type  num_seqs        sum_len  min_len   avg_len  max_len
    mgv_contigs.fna  FASTA   DNA    189,680  8,803,222,510    1,244  46,410.9  553,716

    # Genome completeness
    Complete: n=26,030
    >90% complete: n=53,220
    50-90% complete: n=110,430
    
    $ seqkit  seq -n mgv_contigs.fna | head -n 3
    MGV-GENOME-0364295
    MGV-GENOME-0364296
    MGV-GENOME-0364303

Stats (optional)
    
    cat mgv_contig_info.tsv \
        | csvtk cut -t -f completeness \
        | csvtk plot hist -o completeness.hist.png

    # Complete
    $ cat mgv_contig_info.tsv \
        | csvtk filter2 -t -f '$completeness == 100' \
        | csvtk nrows 
    32577
    
    # >90% complete
    $ cat mgv_contig_info.tsv \
        | csvtk filter2 -t -f '$completeness >= 90' \
        | csvtk nrows 
    78814  # < 79250
    
    # checkv_quality
    $ cat mgv_contig_info.tsv \
        | csvtk cut -t -f checkv_quality \
        | csvtk freq -t -nr | more
    checkv_quality  frequency
    Medium-quality  110430
    High-quality    53220
    Complete        26030
    
    # >90% complete && checkv_quality == High-quality/Complete
    $ cat mgv_contig_info.tsv \
        | csvtk filter2 -t -f '$completeness >= 90 && $checkv_quality != "Medium-quality"' \
        | csvtk nrows 
    78813
    
    
Stats of high-quality genomes (optional)

    # high-quality genomes
    cat mgv_contig_info.tsv \
        | csvtk filter2 -t -f '$completeness >= 90 && $checkv_quality != "Medium-quality"' \
        > mgv.hq.tsv
    
    # number of families
    cat mgv.hq.tsv \
        | csvtk freq -t -f ictv_family -nr \
        | csvtk nrow -t
    19
    
    # number of species
    cat mgv.hq.tsv \
        | csvtk freq -t -f votu_id -nr \
        | csvtk nrow -t
    26779
    
    # baltimore
    cat mgv.hq.tsv \
        | csvtk freq -t -f baltimore -nr \
        | csvtk pretty -t
    baltimore   frequency
    ---------   ---------
    dsDNA       76527
    ssDNA       2215
    NULL        62
    DNA         7
    dsDNA-RT    1
    ssRNA-RT    1
    
    # prophage?
    cat mgv.hq.tsv \
        | csvtk freq -t -f prophage -nr \
        | csvtk pretty -t
    prophage   frequency
    --------   ---------
    No         58366
    Yes        20447
    
    # species with both prophage and lytic phages
    cat mgv.hq.tsv \
        | csvtk freq -t -f votu_id,prophage \
        | csvtk freq -t -f votu_id \
        | csvtk filter2 -t -f '$frequency > 1' \
        | csvtk nrow -t
    2745


Prepare genomes:

    # ID
    cat mgv_contig_info.tsv \
        | csvtk filter2 -t -f '$completeness >= 90 && $checkv_quality != "Medium-quality"' \
        | csvtk cut -t -f contig_id \
        | csvtk del-header \
        > mgv.hq.tsv.id

    # extract sequences
    seqkit grep -f mgv.hq.tsv.id mgv_contigs.fna -o mgv.hq.fasta.gz
    
    # split into files with one genome
    seqkit split2 -s 1 mgv.hq.fasta.gz -O mgv
    
    # rename with 
    find mgv/ -name "*.fasta.gz" \
        | rush -j 20 'mv {} {/}/$(seqkit seq -ni {}).fa.gz'
        
Create taxdump files and `taxid.map` with taxonkit (version >= v0.12.0):

    cat mgv_contig_info.tsv \
            | csvtk cut -t -f ictv_order,ictv_family,ictv_genus,votu_id,contig_id \
        | csvtk del-header \
        > mgv.taxonomy.tsv

    taxonkit create-taxdump mgv.taxonomy.tsv --out-dir mgv-taxdump \
        --force -A 5 -R order,family,genus,species
    
    cp mgv-taxdump/taxid.map .
    
    # name.map
    cat mgv-taxdump/taxid.map \
        | taxonkit lineage --data-dir mgv-taxdump/ -i 2 \
        | cut -f 1,3 \
        > name.map
        
    head -n 5 name.map 
    MGV-GENOME-0364295      Caudovirales;crAss-phage;OTU-61123
    MGV-GENOME-0364296      Caudovirales;crAss-phage;OTU-61123
    MGV-GENOME-0364303      Caudovirales;crAss-phage;OTU-05782
    MGV-GENOME-0364311      Caudovirales;crAss-phage;OTU-01114
    MGV-GENOME-0364312      Caudovirales;crAss-phage;OTU-23935
    
Building database:
    
    # compute k-mers
    #   reference genomes are split into 10 chunks
    #   k = 21
    kmcp compute -I mgv/ -k 21 -n 10 -l 150 -O mgv-k21-n10 --force

    # viral genomes are small:
    #   using small false positive rate: 0.05
    #   still using one hash function: 1
    kmcp index -j 32 -I mgv-k21-n10 -O mgv.kmcp \
        -j 32 -f 0.05 -n 1 -x 100K -8 1M \
        --log mgv.kmcp.log
    
    cp taxid.map name.map mgv.kmcp/

## Building custom databases

Files:

1. Genome files
    - (Gzip-compressed) FASTA/Q format.
    - One genome per file **with the reference identifier in the file name**.
2. TaxId mapping file (for metagenomic profiling)
    - Two-column (reference identifier and TaxId) tab-delimited.
3. [NCBI taxonomy dump files](ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz) (for metagenomic profiling).
  You can use [taxonkit create-taxdump](https://bioinf.shenwei.me/taxonkit/usage/#create-taxdump)
  to create NCBI-style taxdump files for custom genome collections, which also generates a TaxId mapping file.
    - `names.dmp`
    - `nodes.dmp`
    - `merged.dmp` (optional)
    - `delnodes.dmp` (optional)

Tools:

- [kmcp](/kmcp/download) for metagenomic profiling.
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

**Input:**

1. Input plain or gzipped FASTA/Q files can be given via positional arguments or
    the flag `-i/--infile-list` with the list of input files,
2. Or a directory containing sequence files via the flag `-I/--in-dir`,
    with multiple-level sub-directories allowed. A regular expression
    for matching sequencing files is available via the flag `-r/--file-regexp`.
    The default pattern matches files with extension of 
    `fa`, `fasta`, `fna`, `fq`, `fastq`, `fa.gz`, `fasta.gz`,
    `fna.gz`, `fq.gz` or `fastq.gz`.


You may rename the sequence files for convenience using [brename](https://github.com/shenwei356/brename).
because the sequence/genome identifier in the index and search results would be:

1. For the default mode (computing k-mers for the whole file):
        
        the basename of file with common FASTA/Q file extension removed,
        captured via the flag -N/--ref-name-regexp.
        
2. For splitting sequence mode (see details below):
        
        same to 1).
        
3. For computing k-mers for each sequence:
        
        the sequence identifier.

Unwanted sequences like plasmid sequences can be **filtered out** by
the name via regular expression(s) (`-B/--seq-name-filter`).

**How to compute k-mers**:
By default, `kmcp` computes k-mers (sketches) of every file,
you can also use `--by-seq` to compute for every sequence,
where sequence IDs in all input files are better to be distinct.
It also **supports splitting sequences into chunks, this
could increase the specificity in profiling results at the cost
of a slower searching speed**. 

**Splitting sequences**:

1. Sequences can be split into chunks by a chunk size 
    (`-s/--split-size`) or number of chunks (`-n/--split-number`)
    with overlap (`-l/--split-overlap`).
    ***In this mode, the sequences of each genome should be saved in an
    individual file***.
2. When splitting by number of chunks, **all sequences (except for
    these matching any regular expression given by `-B/--seq-name-filter`)
    in a sequence file are concatenated with k-1 Ns before splitting**.
3. Both sequence/reference IDs and chunks indices are saved for later use,
    in form of meta/description data in `.unik` files, and will
    be reported in `kmcp search` results.

**Metadata**:

1. Every outputted `.unik` file contains the sequence/reference ID,
    chunk index, number of chunks, and genome size of reference.
2. When parsing whole sequence files or splitting by number of chunks,
    the **identifier of a reference** is the basename of the input file
    by default. It can also be extracted from the input file name via
    `-N/--ref-name-regexp`, e.g., `^(\w{3}_\d{9}\.\d+)` for **RefSeq assembly accessions**.

Multiple **sizes of k-mers** are supported, but a single k-mer size
is good enough. 

**Supported k-mer (sketches) types**:

1. K-mer:
    - ntHash of k-mer (`-k`)
2. K-mer sketchs (all using ntHash):
    - Scaled MinHash (`-k -D`), previously named Scaled MinHash
    - Minimizer      (`-k -W`), optionally scaling/down-sampling (`-D`)
    - Closed Syncmer (`-k -S`), optionally scaling/down-sampling (`-D`)

**Output**:

1. All outputted `.unik` files are saved in `${outdir}`, with path

        ${outdir}/xxx/yyy/zzz/${infile}-id_${seqID}.unik

    where dirctory tree `/xxx/yyy/zzz/` is built for > 1000 output files.

2. For splitting sequence mode (`--split-size > 0` or `--split-number > 0`),
    output files are:
    
        ${outdir}//xxx/yyy/zzz/${infile}/{seqID}-chunk_${chunkIdx}.unik
        
3. A summary file (`${outdir}/_info.txt`) is generated for later use.
     ***Users need to check if the reference IDs (column `name`) are what
     supposed to be***.

**Performance tips**:

1. Decrease value of `-j/--threads` for data in hard disk drives (HDD) to
    reduce I/O pressure.

**Commands** ([usage](/kmcp/usage/#compute)):

    # compute k-mers
    #   sequence containing "plasmid" in name are ignored,
    #   reference genomes are split into 10 chunks,
    #   k = 21
    kmcp compute --in-dir refs/ \
        --kmer 21 \
        --split-number 10 \
        --split-overlap 150 \
        --seq-name-filter plasmid \
        --ref-name-regexp '(.+).fasta.gz' \
        --out-dir refs-k21-n10

    # demo output
    13:11:13.397 [INFO] kmcp v0.9.0
    13:11:13.397 [INFO]   https://github.com/shenwei356/kmcp
    13:11:13.397 [INFO] 
    13:11:13.397 [INFO] checking input files ...
    13:11:13.398 [INFO]   9 input file(s) given
    13:11:13.398 [INFO] 
    13:11:13.398 [INFO] -------------------- [main parameters] --------------------
    13:11:13.398 [INFO] input and output:
    13:11:13.398 [INFO]   input directory: refs/
    13:11:13.398 [INFO]     regular expression of input files: (?i)\.(f[aq](st[aq])?|fna)(.gz)?$
    13:11:13.398 [INFO]     *regular expression for extracting reference name from file name: (?i)(.+).fasta.gz
    13:11:13.398 [INFO]     *regular expressions for filtering out sequences: [plasmid]
    13:11:13.398 [INFO]   output directory: refs-k21-n10
    13:11:13.398 [INFO] 
    13:11:13.398 [INFO] sequences splitting: true
    13:11:13.398 [INFO]   split parts: 10, overlap: 100 bp
    13:11:13.398 [INFO] 
    13:11:13.398 [INFO] k-mer (sketches) computing:
    13:11:13.398 [INFO]   k-mer size(s): 21
    13:11:13.398 [INFO]   circular genome: false
    13:11:13.398 [INFO]   saving exact number of k-mers: true
    13:11:13.398 [INFO] 
    13:11:13.398 [INFO] -------------------- [main parameters] --------------------
    13:11:13.398 [INFO] 
    13:11:13.398 [INFO] computing ...
    processed files:  9 / 9 [======================================] ETA: 0s. done
    13:11:13.845 [INFO] 
    13:11:13.845 [INFO] elapsed time: 453.870411ms
    13:11:13.845 [INFO]

A summary file (`_info.txt`) is generated for later use.
***Users need to check if the reference IDs (column `name`) are what supposed to be***.

|#path                                                                |name       |chunkIdx|idxNum|genomeSize|kmers |
|:--------------------------------------------------------------------|:----------|:------|:-----|:---------|:-----|
|refs-k21-n10/672/060/672/NC_010655.1.fasta.gz/NC_010655.1-chunk_0.unik|NC_010655.1|0      |10    |2664102   |264247|
|refs-k21-n10/672/060/672/NC_010655.1.fasta.gz/NC_010655.1-chunk_1.unik|NC_010655.1|1      |10    |2664102   |266237|
|refs-k21-n10/595/162/595/NC_012971.2.fasta.gz/NC_012971.2-chunk_0.unik|NC_012971.2|0      |10    |4558953   |450494|
|refs-k21-n10/292/039/292/NC_000913.3.fasta.gz/NC_000913.3-chunk_0.unik|NC_000913.3|0      |10    |4641652   |459277|
|refs-k21-n10/934/859/934/NC_013654.1.fasta.gz/NC_013654.1-chunk_0.unik|NC_013654.1|0      |10    |4717338   |470575|

Meta data in the `.unik` file can be showed using `kmcp utils unik-info`:

    kmcp utils unik-info refs-k21-n10/072/380/072/NZ_CP028116.1.fasta.gz/NZ_CP028116.1-chunk_0.unik -a


### Step 2. Building databases

KMCP builds index for k-mers (sketches) with a modified Compact Bit-sliced
Signature Index ([COBS](https://arxiv.org/abs/1905.09624)). 
We completely rewrite the algorithms, data structure, and file format,
and have improved the indexing and searching speed
(check [benchmark](/kmcp/benchmark/searching)).

**Input**:

- The output directory generated by `kmcp compute`.

**Database size and searching accuracy**:

0. **Use `--dry-run` to adjust parameters and check the final number of 
    index files (#index-files) and the total file size**.
1. `-f/--false-positive-rate`: the default value `0.3` is enough for a
    query with tens of k-mers (see BIGSI/COBS paper).
    Small values could largely increase the size of the database.
2. `-n/--num-hash`: large values could reduce the database size,
    at the cost of a slower searching speed. Values <=4 are recommended.
3. The value of block size `-b/--block-size` is better to be multiple of 64.
    The default value is: 
    
        (#unikFiles/#threads + 7) / 8 * 8

4. Use flag `-x/--block-sizeX-kmers-t`, `-8/--block-size8-kmers-t`,
    and `-1/--block-size1-kmers-t` to separately create indexes for
    inputs with huge number of k-mers, for precise control of
    database size.

**Taxonomy data**:

1. No taxonomy data are included in the database.
2. Taxonomy information are only needed in `kmcp profile`.

**Performance tips**:

1. The umber of blocks (`.uniki` files) is better to be smaller than or equal
   to the number of CPU cores for faster searching speed. 
   **We can set `-j/--threads` to control the blocks number**.
   When more threads (>= 1.3 * #blocks) are given, extra workers are
   automatically created.
2. `#threads` files are simultaneously opened, and the max number
    of opened files is limited by the flag `-F/--max-open-files`.
    You may use a small value of `-F/--max-open-files` for 
    hard disk drive storages.
3. When the database is used in a new computer with more CPU cores,
   `kmcp search` could automatically scale to utilize as many cores
   as possible.

**Commands** ([usage](/kmcp/usage/#compute)):

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
    22:44:17.710 [INFO] kmcp v0.9.0
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

Where `num-sigs` is the size of the bloom filters, and `num-names` is the number of genome (chunks).

What's next? Check the [tutorials](/kmcp/tutorial).

## Merging GTDB and NCBI taxonomy

By default, we use NCBI taxonomy database for reference genomes from GTDB and Refseq.
Alternately, we can also use the GTDB taxonomy.

The idea is to export lineages from both GTDB and NCBI using [taxonkit reformat](https://bioinf.shenwei.me/taxonkit/usage/#reformat), and then create taxdump files from them with [taxonkit create-taxdump](https://bioinf.shenwei.me/taxonkit/usage/#create-taxdump).

1. Mapping genome assembly accessions to GTDB lineages using existing [GTDB-taxdump](https://github.com/shenwei356/gtdb-taxdump),
   which provids GTDB taxonomy taxdump files for each GTDB version and a `taxid.map` file mapping accessions to TaxIds of taxa at the subspecies rank.
   Here, we only output lineages down to the species rank and leave name of taxa at the subspecies/strain rank empty.

        taxonkit reformat \
            --data-dir gtdb-taxdump/R214/ \
            --taxid-field 2 \
            --format "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t" \
            gtdb-taxdump/R214/taxid.map \
        | csvtk cut -H -t -f -2 -o gtdb-r214.tsv


1. Exporting taxonomic lineages of viral and fungal taxa with rank equal to or lower than species from NCBI taxdump.
   For taxa whose rank is "no rank" below the species, we treat them as tax of strain rank (`--pseudo-strain`, taxonkit v0.14.1 or later versionsneeded).

        # for viral data
        taxonkit reformat \
            --data-dir taxdump/ \
            --taxid-field 2 \
            --pseudo-strain \
            --format "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}" \
            taxid.map \
        | csvtk cut -H -t -f -2 -o genbank-viral.tsv

        # for fungal data
        taxonkit reformat \
            --data-dir taxdump/ \
            --taxid-field 2 \
            --pseudo-strain \
            --format "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}" \
            taxid.map \
        | csvtk cut -H -t -f -2 -o refseq-fungi.tsv


1. Creating taxdump from lineages above.

        # move the lineage files to current directory
        # cp gtdb/gtdb-r214.tsv genbank-viral/genbank-viral.tsv refseq-fungi/refseq-fungi.tsv .

        cat gtdb-r214.tsv genbank-viral.tsv refseq-fungi.tsv \
            | taxonkit create-taxdump \
                --field-accession 1 \
                -R "superkingdom,phylum,class,order,family,genus,species,strain" \
                -O taxdump.gtdb+ncbi

1. NCBI taxonomy of Viruses changed rapidly, of which some TaxIds might be deleted.
   You may added them to the the `taxid.map` file for compatibility.

        cat genbank-viral.tsv refseq-fungi.tsv gtdb-r214.tsv \
            | csvtk grep -Ht -f 2 -r -p "^$" \
            | cut -f 1 \
            | awk '{print $1"\t0"}' \
            >> taxdump.gtdb+ncbi/taxid.map

Some tests:

    # SARS-COV-2 in NCBI taxonomy
    $ echo 2697049 \
        | taxonkit lineage -t --data-dir taxdump/ \
        | csvtk cut -Ht -f 3 \
        | csvtk unfold -Ht -f 1 -s ";" \
        | taxonkit lineage -r -n -L --data-dir taxdump/ \
        | csvtk cut -Ht -f 1,3,2 \
        | csvtk pretty -Ht
    10239     superkingdom   Viruses
    2559587   clade          Riboviria
    2732396   kingdom        Orthornavirae
    2732408   phylum         Pisuviricota
    2732506   class          Pisoniviricetes
    76804     order          Nidovirales
    2499399   suborder       Cornidovirineae
    11118     family         Coronaviridae
    2501931   subfamily      Orthocoronavirinae
    694002    genus          Betacoronavirus
    2509511   subgenus       Sarbecovirus
    694009    species        Severe acute respiratory syndrome-related coronavirus
    2697049   no rank        Severe acute respiratory syndrome coronavirus 2

    $ echo "Severe acute respiratory syndrome coronavirus 2" | taxonkit name2taxid --data-dir taxdump.gtdb+ncbi/
    Severe acute respiratory syndrome coronavirus 2 216305222

    $ echo 216305222 \
        | taxonkit lineage -t --data-dir taxdump.gtdb+ncbi/ \
        | csvtk cut -Ht -f 3 \
        | csvtk unfold -Ht -f 1 -s ";" \
        | taxonkit lineage -r -n -L --data-dir taxdump.gtdb+ncbi/ \
        | csvtk cut -Ht -f 1,3,2 \
        | csvtk pretty -Ht
    1287770734   superkingdom   Viruses
    1506901452   phylum         Pisuviricota
    1091693597   class          Pisoniviricetes
    37745009     order          Nidovirales
    738421640    family         Coronaviridae
    906833049    genus          Betacoronavirus
    1015862491   species        Severe acute respiratory syndrome-related coronavirus
    216305222    strain         Severe acute respiratory syndrome coronavirus 2



    $ echo "Escherichia coli"  | taxonkit name2taxid --data-dir taxdump.gtdb+ncbi/
    Escherichia coli        1945799576

    $ echo 1945799576 \
        | taxonkit lineage -t --data-dir taxdump.gtdb+ncbi/ \
        | csvtk cut -Ht -f 3 \
        | csvtk unfold -Ht -f 1 -s ";" \
        | taxonkit lineage -r -n -L --data-dir taxdump.gtdb+ncbi/ \
        | csvtk cut -Ht -f 1,3,2 \
        | csvtk pretty -Ht
    609216830    superkingdom   Bacteria
    1641076285   phylum         Proteobacteria
    329474883    class          Gammaproteobacteria
    1012954932   order          Enterobacterales
    87250111     family         Enterobacteriaceae
    1187493883   genus          Escherichia
    1945799576   species        Escherichia coli
