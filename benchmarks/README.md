# Benchmarks

- [Benchmarks on 64 simulated mouse gut short-read datasets from CAMI2 challenge](cami2-mouse-gut)
- [Benchmarks on 25 simulated prokaryotic communities from Sun et al](sim-bact-sun2021)
- [Benchmarks on 16 mock virome communities from Roux et al](mock-virome-roux2016)
- [Benchmarks on 87 metagenomic samples of infected body fluids](real-pathogen-gu2020)
- [Benchmarks on the ZymoBIOMICS Gut Microbiome Standard D6331](mock-hifi-zymo)

## Custom database for Kraken and Bracken

Preparing genomes

    # link to genomes
    $ ls
    genbank-viral
    gtdb
    refseq-fungi
    taxid.genbank-viral.map
    taxid.gtdb.map
    taxid.refseq-fungi.map
    
   
Formating FASTA ID for Kraken.

    time for db in refseq-fungi genbank-viral gtdb; do    
        outdir=$db.formated
        mkdir -p $outdir
        
        find $db/ -name "*.fna.gz" \
            | rush -v outdir=$outdir -v map=taxid.$db.map \
                'taxid=$(grep {%..} {map} | cut -f 2); \
                seqkit seq -i {} | seqkit replace -p "\$" -r "|kraken:taxid|$taxid" -o {outdir}/{%.}'
    done

Building Kraken database. https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown

    dbname=kmcp
    
    # link to taxonomy dump files
    ln -s ~/ws/data/kmcp/genbank-viral/taxdump $dbname/taxonomy
            
    
    memusg -t -s "find *.formated -name '*.fna' \
        | rush -j 40 -v db=$dbname \
            'kraken2-build --add-to-library {} --db {db}' " > kraken.add.log 2>&1
            
    memusg -t -s "kraken2-build --build --db $dbname --threads 40" > kraken.build.log 2>&1
 
Building Bracken database. https://github.com/jenniferlu717/Bracken, recommend https://ccb.jhu.edu/software/bracken/index.shtml?t=manual#step1 .

    find *.formated -name '*.fna' | rush -j 1 'cat {}' > all.fasta
    
    # very important, necessary
    export KRAKEN2_DB_PATH=$(pwd)
    
    dbname=kmcp
    
    memusg -t -s "kraken2 --db=${dbname} --threads=40 all.fasta > ${dbname}/database.kraken" \
        > bracken.build_1a.log 2>&1
    
    for READ_LEN in 50 75 100 125 150; do
        memusg -t -s "kmer2read_distr --seqid2taxid ${dbname}/seqid2taxid.map --taxonomy ${dbname}/taxonomy \
            --kraken ${dbname}/database.kraken --output ${dbname}/database${READ_LEN}mers.kraken -k 35 -l ${READ_LEN} -t 40; \
            generate_kmer_distribution.py -i ${dbname}/database${READ_LEN}mers.kraken -o ${dbname}/database${READ_LEN}mers.kmer_distrib" \
            > bracken.build_1bc.${READ_LEN}.log 2>&1
    done

The path of database is:

    ~/ws/db/kraken/kmcp/kmcp
    
## Custom databases for Centrifuge

https://github.com/DaehwanKimLab/centrifuge/blob/master/MANUAL.markdown

We use the files above to create centrifuge index.

    # reformat sequence header
    seqkit replace -p '\|.+' all.fasta -o all.cf.fasta
    
    # reformat seq2taxid.map
    csvtk replace -Ht -p '\|.+' kmcp/seqid2taxid.map -o seqid2taxid.cf.map
    
    # elapsed time: 1.0days 20h:35m:33s
    # peak rss: 522.41 GB
    memusg -t -s "centrifuge-build --threads 40 --conversion-table seqid2taxid.cf.map \
        --taxonomy-tree kmcp/taxonomy/nodes.dmp \
        --name-table kmcp/taxonomy/names.dmp \
        all.cf.fasta \
        kmcp " > centrifuge.build.log 2>&1

    # generated files
    -rw-rw-r-- 1 shenwei shenwei  56G Jan  1 15:51 kmcp.1.cf
    -rw-rw-r-- 1 shenwei shenwei  42G Jan  1 15:50 kmcp.2.cf
    -rw-rw-r-- 1 shenwei shenwei 147M Dec 30 19:44 kmcp.3.cf
    -rw-rw-r-- 1 shenwei shenwei  70M Jan  1 15:50 kmcp.4.cf

    # move out of the kraken directory:
    mkdir -p ~/ws/db/centrifuge
    mv kmcp.*.cf ~/ws/db/centrifuge
    
The path of database is:

    ~/ws/db/centrifuge/kmcp

## Custom database for Metamaps

https://github.com/DiltheyLab/MetaMaps/issues/5#issuecomment-414301437

Preparing files:

    genbank-viral/  gtdb/  refseq-fungi/  taxdump/ taxid.map
    
    $ tree refseq-fungi/ | head -n 3
    refseq-fungi/
    ├── GCF_000001985.1.fna
    ├── GCF_000002495.2.fna
    
File list

    fd fn?a$ genbank-viral/ gtdb/ refseq-fungi/ > files.tsv
    
    csvtk join -Ht -f '2;1' <(csvtk mutate -Ht -p '/(\w{3}_\d{9}\.\d+)' files.tsv) taxid.map \
        | csvtk cut -Ht -f 3,1 \
        > taxid2file.tsv
        
Build it!

    memusg -t -s "perl $(which combineAndAnnotateReferences.pl) \
        --inputFileList taxid2file.tsv \
        --outputFile kmcp.fa \
        --taxonomyInDirectory taxdump/ \
        --taxonomyOutDirectory kmcp-taxdump " > metamaps.a.log 2>&1
        
        
     memusg -t -s "perl $(which buildDB.pl) --DB kmcp \
        --FASTAs kmcp.fa \
        --taxonomy kmcp-taxdump " > metamaps.b.log 2>&1
