# Benchmarks

- [Benchmarks on 64 simulated mouse gut short-read datasets from CAMI2 challenge](cami2-mouse-gut)
- [Benchmarks on 25 simulated bacterial communities from Sun et al](sim-bact-sun2021)
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
    
    for READ_LEN in 75 100 125 150; do
        memusg -t -s "kmer2read_distr --seqid2taxid ${dbname}/seqid2taxid.map --taxonomy ${dbname}/taxonomy \
            --kraken ${dbname}/database.kraken --output ${dbname}/database${READ_LEN}mers.kraken -k 35 -l ${READ_LEN} -t 40; \
            generate_kmer_distribution.py -i ${dbname}/database${READ_LEN}mers.kraken -o ${dbname}/database${READ_LEN}mers.kmer_distrib" \
            > bracken.build_1bc.${READ_LEN}.log 2>&1
    done

## Custom databases fro Centrifuge

 https://github.com/DaehwanKimLab/centrifuge/blob/master/MANUAL.markdown


    
    memusg -t -s "centrifuge-build --conversion-table kmcp/seqid2taxid.map \
        --taxonomy-tree kmcp/taxonomy/nodes.dmp \
        --name-table kmcp/taxonomy/names.dmp \
        all.fasta \
        kmcp.cf " > centrifuge.build.log 2>&1
