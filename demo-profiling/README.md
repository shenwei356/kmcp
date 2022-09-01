# Demo of taxonomic profiling

## Data

### Reference genomes

We choose 13 bacterial genomes to make mock metagenomic communities.

Taxonomy information:

    id                superkingdom   phylum           class                 order              family               genus            species
    ---------------   ------------   --------------   -------------------   ----------------   ------------------   --------------   --------------------------
    GCF_003697165.2   Bacteria       Proteobacteria   Gammaproteobacteria   Enterobacterales   Enterobacteriaceae   Escherichia      Escherichia coli
    GCF_002949675.1   Bacteria       Proteobacteria   Gammaproteobacteria   Enterobacterales   Enterobacteriaceae   Shigella         Shigella dysenteriae
    GCF_002950215.1   Bacteria       Proteobacteria   Gammaproteobacteria   Enterobacterales   Enterobacteriaceae   Shigella         Shigella flexneri
    GCF_000742135.1   Bacteria       Proteobacteria   Gammaproteobacteria   Enterobacterales   Enterobacteriaceae   Klebsiella       Klebsiella pneumoniae
    GCF_000006945.2   Bacteria       Proteobacteria   Gammaproteobacteria   Enterobacterales   Enterobacteriaceae   Salmonella       Salmonella enterica
    
    GCF_001544255.1   Bacteria       Firmicutes       Bacilli               Lactobacillales    Enterococcaceae      Enterococcus     Enterococcus faecium
    GCF_000392875.1   Bacteria       Firmicutes       Bacilli               Lactobacillales    Enterococcaceae      Enterococcus     Enterococcus faecalis
    
    GCF_001027105.1   Bacteria       Firmicutes       Bacilli               Bacillales         Staphylococcaceae    Staphylococcus   Staphylococcus aureus
    GCF_006742205.1   Bacteria       Firmicutes       Bacilli               Bacillales         Staphylococcaceae    Staphylococcus   Staphylococcus epidermidis
    
    GCF_001096185.1   Bacteria       Firmicutes       Bacilli               Lactobacillales    Streptococcaceae     Streptococcus    Streptococcus pneumoniae
    GCF_000148585.2   Bacteria       Firmicutes       Bacilli               Lactobacillales    Streptococcaceae     Streptococcus    Streptococcus mitis
    
    GCF_009759685.1   Bacteria       Proteobacteria   Gammaproteobacteria   Moraxellales       Moraxellaceae        Acinetobacter    Acinetobacter baumannii
    
    GCF_000017205.1   Bacteria       Proteobacteria   Gammaproteobacteria   Pseudomonadales    Pseudomonadaceae     Pseudomonas      Pseudomonas aeruginosa

Genome details:

    # name.map
    ls refs/*.gz  | rush 'echo -ne "{%..}\t"; seqkit head -n 1 {} | seqkit seq -n | cut -d " " -f 2-' > name.map
    
    # genome_size.tsv
    seqkit stats -T refs/*.gz \
        | csvtk cut -t -f file,sum_len \
        | csvtk rename -t -f 1,2 -n id,genome_size \
        | csvtk replace -p '.+/(.+)\.fa\.gz' -r '$1' \
        > genome_size.tsv

    
    csvtk join -t -f id \
            <(seqkit stats -j 10 refs/*.gz -T -b \
                | csvtk mutate -t -n id -p "(.+)\.fa") \
            <(csvtk add-header -t -n id,name name.map) \
        | csvtk cut -t -f file,num_seqs,sum_len,name \
        | csvtk sort -t -k name \
        | csvtk csv2md -t
        
|file                 |num_seqs|sum_len|name                                                                              |
|:--------------------|:-------|:------|:---------------------------------------------------------------------------------|
|GCF_009759685.1.fa.gz|2       |3990388|Acinetobacter baumannii strain ATCC 19606 chromosome, complete genome             |
|GCF_000392875.1.fa.gz|3       |2881400|Enterococcus faecalis ATCC 19433 acAqW-supercont1.1, whole genome shotgun sequence|
|GCF_001544255.1.fa.gz|38      |2484851|Enterococcus faecium NBRC 100486, whole genome shotgun sequence                   |
|GCF_003697165.2.fa.gz|2       |5034834|Escherichia coli DSM 30083 = JCM 1649 = ATCC 11775 chromosome, complete genome    |
|GCF_000742135.1.fa.gz|5       |5545784|Klebsiella pneumoniae strain ATCC 13883 scaffold1, whole genome shotgun sequence  |
|GCF_000017205.1.fa.gz|1       |6588339|Pseudomonas aeruginosa PA7, complete genome                                       |
|GCF_000006945.2.fa.gz|2       |4951383|Salmonella enterica subsp. enterica serovar Typhimurium str. LT2, complete genome |
|GCF_002949675.1.fa.gz|2       |4578459|Shigella dysenteriae strain ATCC 13313 chromosome, complete genome                |
|GCF_002950215.1.fa.gz|3       |4938295|Shigella flexneri 2a strain ATCC 29903 chromosome, complete genome                |
|GCF_001027105.1.fa.gz|2       |2782562|Staphylococcus aureus subsp. aureus DSM 20231 chromosome, complete genome         |
|GCF_006742205.1.fa.gz|2       |2427041|Staphylococcus epidermidis NBRC 100911 DNA, complete genome                       |
|GCF_000148585.2.fa.gz|1       |1868883|Streptococcus mitis NCTC 12261 chromosome, complete genome                        |
|GCF_001096185.1.fa.gz|24      |2117177|Streptococcus pneumoniae strain SMRU824, whole genome shotgun sequence            |

### Taxonomy data

Please download and uncompress [taxdump.tar.gz](ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz),
and then copy `names.dmp`, `nodes.dmp`, `delnodes.dmp` and `merged.dmp` to directory `taxdump`.

Or create custom taxdump files with `taxonomy.tsv` using [taxonkit create-taxdump](https://bioinf.shenwei.me/taxonkit/usage/#create-taxdump):

    taxonkit create-taxdump -A 1 taxonomy.tsv -O taxdump-custom/


### Mock community 

Generating a mock dataset using https://github.com/iqbal-lab-org/simutator, 
every 2kb where a region of length 1.5kb gets 30 snps, 2 insertion and 4 deletions up to length 10bp.

Designed relative depths:

    $ cat depth.tsv 
    GCF_003697165.2.fa      1
    GCF_002949675.1.fa      1
    GCF_002950215.1.fa      1
    GCF_000742135.1.fa      1
    GCF_000006945.2.fa      1
    GCF_000392875.1.fa      0.1
    GCF_001544255.1.fa      0.1
    GCF_001027105.1.fa      0.05
    GCF_006742205.1.fa      0.05
    GCF_000148585.2.fa      0.01
    GCF_001096185.1.fa      0.01
    GCF_000017205.1.fa      0.005
    GCF_009759685.1.fa      0.005

Steps:
    
    # tools:
    #   - https://github.com/shenwei356/rush

    # unzip all references which are required by 'simutator'
    mkdir -p refs.plain
    ls refs/*.gz | rush 'seqkit seq {} -o refs.plain/{%.}'
    cd refs.plain
    
    # mutate
    ls *.fa | rush 'simutator mutate_fasta --seed 1 --complex 2000:1500:30:2:4:10 {} {}.m.fasta'
    rm *.vcf
    
    # simulate reads. since simutator only supports depth of integer,
    # we'll generate reads of depth 2 first and then perform downsampling.
    ls *m.fasta*.fa | rush 'simutator simulate_reads {} {@^(.+?\.fa)} --read_length 150 --read_depth 2 --fragment_length 350'

    # downsampling
    ls *.fa | grep -v m.fasta \
        | rush 'p=$(grep {%} ../depth.tsv | cut -f 2); \
            for f in {}.*.{1,2}.fq.gz; do seqkit sample -p $p $f -o $f.sample.fq.gz; done'

    # concatenate all sampled reads
    cd ..
    seqkit seq refs.plain/*.1.fq.gz.sample.fq.gz -o mock_1.fastq.gz
    seqkit seq refs.plain/*.2.fq.gz.sample.fq.gz -o mock_2.fastq.gz

    # ----------------------------------    
    # stats

    # stats of every reference
    ls *.fa | grep -v m.fasta | rush 'seqkit seq {}.*.sample.fq.gz | seqkit stats -i {%.}  -T -o {}.stats -a'
    csvtk concat -t *.stats | csvtk sort -t -k num_seqs:nr | csvtk pretty -t
    file                 format   type   num_seqs   sum_len    min_len   avg_len   max_len   Q1     Q2      Q3     sum_gap   N50   Q20(%)   Q30(%)   GC(%)
    ------------------   ------   ----   --------   --------   -------   -------   -------   ----   -----   ----   -------   ---   ------   ------   -----
    GCF_000742135.1.fa   FASTQ    DNA    73496      11024400   150       150.0     150       75.0   150.0   75.0   0         150   97.66    90.66    56.82
    GCF_003697165.2.fa   FASTQ    DNA    66770      10015500   150       150.0     150       75.0   150.0   75.0   0         150   97.65    90.65    50.58
    GCF_000006945.2.fa   FASTQ    DNA    65664      9849600    150       150.0     150       75.0   150.0   75.0   0         150   97.66    90.65    52.24
    GCF_002950215.1.fa   FASTQ    DNA    65490      9823500    150       150.0     150       75.0   150.0   75.0   0         150   97.66    90.65    50.62
    GCF_002949675.1.fa   FASTQ    DNA    60720      9108000    150       150.0     150       75.0   150.0   75.0   0         150   97.66    90.66    50.93
    GCF_000392875.1.fa   FASTQ    DNA    3796       569400     150       150.0     150       75.0   150.0   75.0   0         150   97.64    90.64    37.96
    GCF_001544255.1.fa   FASTQ    DNA    3294       494100     150       150.0     150       75.0   150.0   75.0   0         150   97.65    90.62    38.39
    GCF_001027105.1.fa   FASTQ    DNA    1814       272100     150       150.0     150       75.0   150.0   75.0   0         150   97.65    90.68    33.31
    GCF_006742205.1.fa   FASTQ    DNA    1594       239100     150       150.0     150       75.0   150.0   75.0   0         150   97.66    90.66    32.52
    GCF_000017205.1.fa   FASTQ    DNA    434        65100      150       150.0     150       75.0   150.0   75.0   0         150   97.66    90.71    66.12
    GCF_001096185.1.fa   FASTQ    DNA    284        42600      150       150.0     150       75.0   150.0   75.0   0         150   97.69    90.71    40.32
    GCF_009759685.1.fa   FASTQ    DNA    272        40800      150       150.0     150       75.0   150.0   75.0   0         150   97.71    90.95    39.42
    GCF_000148585.2.fa   FASTQ    DNA    242        36300      150       150.0     150       75.0   150.0   75.0   0         150   97.66    90.64    40.18
      
    # stats
    seqkit stats -a  mock*.fastq.gz
    file              format  type  num_seqs     sum_len  min_len  avg_len  max_len  Q1   Q2  Q3  sum_gap  N50  Q20(%)  Q30(%)  GC(%)
    mock_1.fastq.gz  FASTQ   DNA    171,935  25,790,250      150      150      150  75  150  75        0  150      98   91.66  51.87
    mock_2.fastq.gz  FASTQ   DNA    171,935  25,790,250      150      150      150  75  150  75        0  150   97.31   89.65  51.86

The gold standard (gound-truth) of taxonomic abundance:

    # depth of every genome
    csvtk concat -t refs.plain/*.stats \
        | csvtk cut -t -f file,num_seqs,sum_len \
        | csvtk join -t - name2.map genome_size.tsv \
        | csvtk rename -t -f 1-3 -n 'id,reads,bases' \
        | csvtk mutate2 -t -n depth -e '$bases/${genome_size}' -w 6 \
        | csvtk sort -t -k depth:nr \
        | csvtk cut -t -f id,species,reads,bases,genome_size,depth \
        > mock.depth.tsv
    
    # total coverage
    d=$(csvtk summary -t -f depth:sum mock.depth.tsv -w 6 | csvtk del-header)
    echo $d
    10.597366
    
    # relative abundance (percentage)    
    csvtk -t mutate2 -t -n abundance -e "\$depth/$d*100" -w 6 mock.depth.tsv \
        | tee mock.gs.tsv \
        | csvtk csv2md -t
    
|id             |species                   |reads|bases   |genome_size|depth   |abundance|
|:--------------|:-------------------------|:----|:-------|:----------|:-------|:--------|
|GCF_002949675.1|Shigella dysenteriae      |60720|9108000 |4578459    |1.989316|18.771797|
|GCF_000006945.2|Salmonella enterica       |65664|9849600 |4951383    |1.989262|18.771287|
|GCF_002950215.1|Shigella flexneri         |65490|9823500 |4938295    |1.989249|18.771164|
|GCF_003697165.2|Escherichia coli          |66770|10015500|5034834    |1.989241|18.771089|
|GCF_000742135.1|Klebsiella pneumoniae     |73496|11024400|5545784    |1.987888|18.758322|
|GCF_001544255.1|Enterococcus faecium      |3294 |494100  |2484851    |0.198845|1.876362 |
|GCF_000392875.1|Enterococcus faecalis     |3796 |569400  |2881400    |0.197612|1.864728 |
|GCF_006742205.1|Staphylococcus epidermidis|1594 |239100  |2427041    |0.098515|0.929618 |
|GCF_001027105.1|Staphylococcus aureus     |1814 |272100  |2782562    |0.097788|0.922758 |
|GCF_001096185.1|Streptococcus pneumoniae  |284  |42600   |2117177    |0.020121|0.189868 |
|GCF_000148585.2|Streptococcus mitis       |242  |36300   |1868883    |0.019423|0.183281 |
|GCF_009759685.1|Acinetobacter baumannii   |272  |40800   |3990388    |0.010225|0.096486 |
|GCF_000017205.1|Pseudomonas aeruginosa    |434  |65100   |6588339    |0.009881|0.093240 |
    
## Metagenomic profiling

### Building database

    # computing k-mers
    kmcp compute \
        --in-dir refs/ \
        --ref-name-regexp "^([\w\.\_]+\.\d+)" \
        --seq-name-filter "plasmid" \
        --kmer 21 \
        --split-number 10 \
        --split-overlap 150 \
        --out-dir refs-k21-n10 \
        --force

    # indexing k-mers
    kmcp index \
        --in-dir refs-k21-n10/\
        --num-hash 1 \
        --false-positive-rate 0.3 \
        --out-dir refs-k21-n10.kmcp \
        --force

### Searching and profiling

Searching:
    
    # paired information are not used.
    kmcp search \
        --db-dir refs-k21-n10.kmcp/ \
        --min-query-cov 0.55 \
         mock_1.fastq.gz mock_2.fastq.gz \
        --out-file mock.kmcp.gz \
        --log mock.kmcp.gz.log
    
    # matched reads
    grep "queries matched" mock.kmcp.gz.log
    11:25:33.320 [INFO] 89.4297% (307522/343870) queries matched

Profiling using mode 1 for low coverage data:

    for f in *.kmcp.gz; do
        kmcp profile \
            --taxid-map taxdump-custom/taxid.map \
            --taxdump taxdump-custom/ \
            $f \
            --mode 1 \
            --out-prefix $f.kmcp.profile \
            --metaphlan-report $f.metaphlan.profile \
            --cami-report $f.cami.profile \
            --binning-result $f.binning.gz \
            --log $f.kmcp.profile.log
    done
    

Profiling results:

    cat mock.kmcp.gz.kmcp.profile \
        | csvtk cut -t -f ref,reads,percentage,taxname \
        | csvtk join -t - <(csvtk cut -t -f id,abundance mock.gs.tsv) \
        | csvtk rename -t -f abundance -n ground_truth \
        | csvtk cut -t -f ref,taxname,ground_truth,percentage \
        | csvtk csv2md -t

|ref            |taxname                   |ground_truth|percentage|
|:--------------|:-------------------------|:-----------|:---------|
|GCF_002949675.1|Shigella dysenteriae      |18.771797   |19.403161 |
|GCF_000742135.1|Klebsiella pneumoniae     |18.758322   |18.921227 |
|GCF_002950215.1|Shigella flexneri         |18.771164   |18.813278 |
|GCF_000006945.2|Salmonella enterica       |18.771287   |18.679621 |
|GCF_003697165.2|Escherichia coli          |18.771089   |18.095846 |
|GCF_001544255.1|Enterococcus faecium      |1.876362    |1.884776  |
|GCF_000392875.1|Enterococcus faecalis     |1.864728    |1.809396  |
|GCF_006742205.1|Staphylococcus epidermidis|0.929618    |0.923508  |
|GCF_001027105.1|Staphylococcus aureus     |0.922758    |0.907424  |
|GCF_000148585.2|Streptococcus mitis       |0.183281    |0.187531  |
|GCF_001096185.1|Streptococcus pneumoniae  |0.189868    |0.184554  |
|GCF_009759685.1|Acinetobacter baumannii   |0.096486    |0.096786  |
|GCF_000017205.1|Pseudomonas aeruginosa    |0.093240    |0.092891  |
