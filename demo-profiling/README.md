# Demo of taxonomic profiling

## Data

### Reference genomes

We choose 15 bacterial genomes to make a mock metagenomic community.

Taxonomy information (NCBI Taxonomy):
    
    cat taxid.map \
        | taxonkit reformat -I 2 -f '{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}' \
        | csvtk cut -t -f -2 \
        | csvtk add-header -t -n id,superkingdom,phylum,class,order,family,genus,species \
        > taxonomy.tsv
    
    csvtk pretty -t taxonomy.tsv

    id                superkingdom   phylum           class                 order              family               genus            species
    ---------------   ------------   --------------   -------------------   ----------------   ------------------   --------------   --------------------------
    GCF_003697165.2   Bacteria       Proteobacteria   Gammaproteobacteria   Enterobacterales   Enterobacteriaceae   Escherichia      Escherichia coli
    GCF_002949675.1   Bacteria       Proteobacteria   Gammaproteobacteria   Enterobacterales   Enterobacteriaceae   Shigella         Shigella dysenteriae
    GCF_002950215.1   Bacteria       Proteobacteria   Gammaproteobacteria   Enterobacterales   Enterobacteriaceae   Shigella         Shigella flexneri
    GCF_000742135.1   Bacteria       Proteobacteria   Gammaproteobacteria   Enterobacterales   Enterobacteriaceae   Klebsiella       Klebsiella pneumoniae
    GCF_000006945.2   Bacteria       Proteobacteria   Gammaproteobacteria   Enterobacterales   Enterobacteriaceae   Salmonella       Salmonella enterica
    GCF_001544255.1   Bacteria       Firmicutes       Bacilli               Lactobacillales    Enterococcaceae      Enterococcus     Enterococcus faecium
    GCF_000392875.1   Bacteria       Firmicutes       Bacilli               Lactobacillales    Enterococcaceae      Enterococcus     Enterococcus faecalis
    GCF_001457655.1   Bacteria       Proteobacteria   Gammaproteobacteria   Pasteurellales     Pasteurellaceae      Haemophilus      Haemophilus influenzae
    GCF_900638025.1   Bacteria       Proteobacteria   Gammaproteobacteria   Pasteurellales     Pasteurellaceae      Haemophilus      Haemophilus parainfluenzae
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
|GCF_001457655.1.fa.gz|1       |1890645|Haemophilus influenzae genome assembly NCTC8143, chromosome : 1                   |
|GCF_900638025.1.fa.gz|1       |2062405|Haemophilus parainfluenzae strain NCTC10665 genome assembly, chromosome: 1        |
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

    taxonkit create-taxdump --field-accession 1 taxonomy.tsv -O taxdump-custom/


### Mock community 

Generating a mock dataset using https://github.com/iqbal-lab-org/simutator, 
every 2kb where a region of length 1.5kb gets 30 snps, 2 insertion and 4 deletions up to length 10bp.

Designed relative depths:

    cat depth.tsv \
        | csvtk mutate -Ht \
        | csvtk replace -Ht -f 3 -p '(.+)\.fa' -r '{kv}' -k name2.map 

    GCF_003697165.2.fa      1       Escherichia coli
    GCF_002949675.1.fa      1       Shigella dysenteriae
    GCF_002950215.1.fa      1       Shigella flexneri
    GCF_000742135.1.fa      1       Klebsiella pneumoniae
    GCF_000006945.2.fa      1       Salmonella enterica
    GCF_001544255.1.fa      0.1     Enterococcus faecium
    GCF_000392875.1.fa      0.1     Enterococcus faecalis
    GCF_001457655.1.fa      0.1     Haemophilus influenzae
    GCF_900638025.1.fa      0.1     Haemophilus parainfluenzae
    GCF_001027105.1.fa      0.05    Staphylococcus aureus
    GCF_006742205.1.fa      0.05    Staphylococcus epidermidis
    GCF_000148585.2.fa      0.01    Streptococcus mitis
    GCF_001096185.1.fa      0.01    Streptococcus pneumoniae
    GCF_000017205.1.fa      0.005   Pseudomonas aeruginosa
    GCF_009759685.1.fa      0.005   Acinetobacter baumannii

Steps:
    
    # tools:
    #   - https://github.com/shenwei356/rush
    #   - https://github.com/shenwei356/taxonkit

    # unzip all references which are required by 'simutator'
    mkdir -p mock
    ls refs/*.gz | rush 'seqkit seq {} -o mock/{%.}'
    cd mock
    
    # mutate
    ls *.fa | rush 'simutator mutate_fasta --seed 1 --complex 2000:1500:30:2:4:10 {} {}.m.fasta'
    rm *.vcf
    
    # simulate reads. since simutator only supports depth of integer,
    # we'll generate reads of depth 2 first and then perform downsampling.
    ls *m.fasta*.fa \
        | rush 'simutator simulate_reads {} {@^(.+?\.fa)} \
            --read_length 150 --read_depth 2 --fragment_length 350; \
            rm {}'

    # downsampling
    ls *.fa | grep -v m.fasta \
        | rush 'p=$(grep {%} ../depth.tsv | cut -f 2); \
            for f in {}.*.{1,2}.fq.gz; do seqkit sample -p $p $f -o $f.sample.fq.gz; done'

    # concatenate all sampled reads
    cd ..
    seqkit seq mock/*.1.fq.gz.sample.fq.gz -o mock_1.fastq.gz
    seqkit seq mock/*.2.fq.gz.sample.fq.gz -o mock_2.fastq.gz

    # ----------------------------------    
    # stats

    # stats of every reference
    ls mock/*.fa | grep -v m.fasta | rush 'seqkit seq {}.*.sample.fq.gz | seqkit stats -i {%.}  -T -o {}.stats -a'
    csvtk concat -t mock/*.stats | csvtk sort -t -k num_seqs:nr | csvtk pretty -t
    
    file              format   type   num_seqs   sum_len    min_len   avg_len   max_len   Q1     Q2      Q3     sum_gap   N50   Q20(%)   Q30(%)   GC(%)
    ---------------   ------   ----   --------   --------   -------   -------   -------   ----   -----   ----   -------   ---   ------   ------   -----
    GCF_000742135.1   FASTQ    DNA    73472      11020800   150       150.0     150       75.0   150.0   75.0   0         150   97.66    90.66    56.81
    GCF_003697165.2   FASTQ    DNA    66770      10015500   150       150.0     150       75.0   150.0   75.0   0         150   97.66    90.65    50.63
    GCF_000006945.2   FASTQ    DNA    65664      9849600    150       150.0     150       75.0   150.0   75.0   0         150   97.66    90.65    52.19
    GCF_002950215.1   FASTQ    DNA    65490      9823500    150       150.0     150       75.0   150.0   75.0   0         150   97.65    90.64    50.67
    GCF_002949675.1   FASTQ    DNA    60720      9108000    150       150.0     150       75.0   150.0   75.0   0         150   97.66    90.65    50.95
    GCF_000392875.1   FASTQ    DNA    3798       569700     150       150.0     150       75.0   150.0   75.0   0         150   97.65    90.66    37.96
    GCF_001544255.1   FASTQ    DNA    3294       494100     150       150.0     150       75.0   150.0   75.0   0         150   97.65    90.63    38.60
    GCF_900638025.1   FASTQ    DNA    2744       411600     150       150.0     150       75.0   150.0   75.0   0         150   97.64    90.60    39.70
    GCF_001457655.1   FASTQ    DNA    2492       373800     150       150.0     150       75.0   150.0   75.0   0         150   97.65    90.60    38.57
    GCF_001027105.1   FASTQ    DNA    1814       272100     150       150.0     150       75.0   150.0   75.0   0         150   97.66    90.57    33.61
    GCF_006742205.1   FASTQ    DNA    1594       239100     150       150.0     150       75.0   150.0   75.0   0         150   97.64    90.61    32.85
    GCF_000017205.1   FASTQ    DNA    434        65100      150       150.0     150       75.0   150.0   75.0   0         150   97.66    90.68    65.97
    GCF_001096185.1   FASTQ    DNA    284        42600      150       150.0     150       75.0   150.0   75.0   0         150   97.65    90.72    40.44
    GCF_009759685.1   FASTQ    DNA    272        40800      150       150.0     150       75.0   150.0   75.0   0         150   97.65    90.67    39.23
    GCF_000148585.2   FASTQ    DNA    242        36300      150       150.0     150       75.0   150.0   75.0   0         150   97.65    90.57    40.33
      
    # stats
    seqkit stats -a  mock*.fastq.gz
    file             format  type  num_seqs     sum_len  min_len  avg_len  max_len  Q1   Q2  Q3  sum_gap  N50  Q20(%)  Q30(%)  GC(%)
    mock_1.fastq.gz  FASTQ   DNA    174,542  26,181,300      150      150      150  75  150  75        0  150   98.01   91.67  51.69
    mock_2.fastq.gz  FASTQ   DNA    174,542  26,181,300      150      150      150  75  150  75        0  150   97.31   89.63  51.69

The gold standard (gound-truth) of taxonomic abundance:

    # depth of every genome
    csvtk concat -t mock/*.stats \
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
    10.994104
    
    # relative abundance (percentage)    
    csvtk join -t mock.depth.tsv <(csvtk add-header -t -n id,taxid taxdump-custom/taxid.map) \
        | csvtk mutate2 -t -n abundance -e "\$depth/$d*100" -w 6 \
        | tee mock.gs.tsv \
        | csvtk csv2md -t
    
|id             |species                   |reads|bases   |genome_size|depth   |taxid     |abundance|
|:--------------|:-------------------------|:----|:-------|:----------|:-------|:---------|:--------|
|GCF_002949675.1|Shigella dysenteriae      |60720|9108000 |4578459    |1.989316|524994882 |18.094390|
|GCF_000006945.2|Salmonella enterica       |65664|9849600 |4951383    |1.989262|1678121664|18.093898|
|GCF_002950215.1|Shigella flexneri         |65490|9823500 |4938295    |1.989249|2695851945|18.093780|
|GCF_003697165.2|Escherichia coli          |66770|10015500|5034834    |1.989241|4093283224|18.093707|
|GCF_000742135.1|Klebsiella pneumoniae     |73472|11020800|5545784    |1.987239|3958205156|18.075498|
|GCF_900638025.1|Haemophilus parainfluenzae|2744 |411600  |2062405    |0.199573|1063930303|1.815273 |
|GCF_001544255.1|Enterococcus faecium      |3294 |494100  |2484851    |0.198845|4145431389|1.808651 |
|GCF_000392875.1|Enterococcus faecalis     |3798 |569700  |2881400    |0.197716|3809813362|1.798382 |
|GCF_001457655.1|Haemophilus influenzae    |2492 |373800  |1890645    |0.197710|328800344 |1.798328 |
|GCF_006742205.1|Staphylococcus epidermidis|1594 |239100  |2427041    |0.098515|1920251658|0.896071 |
|GCF_001027105.1|Staphylococcus aureus     |1814 |272100  |2782562    |0.097788|1569132721|0.889459 |
|GCF_001096185.1|Streptococcus pneumoniae  |284  |42600   |2117177    |0.020121|2983929374|0.183016 |
|GCF_000148585.2|Streptococcus mitis       |242  |36300   |1868883    |0.019423|1527235303|0.176667 |
|GCF_009759685.1|Acinetobacter baumannii   |272  |40800   |3990388    |0.010225|72054943  |0.093004 |
|GCF_000017205.1|Pseudomonas aeruginosa    |434  |65100   |6588339    |0.009881|3843752343|0.089875 |

    # CAMI format
    csvtk del-header mock.gs.tsv \
        | taxonkit profile2cami -i 7 -a 8 --data-dir taxdump-custom/ -s 0 \
        > mock.gs.profile
    
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
    
    # clean tmp files
    rm -rf refs-k21-n10

### Searching and profiling

Searching:
    
    # paired information are not used.
    kmcp search \
        --db-dir refs-k21-n10.kmcp/ \
        mock_1.fastq.gz \
        mock_2.fastq.gz \
        --out-file mock.kmcp.gz \
        --log mock.kmcp.gz.log
    
    # matched reads
    grep "queries matched" mock.kmcp.gz.log
    10:18:31.492 [INFO] 88.4713% (308839/349084) queries matched

Profiling using mode 1 for low coverage data:

    for f in *.kmcp.gz; do
        kmcp profile \
            --taxid-map taxdump-custom/taxid.map \
            --taxdump taxdump-custom/ \
            $f \
            --mode 1 \
            --out-prefix $f.kmcp.profile \
            --metaphlan-report $f.metaphlan.profile \
            --sample-id 0 \
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
    |GCF_003697165.2|Escherichia coli          |18.093707   |18.663804 |
    |GCF_002949675.1|Shigella dysenteriae      |18.094390   |18.201855 |
    |GCF_000006945.2|Salmonella enterica       |18.093898   |18.143627 |
    |GCF_000742135.1|Klebsiella pneumoniae     |18.075498   |17.738253 |
    |GCF_002950215.1|Shigella flexneri         |18.093780   |17.728060 |
    |GCF_900638025.1|Haemophilus parainfluenzae|1.815273    |1.809292  |
    |GCF_000392875.1|Enterococcus faecalis     |1.798382    |1.800250  |
    |GCF_001544255.1|Enterococcus faecium      |1.808651    |1.795723  |
    |GCF_001457655.1|Haemophilus influenzae    |1.798328    |1.787560  |
    |GCF_006742205.1|Staphylococcus epidermidis|0.896071    |0.906778  |
    |GCF_001027105.1|Staphylococcus aureus     |0.889459    |0.881014  |
    |GCF_000148585.2|Streptococcus mitis       |0.176667    |0.178996  |
    |GCF_001096185.1|Streptococcus pneumoniae  |0.183016    |0.177453  |
    |GCF_009759685.1|Acinetobacter baumannii   |0.093004    |0.098158  |
    |GCF_000017205.1|Pseudomonas aeruginosa    |0.089875    |0.089177  |
    
Assessing with [OPAL](https://github.com/CAMI-challenge/OPAL):

    opal.py -g mock.gs.profile mock.kmcp.gz.cami.profile -l KMCP -o opal/
