# Detecting contaminated sequences in a bacteria assembly

## Changelog

Tracked changes in [github](https://github.com/shenwei356/kmcp/commits/main/docs/tutorial/detecting-contaminated-seqs/index.md)

- 2023-07-20
    - Recommend using masked GTDB index files.
    - Recommend using GTDB taxdump files rather than NCBI taxdump files.
    - Adjust profiling prameters.
- 2023-07-16
    - First public version.

## Tools

- kmcp: https://github.com/shenwei356/kmcp
- seqkit: https://github.com/shenwei356/seqkit, >= v2.5.0 which has the new command `seqkit merge-slides`.
- taxonkit: https://github.com/shenwei356/taxonkit
- csvtk: https://github.com/shenwei356/csvtk, >=v0.26.0 is needed for the new features of `pretty` subcommand.

## Database

- KMCP database v2023.05: please choose files below to download from [here](https://1drv.ms/f/s!Ag89cZ8NYcqtlhYlyHuaKweYNu93?e=9NOxqg).
    - GTDB index files (need to be uncompressed)
        - Option 1: **(Recommended) with both plasmid and prophage masked** with [genomad](https://github.com/apcamargo/genomad).
            - `gtdb_masked.part_1.kmcp.tar.gz`
            - `gtdb_masked.part_2.kmcp.tar.gz`
        - Option 2: **only plasmid filtered out** according to sequence names.
            - `gtdb.part_1.kmcp.tar.gz`
            - `gtdb.part_2.kmcp.tar.gz`
    - Taxdump files and TaxId mapping file:
        - Option 1: **(Recommended) GTDB+NCBI taxonomy**. Actually, only the GTDB taxonomy is used here, cause we only use the GTDB KMCP database, no viral or fungal databases are used.
            - `taxdump.gtdbR214.1+ncbi.tar.gz`, uncompressed directory is `taxdump.gtdb+ncbi/`
            - `taxid.gtdbR214.1+ncbi.map`
        - Option 2: **NCBI taxonomy**
            - `taxdump.tar.gz`, uncompressed directory is `taxdump/`
            - `taxid.map`

## Hardware

- RAM >= 64GB
- CPUs >= 32 preferred.

## Steps

Generating reads (sliding a 200-bp windows with a step of 50-bp) from the genome,
and performing metagenomoic profiling with them.

    # variables
    input=SAMN02360712.contigs.fa.gz

    taxdump=~/ws/data/kmcp2023/taxdump.gtdb+ncbi/
    taxmap=~/ws/data/kmcp2023/taxid.gtdbR214.1+ncbi.map

    # search against GTDB databases
    # !!! if the KMCP databases are in a network-attached storage disk (NAS),
    # !!! please add the flag "-w" to "kmcp search"
    for db in ~/ws/data/kmcp2023/gtdb_masked.part_{1,2}.kmcp/; do
        seqkit sliding -g -s 50 -W 200 $input \
            | kmcp search -w -d $db -o $input.kmcp@$(basename $db).tsv.gz
    done

    # merge seach results
    kmcp merge -o $input.kmcp.tsv.gz $input.kmcp@*.tsv.gz

    # profiling and outputing read classification results.
    # here the profiling mode 0 with some modification is used.
    kmcp profile -X $taxdump -T $taxmap $input.kmcp.tsv.gz \
        --mode 0 \
        --min-chunks-reads 2 \
        --min-chunks-fraction 0.1 \
        --min-uniq-reads 2 \
        --min-hic-ureads 2 \
        --min-hic-ureads-qcov 0.75  \
        -o $input.kmcp.tsv.gz.k.profile \
        -B $input.kmcp.tsv.gz.binning.gz

    cat $input.kmcp.tsv.gz.k.profile \
        | csvtk cut -t -f ref,percentage,score,chunksFrac,reads,taxid,taxname \
        | csvtk pretty -t

    ref               percentage   score    chunksFrac   reads    taxid        taxname
    ---------------   ----------   ------   ----------   ------   ----------   -------------------------
    GCF_001457615.1   99.404207    100.00   1.00         113705   1696268695   Pseudomonas aeruginosa
    GCF_000017205.1   0.287847     90.83    1.00         343      1859336444   Pseudomonas aeruginosa_A
    GCF_000008985.1   0.062836     100.00   0.70         23       418881019    Francisella tularensis
    GCF_002091755.1   0.036987     83.61    0.20         38       744524063    Pseudomonas nitroreducens
    GCA_003248965.1   0.043636     100.00   0.10         24       1907132177   Brevundimonas sp003248965
    GCF_014323605.1   0.026607     100.00   0.10         9        1727494734   Lactobacillus kimbladii_B
    GCF_002086125.1   0.012896     94.72    0.10         11       2029805517   Mycobacterium arosiense
    GCF_015070855.1   0.072767     91.67    0.10         60       1955376258   Stutzerimonas lopnurensis
    GCF_008369105.1   0.009527     90.83    0.10         6        1075315903   Aquicoccus porphyridii
    GCF_014199795.1   0.030303     86.67    0.10         27       775010010    Xanthomonas_A sp014199795
    GCF_003797765.1   0.012387     85.83    0.10         14       49646394     Ramlibacter sp003797765


Checking contaminated regions. Here we check "inter-genus" contamination.

    # the taxid of the predominant target of the highest percentage
    taxidMain=$(csvtk cut -t -f taxid $input.kmcp.tsv.gz.k.profile | sed 1d | head -n 1)

    # the genus taxid of the target.
    # the placeholders of different ranks:
    #    {k}: superkingdom
    #    {K}: kingdom
    #    {p}: phylum
    #    {c}: class
    #    {o}: order
    #    {f}: family
    #    {g}: genus
    #    {s}: species
    taxidMainGenus=$(echo $taxidMain | taxonkit reformat --data-dir $taxdump -I 1 -t -f '{g}' | cut -f 3)

    # find out sequence IDs belonging to other genera (potential contaminated contigs).
    #   1. filter out reads with taxids equal to and below the genus node of the predominant target.
    #   2. filter out reads with taxids above the genus node of the predominant target.
    zcat $input.kmcp.tsv.gz.binning.gz | grep -v ^@ \
        | csvtk grep -Ht -f 2 -v -P <(taxonkit list --data-dir $taxdump --ids $taxidMainGenus -I "") \
        | csvtk grep -Ht -f 2 -v -P <(csvtk cut -t -f taxpathsn $input.kmcp.tsv.gz.k.profile \
                                        | sed 1d | head -n 1 | sed "s/;/\n/g") \
            -o $input.kmcp.tsv.gz.binning.filtered.tsv

    # merge regions. seqkit >=v2.5.0 is needed.
    seqkit merge-slides $input.kmcp.tsv.gz.binning.filtered.tsv --quiet \
        -o $input.kmcp.tsv.gz.cont.tsv

    # collect taxonomy info for each contig
    csvtk mutate -Ht -p '^(.+)_sliding' $input.kmcp.tsv.gz.binning.filtered.tsv \
        | csvtk freq -Ht -f 3,2 -nr \
        | taxonkit lineage --data-dir $taxdump -i 2 -n -L \
        | awk -F'\t' '{print $1"\t"$2"("$4")"}' \
        | csvtk fold -Ht -f 1 -v 2 -s ',' \
            -o $input.kmcp.tsv.gz.binning.filtered.tsv.taxa

    # sum up regions of each contigs and append taxonomy information
    awk '{print $0"\t"($3-$2)}' $input.kmcp.tsv.gz.cont.tsv \
        | csvtk summary -Ht -g 1 -f 4:count -f 4:sum  -w 0 \
        | csvtk join -Ht - <(seqkit fx2tab -ni -l $input) \
        | awk '{print $0"\t"($3/$4)}' \
        | csvtk join -Ht - $input.kmcp.tsv.gz.binning.filtered.tsv.taxa \
        | csvtk add-header -Ht -n chr,regions,len,contig_len,proportion,taxa \
        | csvtk sort -t -k proportion:nr \
        | tee $input.kmcp.tsv.gz.cont.details2.tsv \
        | csvtk pretty -t -x , -W 40 -S bold

    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━┳━━━━━━┳━━━━━━━━━━━━┳━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    ┃ chr                      ┃ regions ┃ len  ┃ contig_len ┃ proportion ┃ taxa                                     ┃
    ┣━━━━━━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━╋━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
    ┃ SAMN02360712.contig00044 ┃ 1       ┃ 1151 ┃ 1151       ┃ 1          ┃ 418881019(Francisella tularensis)        ┃
    ┣━━━━━━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━╋━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
    ┃ SAMN02360712.contig00024 ┃ 1       ┃ 1300 ┃ 54646      ┃ 0.0237895  ┃ 1907132177(Brevundimonas sp003248965)    ┃
    ┣━━━━━━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━╋━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
    ┃ SAMN02360712.contig00003 ┃ 6       ┃ 4900 ┃ 583201     ┃ 0.00840191 ┃ 1955376258(Stutzerimonas lopnurensis),   ┃
    ┃                          ┃         ┃      ┃            ┃            ┃ 2029805517(Mycobacterium arosiense),     ┃
    ┃                          ┃         ┃      ┃            ┃            ┃ 1075315903(Aquicoccus porphyridii)       ┃
    ┣━━━━━━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━╋━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
    ┃ SAMN02360712.contig00028 ┃ 1       ┃ 200  ┃ 47311      ┃ 0.00422735 ┃ 49646394(Ramlibacter sp003797765)        ┃
    ┣━━━━━━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━╋━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
    ┃ SAMN02360712.contig00002 ┃ 1       ┃ 1500 ┃ 622965     ┃ 0.00240784 ┃ 775010010(Xanthomonas_A sp014199795)     ┃
    ┣━━━━━━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━╋━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
    ┃ SAMN02360712.contig00001 ┃ 2       ┃ 1350 ┃ 670989     ┃ 0.00201196 ┃ 49646394(Ramlibacter sp003797765),       ┃
    ┃                          ┃         ┃      ┃            ┃            ┃ 1727494734(Lactobacillus kimbladii_B)    ┃
    ┣━━━━━━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━╋━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
    ┃ SAMN02360712.contig00012 ┃ 1       ┃ 200  ┃ 164357     ┃ 0.00121686 ┃ 49646394(Ramlibacter sp003797765)        ┃
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━┻━━━━━━┻━━━━━━━━━━━━┻━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛


We can see the whole (proportion: 1) contig `SAMN02360712.contig00044` is from a species belonging to a different genus even a different order (see the table below). So it is likely be a contaminated sequence.

    echo 1696268695 418881019 \
        | sed -E 's/\s+/\n/g' \
        | taxonkit reformat --data-dir $taxdump -I 1 -f '{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}' \
        | csvtk add-header -t -n taxid,kingdom,phylum,class,order,family,genus,species \
        | csvtk csv2md -t

|taxid     |kingdom |phylum        |class              |order          |family          |genus      |species               |
|:---------|:-------|:-------------|:------------------|:--------------|:---------------|:----------|:---------------------|
|1696268695|Bacteria|Pseudomonadota|Gammaproteobacteria|Pseudomonadales|Pseudomonadaceae|Pseudomonas|Pseudomonas aeruginosa|
|418881019 |Bacteria|Pseudomonadota|Gammaproteobacteria|Francisellales |Francisellaceae |Francisella|Francisella tularensis|

## Results of other genomes

SAMEA2437751.contigs.fa.gz

    ref               percentage   score    chunksFrac   reads   taxid        taxname
    ---------------   ----------   ------   ----------   -----   ----------   -----------------------------------
    GCF_900636965.1   96.429166    100.00   1.00         46468   1117865993   Lacticaseibacillus rhamnosus
    GCF_000829055.1   0.549915     100.00   0.60         273     825607537    Lacticaseibacillus casei
    GCF_000829035.1   0.354971     100.00   0.60         180     48517550     Lacticaseibacillus paracasei
    GCF_000160855.1   0.644110     76.11    0.30         206     874525029    Lactobacillus helveticus
    GCF_004010835.2   0.086689     71.39    0.30         30      1671588307   Oenococcus sp004010835
    GCF_003946165.1   0.634638     100.00   0.20         267     421700347    Lacticaseibacillus baoqingensis
    GCF_018314255.1   0.076443     100.00   0.20         31      1789109845   Lentilactobacillus buchneri
    GCF_001435035.1   0.063781     100.00   0.20         25      616476218    Lacticaseibacillus manihotivorans
    GCF_000425885.1   0.407963     98.89    0.20         181     967313734    Schleiferilactobacillus harbinensis
    GCF_009687905.1   0.099333     93.61    0.20         35      1782720633   Lacticaseibacillus zhaodongensis
    GCF_000241055.1   0.343355     92.78    0.20         103     1598366773   Oenococcus kitaharae
    GCF_004123795.1   0.045943     91.67    0.20         16      433079397    Lacticaseibacillus chiayiensis
    GCF_003946675.1   0.024679     89.72    0.10         8       754863007    Lapidilactobacillus gannanensis
    GCF_001434935.1   0.012219     88.61    0.10         4       1458176950   Liquorilactobacillus uvarum
    GCF_001042675.1   0.049260     88.34    0.10         15      688400175    Parascardovia denticolens
    GCF_016908275.1   0.049721     86.11    0.10         13      697696210    Periweissella beninensis
    GCF_001312865.1   0.069118     79.72    0.10         22      1758332344   Lacticaseibacillus thailandensis
    GCF_000146325.1   0.016039     78.33    0.10         5       973442108    Pediococcus acidilactici
    GCF_020180945.1   0.042655     73.33    0.10         21      332970717    Lacticaseibacillus sp0201809

    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━┳━━━━━━━┳━━━━━━━━━━━━┳━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    ┃ chr                      ┃ regions ┃ len   ┃ contig_len ┃ proportion ┃ taxa                                               ┃
    ┣━━━━━━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━╋━━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
    ┃ SAMEA2437751.contig00033 ┃ 11      ┃ 8331  ┃ 19131      ┃ 0.435471   ┃ 967313734(Schleiferilactobacillus harbinensis)     ┃
    ┣━━━━━━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━╋━━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
    ┃ SAMEA2437751.contig00030 ┃ 37      ┃ 13249 ┃ 31099      ┃ 0.426027   ┃ 874525029(Lactobacillus helveticus),               ┃
    ┃                          ┃         ┃       ┃            ┃            ┃ 1598366773(Oenococcus kitaharae),                  ┃
    ┃                          ┃         ┃       ┃            ┃            ┃ 688400175(Parascardovia denticolens),              ┃
    ┃                          ┃         ┃       ┃            ┃            ┃ 1671588307(Oenococcus sp004010835),                ┃
    ┃                          ┃         ┃       ┃            ┃            ┃ 1448562449(Oenococcus)                             ┃
    ┣━━━━━━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━╋━━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
    ┃ SAMEA2437751.contig00038 ┃ 11      ┃ 5450  ┃ 14019      ┃ 0.388758   ┃ 1789109845(Lentilactobacillus buchneri),           ┃
    ┃                          ┃         ┃       ┃            ┃            ┃ 874525029(Lactobacillus helveticus),               ┃
    ┃                          ┃         ┃       ┃            ┃            ┃ 1671588307(Oenococcus sp004010835),                ┃
    ┃                          ┃         ┃       ┃            ┃            ┃ 973442108(Pediococcus acidilactici),               ┃
    ┃                          ┃         ┃       ┃            ┃            ┃ 1458176950(Liquorilactobacillus uvarum),           ┃
    ┃                          ┃         ┃       ┃            ┃            ┃ 1598366773(Oenococcus kitaharae)                   ┃
    ┣━━━━━━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━╋━━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
    ┃ SAMEA2437751.contig00047 ┃ 2       ┃ 950   ┃ 6428       ┃ 0.147791   ┃ 697696210(Periweissella beninensis)                ┃
    ┣━━━━━━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━╋━━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
    ┃ SAMEA2437751.contig00059 ┃ 1       ┃ 93    ┃ 643        ┃ 0.144635   ┃ 1671588307(Oenococcus sp004010835)                 ┃
    ┣━━━━━━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━╋━━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
    ┃ SAMEA2437751.contig00046 ┃ 1       ┃ 900   ┃ 6542       ┃ 0.137573   ┃ 1598366773(Oenococcus kitaharae)                   ┃
    ┣━━━━━━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━╋━━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
    ┃ SAMEA2437751.contig00045 ┃ 1       ┃ 850   ┃ 8566       ┃ 0.0992295  ┃ 1598366773(Oenococcus kitaharae)                   ┃
    ┣━━━━━━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━╋━━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
    ┃ SAMEA2437751.contig00028 ┃ 3       ┃ 700   ┃ 38428      ┃ 0.0182159  ┃ 754863007(Lapidilactobacillus gannanensis)         ┃
    ┣━━━━━━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━╋━━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
    ┃ SAMEA2437751.contig00008 ┃ 2       ┃ 850   ┃ 120682     ┃ 0.0070433  ┃ 967313734(Schleiferilactobacillus harbinensis)     ┃
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━┻━━━━━━━┻━━━━━━━━━━━━┻━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

SAMN10397458.contigs.fa.gz


    ref               percentage   score    chunksFrac   reads   taxid        taxname
    ---------------   ----------   ------   ----------   -----   ----------   --------------------
    GCF_000006945.2   99.330920    100.00   1.00         91824   1678121664   Salmonella enterica
    GCA_900478215.1   0.444068     92.78    1.00         392     1270053497   Salmonella houtenae
    GCF_008692845.1   0.213823     82.78    1.00         161     1647206931   Salmonella arizonae
    GCF_011064845.1   0.011190     75.55    0.10         11      1757895327   Citrobacter freundii

    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━┳━━━━━┳━━━━━━━━━━━━┳━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    ┃ chr                      ┃ regions ┃ len ┃ contig_len ┃ proportion ┃ taxa                             ┃
    ┣━━━━━━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━╋━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
    ┃ SAMN10397458.contig00041 ┃ 1       ┃ 40  ┃ 1740       ┃ 0.0229885  ┃ 1757895327(Citrobacter freundii) ┃
    ┣━━━━━━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━╋━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
    ┃ SAMN10397458.contig00006 ┃ 1       ┃ 800 ┃ 271389     ┃ 0.0029478  ┃ 1757895327(Citrobacter freundii) ┃
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━┻━━━━━┻━━━━━━━━━━━━┻━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
