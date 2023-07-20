# Detecting contaminated sequences in a bacteria assembly

## Changelog

Tracked changes in [github](https://github.com/shenwei356/kmcp/commits/main/docs/tutorial/detecting-contaminated-seqs/index.md)

- 2023-07-20
    - Recommend using masked GTDB index files
    - Recommend using GTDB taxdump files rather than NCBI taxdump files.
- 2023-07-16
    - First public version.

## Tools

- kmcp: https://github.com/shenwei356/kmcp
- seqkit: >= v2.5.0 which has the new command `seqkit merge-slides`.
- taxonkit: https://github.com/shenwei356/taxonkit
- csvtk: https://github.com/shenwei356/csvtk

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
    # here the profiling mode 0 with the highest sensitivity is used.
    kmcp profile -X $taxdump -T $taxmap $input.kmcp.tsv.gz -m 0 \
        -o $input.kmcp.tsv.gz.k.profile -B $input.kmcp.tsv.gz.binning.gz

    cat $input.kmcp.tsv.gz.k.profile \
        | csvtk cut -t -f ref,percentage,score,chunksFrac,reads,taxid,taxname \
        | csvtk pretty -t

    ref               percentage   score    chunksFrac   reads    taxid        taxname
    ---------------   ----------   ------   ----------   ------   ----------   -------------------------
    GCF_001457615.1   99.579056    100.00   1.00         113704   1696268695   Pseudomonas aeruginosa
    GCF_000008985.1   0.062945     100.00   1.00         23       418881019    Francisella tularensis
    GCF_000017205.1   0.289416     90.83    1.00         345      1859336444   Pseudomonas aeruginosa_A
    GCF_002091755.1   0.037182     83.61    0.40         39       744524063    Pseudomonas nitroreducens
    GCF_003797765.1   0.012409     85.83    0.20         14       49646394     Ramlibacter sp003797765
    GCF_009735585.1   0.006522     77.50    0.20         4        182564388    Gramella aestuarii
    GCF_900112375.1   0.012470     66.94    0.20         13       132541362    Pseudomonas citronellolis


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

    # append contig length and taxid(s)
    csvtk join -Ht $input.kmcp.tsv.gz.cont.tsv <(seqkit fx2tab -ni -l $input) \
        | awk '{print $0"\t"($3-$2)"\t"($3-$2)/$4}' \
        | csvtk join -Ht - $input.kmcp.tsv.gz.binning.filtered.tsv.taxa \
        | csvtk add-header -Ht -n chr,begin,end,contig_len,len,proportion,taxa \
        | csvtk sort -t -k proportion:nr \
        | tee $input.kmcp.tsv.gz.cont.details.tsv \
        | csvtk pretty -t

    chr                        begin    end      contig_len   len    proportion    taxa
    ------------------------   ------   ------   ----------   ----   -----------   ---------------------------------------------------------------
    SAMN02360712.contig00044   0        1151     1151         1151   1             418881019(Francisella tularensis)
    SAMN02360712.contig00028   26600    26800    47311        200    0.00422735    49646394(Ramlibacter sp003797765)
    SAMN02360712.contig00012   163700   163900   164357       200    0.00121686    182564388(Gramella aestuarii),49646394(Ramlibacter sp003797765)
    SAMN02360712.contig00012   155750   155950   164357       200    0.00121686    182564388(Gramella aestuarii),49646394(Ramlibacter sp003797765)
    SAMN02360712.contig00001   454500   455250   670989       750    0.00111775    49646394(Ramlibacter sp003797765)
    SAMN02360712.contig00002   362200   362500   622965       300    0.000481568   182564388(Gramella aestuarii)

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
