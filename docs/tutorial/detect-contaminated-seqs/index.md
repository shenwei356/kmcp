## Detecting contaminated sequences in a bacteria genome

## tools

- kmcp: https://github.com/shenwei356/kmcp
- seqkit: >= v2.5.0, https://github.com/shenwei356/seqkit/issues/390#issuecomment-1633495130
- taxonkit: https://github.com/shenwei356/taxonkit
- csvtk: https://github.com/shenwei356/csvtk

## database

- KMCP r214: https://1drv.ms/f/s!Ag89cZ8NYcqtlgEtskAjVVosItE7?e=CdkwhG
    - gtdb.part_1.kmcp.tar.gz, need to be uncompressed.
    - gtdb.part_2.kmcp.tar.gz, need to be uncompressed.
    - taxdump.tar.gz, need to be uncompressed.
    - taxid.map


## hardware

- RAM >= 64GB
- #CPUs >= 32 preferred.

## steps

Generating reads from the genome, and perform metagenomoic profiling with them.

    # variables
    input=SAMN02360712.contigs.fa.gz
    taxdump=~/ws/data/kmcp2023/taxdump/
    taxmap=~/ws/data/kmcp2023/taxid.map

    # search against GTDB databases
    seqkit sliding -s 30 -W 150 $input \
        | kmcp search -d ~/ws/data/kmcp2023/gtdb.part_1.kmcp/ -o $input.kmcp@gtdb.part_1.tsv.gz

    seqkit sliding -s 30 -W 150 $input \
        | kmcp search -d ~/ws/data/kmcp2023/gtdb.part_2.kmcp/ -o $input.kmcp@gtdb.part_2.tsv.gz

    # merge result
    kmcp merge -o $input.kmcp.tsv.gz $input.kmcp@gtdb.part_*.tsv.gz

    # profiling and output reads classification results.
    # here the profiling mode 0 with the highest sensitivity is used.
    kmcp profile -X $taxdump -T $taxmap t.kmcp.tsv.gz -m 0 \
        -o $input.kmcp.tsv.gz.k.profile -B $input.kmcp.tsv.gz.binning.gz

    cat $input.kmcp.tsv.gz.k.profile \
        | csvtk cut -t -f ref,percentage,score,chunksFrac,reads,taxid,taxname \
        | csvtk pretty -t

    ref               percentage   score    chunksFrac   reads    taxid     taxname
    ---------------   ----------   ------   ----------   ------   -------   -------------------------------------------------
    GCF_001457615.1   99.327059    100.00   1.00         189446   287       Pseudomonas aeruginosa
    GCF_000008985.1   0.059509     100.00   1.00         34       177416    Francisella tularensis subsp. tularensis SCHU S4
    GCF_000017205.1   0.490404     91.54    1.00         976      381754    Pseudomonas aeruginosa PA7
    GCF_900112375.1   0.042622     100.00   0.40         93       53408     Pseudomonas citronellolis
    GCF_002091755.1   0.020159     86.16    0.40         37       1215112   Pseudomonas nitroreducens NBRC 12694
    GCF_003797765.1   0.009539     91.54    0.30         20       1882741   Ramlibacter sp. WS9
    GCF_009735585.1   0.015466     100.00   0.20         16       1028746   Christiangramia aestuarii
    GCF_008369105.1   0.006725     90.00    0.20         9        1852029   Aquicoccus porphyridii
    GCF_002086125.1   0.004465     90.00    0.20         8        1265311   Mycobacterium arosiense ATCC BAA-1401 = DSM 45069
    GCF_005877905.1   0.020988     74.61    0.20         42       2571105   Pseudomonas nicosulfuronedens
    GCF_900187975.1   0.003064     70.39    0.20         6        366289    Pseudomonas delhiensis

Check contaminated regions

    # the taxid of the predominant target
    taxidMain=$(csvtk cut -t -f taxid $input.kmcp.tsv.gz.k.profile | sed 1d | head -n 1)
    # the genus taxid of the target
    taxidMainGenus=$(echo taxidMain | taxonkit reformat --data-dir $taxdump -I 1 -t -f '{g}' | cut -f 3)

    # find out sequence IDs belonging to other genera (potential contaminated contigs)
    zcat $input.kmcp.tsv.gz.binning.gz | grep -v ^@ \
        | csvtk grep -Ht -f 2 -v -P <(taxonkit list --data-dir $taxdump --ids 286 -I "") \
        | csvtk grep -Ht -f 2 -v -P <(csvtk cut -t -f taxpathsn $input.kmcp.tsv.gz.k.profile \
                                        | sed 1d | head -n 1 | sed "s/;/\n/g") \
            -o $input.kmcp.tsv.gz.binning.filtered.tsv

    # merge regions
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
        | csvtk add-header -Ht -n chr,begin,end,contig_len,len,frac,taxa \
        | csvtk sort -t -k frac:nr \
        | tee $input.kmcp.tsv.gz.cont.details.tsv \
        | csvtk pretty -t

    chr                        begin    end      contig_len   len    frac          taxa
    ------------------------   ------   ------   ----------   ----   -----------   ------------------------------------------------------------------------------------------
    SAMN02360712.contig00044   0        1140     1151         1140   0.990443      177416(Francisella tularensis subsp. tularensis SCHU S4)
    SAMN02360712.contig00021   3720     3900     66071        180    0.00272434    1028746(Christiangramia aestuarii)
    SAMN02360712.contig00023   45660    45810    55929        150    0.00268197    1265311(Mycobacterium arosiense ATCC BAA-1401 = DSM 45069)
    SAMN02360712.contig00012   163650   163890   164357       240    0.00146024    1028746(Christiangramia aestuarii)
    SAMN02360712.contig00016   11550    11700    106327       150    0.00141074    1028746(Christiangramia aestuarii)
    SAMN02360712.contig00001   454530   455190   670989       660    0.000983623   1882741(Ramlibacter sp. WS9)
    SAMN02360712.contig00008   64890    65130    279605       240    0.000858354   1028746(Christiangramia aestuarii)
    SAMN02360712.contig00010   90000    90150    183802       150    0.000816096   1882741(Ramlibacter sp. WS9)
    SAMN02360712.contig00003   471450   471840   583201       390    0.000668723   1852029(Aquicoccus porphyridii),1265311(Mycobacterium arosiense ATCC BAA-1401 = DSM 45069)
    SAMN02360712.contig00003   36810    37170    583201       360    0.000617283   1852029(Aquicoccus porphyridii),1265311(Mycobacterium arosiense ATCC BAA-1401 = DSM 45069)
    SAMN02360712.contig00002   362220   362490   622965       270    0.000433411   1028746(Christiangramia aestuarii),1852029(Aquicoccus porphyridii)
    SAMN02360712.contig00004   272760   272910   478274       150    0.000313628   1882741(Ramlibacter sp. WS9)
    SAMN02360712.contig00002   495990   496140   622965       150    0.000240784   1028746(Christiangramia aestuarii),1852029(Aquicoccus porphyridii)
