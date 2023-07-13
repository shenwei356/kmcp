## Detecting contaminated sequences in a bacteria genome

## tools

- kmcp: https://github.com/shenwei356/kmcp
- seqkit: >= v2.5.0, https://github.com/shenwei356/seqkit/issues/390#issuecomment-1633495130
- taxonkit: https://github.com/shenwei356/taxonkit
- csvtk: https://github.com/shenwei356/csvtk

## database

- GTDB r214: https://1drv.ms/f/s!Ag89cZ8NYcqtlgEtskAjVVosItE7?e=CdkwhG
    - gtdb.part_1.kmcp.tar.gz, need to be uncompressed.
    - gtdb.part_2.kmcp.tar.gz, need to be uncompressed.
    - taxdump.tar.gz, need to be uncompressed.
    - taxid.map


## hardware

- RAM >= 64GB
- #CPUs >= 32 preferred.

## steps

Generating reads (sliding a 200-bp windows with a step of 50-bp) from the genome,
and performing metagenomoic profiling with them.

    # variables
    input=SAMN02360712.contigs.fa.gz
    taxdump=~/ws/data/kmcp2023/taxdump/
    taxmap=~/ws/data/kmcp2023/taxid.map

    # search against GTDB databases
    # !!! if the KMCP databases are in a network-attached storage disk (NAS),
    # !!! please add the flag "-w" to kmcp
    seqkit sliding -g -s 50 -W 200 $input \
        | kmcp search -d ~/ws/data/kmcp2023/gtdb.part_1.kmcp/ -o $input.kmcp@gtdb.part_1.tsv.gz

    seqkit sliding -g -s 50 -W 200 $input \
        | kmcp search -d ~/ws/data/kmcp2023/gtdb.part_2.kmcp/ -o $input.kmcp@gtdb.part_2.tsv.gz

    # merge seach results
    kmcp merge -o $input.kmcp.tsv.gz $input.kmcp@gtdb.part_*.tsv.gz

    # profiling and outputing read classification results.
    # here the profiling mode 0 with the highest sensitivity is used.
    kmcp profile -X $taxdump -T $taxmap $input.kmcp.tsv.gz -m 0 \
        -o $input.kmcp.tsv.gz.k.profile -B $input.kmcp.tsv.gz.binning.gz

    cat $input.kmcp.tsv.gz.k.profile \
        | csvtk cut -t -f ref,percentage,score,chunksFrac,reads,taxid,taxname \
        | csvtk pretty -t

    ref               percentage   score    chunksFrac   reads    taxid     taxname
    ---------------   ----------   ------   ----------   ------   -------   ------------------------------------------------
    GCF_001457615.1   99.405509    100.00   1.00         114023   287       Pseudomonas aeruginosa
    GCF_000008985.1   0.062663     100.00   1.00         23       177416    Francisella tularensis subsp. tularensis SCHU S4
    GCF_000017205.1   0.453168     90.56    1.00         542      381754    Pseudomonas aeruginosa PA7
    GCF_900112375.1   0.037236     100.00   0.40         49       53408     Pseudomonas citronellolis
    GCF_002091755.1   0.026943     86.39    0.20         30       1215112   Pseudomonas nitroreducens NBRC 12694
    GCF_009735585.1   0.014481     85.56    0.20         9        1028746   Christiangramia aestuarii

Checking contaminated regions

    # the taxid of the predominant target
    taxidMain=$(csvtk cut -t -f taxid $input.kmcp.tsv.gz.k.profile | sed 1d | head -n 1)
    # the genus taxid of the target
    taxidMainGenus=$(echo $taxidMain | taxonkit reformat --data-dir $taxdump -I 1 -t -f '{g}' | cut -f 3)

    # find out sequence IDs belonging to other genera (potential contaminated contigs)
    zcat $input.kmcp.tsv.gz.binning.gz | grep -v ^@ \
        | csvtk grep -Ht -f 2 -v -P <(taxonkit list --data-dir $taxdump --ids $taxidMainGenus -I "") \
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
    ------------------------   ------   ------   ----------   ----   -----------   --------------------------------------------------------
    SAMN02360712.contig00044   0        1151     1151         1151   1             177416(Francisella tularensis subsp. tularensis SCHU S4)
    SAMN02360712.contig00012   163600   163900   164357       300    0.00182529    1028746(Christiangramia aestuarii)
    SAMN02360712.contig00008   64850    65150    279605       300    0.00107294    1028746(Christiangramia aestuarii)
    SAMN02360712.contig00002   362200   362500   622965       300    0.000481568   1028746(Christiangramia aestuarii)
