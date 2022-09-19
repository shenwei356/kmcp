# Profiling accuracies with different genome chunk numbers

Here we use the same reference genomes and simulated community to analyze
how would the number of genome chunks affects the profiling accuracy.

## Building databases

Reference genomes are split into 1, 5, 10, and 20 chunks respectively.

    for n in 1 5 10 20; do
        # computing k-mers
        kmcp compute \
            --in-dir refs/ \
            --ref-name-regexp "^([\w\.\_]+\.\d+)" \
            --seq-name-filter "plasmid" \
            --kmer 21 \
            --split-number $n \
            --split-overlap 150 \
            --out-dir refs-k21-n$n \
            --force

        # indexing k-mers
        kmcp index \
            --in-dir refs-k21-n$n/\
            --num-hash 1 \
            --false-positive-rate 0.3 \
            --out-dir refs-k21-n$n.kmcp \
            --force
        
        # clean tmp files
        rm -rf refs-k21-n$n
    done

## Searching

    for db in *.kmcp; do
        # paired information are not used.
        kmcp search \
            --db-dir $db \
            mock_1.fastq.gz \
            mock_2.fastq.gz \
            --out-file mock.kmcp@$db.gz
    done

## False discovery rates of read mapping

FDR = FP / (FP + TP)
 
    # read IDs of each reference assembly accession
    ls mock/*.fa | rush 'seqkit seq -n -i {}.*.sample.fq.gz -o {}.ids'

    for f in *.kmcp.gz; do
        cut -f 1 taxid.map | while read id; do
            tp=$(csvtk grep -Ht -f 6 -p $id $f \
                | csvtk uniq -Ht -f 1 \
                | csvtk grep -Ht -f 1 -P mock/$id.fa.ids \
                | csvtk nrow -Ht)
            fp=$(csvtk grep -Ht -f 6 -p $id $f \
                | csvtk uniq -Ht -f 1 \
                | csvtk grep -Ht -f 1 -P mock/$id.fa.ids -v \
                | csvtk nrow -Ht)
            fdr=$(echo "scale=6; $fp/($fp+$tp)" | bc)
            echo -e "$f\t$id\t$tp\t$fp\t$fdr"
        done
    done \
        | csvtk replace -Ht -p '.+n(\d+).+' -r '$1' \
        | csvtk add-header -t -n chunks,ref,tp,fp,fdr \
        > fdr.tsv

    cat fdr.tsv \
        | csvtk sort -t -k chunks:n -k ref \
        | csvtk fold -t -f ref -v fdr \
        | csvtk sep -t -R -f fdr -s '; ' -n '1,5,10,20' \
        | csvtk join -t - <(csvtk cut -t -f id,species,reads mock.gs.tsv) \
        | tee fdr2.tsv \
        | csvtk sort -t -k 2:n \
        | csvtk rename2 -t -f 2-5 -p '^' -r 'chunks=' \
        | csvtk csv2md -t
        
|ref            |chunks=1|chunks=5|chunks=10|chunks=20|species                   |reads|
|:--------------|:-------|:-------|:--------|:--------|:-------------------------|:----|
|GCF_000392875.1|.009667 |.011386 |.010841  |.008115  |Enterococcus faecalis     |3798 |
|GCF_001544255.1|.013689 |.015161 |.012991  |.015528  |Enterococcus faecium      |3294 |
|GCF_006742205.1|.021379 |.038841 |.038208  |.034098  |Staphylococcus epidermidis|1594 |
|GCF_001027105.1|.022370 |.025686 |.027027  |.026978  |Staphylococcus aureus     |1814 |
|GCF_000006945.2|.052860 |.067012 |.066121  |.066822  |Salmonella enterica       |65664|
|GCF_000742135.1|.075247 |.077225 |.076650  |.081913  |Klebsiella pneumoniae     |73472|
|GCF_900638025.1|.098618 |.163709 |.167558  |.168662  |Haemophilus parainfluenzae|2744 |
|GCF_001457655.1|.122909 |.192624 |.196218  |.202635  |Haemophilus influenzae    |2492 |
|GCF_000148585.2|.159183 |.258389 |.250000  |.275080  |Streptococcus mitis       |242  |
|GCF_001096185.1|.168384 |.181229 |.189102  |.200636  |Streptococcus pneumoniae  |284  |
|GCF_000017205.1|.348837 |.388888 |.413946  |.414201  |Pseudomonas aeruginosa    |434  |
|GCF_003697165.2|.465953 |.490760 |.490679  |.483973  |Escherichia coli          |66770|
|GCF_002950215.1|.488692 |.515451 |.517123  |.521609  |Shigella flexneri         |65490|
|GCF_002949675.1|.492149 |.525276 |.526293  |.527906  |Shigella dysenteriae      |60720|
|GCF_009759685.1|.504366 |.581548 |.608895  |.601265  |Acinetobacter baumannii   |272  |


## Profiling accuracies

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
            --log $f.kmcp.profile.log
    done

    for f in *.kmcp.profile; do
        chunks=$(echo $f | sed 's/mock.kmcp@refs-k21-n//' | sed 's/.kmcp.gz.kmcp.profile//')
        cat $f \
            | csvtk cut -t -f ref,reads,percentage,taxname \
            | csvtk join -t - <(csvtk cut -t -f id,abundance mock.gs.tsv) \
            | csvtk rename -t -f abundance -n ground_truth \
            | csvtk cut -t -f ref,taxname,ground_truth,percentage \
            | csvtk rename -t -f percentage -n $chunks \
            > $f.tsv;
    done
    
    csvtk join -t -f ref,taxname,ground_truth --infile-list <(ls *.kmcp.profile.tsv | sort -V) \
        | tee abundance.tsv \
        | csvtk rename2 -t -f 4-7 -p '^' -r 'chunks=' \
        | csvtk csv2md -t
    
|ref            |taxname                   |ground_truth|chunks=1 |chunks=5 |chunks=10|chunks=20|
|:--------------|:-------------------------|:-----------|:--------|:--------|:--------|:--------|
|GCF_003697165.2|Escherichia coli          |18.093707   |19.725027|18.873090|18.663804|17.478403|
|GCF_000742135.1|Klebsiella pneumoniae     |18.075498   |18.849011|17.769042|17.738253|17.879812|
|GCF_000006945.2|Salmonella enterica       |18.093898   |18.039382|18.188425|18.143627|18.078030|
|GCF_002950215.1|Shigella flexneri         |18.093780   |17.440064|17.466549|17.728060|18.495286|
|GCF_002949675.1|Shigella dysenteriae      |18.094390   |16.190428|18.135503|18.201855|18.597809|
|GCF_000392875.1|Enterococcus faecalis     |1.798382    |1.909217 |1.805773 |1.800250 |1.708806 |
|GCF_001544255.1|Enterococcus faecium      |1.808651    |1.889844 |1.812982 |1.795723 |1.840437 |
|GCF_900638025.1|Haemophilus parainfluenzae|1.815273    |1.800269 |1.816523 |1.809292 |1.812660 |
|GCF_001457655.1|Haemophilus influenzae    |1.798328    |1.757202 |1.791787 |1.787560 |1.772003 |
|GCF_006742205.1|Staphylococcus epidermidis|0.896071    |0.929816 |0.908313 |0.906778 |0.913170 |
|GCF_001027105.1|Staphylococcus aureus     |0.889459    |0.927290 |0.887062 |0.881014 |0.879919 |
|GCF_001096185.1|Streptococcus pneumoniae  |0.183016    |0.185368 |0.176324 |0.177453 |0.175524 |
|GCF_000148585.2|Streptococcus mitis       |0.176667    |0.170479 |0.181657 |0.178996 |0.182886 |
|GCF_000017205.1|Pseudomonas aeruginosa    |0.089875    |0.094310 |0.089483 |0.089177 |0.089206 |
|GCF_009759685.1|Acinetobacter baumannii   |0.093004    |0.092294 |0.097489 |0.098158 |0.096050 |

Assessing with [OPAL](https://github.com/CAMI-challenge/OPAL):

    opal.py -g mock.gs.profile -o opal2/ \
        mock.kmcp@refs-k21-n1.kmcp.gz.cami.profile \
        mock.kmcp@refs-k21-n5.kmcp.gz.cami.profile \
        mock.kmcp@refs-k21-n10.kmcp.gz.cami.profile \
        mock.kmcp@refs-k21-n20.kmcp.gz.cami.profile \
        -l 1,5,10,20 
        
    cat opal2/results.tsv \
        | csvtk grep -t -f tool -p 'Gold standard' -v \
        | csvtk grep -t -f rank -p na -p species \
        | csvtk grep -t -f metric -p 'L1 norm error' -p 'Weighted UniFrac error' \
        | csvtk sort -t -k metric -k tool:N \
        | csvtk rename -t -f tool -n chunks \
        | tee abundance2.tsv \
        | csvtk cut -t -f -sample \
        | csvtk csv2md -t

|chunks|rank   |metric                |value                |
|:-----|:------|:---------------------|:--------------------|
|1     |species|L1 norm error         |0.053504460000000024 |
|5     |species|L1 norm error         |0.018994210000000008 |
|10    |species|L1 norm error         |0.014946969999999966 |
|20    |species|L1 norm error         |0.019261479999999997 |
|1     |na     |Weighted UniFrac error|0.023275701493843635 |
|5     |na     |Weighted UniFrac error|0.0062703427346909885|
|10    |na     |Weighted UniFrac error|0.004712290601693156 |
|20    |na     |Weighted UniFrac error|0.007459834023109471 |
