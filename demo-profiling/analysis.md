### Analysis of ambiguous and unambiguous reads
    
    # clean
    zcat mock2.kmcp.gz  | csvtk rename -C$ -t -f 1 -n query | csvtk uniq -t -f query,target -o t.gz
    
    # frequency of the number of hits
    csvtk freq -t -f query -nr t.gz -o t.gz.freq1
    
    # ID of unambiguous reads
    csvtk filter2 -t -f '$frequency==1' t.gz.freq1 | csvtk cut -t -f query -o t.gz.uniq.reads
        
    # total number of reads
    csvtk nrow t.gz.freq1
    307276
    
    # number of unambiguous reads
    csvtk nrow t.gz.uniq.reads
    195601
    
    # distribution of target number
    csvtk plot hist -t -f frequency t.gz.freq1 -o t.gz.freq1.hist.png \
        --title "distribution of target number" --xlab "target number" --ylab "#reads"
    
The distribution of target number:

![](t.gz.freq1.hist.png)
    
    # pair of query and target
    csvtk freq -t -f query,target  t.gz -o t.gz.pair
    
The number of unambiguous reads:

    csvtk grep -t -f query -P t.gz.uniq.reads t.gz.pair \
        | csvtk freq -t -f target -nr \
        | csvtk join -t -k - <(csvtk add-header -t -n target,name name2.map) genome-size.tsv \
        | csvtk mutate2 -t -n fraction -e '$frequency / ${genome-size}' -w 6 \
        | csvtk csv2md -t
        
|target         |frequency|name                      |genome-size|fraction|
|:--------------|:--------|:-------------------------|:----------|:-------|
|GCF_000742135.1|65908    |Klebsiella pneumoniae     |5545784    |0.011884|
|GCF_000006945.2|57244    |Salmonella enterica       |4951383    |0.011561|
|GCF_003697165.2|26412    |Escherichia coli          |5034834    |0.005246|
|GCF_002950215.1|18303    |Shigella flexneri         |4938295    |0.003706|
|GCF_002949675.1|17294    |Shigella dysenteriae      |4578459    |0.003777|
|GCF_000392875.1|3322     |Enterococcus faecalis     |2881400    |0.001153|
|GCF_001544255.1|3012     |Enterococcus faecium      |2484851    |0.001212|
|GCF_001027105.1|1612     |Staphylococcus aureus     |2782562    |0.000579|
|GCF_006742205.1|1425     |Staphylococcus epidermidis|2427041    |0.000587|
|GCF_000017205.1|395      |Pseudomonas aeruginosa    |6588339    |0.000060|
|GCF_009759685.1|256      |Acinetobacter baumannii   |3990388    |0.000064|
|GCF_001096185.1|221      |Streptococcus pneumoniae  |2117177    |0.000104|
|GCF_000148585.2|197      |Streptococcus mitis       |1868883    |0.000105|

The number of ambiguous reads:

    csvtk grep -t -f query -v -P t.gz.uniq.reads t.gz.pair \
        | csvtk replace -t -f target -p '(.+)' -k name2.map -r '{kv}' \
        | csvtk sort -t -k query -k target \
        | csvtk fold -t -f query -v target \
        | csvtk freq -t -f target -nr \
        | csvtk filter2 -t -f '$frequency >= 200' \
        | csvtk csv2md -t
        
|target                                                                                               |frequency|
|:----------------------------------------------------------------------------------------------------|:--------|
|Escherichia coli; Shigella dysenteriae; Shigella flexneri                                            |57564    |
|Shigella dysenteriae; Shigella flexneri                                                              |22195    |
|Escherichia coli; Shigella flexneri                                                                  |12888    |
|Escherichia coli; Shigella dysenteriae                                                               |9224     |
|Escherichia coli; Klebsiella pneumoniae; Salmonella enterica; Shigella dysenteriae; Shigella flexneri|2610     |
|Escherichia coli; Salmonella enterica; Shigella dysenteriae; Shigella flexneri                       |1489     |
|Escherichia coli; Klebsiella pneumoniae; Shigella dysenteriae; Shigella flexneri                     |1404     |
|Klebsiella pneumoniae; Shigella dysenteriae; Shigella flexneri                                       |821      |
|Escherichia coli; Klebsiella pneumoniae; Shigella dysenteriae                                        |666      |
|Escherichia coli; Klebsiella pneumoniae                                                              |602      |
|Klebsiella pneumoniae; Salmonella enterica                                                           |284      |
|Klebsiella pneumoniae; Shigella dysenteriae                                                          |203      |
