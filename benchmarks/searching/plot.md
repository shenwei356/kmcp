## Preparing data

KMCP vs COBS

    cat bench.kmcp-cobs.short.tsv \
        | csvtk cut -t -f 1,2,4 \
        | csvtk rename -t  -f 2,3 -n COBS,KMCP \
        | csvtk gather -t -k app -v time -f 2,3 \
        | csvtk mutate2 -t -n group -e '"1M 150-bp reads"' \
        | csvtk tab2csv \
        > bench.kmcp-cobs.short.tsv.time.csv

    cat bench.kmcp-cobs.long.tsv \
        | csvtk cut -t -f 1,2,4 \
        | csvtk rename -t  -f 2,3 -n COBS,KMCP \
        | csvtk gather -t -k app -v time -f 2,3 \
        | csvtk mutate2 -t -n group -e '"One genome"' \
        | csvtk tab2csv \
        > bench.kmcp-cobs.long.tsv.time.csv
        
    csvtk concat bench.kmcp-cobs.long.tsv.time.csv bench.kmcp-cobs.short.tsv.time.csv \
        > bench.kmcp-cobs.csv
        
    cat bench.kmcp-cobs.csv \
        | csvtk summary -g group,app -f time:mean \
        | csvtk sort -k group \
        | csvtk cut -f group,app,time:mean \
        | csvtk csv2md 

|group          |app |time:mean|
|:--------------|:---|:--------|
|1M 150-bp reads|KMCP|62.22    |
|1M 150-bp reads|COBS|749.32   |
|One genome     |KMCP|13.32    |
|One genome     |COBS|25.70    |

KMCP vs Mash and Sourmash

    cat bench.kmcp-mash-sourmash.thread1.tsv \
        | csvtk cut -t -f 1,2,4,6 \
        | csvtk rename -t  -f 2,3,4 -n Mash,Sourmash,KMCP \
        | csvtk gather -t -k app -v time -f 2,3,4 \
        | csvtk mutate2 -t -n group -e '"Threads=1"' \
        | csvtk tab2csv \
        > bench.kmcp-mash-sourmash.thread1.tsv.csv
        
    cat bench.kmcp-mash-sourmash.thread8.tsv \
        | csvtk cut -t -f 1,2,4,6 \
        | csvtk rename -t  -f 2,3,4 -n Mash,Sourmash,KMCP \
        | csvtk gather -t -k app -v time -f 2,3,4 \
        | csvtk mutate2 -t -n group -e '"Threads=8"' \
        | csvtk tab2csv \
        > bench.kmcp-mash-sourmash.thread8.tsv.csv
        
    csvtk concat bench.kmcp-mash-sourmash.thread1.tsv.csv bench.kmcp-mash-sourmash.thread8.tsv.csv\
        > bench.kmcp-mash-sourmash.csv

    cat bench.kmcp-mash-sourmash.csv \
        | csvtk summary -g group,app -f time:mean \
        | csvtk sort -k group \
        | csvtk cut -f group,app,time:mean \
        | csvtk csv2md 

|group    |app     |time:mean|
|:--------|:-------|:--------|
|Threads=1|Sourmash|3.48     |
|Threads=1|Mash    |6.46     |
|Threads=1|KMCP    |1.11     |
|Threads=8|Sourmash|3.48     |
|Threads=8|Mash    |4.20     |
|Threads=8|KMCP    |0.60     |
    
## Plotting

    ./plot.R
    
![](bench.searching.jpg)

    
