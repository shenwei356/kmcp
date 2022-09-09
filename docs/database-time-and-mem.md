# Building KMCP database

The steps below is same to theses in [database.md](database.md), but `memusg` is used to

## GTDB 
  
    input=gtdb
    
    # compute k-mers
    #   sequence containing "plasmid" in name are ignored,
    #   reference genomes are split into 10 chunks with 100bp overlap
    #   k = 21
    #
    # elapsed time: 10m:34s
    # peak rss: 3.87 GB
    # file size: 978.37 GB
    #
    memusg -t -s "kmcp compute -I $input -O gtdb-r202-k21-n10 -k 21 -n 10 -l 150 -B plasmid \
        --log gtdb-r202-k21-n10.log -j 32 --force" > gtdb.kmcp.s1.log2 2>&1

    # build database
    #   number of index files: 32, for server with >= 32 CPU cores
    #   bloom filter parameter:
    #     number of hash function: 1
    #     false positive rate: 0.3
    #
    # elapsed time: 11m:48s
    # peak rss: 13.95 GB
    # file size: 58.03 GB
    #
    memusg -t -s "kmcp index -j 32 -I gtdb-r202-k21-n10 -O gtdb.kmcp -n 1 -f 0.3 \
        --log gtdb.kmcp.log --force" > gtdb.kmcp.s2.log2 2>&1
    
    # cp taxid and name mapping file to database directory
    cp taxid.map name.map gtdb.kmcp/
    
## RefSeq fungi

    name=fungi
    
    input=files.renamed    
    
    # elapsed time: 1m:02s
    # peak rss: 11.72 GB
    # file size: 70.52 GB
    #
    memusg -t -s "kmcp compute -I $input -O refseq-$name-k21-n10 \
        -k 21 --seq-name-filter plasmid \
        --split-number 10 --split-overlap 150 \
        --log refseq-$name-k21-n10.log -j 32 --force" > refseq-fungi.kmcp.s1.log2 2>&1
      
    # elapsed time: 52.204s
    # peak rss: 1.19 GB
    # file size: 4.18 GB
    memusg -t -s "kmcp index -I refseq-$name-k21-n10/ -O refseq-fungi.kmcp \
        -j 32 -f 0.3 -n 1 \
        --log refseq-fungi.kmcp.log --force" > refseq-fungi.kmcp.s2.log2 2>&1
    
    # cp taxid and name mapping file to database directory
    cp taxid.map name.map refseq-fungi.kmcp/

## Genbank viral

    name=viral
    
    input=files.renamed.slim
        
    # elapsed time: 21.051s
    # peak rss: 3.53 GB
    # file size: 9.16 GB
    # 
    memusg -t -s "kmcp compute -I $input -O genbank-$name-k21-n10 \
        -k 21 --seq-name-filter plasmid \
        --split-number 10 --split-overlap 150 \
        --log genbank-$name-k21-n10.log -j 32 --force" > genbank-viral.kmcp.s1.log2 2>&1
    
    # viral genomes are small:
    #   using small false positive rate: 0.05
    #   still using one hash function: 1
    #
    # elapsed time: 31.828s
    # peak rss: 3.36 GB
    # file size: 4.72 GB
    #
    memusg -t -s "kmcp index -I genbank-$name-k21-n10/ -O genbank-viral.kmcp \
        -j 32 -f 0.05 -n 1 -x 100K -8 1M \
        --log genbank-viral.kmcp.log --force" > genbank-viral.kmcp.s2.log2 2>&1
    
    # cp taxid and name mapping file to database directory
    cp taxid.map name.map genbank-viral.kmcp/
