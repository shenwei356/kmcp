
## Database

Prebuilt databases

- DB for bacteria: [refseq-cami2-k21-n10.db.tar.gz](http://app.shenwei.me/data/tmp/refseq-cami2-k21-n10.db.tar.gz)
- DB for virus: [refseq-cami2-virus-k21-n5.db.tar.gz](http://app.shenwei.me/data/tmp/refseq-cami2-virus-k21-n5.db.tar.gz)
- TaxId mapping: [taxid.map](http://app.shenwei.me/data/tmp/taxid-virus.map), and [taxid-virus.map](http://app.shenwei.me/data/tmp/taxid-virus.map)
- [taxdump.tar.gz](http://app.shenwei.me/data/tmp/taxdump.tar.gz)

## Steps

    # search on multiple database
    file=sample_0.fq.gz
    sampleId=0
    for db in refseq-cami2-k21-n10.db refseq-cami2-virus-k21-n5.db; do
        kmcp search -j 40 -d $db $file -o $file.kmcp@$db.tsv.gz
    done
    
    # merge search result
    kmcp merge $file.kmcp@*.tsv.gz -o $file.kmcp.tsv.gz
    
    # profile
    kmcp profile -T taxid.map -T taxid-virus.map -X taxdump/ \
        $file.kmcp.tsv.gz \
        -o $file.kmcp.tsv.gz.k.profile \
        -C $file.kmcp.tsv.gz.cami.profile -s $sampleId        
    
