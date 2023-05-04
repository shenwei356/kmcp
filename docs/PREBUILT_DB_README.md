# KMCP databases v2023.05

Source code: https://github.com/shenwei356/kmcp
Documents  : https://bioinf.shenwei.me/kmcp/database/

## Data source

Access date: 2023-05-03

|DB                  |source      |#NCBI-species| #assemblies|db-parameters                     |size     |
|:-------------------|:-----------|:------------|:-----------|:---------------------------------|:--------|
|Bacteria and Archaea|GTDB r214   |34395+       |85205       |k=21, chunks=10; fpr=0.3, hashes=1|50+49 GB |
|Fungi               |Refseq r217 |491          |496         |k=21, chunks=10; fpr=0.3, hashes=1|5.5 GB   |
|Viruses             |GenBank r255|26680        |33479       |k=21, chunks=5; fpr=0.05, hashes=1|5.9 GB   |

Taxdump files:

- NCBI
    - file: taxdump.tar.gz
    - version: taxdmp_2023-05-01 https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2023-05-01.zip
- GTDB+NCBI
    - file: taxdump.gtdb+ncbi.tar.gz
    - version: GTDB r214 + NCBI taxdmp_2023-05-01

## Files

Please download files according to your purposes:

- For metagenomic profiling: `metagenomic-profiling`
- For genome similarity estimation: `genome-search`
- For developers: `genomes`

File tree:

    v2023.05/
    ├── genomes                                            Genomes used to build the databases
    │   ├── genbank-viral.tar
    │   ├── genbank-viral.tar.md5.txt
    │   ├── genbank-viral.taxid.map.stats.tsv              Complete lineages (NCBI Taxonomy) of all TaxIds and the number of genomes
    │   ├── gtdb.taxid.map.stats.tsv
    │   ├── gtdb.txt
    │   ├── refseq-fungi.tar
    │   ├── refseq-fungi.tar.md5.txt
    │   └── refseq-fungi.taxid.map.stats.tsv
    ├── genome-search                                      KMCP databases for genome similarity estimation
    │   ├── genbank-viral.minhash.kmcp.tar.gz
    │   ├── genbank-viral.minhash.kmcp.tar.gz.md5.txt
    │   ├── gtdb.minhash.kmcp.tar.gz
    │   ├── gtdb.minhash.kmcp.tar.gz.md5.txt
    │   ├── name.map                                       Mapping assembly accessions to genome names (may be the name of the larget contig)
    │   ├── refseq-fungi.minhash.kmcp.tar.gz
    │   └── refseq-fungi.minhash.kmcp.tar.gz.md5.txt
    └── metagenomic-profiling                              KMCP databases for metagenomic profiling
        ├── genbank-viral.kmcp.tar.gz
        ├── genbank-viral.kmcp.tar.gz.md5.txt
        ├── gtdb.part_1.kmcp.tar.gz                        GTDB representative genomes are split into 2 parts before building the database
        ├── gtdb.part_1.kmcp.tar.gz.md5.txt
        ├── gtdb.part_2.kmcp.tar.gz
        ├── gtdb.part_2.kmcp.tar.gz.md5.txt
        ├── refseq-fungi.kmcp.tar.gz
        ├── refseq-fungi.kmcp.tar.gz.md5.txt
        ├── taxdump.gtdb+ncbi.tar.gz                       Taxdump files (Combining GTDB and NCBI Taxonomy)
        ├── taxdump.gtdb+ncbi.tar.gz.md5.txt
        ├── taxdump.tar.gz                                 Taxdump files (NCBI Taxonomy, 2023-05-01)
        ├── taxdump.tar.gz.md5.txt
        ├── taxid.gtdb+ncbi.map                            Mapping assembly accessions to TaxId (Combining GTDB and NCBI Taxonomy)
        └── taxid.map                                      Mapping assembly accessions to TaxId (NCBI Taxonomy)
