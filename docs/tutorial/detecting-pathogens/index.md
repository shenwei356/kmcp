# Detecting specific pathogens in sequencing data

KMCP v0.9.3 or later versions is needed, which fixed a bug in chunk computation when splitting circular genomes.

## Data

1. Reference genomes (taken HBV virus for example).

        wget https://hbvdb.lyon.inserm.fr/data/nucleic/fasta/all_Genomes.fas

        # removing duplicated sequences
        seqkit rmdup -s all_Genomes.fas -o refs.fa

        seqkit stats refs.fa
        file    format  type  num_seqs     sum_len  min_len  avg_len  max_len
        refs.fa  FASTA   DNA      9,248  29,673,289    3,182  3,208.6    3,275


        # --- only for HBV ref genomes -----

        # only keep at most 10 genomes for each genotype
        seqkit fx2tab -i refs.fa \
            | csvtk mutate -Ht -p '_([A-Z\-]+)$' \
            | csvtk uniq -Ht -f 4 -n 10 \
            | seqkit tab2fx -o refs.slim.fa

        file          format  type  num_seqs  sum_len  min_len  avg_len  max_len
        refs.slim.fa  FASTA   DNA        221  709,855    3,182    3,212    3,257

2. NCBI taxdump files, only needed for multiple reference genomes.

        wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

        mkdir -p taxdump
        tar -zxvf taxdump.tar.gz -C taxdump/

## Preprocessing reference genomes

Saving genomes in separated files, with the genome ID as the file name.

    # 1. split into separated files, each with one genome
    seqkit split2 -s 1 refs.slim.fa -O refs --force

    # 2. rename files
    # rush: https://github.com/shenwei356/rush
    find refs/ -name "*.fa" \
        | rush 'id=$(seqkit seq -ni {} | head -n 1); mv {} "refs/$id.fa"'

    tree refs/ | head -n 8
    refs/
    ├── gnl|hbvnuc|AB010290_FT00000_P-B.fa
    ├── gnl|hbvnuc|AB014392_FT00000_P-C.fa
    ├── gnl|hbvnuc|AB033558_FT00000_P-RF-DE.fa
    ├── gnl|hbvnuc|AB059660_FT00000_P-H.fa
    ├── gnl|hbvnuc|AB059661_FT00000_P-H.fa
    ├── gnl|hbvnuc|AB064310_FT00000_P-G.fa
    ├── gnl|hbvnuc|AB064311_FT00000_P-G.fa

## Creating a KMCP database

Spliting each genome into 10 chunks with 150-bp overlaps, and computing k-mers.

    # the default k-mer size is 21.
    # but for many highly similar reference genomes, you may need a larger k value.
    kmcp compute --circular -k 31 -n 10 -l 150 -I refs/ -O refs.tmp --force

Indexing the k-mers, using a small false-positive rate for small genomes.

    # you can add the flag --dry-run to see the size of database in advance
    kmcp index -f 0.001 -I refs.tmp/ -O refs.kmcp --force

Creating the taxid mapping file (only needed for multiple reference genomes).

    # 10407 is the taxid of Hepatitis B virus species
    cut -f 1 refs.kmcp/R001/__name_mapping.tsv \
        | awk '{print $0"\t10407"}' \
        > taxid.map

    cp taxid.map refs.kmcp/

## Searching reads against the KMCP database

    kmcp search -d refs.kmcp/ sample_1.fq.gz sample_2.fq.gz -o sample.kmcp.tsv.gz

## Profiling

- For the database with a single reference genome.

        # the preset profiling mode 1 with a higher sensitivity is used here.
        kmcp profile --level strain -m 1 \
            sample.kmcp.tsv.gz -o sample.k.profile -C sample.c.profile -M sample.m.profile

- For the database with multiple reference genomes.

        # the preset profiling mode 1 with a higher sensitivity is used here.
        kmcp profile -X taxdump/ -T taxid.map -m 1 \
            sample.kmcp.tsv.gz -o sample.k.profile -C sample.c.profile -M sample.m.profile
