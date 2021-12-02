# Frequently Asked Questions 

## How can I run KMCP on a computer without enough main memory?

By default, `kmcp search` loads the whole database into main memory (RAM) for fast searching.
Optionally, the flag `--low-mem` can be set to avoid loading the whole database,
while it's much slower, >10X slower on SSD and should be much slower on HDD disks.

Another way is dividing the reference genomes into several parts
and building smaller databases for all parts, so that the biggest
database can be loaded into RAM. After performing database searching,
search results on all small databases can be merged with `kmcp merge`
for downstream analysis.

Buiding small databases can also accelerate searching on a *computer cluster*,
where every node searches a part of the database.

## What k-mer size should I use to build the database?

Multiple k-mer sizes are supported, but one value is good enough.

Bigger k-mer sizes bring high specificity in cost of decrease
of sensitivity. `k = 21` is recommended for metagenomic profiling.

## How to add new genomes to the database?

KMCP builds database very fast,
you can eigher rebuilt the database after adding new genomes,
or create a separate database with the new genomes,
search against these databases, and merge the results.

## unexpected EOF error

Some files could corrupt during downloading, we recommend checking
sequence file integrity using seqkit (`gzip -t` failed for some files in
my tests).

1. List corrupted files

        # corrupted files
        find $genomes -name "*.gz" \
            | rush 'seqkit seq -w 0 {} > /dev/null; if [ $? -ne 0 ]; then echo {}; fi' \
            > failed.txt

        # empty files
        find $genomes -name "*.gz" -size 0 >> failed.txt
    
2. Delete these files:

        cat failed.txt | rush '/bin/rm {}'

3. Redownload these files. For example, from `rush` cache file `download.rush` (list of succeed commands), where the commands are in one line.

        grep -f failed.txt download.rush \
            | sed 's/__CMD__//g' \
            | rush '{}'

    For `genome_updater`, URLs of genomes can be found in files like `2021-09-30_13-32-30_url_downloaded.txt`, 
    you can extract URLs using `grep -f failed.txt -v *url_downloaded.txt` or some other ways,
    and batch redownload them using `parallel`.
