# Frequently Asked Questions 

## General

### How can I run KMCP on a computer without enough main memory?

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

## Database building

### What k-mer size should I use to build the database?

Multiple k-mer sizes are supported, but one value is good enough.

Bigger k-mer sizes bring high specificity in cost of decrease
of sensitivity. `k = 21` is recommended for metagenomic profiling.

### How to add new genomes to the database?

KMCP builds database very fast,
you can either rebuilt the database after adding new genomes,
or create a separate database with the new genomes,
search against these databases, and [merge](/usage/#merge) the results.

### Unexpected EOF error

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

## Searching

### Why are the CPU usages are very low, not 100%?

Please check the log (in terminal not the log file via the option `--log`),
it should report current searching speed, like:

    21:07:42.949 [INFO] reading sequence file: t_t3.fa.gz
    processed queries: 1253376, speed: 8.127 million queries per minute

If the speed is very slow. 

- Are you running KMCP in a computer cluster (HPC)? 
    - If yes, please switch on `-w/--load-whole-db`.
      Because the default database loading mode would be very slow for network-attached storage (NAS).

- Are the reference genomes are highly similar? E.g., tens of thousands of genomes of a same species?
    - If yes, check the search result to see if there are thousands of matches for a read.
      You may choose another graph-based sequence searching tool.

### Can I run multiple KMCP processes in a machine?

Yes you can. But note that KMCP search is CPU- and RAM-intense. So please to **limit the number of CPUs cores to use for each process** with the falg `-j/--threads`, of which the default value is the available CPUs cores of the machine.

## Profiling

### How to tune some options when using preset profiling modes?

Sorry it's not supported due to the limitation of the command-line argument parsers.
You need explicitly set all relevant options of the mode.