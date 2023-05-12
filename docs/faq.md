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

### Can I create a database with custom genome collections for profiling?

Yes, you can use [taxonkit create-taxdump](https://bioinf.shenwei.me/taxonkit/usage/#create-taxdump)
to create NCBI-style taxdump files for profiling, which also generates a `taxid.map` file.

### What k-mer size should I use to build the database?

Multiple k-mer sizes are supported, but one value is good enough.

Bigger k-mer sizes bring high specificity at the cost of decrease
of sensitivity. `k = 21` is recommended for metagenomic profiling.

### How to add new genomes to the database?

KMCP builds database very fast,
you can either rebuilt the database after adding new genomes,
or create a separate database with the new genomes,
search against these databases, and [merge](/kmcp/usage/#merge) the results.

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
      
- Are you searching with metatranscriptomics data?
    - If yes, the search results would show that a huge number of reads from 16 rRNA genes have thousands of matches,
      therefore, writing results slow down the search.
      So these reads should be filtered out before the search using tools like https://github.com/hzi-bifo/RiboDetector.

### Can I run multiple KMCP processes in a machine?

Yes you can. But note that KMCP search is CPU- and RAM-intense. So please to **limit the number of CPUs cores to use for each process** with the falg `-j/--threads`, of which the default value is the available CPUs cores of the machine.

## Profiling

### Where is the taxid.map?

Each prebuilt database contains a `taxid.map` file in its directory.
You can concatenate them into a big one:

    $ cat gtdb.kmcp/taxid.map refseq-viral.kmcp/taxid.map refseq-fungi.kmcp/taxid.map > taxid.map
    
    $ head -n 5 taxid.map
    GCA_004025655.1 10243
    GCA_004025395.1 10243
    GCA_004025355.1 10243
    GCA_003971405.1 10243
    GCA_003971385.1 10243
    
Or set the options `-T/--taxid-map` multiple times:

    kmcp profile -T gtdb.kmcp/taxid.map -T refseq-viral.kmcp/taxid.map -T refseq-fungi.kmcp/taxid.map ...

For other custom genome collections, you can use
[taxonkit create-taxdump](https://bioinf.shenwei.me/taxonkit/usage/#create-taxdump)
to create NCBI-style taxdump files for custom taxonomy, e.g.,
[GTDB](https://github.com/shenwei356/gtdb-taxdump) and
[ICTV](https://github.com/shenwei356/ictv-taxdump), which also generates a `taxid.map` file.

### Unknown taxid?

> 19:54:54.632 [ERRO] unknown taxid for NZ_CP028116.1, please check taxid mapping file(s)

If the `kmcp profile` reports this, you may need to check if the taxid mapping file contain all the reference IDs.
And make sure the reference IDs match these in the database, the later ones are listed in: 

    $ head -n 5 $kmcp_db_dir/R001/__name_mapping.tsv
    NC_013654.1     NC_013654.1
    NC_000913.3     NC_000913.3
    NC_010655.1     NC_010655.1
    NC_012971.2     NC_012971.2
    NC_011750.1     NC_011750.1

There's another case: you used `--name-map` in `kmcp search`.
Please don't do this if you will use the search result for metagenomic
profiling which needs the original reference IDs. We add a note now:

```
-N, --name-map strings  ► Tabular two-column file(s) mapping reference IDs to user-defined
                        values. Don't use this if you will use the result for metagenomic
                        profiling which needs the original reference IDs.

```

### How to tune parts of options when using preset profiling modes?

<s>Sorry it's not supported due to the limitation of the command-line argument parsers.
You need explicitly set all relevant options of the mode.</s>

It's available since v0.8.2.

### Are the elements in the bloom filters uniformly distributed?

Some one asked me this when I were giving a talk of KMCP.
The answers is yes. I created a new command [`kmcp utils index-density`](https://bioinf.shenwei.me/kmcp/usage/#index-density) to plot the density of a index file.

Here's are some example outputs. The width is the number of bins and the height is the number of names (references or reference chunks).

1. v2023.05-genbank-viral-_block001.uniki (FPR of bloom filter: 0.3), only a part of image is shown.

    |file           |k  |canonical|num-hashes|num-sigs|num-names|
    |:--------------|:--|:--------|:---------|:-------|:--------|
    |_block001.uniki|21 |true     |1         |6395    |10400    |

    ![](v2023.05-genbank-viral-_block001.uniki.part.jpg)

1. v2023.05-refseq-fungi-_block001.uniki (FPR of bloom filter: 0.3)

    |file           |k  |canonical|num-hashes|num-sigs|num-names|
    |:--------------|:--|:--------|:---------|:-------|:--------|
    |_block001.uniki|21 |true     |1         |2089979 |160      |

    ![](v2023.05-refseq-fungi-_block001.uniki.jpg)

2. v2023.05-refseq-fungi-_block002.uniki (FPR of bloom filter: 0.3)

    |file           |k  |canonical|num-hashes|num-sigs|num-names|
    |:--------------|:--|:--------|:---------|:-------|:--------|
    |_block002.uniki|21 |true     |1         |2599648 |160      |

    ![](v2023.05-refseq-fungi-_block002.uniki.jpg)
