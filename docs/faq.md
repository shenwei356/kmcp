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

## What k-mer size should I use to build the database?

Multiple k-mer sizes are supported, but one value is good enough.

Bigger k-mer sizes bring high specificity in cost of decrease
of sensitivity. `k = 21` is recommended for metagenomic profiling.

## How to add new genomes to the database?
