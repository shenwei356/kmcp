# Changelog

### v0.7.0 - 2021-09

- `search`:
    - better log.

### v0.6.0 - 2021-08-13

- new command:
    - `merge`: merge search results from multiple databases.
- `compute`:
    - fix splitting very short genomes.
    - remove flag `-e/--exact-number`, making it default.
- `index`:
    - do not roundup sizes of indexes. The searching speed is not
      affected and even faster due to optimization of `search` command.
    - use three k-mers thresholds to control index file size.
    - better control of cocurrency number and better progress bar.
    - do not support RAMBO index anymore.
- `search`:
    - 1.37X speedup, and faster for database with two or more hash functions.
    - new flag `-S/--do-not-sort`.
- `profile`:
    - fix a nil pointer bug when no taxid mapping data given.
    - fix number of ureads.
    - new flag `-m/--keep-main-matches` and `--max-score-gap`

### v0.5.0 - 2021-06-24

- `compute`:
    - support multiple sizes of k-mer.
    - fix bug of `--by-seq`.
    - more log.
- `index`:
    - default block size is computed by `-j/--threads` instead of number of CPUs.
- `search`:
    - show real-time processing speed.
    - new flag `-g/--query-whole-file`.
    - new flag `-u/--kmer-dedup-threshold`.
    - new flag `-m/--min-query-len`.
    - increase speed for database with mulitple hashes. 
- `profile`:
    - better decision of the existence of a reference.
    - new flag `-B/--binning-result` for output reads binning result.
    - new flag `-m/--norm-abund`.
    
### v0.4.0 - 2021-04-08

- new command:
    - `profile` for generating taxonomic profile from search result.
- `compute`:
    - new flag `-B/--seq-name-filter` for filtering out unwanted sequences like plasmid.
    - new flag `-N/--ref-name-regexp` for extracting reference name from sequence file.
- `search`:
    - change default threshold value.
    - new flag `-n/--keep-top-scores` for keeping matches with the top N score.

### v0.3.0 - 2021-03-16

- use `--quiet` to replace `--verbose`, making printing log info default.
- `search`: 
    - fix computing intersetion between repeats.
    - fix closing mmap on Windows.
    - change output format and add Jaccard Index.
    - speedup by parallelizing name mapping and database closing.
    - flush result immediately.
    - keep the output order by default
- `compute`: change default file regexp for matching `.fna` files.
- `autocomplete`: support bash, zsh, fish, powershell.

### v0.2.1 - 2020-12-31

- `index`: reduce memory occupation.
  
### v0.2.0 - 2020-12-30

- Add support of RAMBO like indexing.
- Limit to only one input database.
- Change output format.

### v0.1.0 - 2020-xx-xx

- First release with basic function.
