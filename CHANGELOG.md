# Changelog

## v0.4.1 - 2021-04

- `search`:
    - show real-time processing speed.
    - new flag `-g/--query-whole-file`.
- `compute`:
    - fix bug of `--by-seq`.
    - more log.

## v0.4.0 - 2021-04-08

- new command:
    - `profile` for generating taxonomic profile from search result.
- `compute`:
    - new flag `-B/--seq-name-filter` for filtering out unwanted sequences like plasmid.
    - new flag `-N/--ref-name-regexp` for extracting reference name from sequence file.
- `search`:
    - change default threshold value.
    - new flag `-n/--keep-top-scores` for keeping matches with the top N score.

## v0.3.0 - 2021-03-16

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

## v0.2.1 - 2020-12-31

- `index`: reduce memory occupation.
  
## v0.2.0 - 2020-12-30

- Add support of RAMBO like indexing.
- Limit to only one input database.
- Change output format.

## v0.1.0 - 2020-xx-xx

- First release with basic function.
