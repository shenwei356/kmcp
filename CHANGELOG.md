# Changelog

## v0.3.0 - 2021-01-xx

- use `--quiet` to replace `--quiet`, making printing log info default.
- `search`: 
  - fix computing intersetion between repeats.
  - fix closing mmap on Windows.
  - change output format and add Jaccard Index.  
  - speedup by parallelizing name mapping and database closing.
  - flush result immediately.
- `compute`: change default file regexp for matching `.fna` files.

## v0.2.1 - 2020-12-31

- `index`: reduce memory occupation.
  
## v0.2.0 - 2020-12-30

- Add support of RAMBO like indexing$$
- Limit to only one input database.
- Change output format.

## v0.1.0 - 2020-xx-xx

- First release with basic function.
