# metadata

Files that describe the raw data. They are small, written by hand and committed.

## sources.tsv

Where each raw file comes from. Columns:

| Column | Meaning |
|--------|---------|
| `name` | short label for the file |
| `url` | download link (HTTPS) |
| `md5` | checksum published by the provider, if any; `NA` when the provider publishes none (true for these GEO and GO files), in which case `checksums.sha256` is the reference |
| `unpack` | `none` = keep the file as downloaded (the code reads the `.gz` files directly) |
| `dest` | folder in this project where the file goes |

Sources (checked 2026-11-27):

- NCBI GEO series **GSE96583**, batch 2 (Kang et al. 2018, *Nat Biotechnol* 36:89–94,
  <https://doi.org/10.1038/nbt.4042>): samples GSM2560248 and GSM2560249, plus the
  batch-2 gene list and cell metadata. GEO places no restrictions on use; submitters may
  hold rights in their data. Not redistributed here.
- Gene Ontology release **2026-08-05** (<https://release.geneontology.org/2026-08-05/>):
  `HUMAN-uniprot.gaf.gz` (annotations) and `go-basic.obo` (ontology). CC BY 4.0.
  A dated release is used, not `current`, so everyone gets the same files.

`scripts/download.sh` downloads every file in this table, makes it read-only and checks it.

## checksums.sha256

SHA-256 checksum of every raw file. Check with
`shasum -a 256 -c metadata/checksums.sha256` (or `sha256sum -c` on Linux);
`scripts/download.sh` does this at the end.

## samples.tsv

One row per GEO sample of batch 2.

| Column | Meaning |
|--------|---------|
| `sample` | GEO sample accession (GSM…) |
| `series` | GEO series accession |
| `batch` | batch in the original study (batch 2 = 8 donors with lupus, 10x Chromium) |
| `condition` | `ctrl` = unstimulated PBMCs; `stim` = stimulated with interferon-β |
| `used_in_analysis` | `yes` if cells from this sample enter the analysis (only `ctrl` cells are used; the `stim` files are read but filtered out) |
| `matrix` | count matrix (Matrix Market, genes × cell barcodes), relative to the project folder |
| `barcodes` | cell barcodes for the matrix columns |

## Cell metadata (`GSE96583_batch2.total.tsne.df.tsv.gz`)

One row per cell barcode. The columns used by the code (named under `dataset: metadata_columns` in `config/config.yaml`):

| Column | Meaning |
|--------|---------|
| `ind` | donor ID (8 donors) |
| `stim` | condition (`ctrl` or `stim`) |
| `cell` | cell-type label from the original study (we use `CD14+ Monocytes`) |
| `multiplets` | `singlet` or doublet call from demuxlet (only singlets are used) |
