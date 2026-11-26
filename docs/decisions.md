# Decisions

Choices made in this project and why, newest last. Settings themselves live in
`scripts/p36/config.py` (later `config/config.yaml`).

## Data

- **2026-11-27. One dataset: GEO GSE96583, batch 2.** Single-cell data from 8 donors in one
  10x experiment, so cells can be resampled within donors and edge replication tested in
  independent donors. Gene Ontology is used as reference labels, not as a second dataset.
- **2026-11-27. Download only the batch-2 sample files, not `GSE96583_RAW.tar`.** The tar
  holds every batch; the two batch-2 samples (GSM2560248, GSM2560249) and the two batch-2
  metadata files are all the analysis reads.
- **2026-11-27. Keep the stimulated sample's files (GSM2560249).** Only control cells are
  analysed, but the loading code reads both samples and filters by condition, so the
  files are kept for now (28 MB).
- **2026-11-27. Pin the Gene Ontology release to 2026-08-05** instead of `current`. The
  current files change every month (and their names and folders changed: the human file
  is now `annotations/gaf/HUMAN-uniprot.gaf.gz`), which would change the labels and the
  results.
- **2026-11-27. Own SHA-256 checksums.** GEO and GO publish no checksums for these files,
  so `metadata/checksums.sha256` is the reference that anyone can check against.
- **2026-11-27. Cell type: CD14+ monocytes, control condition, singlets only:** 2,785 cells,
  126–779 per donor.
