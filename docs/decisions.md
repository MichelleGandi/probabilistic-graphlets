# Decisions

Choices made in this project and why, newest last. Settings themselves live in
`config/config.yaml`.

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

## Software environment

- **2026-11-27. Platforms: linux-64, osx-arm64, osx-64; not win-64.** Graphlet counting uses
  `orbit-count` (ORCA), which PyPI ships as built wheels for Linux and macOS only, with
  no Windows wheel and no source package, so the environment cannot be built on Windows.
  Windows users can run the project with the container image (lesson 08) or WSL.
- **2026-11-27. `orbit-count` comes from PyPI** (`[pypi-dependencies]` in `pixi.toml`): it
  is not packaged on conda-forge or bioconda. Its source is licensed GPL-3.0; this project
  imports it but does not copy or distribute its code.
- **2026-11-27. One task per analysis step, chained with `depends-on`;** `pixi run
  reproduce` runs them all from the download onwards.
- **2026-11-27. No separate test mode.** The full analysis on the real data took 2 min 34 s
  on a MacBook Air (10 cores) including the checksum check, so the GitHub *Reproduce*
  check and a partner can run the real analysis rather than a synthetic one.

## Method

- **2026-11-27. Pairs never co-expressed in any bootstrap get edge probability 0**
  (`networks: zero_if_never_stable` in `config/config.yaml`). Isotonic calibration maps stability 0
  to the background replication rate (2–9 %), so in the first run every one of the
  1,263,255 gene pairs had p ≥ 0.01 and the probabilistic network was complete
  (density 1.0). A pair that never passes the co-expression threshold in any of the 100
  resamples has no evidence for an edge, so its probability is set to 0. Results before
  this change: median AUPR 0.056 (probabilistic GDV) vs 0.057 (binary, matched density).
  After it: 390,929 candidate edges (density 0.31, 21,021 expected edges); median AUPR
  0.076 (probabilistic GDV) vs 0.061 (binary, matched density), 0.054 (expected degree)
  and 0.084 (expression only). Because this change was made after seeing the first
  results, both runs are reported.

## Automation

- **2026-11-27. All analysis settings in `config/config.yaml`.** They used to be Python
  constants in `scripts/p36/config.py`; that module now only reads the YAML file and works
  out the folder paths, so the scripts are unchanged. Folder paths stay in code because
  they are fixed by the project layout, not analysis choices.
- **2026-11-27. `curl` and `coreutils` (for `sha256sum`, `md5sum`) come from the Pixi
  environment.** `scripts/download.sh` first used the system's `shasum`, which macOS and
  GitHub's runner have but the container image (plain Ubuntu) does not. With these tools in
  `pixi.toml`, every tool the analysis uses comes from the environment, and the download
  script works the same on macOS, Linux and in the container.

## Container

- **2026-11-27. Container image for linux/amd64, built from `pixi.lock`** with the
  template's `containers/Dockerfile`. The image holds the software only; the project is
  mounted at `/work`. Run on an Apple Silicon Mac under emulation, the analysis gave
  results identical to the native macOS run (median AUPR 0.076, 0.061, 0.084), so the
  results do not depend on the operating system or processor.
