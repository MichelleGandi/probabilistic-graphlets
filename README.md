# Calibrated Probabilistic Graphlets

Do graphlet features from a **calibrated probabilistic** gene co-expression network
predict gene function better than features from **thresholded (binary)** networks made
from the same data?

I use single-cell RNA-seq of CD14+ monocytes from 8 donors (Kang et al. 2018, NCBI GEO
GSE96583). For each pair of genes I estimate the probability of a co-expression edge by
resampling cells within donors and calibrating it against replication in independent
donors. I then compute probabilistic graphlet degree vectors (Doria-Belenguer et al.
2020) and compare them with binary graphlet features, an expected-degree baseline and an
expression baseline for predicting Gene Ontology biological process terms.

Author: Michelle Gandidzanwa ([ORCID 0009-0004-1593-1446](https://orcid.org/0009-0004-1593-1446))

## Related Publications

This method builds on Milenković and Pržulj (2008), *Cancer Informatics* 6:257–273
(graphlet degree signatures) and Doria-Belenguer et al. (2020), *Bioinformatics*
36(Suppl 2):i804–i812, <https://doi.org/10.1093/bioinformatics/btaa812> (probabilistic
graphlets). The data are from Kang et al. (2018), *Nature Biotechnology* 36:89–94,
<https://doi.org/10.1038/nbt.4042>.

## Reproduce the results

Requirements: [Pixi](https://pixi.prefix.dev/latest/) and an internet connection, on
Linux or macOS. Windows is not supported because the graphlet counter `orbit-count` has
no Windows build (see `docs/decisions.md`) on Windows, use WSL2.

```bash
git clone https://github.com/MichelleGandi/probabilistic-graphlets.git
cd probabilistic-graphlets
pixi run reproduce
```

`pixi run reproduce` installs the exact software versions in `pixi.lock`, downloads the
raw data and checks its checksums, runs every analysis step and renders the report. All
settings are in `config/config.yaml`. A GitHub Action
(`.github/workflows/reproduce.yml`) runs the same command on a fresh machine on every
push.

Other entry points: `pixi task list` lists every task, one per analysis step (for
example `pixi run preprocess`), `pixi run test` checks the graphlet counting against
exact enumeration.

### Container image

The software environment can also be run as a Docker image built from `pixi.lock`
(`containers/Dockerfile`), for example on Windows. From the project folder:

```bash
docker build --platform linux/amd64 -f containers/Dockerfile \
  --build-arg ENV_NAME=default -t probabilistic-graphlets:0.1.0 .
docker run --rm --platform linux/amd64 -v "$PWD":/work -w /work probabilistic-graphlets:0.1.0 \
  bash -c 'bash scripts/download.sh && for s in 01_verify_data 02_preprocess 03_edge_confidence 04_calibrate 05_networks 06_graphlets 07_go_labels 08_predict 09_results; do python scripts/$s.py || exit 1; done'
```

The image contains the software only, the project folder is in `/work` and the
steps run in the same order as `pixi run reproduce` (without the report). Image built on
2026-11-27 (linux/amd64): ID
`sha256:e6f7afdfad84729e14f1e924d4643dfe6e87649abf391763aa98d68a93463138`. Run in the
container, the analysis gave the same results.

## Hardware

Measured on a MacBook Air (10 cores): about 1.2 GB of
disk for the software environment, 99 MB for the raw data and 2 MB for the results and a
full run took about 2 minute with at most about 0.9 GB of memory. On GitHub's standard Linux runner it takes about 3 minutes.

## Outputs

| Path | Content |
|------|---------|
| `reports/report.html` | The report: question, data, methods, calibration, networks, prediction results, discussion, software versions |
| `results/tables/` | Every table, numbered by the step that writes it |
| `results/figures/` | Calibration, network size and prediction figures |
| `data/processed/` | Intermediate files (selected cells, edge scores, networks, graphlet counts) |


## Data

| Data | Source | Size | License / terms |
|------|--------|------|-----------------|
| Single-cell RNA-seq of PBMCs, batch 2: 8 donors with lupus, control and interferon-β-stimulated (10x Chromium) | NCBI GEO [GSE96583](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583) (samples GSM2560248, GSM2560249); Kang et al. 2018 | 56 MB (6 files) | No restrictions on use from GEO, submitters may hold rights in their data. Cite Kang et al. 2018. |
| Gene Ontology annotations and ontology, release 2026-08-05 | <https://release.geneontology.org/2026-08-05/> | 43 MB (2 files) | CC BY 4.0 |

`bash scripts/download.sh` (run by `pixi run reproduce`) downloads every file listed in
`metadata/sources.tsv` into `data/raw/`, makes the files read-only and checks them
against `metadata/checksums.sha256`. The analysis uses the 2,785 unstimulated singlet
CD14+ monocytes (126–779 per donor). Every metadata file is explained in
`metadata/README.md`.

### Data availability statement

The single-cell RNA-seq data analysed here are publicly available from NCBI GEO under
accession GSE96583 (Kang et al. 2018), this project uses samples GSM2560248 and
GSM2560249 and the batch-2 gene and cell metadata files. Gene Ontology annotations and
ontology are from the Gene Ontology release of 2026-08-05 (CC BY 4.0). The exact URLs and
SHA-256 checksums of all raw files are in `metadata/sources.tsv` and
`metadata/checksums.sha256` and `scripts/download.sh` retrieves them. The data are not
stored in this repository. All derived results are regenerated from the raw data by the
code in this repository.

### Responsible use of human data (CARE)

The expression data come from human donors and were deposited by the original study in
de-identified, processed form for research. This project uses only those processed
counts and the published cell annotations, reports results at the level of genes
and cell type and makes no attempt to identify any donor. Anyone reusing the data should get it from GEO and cite the original study.

## Repository contents

See `docs/structure.md` (each folder also has its own README). Analysis decisions and
their reasons are in `docs/decisions.md`, the exact software versions used are in
`docs/software-versions.txt` and `pixi.lock`.

## Cite

See `CITATION.cff` (GitHub: "Cite this repository").

## License

- Code: MIT (`LICENSE`).
- Report text, figures and documentation: CC BY 4.0.
- Data: not covered by these licenses, they keep their providers' terms (see Data).
- `reports/vancouver.csl` (citation style, from the Citation Style Language project):
  CC BY-SA 3.0.
- `orbit-count`, installed by Pixi from PyPI: GPL-3.0 (not part of this repository).

## Contact

Questions and problems: open an issue on GitHub or email
Michelle.Gandidzanwa@mbzuai.ac.ae. Contributions: see `CONTRIBUTING.md`. Everyone
taking part follows `CODE_OF_CONDUCT.md`.

---
Project layout generated from the Reproducible and Open Science project template
(<https://github.com/khan-lab/cbio7501-f26/tree/main/templates/cookiecutter-repro-project>).
