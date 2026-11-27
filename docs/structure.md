# Project structure

| Path | Purpose |
| ------ | --------- |
| `data/raw/` | Downloaded data, exactly as received (by `scripts/download.sh`) |
| `data/processed/` | Data derived from the raw data |
| `metadata/` | Sample sheet, data sources, checksums, data dictionary |
| `scripts/` | Code to run the analysis (Python, shell); `scripts/p36/` holds shared functions |
| `tests/` | Check of the graphlet counting (`pixi run test`) |
| `config/` | `config.yaml`: every analysis setting |
| `workflow/` | Snakemake and Nextflow pipelines |
| `notebooks/` | Jupyter notebooks for exploration |
| `reports/` | The Quarto report (`report.qmd`), references and citation style |
| `results/` | Tables and figures made by the code |
| `docs/` | Decisions, structure, software versions |
| `envs/` | Conda `environment.yml` exports, for people without Pixi |
| `containers/` | Dockerfile and Apptainer recipe |

Each folder has its own `README.md` with more detail.
