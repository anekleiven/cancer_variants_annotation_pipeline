## Variants — annotation_pipeline

This repository contains a small annotation pipeline for cancer variants. It performs hotspot annotation, maps UniProt identifiers, annotates protein domains (Pfam), and functional sites.

This README explains the repository layout, required input data, how to run the scripts and expected outputs.

## Repository layout

- `scripts/` — pipeline scripts (Python and R)
  - `01_hotspots_long.py`
  - `02_annotate_hotspots.py`
  - `03_add_uniprot_mapping.R`
  - `04_annotate_protein_domains.py`
  - `05_annotate_functional_sites.py`

## Quick summary

The pipeline takes variant lists (expected as TSV/MAF-like files with gene and protein coordinates) and annotates them using hotspots, Pfam domains and curated functional sites. It does not download external databases automatically — you must provide the necessary reference data files in `data/` before running.

## Data requirements (you must provide)

Place the following files in the `data/` directory. Filenames used in the repo are shown; adapt the scripts or rename your files to match these names if needed.

- `hotspots.txt` — hotspot definitions for proteins
- `Pfam-A.regions.tsv.gz` — Pfam domains and regions (download from Pfam)
- `Pfam-A.clans.tsv` — Pfam clans and family data (download from Pfam).

Where to obtain common resources:
- Cancer hotspots: https://www.cancerhotspots.org/#/home
- Pfam (family definitions, clans): https://pfam.xfam.org/
- UniProt mappings and accession lists: https://www.uniprot.org/

Notes:
- The scripts currently expect these files to be present in `data/`. If your filenames differ, either rename them or update the script arguments.

## Dependencies / Environment

- Python: 3.8+ recommended. The Python scripts likely use common libraries such as `pandas` — check the top of each script for imports. 
- R: R 4.x if you run the R script (`03_add_uniprot_mapping.R`). Ensure required R packages used in the script are installed.

## How to run (high-level)

Run the pipeline scripts in order. The scripts are located in `scripts/` and may accept command-line arguments — inspect each script's top-level docstring or `--help` if implemented.

Example (run sequentially):

```powershell
python scripts\01_hotspots_long.py
python scripts\02_annotate_hotspots.py
Rscript scripts\03_add_uniprot_mapping.R
python scripts\04_annotate_protein_domains.py
python scripts\05_annotate_functional_sites.py
```

Notes:
- These are example invocations. Some scripts might require input paths or flags. Inspect each script for usage details and required inputs.

## Expected outputs (examples)

- `output/annotated_with_hotspots.maf` — annotated MAF (example)
- `output/variants_with_pfam_domains.tsv` — variants annotated with Pfam domains

The exact output filenames are created by the scripts. 

## Contact / author

Contact me by e-mail (ane.kleiven@gmail.com) if you have any questions about how to use the pipeline. 

