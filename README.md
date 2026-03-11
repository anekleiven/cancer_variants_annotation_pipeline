# Cancer Variants Annotation Pipeline 👩🏽‍🔬

Annotation pipeline that annotates somatic cancer variants with classification evidence such as cancer hotspots, protein domains, functional sites, pathogenic germline variant distances and MAVEs 🧬

## Requirements 💻
- Python 3.10+
- R 4.2+

## Setup Instructions 🔧

1. **Create Virtual Environment:**
`python -m venv .venv`
`. .venv/bin/activate`

2. **Install Python Requirements:**
`pip install -r requirements.txt`

3. **Install R requirements:**
`Rscript install_deps.R`

## Pipeline Order 🤓🤓

`01A_hotspots_long.py`: transform hotspots file to long format (one row per variant) 

`01B_annotate_hotspots.py`: annotate somatic cancer variants with cancer hotspots 

`02_add_uniprot_mapping.R`: annotate somatic cancer variants with UniProt accession numbers 

`03_annotate_protein_domains.py`: annotate somatic cancer variants with protein domains 

`04_annotate_functional_sites.py`: annotate somatic cancer variants with functional sites 

`05_annotate_germline_proximity.py`: annotate somatic cancer variants with germline variant distances 

`06_mave_tsv_to_vcf.py`: convert tsv file to vcf file for assembly liftover and MAVE annotation 

`07_annotate_maves.py`: annotate somatic cancer variants with MAVE data 

