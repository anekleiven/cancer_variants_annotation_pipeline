# Cancer Variants Annotation Pipeline 

## Setup Instructions 

1. **Create Virtual Environment:**
`python -m venv .venv`
`. .venv/bin/activate`

2. **Install Python Requirements:**
`pip install -r requirements.txt`

3. **Install R requirements:**
`Rscript install_deps.R`

## Pipeline Order
`01_hotspots_long.py`
`02_annotate_hotspots.py`
`03_add_uniprot_mapping.R`
`04_annotate_protein_domains.py`
`05_annotate_functional_sites.py`
`06_annotate_germline_proximity.py`

