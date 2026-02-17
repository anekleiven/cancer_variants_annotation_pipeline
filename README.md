# Cancer Variants Annotation Pipeline 

Annotation pipeline that annotates somatic cancer variants with classification evidence such as cancer hotspots, protein domains, functional sites, pathogenic germline variant distances and MAVEs. 

## Requirements
- Python 3.10+
- R 4.2+

## Setup Instructions 

1. **Create Virtual Environment:**
`python -m venv .venv`
`. .venv/bin/activate`

2. **Install Python Requirements:**
`pip install -r requirements.txt`

3. **Install R requirements:**
`Rscript install_deps.R`

## Pipeline Order
`01_hotspots_long.py` - converts hotspot data to long format 

`02_annotate_hotspots.py` - annotates cancer variants with cancer hotspots 

`03_add_uniprot_mapping.R` - annotates cancer variants with uniprot accession numbers 

`04_annotate_protein_domains.py` - annotates cancer variants with protein domains 

`05_annotate_functional_sites.py` - annotates cancer variants with functional sites 

`06_annotate_germline_proximity.py` - annotates cancer variants with germline variant distances 

`07_mave_tsv_to_vcf.py` - converts tsv file to vcf file for assembly liftover and MAVE annotation 

`08_annotate_maves.py` - annotates original variant file with MAVE data 

