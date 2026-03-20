# Cancer Variants Annotation Pipeline 🧬🧬

Annotation pipeline for somatic cancer variants. 
Annotates variants with classification evidence such as:
Cancer hotspots, protein domains, functional sites, pathogenic germline variant distances and MAVEs 

## Requirements 💻
- Python 3.10+
- R 4.2+

## Setup Instructions 🔧

1. **Create Virtual Environment:**
`python -m venv .venv`

2. **Activate Virtual Environment:**
`. .venv/bin/activate`

3. **Install Python Requirements:**
`pip install -r requirements.txt`

4. **Install R requirements:**
`Rscript install_deps.R`

## Pipeline Order 

`01A_hotspots_long.py`: Transforms hotspot.txt to long format (one row per variant). 

`01B_annotate_hotspots.py`: Annotates somatic cancer variants with hotspot data from cancerhotspots.org. 

`02_add_uniprot_mapping.R`: Annotates somatic cancer variants with UniProt accession numbers using the geneOncoX R package. 

`03_annotate_protein_domains.py`: Annotates somatic cancer variants with protein domains from Pfam. 

`04_annotate_functional_sites.py`: Annotates somatic cancer variants with functional sites from UniProt. 

`05_annotate_germline_proximity.py`: Annotates somatic cancer variants with the distance to pathogenic germline variants from ClinVar. 

`06_mave_tsv_to_vcf.py`: Converts .tsv files to .vcf files for assembly liftover and MAVE annotation. 

`07_annotate_maves.py`: Annotates somatic cancer variants with functional data from MaveDB. 

