# Cancer Variants Annotation Pipeline 🧬🧬

Annotation pipeline for somatic cancer variants. 
Annotates variants with classification evidence such as:
Cancer hotspots, protein domains, functional sites, pathogenic germline variant distances and MAVEs. 

Required files and downloads are described in each script. 

## Data Flow 📋
The pipeline expects a variant file in `.tsv` format as starting input. 
Each script (01-07) appends specific annotations to the dataset, resulting in a 
final master file used for downstream analysis.

## Requirements 💻
- Python 3.10+
- R 4.2+

## External dependencies 
- CrossMap (v.0.7.3): Used for GRCh37 to GRCh38 liftover.
- Docker: Required to run the Ensembl VEP container. 
- Ensembl VEP (v.115): With the MaveDB Plugin and GRCh38 cache. 

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

`---Manual step: Assembly Liftover---`: Run CrossMap in terminal to convert coordinates from GRCh37 to GRCh38. 

`06_mave_tsv_to_vcf.py`: Converts .tsv files to .vcf files for assembly liftover and MAVE annotation. 

`---Manual step: VEP annotation---`: Run Ensembl VEP via Docker to annotate .VCF with MaveDB scores. 

`07_annotate_maves.py`: Annotates somatic cancer variants with functional data from MaveDB. 

