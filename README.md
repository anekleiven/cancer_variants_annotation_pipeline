# Cancer Variants Annotation Pipeline 🧬🧬

Annotation pipeline for somatic cancer variants. 

Annotates variants with classification evidence:
- Cancer hotspots from cancerhotspots.org
- Protein domains from Pfam
- Functional sites from UniProt
- Pathogenic germline variant distances from ClinVar 
- Functional data from MaveDB

## Data Flow 📋
The pipeline expects a variant file in `.tsv` format as starting input. 
Each script (01-07) appends specific annotations to the dataset, resulting in a 
final master file used for downstream analysis.

## Requirements 💻
- Python 3.10+
- R 4.2+

## External Data Requirements 
Before running the pipeline, the following files must be downloaded manually:
| File | Source | Used in Script |
| :--- | :--- | :--- |
| `hotspots.txt` | [cancerhotspots.org](https://www.cancerhotspots.org) | `01A_hotspots_long.py` |
| `Pfam-A.regions.tsv.gz` | [EMBL-EPI FTP](https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.2/) | `03_annotate_protein_domains.py` |
| `Pfam-A.clans.tsv.gz` | [EMBL-EPI FTP](https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.2/) | `03_annotate_protein_domains.py` |
| `variant_summary.txt.gz` | [NCBI ClinVar FTP](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/) | `05_annotate_germline_proximity.py` |


## External Dependencies 
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

