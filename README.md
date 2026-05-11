# Cancer Variants Annotation Pipeline 🧬🧬

Annotation pipeline for somatic cancer variants. 

Annotates variants with classification evidence:
- Cancer hotspots from cancerhotspots.org
- Protein domains from Pfam
- Functional sites from UniProt
- Pathogenic germline variant distances from ClinVar 
- Functional data from MaveDB
- Tumor suppressor gene and oncogene annotations from the Network of Cancer Genes (NCG) 


## Requirements 💻
- Python 3.10+
- R 4.2+

## Setup Instructions 🔧

1. **Clone the repository:**
```bash
   git clone https://github.com/anekleiven/cancer_variants_annotation_pipeline.git
   cd cancer_variants_annotation_pipeline
```
2. **Create Virtual Environment:**
`python -m venv .venv`

3. **Activate Virtual Environment:**
`. .venv/bin/activate`

4. **Install Python Requirements:**
`pip install -r requirements.txt`

5. **Install R requirements:**
`Rscript install_deps.R`


## Prerequisites 📋

This pipeline is designed to annotate variants that have already been labeled for oncogenicity. 
Before running the scripts in this repository, please complete the processing steps described here:

[genie_oncokb_processing_scripts](https://github.com/anekleiven/genie_oncokb_processing_scripts)


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

## Data Flow 📋
The pipeline expects a variant file in `.tsv` format as starting input. 
Each script (01-08) appends specific annotations to the dataset, resulting in a 
final master file used for downstream analysis.


## Pipeline Workflow 🚀

The pipeline consists of a series of scripts. Each step appends specific annotations to the dataset, sequentially building the final variant file.

| Step | Script / Action | Description |
| :--- | :--- | :--- |
| **01A** | `01A_hotspots_long.py` | Transforms raw `hotspot.txt` into a long-format TSV. |
| **01B** | `01B_annotate_hotspots.py` | Annotates variants with hotspot data from *cancerhotspots.org*. |
| **02** | `02_add_uniprot_mapping.R` | Maps variants to UniProt accession numbers using `geneOncoX`. |
| **03** | `03_annotate_protein_domains.py` | Annotates variants with protein domain information from *Pfam*. |
| **04** | `04_annotate_functional_sites.py` | Annotates variants with specific functional sites from *UniProt*. |
| **05** | `05_annotate_germline_proximity.py` | Calculates distance to pathogenic germline variants from *ClinVar*. |
| **06** | `06_mave_tsv_to_vcf.py` | Converts TSV to VCF format for liftover and VEP processing. |
| **Manual**| *Assembly Liftover* | Uses **CrossMap** to convert coordinates from GRCh37 to GRCh38. |
| **Manual**| *VEP Annotation* | Uses **Ensembl VEP** (via Docker) to annotate with *MaveDB* scores. |
| **07** | `07_annotate_maves.py` | Integrates functional data from *MaveDB* into the full variant dataset. |
| **08** | `08_annotate_tsg_og.R` | Adds TSG/Oncogene annotations from *NCG* using `geneOncoX`. |

> **Note:** The manual steps (Liftover and VEP) are necessary because MAVE data integration requires GRCh38 coordinates and the specialized VEP plugin.


## Recommended Sources 🛜

- AACR Project GENIE: https://www.aacr.org/professionals/research/aacr-project-genie/
- OncoKB: https://www.oncokb.org
- geneOncoX: https://sigven.github.io/geneOncoX/index.html


