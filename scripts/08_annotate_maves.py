""" 
====================================================================
Annotate somatic variant data with MAVEs
====================================================================

Script: 08_annotate_maves.py 
Author: Ane Kleiven 

Major outputs: 
  1. Extract MAVE score, gene and HGVSp from the .vcf file
  2. Load somatic variant data 
  3. Merge MAVE data with the original variant data file (.tsv) 


"""

# import libraries 
import argparse
import pandas as pd 
import re 
from pathlib import Path

# create argparse function for user input file paths 
def getargs(): 
    parser = argparse.ArgumentParser(
        description="Merge MAVE annotated data with the original somatic variant file (.tsv)"
    ) 

    parser.add_argument(
        "--mave", 
        type=Path, 
        default="output/somatic_variants_with_maves.vcf",
        help="Path to the MAVE annotated input file (e.g. output/somatic_variants_with_maves.vcf)"
    )

    parser.add_argument(
        "--input", 
        type=Path, 
        default="output/variants_with_germline_proximity.tsv",
        help="Path to the somatic variant input file (e.g. output/variants_with_germline_proximity.tsv)"
    )

    parser.add_argument(
        "--output",
        type=Path,
        default="output/variants_with_maves.tsv",
        help="Path to the annotated output variant file (e.g. output/variants_with_maves.tsv)"
    )

    return parser.parse_args() 


def strip_prefix(hgvsp):
    """Remove transcript prefix: ENSP000123:p.Trp690Ter -> p.Trp690Ter."""
    if hgvsp and ":" in hgvsp:
        return hgvsp.split(":")[-1]
    return hgvsp


# function to extract MAVE information and convert to df 
def vcf_to_dataframe(vcf_path): 
    """
    Analyze VEP VCF file and return a df with one row per CSQ entry, containing: 
    Hugo_Symbol, HGVSp, MaveDB_score, MaveDB_urn, MaveDB_nt, MAVEdb_pro. 
    Only rows with MAVEdb_score are kept. 
    """

    rows = [] 

    with open(vcf_path, "r") as file: 
        for line in file:
            # extract fields in the csq header
            if line.startswith("##INFO=<ID=CSQ"):   
                # search for the word format followed by any characters that aren't a double quote
                match = re.search(r"Format: ([^\"]+)\"", line) 
                # remove whitespace, split by "|" and save as a list 
                csq_fields = match.group(1).strip().split("|") 
                continue 
            # skip remaining lines
            if line.startswith("#"):   
                continue
            # raise value error if header is not found 
            if csq_fields is None:
                raise ValueError("CSQ FORMAT header line not found.") 
        
            # split lines in vcf files by tab and take the INFO column
            info = line.split("\t")[7] 
            # searches for CSQ followed by any characters that arent tab or semicolon
            csq_match = re.search(r"CSQ=([^\t;]+)", info) 
            # lines without CSQ annotation are skipped
            if not csq_match:
                continue 

            # CSQ entries are separated by commas, but field values (one field in one CSQ entry) can also include commas. 
            # Solution: split on commas, then rejoin fragments that belong to the same entry. 
            # Based on the expected number of pipes (|) 
            n_pipes = len(csq_fields) - 1 
            raw_entries = csq_match.group(1).split(",")
            entries = []
            buffer = "" 
            # applies fragments with enough pipes to entries 
            # fragments with fewer pipes are kept in the buffer and glued to the next comma-separated piece. 
            for fragment in raw_entries: 
                buffer = f"{buffer},{fragment}" if buffer else fragment  
                if buffer.count("|") >= n_pipes: 
                    entries.append(buffer) 
                    buffer = "" 
            # map column names to data 
            for entry in entries: 
                f = dict(zip(csq_fields, entry.split("|"))) 

                mave_score_raw = f.get("MaveDB_score", "").strip()
                mave_score = mave_score_raw.split("&")[0] if mave_score_raw else None

                if not mave_score or mave_score == "NA": 
                    continue # skip entries without mave score 
                
                symbol = f.get("SYMBOL", "").strip() 
                if not symbol:
                    continue 


                hgvsp  = strip_prefix(f.get("HGVSp", "").strip())
                mave_pro = f.get("MaveDB_pro", "").strip() 

                #check if hgvsp is missing/NA and we have mave protein data 
                if (not hgvsp or hgvsp == "NA") and mave_pro:
                    hgvsp = mave_pro.split("&")[0].strip() 
                
                if not hgvsp or hgvsp == "NA":
                    continue 

                rows.append({
                    "Hugo_Symbol":  symbol,
                    "HGVSp":        hgvsp,
                    "MaveDB_score": mave_score,
                    "MaveDB_urn":   f.get("MaveDB_urn",  "").strip() or None,
                    "MaveDB_nt":    f.get("MaveDB_nt",   "").strip() or None,
                    "MaveDB_pro":   f.get("MaveDB_pro",  "").strip() or None,
                })

    mave_df = pd.DataFrame(rows).drop_duplicates(subset=["Hugo_Symbol", "HGVSp"])
    print(f"VCF: {len(mave_df):,} unique entries with a MAVE score")
    return mave_df


def main(somatic_variants, vcf_path, out_path):
    
    print(f"Parsing VCF: {vcf_path}")
    mave_df = vcf_to_dataframe(vcf_path)

    if not mave_df.empty:
        mave_df["HGVSp"] = mave_df["HGVSp"].str.strip() 

    print(f"Loading somatic variant file: {somatic_variants}")
    somatic_variants = pd.read_csv(somatic_variants, sep="\t", low_memory=False)
    print(f"Somatic variant file: {len(somatic_variants):,} variants")

    if "HGVSp" in somatic_variants.columns: 
        somatic_variants["HGVSp"] = somatic_variants["HGVSp"].str.strip() 

    print("Joining on Hugo_Symbol + HGVSp...")
    merged_variants = somatic_variants.merge(mave_df, on=["Hugo_Symbol", "HGVSp"], how="left")

    n_scored = merged_variants["MaveDB_score"].notna().sum()
    print(f"Variants with MaveDB score: {n_scored:,} / {len(merged_variants):,} ({100*n_scored/len(merged_variants):.1f}%)")

    merged_variants.to_csv(out_path, sep="\t", index=False)
    print(f"Written: {out_path}")


if __name__ == "__main__":
    args = getargs() 
    main(args.input, args.mave, args.output)

