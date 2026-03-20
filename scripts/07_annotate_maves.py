""" 
====================================================================
Annotate variant data with MAVEs
====================================================================

Script: 07_annotate_maves.py 
Author: Ane Kleiven 

Major outputs: 
  1. Extract MAVE data, gene and HGVSp from the .vcf file
  2. Load variant data (.tsv file)
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
        description="Merge MAVE annotated data with .tsv variant file"
    ) 

    parser.add_argument(
        "--mave_somatic", 
        type=Path, 
        default="output/somatic_variants_with_maves.vcf",
        help="Path to the MAVE annotated input file with somatic variants."
    )

    parser.add_argument(
        "--mave_neutral_clinvar", 
        type=Path, 
        default="output/neutral_clinvar_with_maves.vcf",
        help="Path to the MAVE annotated input file with neutral germline ClinVar variants."
    )

    parser.add_argument(
        "--somatic_variants", 
        type=Path, 
        default="output/variants_with_germline_proximity.tsv",
        help="Path to the somatic variant file."
    )

    parser.add_argument(
        "--neutral_clinvar", 
        type=Path,
        default="output/neutral_clinvar_filtered.tsv", 
        help="Path to the ClinVar variant file containing neutral germline variants."

    )

    parser.add_argument(
        "--output_somatic",
        type=Path,
        default="output/variants_with_maves.tsv",
        help="Path to the annotated output variant file (e.g. output/variants_with_maves.tsv)"
    )

    parser.add_argument(
        "--output_neutral_clinvar",
        type=Path, 
        default="output/clinvar_with_maves.tsv",
        help="Path to the annotated .tsv file from ClinVar."
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
    Analyze VEP VCF file and return a df with one row per CSQ entry per experiment,
    containing: Hugo_Symbol, HGVSp, MaveDB_score, MaveDB_urn, MaveDB_nt, MaveDB_pro.
    Only rows with a MaveDB_score are kept.
    """

    rows = [] 
    csq_fields = None

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
                mave_scores = mave_score_raw.split("&") if mave_score_raw else []

                # skip if all scores are missing/NA
                if not mave_scores or all(s.strip() in ("NA", "") for s in mave_scores):
                    continue

                symbol = f.get("SYMBOL", "").strip() 
                if not symbol:
                    continue 

                hgvsp = strip_prefix(f.get("HGVSp", "").strip())
                mave_pro = f.get("MaveDB_pro", "").strip() 

                # check if hgvsp is missing/NA and we have mave protein data 
                if (not hgvsp or hgvsp == "NA") and mave_pro:
                    hgvsp = mave_pro.split("&")[0].strip() 
                
                if not hgvsp or hgvsp == "NA":
                    continue 

                mave_urns = f.get("MaveDB_urn", "").strip().split("&")
                mave_nts  = f.get("MaveDB_nt",  "").strip().split("&")
                mave_pros = f.get("MaveDB_pro",  "").strip().split("&")

                # one row per experiment (score/urn pair)
                for i, score in enumerate(mave_scores):
                    score = score.strip()
                    if not score or score == "NA":
                        continue
                    rows.append({
                        "Hugo_Symbol":  symbol,
                        "HGVSp":        hgvsp,
                        "MaveDB_score": score,
                        "MaveDB_urn":   mave_urns[i].strip() if i < len(mave_urns) else None,
                        "MaveDB_nt":    mave_nts[i].strip()  if i < len(mave_nts)  else None,
                        "MaveDB_pro":   mave_pros[i].strip() if i < len(mave_pros) else None,
                    })

    mave_df = pd.DataFrame(rows)
    mave_df["MaveDB_score"] = pd.to_numeric(mave_df["MaveDB_score"], errors="coerce")
    mave_df = mave_df.dropna(subset=["MaveDB_score"])
    mave_df = mave_df.drop_duplicates(subset=["Hugo_Symbol", "HGVSp", "MaveDB_urn"])

    # split urn structure into components; experiment_set, experiment, score_set 
    mave_df["experiment_set"] = mave_df["MaveDB_urn"].str.extract(r"urn:mavedb:(\d+)")
    mave_df["experiment"]     = mave_df["MaveDB_urn"].str.extract(r"urn:mavedb:\d+-([\w]+)-")
    mave_df["score_set"]      = mave_df["MaveDB_urn"].str.extract(r"urn:mavedb:\d+-\w+-(\d+)")

    print(f"VCF: {len(mave_df):,} unique entries with a MAVE score (one row per experiment)")

    return mave_df


def process_dataset(vcf_path, tsv_path, out_path, dataset_name):
    """
    Help function to process MAVE datasets
    
    """
    print(f"\n--Processing {dataset_name}--")

    if not vcf_path.exists() or not tsv_path.exists():
        print(f"Skipping {dataset_name}: Input files not found.")
        return

    print(f"Parsing VCF: {vcf_path}")
    mave_df = vcf_to_dataframe(vcf_path)
    if not mave_df.empty:
        mave_df["HGVSp"] = mave_df["HGVSp"].str.strip()

    print(f"Loading TSV: {tsv_path}")
    df = pd.read_csv(tsv_path, sep="\t", low_memory=False)

    # Handle GeneSymbol/Hugo_Symbol column
    if "GeneSymbol" in df.columns and "Hugo_Symbol" not in df.columns:
        print("Renaming 'GeneSymbol' to 'Hugo_Symbol' before merging..")
        df = df.rename(columns={"GeneSymbol": "Hugo_Symbol"})

    if "HGVSp" in df.columns:
        df["HGVSp"] = df["HGVSp"].str.strip()

    print(f"Joining {dataset_name} on Hugo_Symbol + HGVSp...")
    merged = df.merge(mave_df, on=["Hugo_Symbol", "HGVSp"], how="left")

    n_scored = merged["MaveDB_score"].notna().sum()
    print(f"Variants with MaveDB score: {n_scored:,} / {len(merged):,} ({100*n_scored/len(merged):.1f}%)")

    # Save deduplicated file (one row per variant, largest score)
    unique_out = merged.sort_values("MaveDB_score", ascending=False, na_position="last") \
    .drop_duplicates(subset=["Hugo_Symbol", "HGVSp"])
    unique_out.to_csv(out_path, sep="\t", index=False)
    print(f"Written (deduplicated): {out_path}")

    # Save expanded file (for MAVE-analysis) 
    expanded_path = out_path.with_name(out_path.stem + "_expanded" + out_path.suffix)
    merged.to_csv(expanded_path, sep="\t", index=False)
    print(f"Written (expanded): {expanded_path}\n.")


def main():
    args = getargs()

    # 1. Somatic variants 
    process_dataset(
        args.mave_somatic, 
        args.somatic_variants, 
        args.output_somatic, 
        "Somatic"
    )

    # 2. ClinVar neutral variants 
    process_dataset(
        args.mave_neutral_clinvar, 
        args.neutral_clinvar, 
        args.output_neutral_clinvar, 
        "ClinVar"
    )

if __name__ == "__main__":
    main()