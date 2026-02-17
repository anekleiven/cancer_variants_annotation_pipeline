import pandas as pd

df = pd.read_csv("output/variants_with_maves.tsv", sep="\t")

# Filter for rows where the score is NOT null
matches = df[df["MaveDB_score"].notna()]

# Display the key columns for the first 10 matches
print(matches[["Hugo_Symbol", "HGVSp", "MaveDB_score", "MaveDB_pro"]].head(10))