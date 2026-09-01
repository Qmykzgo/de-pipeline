import urllib.request
import json
import pandas as pd
from pathlib import Path

data_dir = Path("/home/yer_kanat/Downloads/lol/depipeline1/data")
data_dir.mkdir(parents=True, exist_ok=True)

print("Fetching NCBI GEO GSE52778 (Airway Smooth Muscle - Dexamethasone Treatment)...")

counts_url = "https://raw.githubusercontent.com/bioconnector/workshops/master/data/airway_scaledcounts.csv"
metadata_url = "https://raw.githubusercontent.com/bioconnector/workshops/master/data/airway_metadata.csv"

# Download counts
req_c = urllib.request.Request(counts_url, headers={"User-Agent": "Mozilla/5.0"})
with urllib.request.urlopen(req_c, timeout=30) as r:
    counts_df = pd.read_csv(r)
print("Raw Counts loaded:", counts_df.shape)

# Download metadata
req_m = urllib.request.Request(metadata_url, headers={"User-Agent": "Mozilla/5.0"})
with urllib.request.urlopen(req_m, timeout=30) as r:
    meta_df = pd.read_csv(r)
print("Raw Metadata loaded:", meta_df.shape)

# Format metadata
meta_df = meta_df.rename(columns={"id": "sample_id", "dex": "condition", "celltype": "batch"})

# Structure counts
gene_col = counts_df.columns[0]
counts_df = counts_df.rename(columns={gene_col: "gene_id"})

# Query gene symbols in batches of 1000 using MyGene.info (standard open bioinformatics API)
gene_ids = counts_df["gene_id"].tolist()
print(f"Mapping {len(gene_ids)} Ensembl IDs to Gene Symbols via MyGene.info batch API...")

id_to_symbol = {}
batch_size = 1000
for i in range(0, len(gene_ids), batch_size):
    batch = gene_ids[i:i + batch_size]
    post_data = urllib.parse.urlencode({"q": ",".join(batch), "scopes": "ensembl.gene", "fields": "symbol"}).encode("utf-8")
    req_sym = urllib.request.Request("https://mygene.info/v3/query", data=post_data, headers={"User-Agent": "Mozilla/5.0"})
    try:
        with urllib.request.urlopen(req_sym, timeout=20) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            for item in data:
                q = item.get("query")
                sym = item.get("symbol")
                if q and sym:
                    id_to_symbol[q] = sym
    except Exception as e:
        print(f"  Batch {i} map error: {e}")
        pass

print(f"Mapped {len(id_to_symbol)} gene symbols successfully.")
counts_df["gene"] = counts_df["gene_id"].map(id_to_symbol).fillna(counts_df["gene_id"])

# Group by gene symbol and take max to avoid duplicate symbols
samples = [s for s in meta_df["sample_id"] if s in counts_df.columns]
counts_clean = counts_df[["gene"] + samples].groupby("gene").max().reset_index()

# Round to integer counts
for s in samples:
    counts_clean[s] = counts_clean[s].round().astype(int)

out_counts = data_dir / "ncbi_airway_counts.csv"
out_meta = data_dir / "ncbi_airway_metadata.csv"

counts_clean.to_csv(out_counts, index=False)
meta_df[["sample_id", "condition", "batch", "geo_id"]].to_csv(out_meta, index=False)

print("\n=== NCBI GSE52778 Airway Dataset Ready ===")
print(f"Saved counts to: {out_counts} ({counts_clean.shape[0]} genes × {len(samples)} samples)")
print(f"Saved metadata to: {out_meta}")
print("\nSample Sheet:")
print(meta_df[["sample_id", "condition", "batch", "geo_id"]].to_string(index=False))
print("\nTop 5 Count Matrix rows:")
print(counts_clean.head().to_string(index=False))
