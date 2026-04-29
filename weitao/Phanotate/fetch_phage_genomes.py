"""
fetch_phage_genomes.py

Reads phage accession IDs from the second column of phage_host_matrix_with_ids.csv,
then downloads each complete genome from NCBI in FASTA format using the Entrez API.
Output .fasta files are ready to use directly as input for PHANOTATE.

Requirements:
    pip install biopython pandas requests

Usage:
    python fetch_phage_genomes.py

Output:
    - One .fasta file per phage, saved in ./phage_genomes/
    - A summary CSV (download_summary.csv) listing success/failure for each accession
"""

import os
import time
import pandas as pd
from Bio import Entrez, SeqIO

# ── CONFIGURATION ─────────────────────────────────────────────────────────────
CSV_FILE      = "phage_host_matrix_with_ids.csv"   # must be in the same folder
OUTPUT_DIR    = "phage_genomes"                     # folder where .fasta files are saved
EMAIL         = "wke29@students.claremontmckenna.edu"            # required by NCBI – change this!
BATCH_SIZE    = 10    # number of records to request per API call
SLEEP_BETWEEN = 0.4  # seconds between calls (NCBI allows ~3 req/s without an API key)
# ─────────────────────────────────────────────────────────────────────────────

Entrez.email = EMAIL
os.makedirs(OUTPUT_DIR, exist_ok=True)


def read_accession_ids(csv_path: str) -> list[str]:
    """Return a deduplicated list of accession IDs from the second column."""
    df = pd.read_csv(csv_path)
    col = df.iloc[:, 1]                        # second column by position
    ids = col.dropna().astype(str).str.strip()
    ids = ids[ids != ""].unique().tolist()
    print(f"[INFO] Found {len(ids)} unique accession IDs in '{csv_path}'")
    return ids


def fetch_genomes(accession_ids: list[str]) -> pd.DataFrame:
    """
    Download FASTA records for each accession ID (ready for PHANOTATE input).
    Returns a summary DataFrame with columns: accession, status, filename, note.
    """
    summary_rows = []

    for i in range(0, len(accession_ids), BATCH_SIZE):
        batch = accession_ids[i : i + BATCH_SIZE]
        id_str = ",".join(batch)
        print(f"\n[FETCH] Batch {i // BATCH_SIZE + 1}: {batch}")

        try:
            # efetch retrieves full GenBank records
            handle = Entrez.efetch(
                db="nucleotide",
                id=id_str,
                rettype="fasta",    # nucleotide FASTA — direct input for PHANOTATE
                retmode="text",
            )
            records = list(SeqIO.parse(handle, "fasta"))
            handle.close()

            fetched_ids = {rec.id.split(".")[0]: rec for rec in records}

            for acc in batch:
                acc_base = acc.split(".")[0]           # strip version number for matching
                matched = fetched_ids.get(acc_base)

                if matched:
                    filename = os.path.join(OUTPUT_DIR, f"{acc}.fasta")
                    with open(filename, "w") as out_handle:
                        SeqIO.write(matched, out_handle, "fasta")
                    seq_len = len(matched.seq)
                    print(f"  [OK]  {acc} → {filename}  ({seq_len:,} bp)")
                    summary_rows.append(
                        {"accession": acc, "status": "success",
                         "filename": filename, "note": f"{seq_len} bp"}
                    )
                else:
                    print(f"  [WARN] {acc} — record not found in batch response")
                    summary_rows.append(
                        {"accession": acc, "status": "not_found",
                         "filename": "", "note": "Not returned by NCBI"}
                    )

        except Exception as e:
            print(f"  [ERROR] Batch failed: {e}")
            for acc in batch:
                summary_rows.append(
                    {"accession": acc, "status": "error",
                     "filename": "", "note": str(e)}
                )

        time.sleep(SLEEP_BETWEEN)   # be polite to NCBI

    return pd.DataFrame(summary_rows)


def main():
    accession_ids = read_accession_ids(CSV_FILE)

    summary_df = fetch_genomes(accession_ids)

    summary_path = "download_summary.csv"
    summary_df.to_csv(summary_path, index=False)

    total    = len(summary_df)
    success  = (summary_df["status"] == "success").sum()
    failed   = total - success

    print("\n" + "=" * 50)
    print(f"DONE — {success}/{total} genomes downloaded successfully.")
    if failed:
        print(f"       {failed} accessions had errors or were not found.")
    print(f"FASTA files saved in    : ./{OUTPUT_DIR}/")
    print(f"Summary saved to        : {summary_path}")


if __name__ == "__main__":
    main()