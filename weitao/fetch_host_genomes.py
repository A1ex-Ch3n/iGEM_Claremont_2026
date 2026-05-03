"""
fetch_host_genomes.py

Reads host accession IDs from the fourth column of phage_host_matrix_with_ids.csv,
then downloads for each host:
  1. Complete genome in FASTA format  → ./host_genomes/
  2. CDS (coding sequences, protein)  → ./host_cds/

Rows with "No Complete Genome Found" in the Host_Accession column are skipped.

Requirements:
    pip install biopython pandas requests

Usage:
    python3 fetch_host_genomes.py
    -- or in ipython --
    import ssl; ssl._create_default_https_context = ssl._create_unverified_context
    %run fetch_host_genomes.py

Output:
    - One .fasta file per host, saved in ./host_genomes/
    - One .faa file per host (CDS amino acid), saved in ./host_cds/
    - A summary CSV (host_download_summary.csv) listing success/failure for each accession
"""

import os
import time
import pandas as pd
from Bio import Entrez, SeqIO

# ── CONFIGURATION ─────────────────────────────────────────────────────────────
CSV_FILE        = "phage_host_matrix_with_ids.csv"  # must be in the same folder
GENOME_DIR      = "host_genomes"                    # complete genome FASTA files
CDS_DIR         = "host_cds"                        # CDS / protein FASTA files
EMAIL           = "wke29@students.claremontmckenna.edu"
BATCH_SIZE      = 5     # smaller batches — host genomes are much larger than phage
SLEEP_BETWEEN   = 0.5  # seconds between calls
SKIP_VALUE      = "No Complete Genome Found"
# ─────────────────────────────────────────────────────────────────────────────

Entrez.email = EMAIL
os.makedirs(GENOME_DIR, exist_ok=True)
os.makedirs(CDS_DIR, exist_ok=True)


def read_host_accession_ids(csv_path: str) -> list[str]:
    """Return a deduplicated list of host accession IDs from the fourth column,
    skipping rows marked as 'No Complete Genome Found'."""
    df = pd.read_csv(csv_path)
    col = df.iloc[:, 3]                         # fourth column (Host_Accession)
    ids = col.dropna().astype(str).str.strip()
    ids = ids[ids != ""]
    ids = ids[ids != SKIP_VALUE]
    ids = ids.unique().tolist()
    print(f"[INFO] Found {len(ids)} unique host accession IDs in '{csv_path}'")
    return ids


def fetch_genomes(accession_ids: list[str]) -> pd.DataFrame:
    """
    For each host accession ID, download:
      - Complete genome FASTA  (.fasta) → GENOME_DIR
      - CDS protein FASTA      (.faa)   → CDS_DIR
    Returns a summary DataFrame.
    """
    summary_rows = []

    for i in range(0, len(accession_ids), BATCH_SIZE):
        batch = accession_ids[i : i + BATCH_SIZE]
        id_str = ",".join(batch)
        print(f"\n[FETCH] Batch {i // BATCH_SIZE + 1}: {batch}")

        # ── 1. Download complete genome FASTA ────────────────────────────────
        try:
            handle = Entrez.efetch(
                db="nucleotide",
                id=id_str,
                rettype="fasta",
                retmode="text",
            )
            records = list(SeqIO.parse(handle, "fasta"))
            handle.close()

            fetched_ids = {rec.id.split(".")[0]: rec for rec in records}

            for acc in batch:
                acc_base = acc.split(".")[0]
                matched = fetched_ids.get(acc_base)

                if matched:
                    filename = os.path.join(GENOME_DIR, f"{acc}.fasta")
                    with open(filename, "w") as out_handle:
                        SeqIO.write(matched, out_handle, "fasta")
                    seq_len = len(matched.seq)
                    print(f"  [GENOME OK]  {acc} → {filename}  ({seq_len:,} bp)")
                    genome_status = "success"
                    genome_note   = f"{seq_len} bp"
                else:
                    print(f"  [GENOME WARN] {acc} — not found in batch response")
                    genome_status = "not_found"
                    genome_note   = "Not returned by NCBI"

                summary_rows.append({
                    "accession":     acc,
                    "genome_status": genome_status,
                    "genome_file":   filename if genome_status == "success" else "",
                    "genome_note":   genome_note,
                })

        except Exception as e:
            print(f"  [GENOME ERROR] Batch failed: {e}")
            for acc in batch:
                summary_rows.append({
                    "accession":     acc,
                    "genome_status": "error",
                    "genome_file":   "",
                    "genome_note":   str(e),
                })

        time.sleep(SLEEP_BETWEEN)

        # ── 2. Download CDS (protein FASTA) ──────────────────────────────────
        try:
            handle = Entrez.efetch(
                db="nucleotide",
                id=id_str,
                rettype="fasta_cds_aa",   # protein sequences of all CDS features
                retmode="text",
            )
            raw_cds = handle.read()
            handle.close()

            for acc in batch:
                acc_base = acc.split(".")[0]
                # Filter CDS records belonging to this accession
                cds_records = [
                    r for r in SeqIO.parse(
                        __import__("io").StringIO(raw_cds), "fasta"
                    )
                    if acc_base in r.id or acc_base in r.description
                ]

                if cds_records:
                    cds_filename = os.path.join(CDS_DIR, f"{acc}.faa")
                    with open(cds_filename, "w") as out_handle:
                        SeqIO.write(cds_records, out_handle, "fasta")
                    print(f"  [CDS OK]     {acc} → {cds_filename}  ({len(cds_records)} CDS)")
                    # Update the matching summary row
                    for row in summary_rows:
                        if row["accession"] == acc:
                            row["cds_status"] = "success"
                            row["cds_file"]   = cds_filename
                            row["cds_note"]   = f"{len(cds_records)} CDS"
                else:
                    print(f"  [CDS WARN]   {acc} — no CDS records found")
                    for row in summary_rows:
                        if row["accession"] == acc:
                            row["cds_status"] = "not_found"
                            row["cds_file"]   = ""
                            row["cds_note"]   = "No CDS returned"

        except Exception as e:
            print(f"  [CDS ERROR] Batch failed: {e}")
            for acc in batch:
                for row in summary_rows:
                    if row["accession"] == acc and "cds_status" not in row:
                        row["cds_status"] = "error"
                        row["cds_file"]   = ""
                        row["cds_note"]   = str(e)

        time.sleep(SLEEP_BETWEEN)

    return pd.DataFrame(summary_rows)


def main():
    accession_ids = read_host_accession_ids(CSV_FILE)

    summary_df = fetch_genomes(accession_ids)

    summary_path = "host_download_summary.csv"
    summary_df.to_csv(summary_path, index=False)

    total          = len(summary_df)
    genome_success = (summary_df["genome_status"] == "success").sum()
    cds_success    = (summary_df["cds_status"]    == "success").sum()

    print("\n" + "=" * 50)
    print(f"DONE")
    print(f"  Genomes : {genome_success}/{total} downloaded successfully → ./{GENOME_DIR}/")
    print(f"  CDS     : {cds_success}/{total} downloaded successfully → ./{CDS_DIR}/")
    print(f"  Summary : {summary_path}")


if __name__ == "__main__":
    main()