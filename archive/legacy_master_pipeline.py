import os
import pandas as pd
import pyrodigal
from Bio import Entrez, SeqIO
from pathlib import Path

# ==========================================
# PART 1: Configuration
# ==========================================
Entrez.email = "your_email@example.com"      # NCBI requires an email for tracking
INPUT_CSV = "test_virus.csv"                  # Your metadata file from Step 1
BASE_DATA_DIR = Path("./ncbi_dataset/data")  # Root directory for all data [cite: 67]

# ==========================================
# PART 2: Core Functions
# ==========================================

def run_integrated_pipeline():
    # 1. Create the base data directory [cite: 48]
    BASE_DATA_DIR.mkdir(parents=True, exist_ok=True)
    
    # 2. Read the CSV to get all Accession IDs
    if not os.path.exists(INPUT_CSV):
        print(f"❌ Error: {INPUT_CSV} not found in the current folder.")
        return
        
    df = pd.read_csv(INPUT_CSV)
    accessions = df['Accession'].unique()
    print(f"🚀 Starting pipeline for {len(accessions)} sequences...")

    # Initialize GeneFinders for Pyrodigal [cite: 45-46]
    single_finder = pyrodigal.GeneFinder(meta=False)
    meta_finder = pyrodigal.GeneFinder(meta=True)

    for acc in accessions:
        # Create a dedicated folder for each phage [cite: 48]
        folder = BASE_DATA_DIR / acc
        folder.mkdir(exist_ok=True)
        
        fna_path = folder / f"{acc}.fna"
        faa_path = folder / "proteins.faa"

        # --- A. DOWNLOAD STAGE (Replaces manual Terminal commands) ---
        if not fna_path.exists():
            print(f"📥 [DOWNLOADING] {acc}...")
            try:
                handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
                with open(fna_path, "w") as f:
                    f.write(handle.read())
                handle.close()
            except Exception as e:
                print(f"❌ {acc} Download failed: {e}")
                continue
        else:
            print(f"✅ [SKIP DOWNLOAD] {acc} already exists.")

        # --- B. ANNOTATION STAGE (Prodigal Logic) [cite: 35-37] ---
        if not faa_path.exists():
            print(f"🧬 [ANNOTATING] {acc}...")
            try:
                for record in SeqIO.parse(fna_path, "fasta"):
                    dna_seq = str(record.seq)
                    
                    # Length check to avoid training errors 
                    if len(dna_seq) >= 20000:
                        single_finder.train(dna_seq) # Self-training for long sequences [cite: 60]
                        genes = single_finder.find_genes(dna_seq)
                    else:
                        genes = meta_finder.find_genes(dna_seq) # Meta-mode for short phages
                    
                    # Write protein translations to .faa file [cite: 62-63]
                    with open(faa_path, "w") as out:
                        genes.write_translations(out, sequence_id=record.id)
                print(f"✨ [SUCCESS] {acc} proteins generated.")
            except Exception as e:
                print(f"❌ {acc} Annotation error: {e}")
        else:
            # Skip if already processed to save time 
            print(f"✅ [SKIP ANNOTATION] {acc} proteins already exist.")

    print("\n🎉 All tasks completed successfully!")

if __name__ == "__main__":
    run_integrated_pipeline()