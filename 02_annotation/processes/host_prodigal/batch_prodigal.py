import os
import pyrodigal
from Bio import SeqIO
from pathlib import Path

def run_pure_python_prodigal(data_dir):
    """
    Runs Prodigal natively in pure Python with self-training.
    """
    base_path = Path(data_dir)
    fna_files = list(base_path.rglob("*.fna"))
    
    if not fna_files:
        print(f"❌ No .fna files found. Path checked: {data_dir}")
        return

    print(f"🔍 Found {len(fna_files)} genome(s). Starting annotation...\n")
    print("-" * 40)

    # Use GeneFinder for pyrodigal 3.0+
    gene_finder = pyrodigal.GeneFinder(meta=False)

    for fna_path in fna_files:
        parent_dir = fna_path.parent
        output_faa = parent_dir / "proteins.faa"
        output_gff = parent_dir / "genes.gff"
        
        if output_faa.exists():
            print(f"✅ [SKIPPED] {parent_dir.name} already processed.")
            continue
            
        print(f"⏳ [PROCESSING] {fna_path.name} ...")
        
        with open(output_faa, "w") as faa_out, open(output_gff, "w") as gff_out:
            for record in SeqIO.parse(fna_path, "fasta"):
                dna_seq = str(record.seq)
                
                # NEW: Train the model on the sequence first
                gene_finder.train(dna_seq)
                
                # Predict genes
                genes = gene_finder.find_genes(dna_seq)
                
                # Write outputs
                genes.write_gff(gff_out, sequence_id=record.id)
                genes.write_translations(faa_out, sequence_id=record.id)
                
        print(f"✨ [SUCCESS] Finished {parent_dir.name}!")

if __name__ == "__main__":
    # Ensure this is the folder path
    TARGET_DIRECTORY = "./ncbi_dataset"
    run_pure_python_prodigal(TARGET_DIRECTORY)