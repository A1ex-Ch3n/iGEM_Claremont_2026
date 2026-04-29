import os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
phanotate_path = "PHANOTATE/phanotate.py"
input_dir = "genomes"
output_dir = "output"

os.makedirs(output_dir, exist_ok=True)

for file in os.listdir(input_dir):

    if file.endswith((".fna", ".fa", ".fasta")):

        genome_path = os.path.join(input_dir, file)
        base = os.path.splitext(file)[0]

        coord_file = os.path.join(output_dir, base + ".txt")
        faa_file = os.path.join(output_dir, base + ".faa")

        # Run PHANOTATE
        with open(coord_file, "w") as out:
            subprocess.run(
                ["python3", phanotate_path, genome_path],
                stdout=out
            )

        # Load genome sequence
        record = SeqIO.read(genome_path, "fasta")
        genome_seq = record.seq

        proteins = []

        # Parse PHANOTATE output
        with open(coord_file) as f:
            for i, line in enumerate(f):

                if line.startswith("#"):
                    continue

                parts = line.strip().split()

                start = int(parts[0])
                end = int(parts[1])