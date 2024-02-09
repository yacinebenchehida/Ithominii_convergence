import os
import subprocess
from Bio import SeqIO

def slice_sequences(fasta_file, window_size):
    # Open the FASTA file and read sequences
    with open(fasta_file, "r") as handle:
        records = list(SeqIO.parse(handle, "fasta"))

    # Get the length of the sequences
    seq_length = len(records[0].seq)

    # Slice sequences into windows
    windows = []
    for start in range(0, seq_length, window_size):
        end = min(start + window_size, seq_length)
        window = [record[start:end] for record in records]
        windows.append(window)

    return windows

def run_raxml_on_window(window_seqs, output_dir):
    # Create a temporary FASTA file for the window sequences
    temp_fasta = os.path.join(output_dir, "window.fasta")
    with open(temp_fasta, "w") as temp_handle:
        for i, seq in enumerate(window_seqs):
            temp_handle.write(f">Seq_{i}\n{seq}\n")

    # Run RAxML on the window sequences
    output_tree = os.path.join(output_dir, "window_tree.tree")
    subprocess.run(["raxmlHPC-PTHREADS", "-m", "GTRGAMMA", "-s", temp_fasta, "-n", "window", "-T", "4", "-p", "12345", "-w", output_dir])

    # Clean up temporary files
    os.remove(temp_fasta)

def main():
    # Input parameters
    fasta_file = "input.fasta"
    window_size = 1000
    output_dir = "output"

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Slice sequences into windows
    windows = slice_sequences(fasta_file, window_size)

    # Run RAxML on each window
    for i, window_seqs in enumerate(windows):
        window_output_dir = os.path.join(output_dir, f"window_{i}")
        os.makedirs(window_output_dir, exist_ok=True)
        run_raxml_on_window(window_seqs, window_output_dir)

if __name__ == "__main__":
    main()
