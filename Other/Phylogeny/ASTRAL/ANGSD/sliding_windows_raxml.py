#!/usr/bin/env python3

import sys  # Importing the sys module for system-specific parameters and functions
import os  # Importing the os module for operating system functionalities
import subprocess  # Importing the subprocess module to spawn new processes, such as running external commands
from Bio import SeqIO  # Importing SeqIO from Bio module for reading sequence files in various formats

# Function to slice sequences into windows
def slice_sequences(fasta_file, window_size):
	# Open the FASTA file and read sequences
	with open(fasta_file, "r") as handle:  # Open the FASTA file in read mode using a context manager
		records = list(SeqIO.parse(handle, "fasta"))  # Parse the sequences using Bio.SeqIO and store them as a list

	# Get the length of the sequences
	seq_length = len(records[0].seq)  # Get the length of the first sequence

	# Slice sequences into windows
	windows = []  # Initialize an empty list to store the windows
	for start in range(0, seq_length, window_size):  # Iterate over the sequence length in steps of window_size
		end = min(start + window_size, seq_length)  # Calculate the end index of the current window
		window = [(record.id, record.seq[start:end]) for record in records]  # Extract sequences for the window
		windows.append(window)  # Append the window to the list of windows

	return windows  # Return the list of windows

# Function to write windows to separate FASTA files
def write_fasta(windows, output_dir):
	# Create output directory if it doesn't exist
	os.makedirs(output_dir, exist_ok=True)  # Create the output directory if it doesn't exist

	# Write each window to a separate FASTA file
	for i, window in enumerate(windows):  # Iterate over the windows
		output_file = os.path.join(output_dir, f"window_{i}.fasta")  # Define the output file path
		with open(output_file, "w") as f:  # Open the output file in write mode
			for id, seq in window:  # Iterate over the sequences in the window
				f.write(f">{id}\n{seq}\n")  # Write the sequence ID and sequence to the file

# Function to check if a window has a high N percentage
def has_high_n_percentage(window_seqs):
	total_bases = sum(seq.count("N") for id, seq in window_seqs)  # Count the total number of 'N' bases in the window
	total_length = sum(len(seq) for id, seq in window_seqs)  # Calculate the total length of the sequences in the window
	n_percentage = total_bases / total_length  # Calculate the percentage of 'N' bases
	return n_percentage > 0.2  # Return True if the percentage of 'N' bases is greater than 0.5, False otherwise

# Function to run RAxML on a window
def run_raxml_on_window(window_file, output_dir, window_index):
	# Check if the window file exists
	if not os.path.exists(window_file):
		return  # Return if the window file does not exist

	# Run RAxML on the window sequences
	subprocess.run(["raxmlHPC", "-m", "GTRCAT", "-s", window_file, "-n", f"window_{window_index}", "-p", "12345", "-T", "4", -w", os.path.abspath(output_dir)])

	# Remove unnecessary files for this window
	for file in os.listdir(output_dir):  # Iterate over files in the output directory
		if file.startswith(f"RAxML_parsimonyTree.window_{window_index}") or \
		   file.startswith(f"RAxML_result.window_{window_index}") or \
		   file.startswith(f"RAxML_info.window_{window_index}") or \
		   file.startswith(f"RAxML_log.window_{window_index}") or \
		   file.startswith(f"window_{window_index}"):
			os.remove(os.path.join(output_dir, file))  # Remove unnecessary files

# Main function
def main():
	# Check if the correct number of arguments is provided
	if len(sys.argv) != 3:
		print("Usage: python script.py input.fasta output_directory")
		sys.exit(1)  # Exit the script with error if the number of arguments is incorrect

	# Get input arguments
	fasta_file = sys.argv[1]  # Get the input FASTA file
	output_dir = sys.argv[2]  # Get the output directory
	window_size = 5000  # Define the window size

	# Slice sequences into windows
	windows = slice_sequences(fasta_file, window_size)
	
	# Filter out windows with high N percentage
	removed_windows = []
	filtered_windows = []
	for i, window in enumerate(windows):
		if has_high_n_percentage(window):
			removed_windows.append(i)
			print(f"Window_{i}.fasta has high N percentage and will be removed")
		else:
			filtered_windows.append(window)

	# Write all filtered windows to the specified output directory
	write_fasta(filtered_windows, output_dir)

	# Remove the fasta files with high N percentage from the output directory
	for i in removed_windows:
		fasta_file_path = os.path.join(output_dir, f"window_{i}.fasta")
		if os.path.exists(fasta_file_path):
			os.remove(fasta_file_path)

	# Run RAxML on the remaining windows
	for i, window in enumerate(filtered_windows):
		window_file = os.path.join(output_dir, f"window_{i}.fasta")
		run_raxml_on_window(window_file, output_dir, i)

# Entry point of the script
if __name__ == "__main__":
	main()
