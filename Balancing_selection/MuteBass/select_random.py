import pandas as pd
import random
import argparse
import os

def process_species(input_file, output_dir):
    # Read the input file
    data = pd.read_csv(input_file, sep="\t", header=None, names=["Sample", "Species"])

    # Process each unique species
    for species in data["Species"].unique():
        species_data = data[data["Species"] == species]

        # Select 5 random rows if more than 5 samples exist, otherwise select all
        selected_rows = species_data if len(species_data) <= 5 else species_data.sample(n=5, random_state=1)

        # Write the output to a species-named file in the specified output directory
        output_file = os.path.join(output_dir, f"{species}.txt")
        selected_rows.to_csv(output_file, sep="\t", index=False, header=False)

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Process species data and write selected samples to output files.")
    parser.add_argument("input_file", help="Path to the input file")
    parser.add_argument("output_dir", help="Path to the output directory")

    # Parse arguments
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Call the function with provided arguments
    process_species(args.input_file, args.output_dir)
