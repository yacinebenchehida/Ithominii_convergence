import argparse

def classify_sequence(sequence):
    # Define the set of valid DNA letters
    dna_letters = set("ACGT")

    # Define the set of valid protein letters
    protein_letters = set("ACDEFGHIKLMNPQRSTVWY")

    # Check if the sequence contains only DNA letters
    if set(sequence.upper()) <= dna_letters:
        return "DNA"
    # Check if the sequence contains only protein letters
    elif set(sequence.upper()) <= protein_letters:
        return "Protein"
    else:
        return "Unknown"

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Classify a sequence in a FASTA file as DNA or Protein.')
    parser.add_argument('fasta_file', help='Path to the input FASTA file')
    args = parser.parse_args()

    # Read the FASTA file
    with open(args.fasta_file, 'r') as file:
        lines = file.readlines()

    # Extract sequence from the FASTA file (ignoring header lines)
    sequence = ''.join(line.strip() for line in lines[1:])

    # Classify the sequence
    classification = classify_sequence(sequence)

    print(classification)

if __name__ == "__main__":
    main()
