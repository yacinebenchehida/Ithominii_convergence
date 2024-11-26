import argparse

def remove_lines(input_file, remove_file, output_file):
    # Read the lines to remove
    with open(remove_file, "r") as rf:
        lines_to_remove = set(int(line.strip()) for line in rf)

    # Process the input file and write the filtered lines
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line_number, line in enumerate(infile, start=1):
            if line_number not in lines_to_remove:
                outfile.write(line)

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Remove specified lines from a file.")
    parser.add_argument("input_file", help="Path to the input file")
    parser.add_argument("remove_file", help="Path to the file with line numbers to remove")
    parser.add_argument("output_file", help="Path to save the output file")
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Call the function with the provided arguments
    remove_lines(args.input_file, args.remove_file, args.output_file)
