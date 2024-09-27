import sys

def read_positions(file_path):
    """Read positions from a file into a list."""
    with open(file_path, 'r') as file:
        return [int(line.strip()) for line in file.readlines()]

def find_line_for_position(positions, x):
    """Find the index of the closest start position less than or equal to x."""
    valid_positions = [(index, pos) for index, pos in enumerate(positions) if pos <= x]
    return max(valid_positions, key=lambda item: item[1], default=None)

def extract_tree_from_file(tree_file, line_number):
    """Extract the tree from the Newick file at the given line number."""
    with open(tree_file, 'r') as file:
        lines = file.readlines()
        if line_number+1 < len(lines):
            return lines[line_number+1].strip()  # Return the Newick tree on that line
    return None

def save_tree_to_file(tree, prefix, position):
    """Save the extracted tree to a file with the given prefix and position."""
    output_file_name = f"{prefix}_{position}.newick"
    with open(output_file_name, 'w') as file:
        file.write(tree)

if __name__ == "__main__":
    # External inputs
    position_file = sys.argv[1]  # Path to positions file
    tree_file = sys.argv[2]       # Path to the Newick tree file
    x = int(sys.argv[3])          # Position of interest (as an integer)
    prefix = sys.argv[4]          # Prefix for the output file
    
    # Read positions
    positions = read_positions(position_file)
    
    # Find the line corresponding to the closest start position
    closest_position = find_line_for_position(positions, x)
    
    if closest_position is not None:
        line_index, _ = closest_position
        
        # Extract the tree corresponding to that line number
        tree = extract_tree_from_file(tree_file, line_index)
        
        if tree:
            save_tree_to_file(tree, prefix, x)  # Save the extracted tree to a file
