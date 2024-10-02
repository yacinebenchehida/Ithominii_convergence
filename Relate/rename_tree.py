import sys

# Function to generate intervals between numbers
def generate_intervals(file_path):
    with open(file_path, 'r') as file:
        numbers = [int(line.strip()) for line in file]

    for i in range(len(numbers) - 1):
        start = numbers[i]
        end = numbers[i + 1] - 1
        print(f"{start} {end}")

    # For the last number, add 1000
    last_number = numbers[-1]
    print(f"{last_number} {last_number + 1000}")

# Main part
if __name__ == "__main__":
    file_path = sys.argv[1]
    generate_intervals(file_path)
(base) [ybc502@login2[viking2] Scripts]$ cat rename_tree.py
from Bio import Phylo
import sys
import re  # Import the regex module

def rename_tree_tips(newick_file, group_file, output_file):
    # Read the group file
    with open(group_file, 'r') as f:
        groups = f.read().splitlines()

    # Read the Newick trees
    trees = Phylo.parse(newick_file, 'newick')  # Use parse to handle multiple trees

    # Dictionary to store tip replacements
    tip_rename_map = {}
    group_mapping = []  # To store the mapping for the group.txt file

    # Loop over pairs of tips in the tree and group names
    for idx, group_name in enumerate(groups):
        # The two tips for each group are idx*2 and idx*2 + 1
        new_name_1 = f"{group_name}_1"
        new_name_2 = f"{group_name}_2"
        tip_rename_map[str(idx*2)] = new_name_1
        tip_rename_map[str(idx*2 + 1)] = new_name_2

        # Extract the desired subspecies using regex
        match = re.search(r'_([^_]+)_(\d)$', new_name_1)  # Extract last part before _1/_2
        if match:
            associated_value = match.group(1)  # The subspecies name before the last underscore
        else:
            associated_value = "unknown"  # Default if no match found

        # Store the mapping for the group.txt file
        group_mapping.append((new_name_1, associated_value))
        group_mapping.append((new_name_2, associated_value))

    # Initialize a list to hold renamed trees
    renamed_trees = []

    # Iterate over each tree in the Newick file
    for tree in trees:
        # Rename the tree tips based on the map
        for clade in tree.find_clades():
            if clade.name in tip_rename_map:
                clade.name = tip_rename_map[clade.name]

        # Add the renamed tree to the list
        renamed_trees.append(tree)

    # Write the renamed trees to the output file
    Phylo.write(renamed_trees, output_file, 'newick')

    # Write the group mapping to group.txt
    with open('group.txt', 'w') as g:
        for new_name, original_name in group_mapping:
            g.write(f"{new_name}\t{original_name}\n")

if __name__ == "__main__":
    # External inputs using sys.argv
    if len(sys.argv) != 4:
        print("Usage: python3 ./rename_tree.py <newick_file> <group_file> <output_file>")
        sys.exit(1)

    newick_file = sys.argv[1]    # Input newick file
    group_file = sys.argv[2]     # Input group file
    output_file = sys.argv[3]    # Output renamed tree file

    # Call the function
    rename_tree_tips(newick_file, group_file, output_file)
