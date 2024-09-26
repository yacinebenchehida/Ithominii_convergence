from Bio import Phylo
import sys

def rename_tree_tips(newick_file, group_file, output_file):
    # Read the group file
    with open(group_file, 'r') as f:
        groups = f.read().splitlines()

    # Read the Newick trees
    trees = Phylo.parse(newick_file, 'newick')  # Use parse to handle multiple trees

    # Dictionary to store tip replacements
    tip_rename_map = {}

    # Loop over pairs of tips in the tree and group names
    for idx, group_name in enumerate(groups):
        # The two tips for each group are idx*2 and idx*2 + 1
        tip_rename_map[str(idx*2)] = group_name
        tip_rename_map[str(idx*2 + 1)] = group_name

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
