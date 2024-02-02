#!/usr/bin/python3

import sys

def extract_columns(column_indices, data_file):
    # Read column indices from numbers.txt and add 1 to each index
    with open(column_indices, 'r') as numbers_file:
        selected_columns = list(map(lambda x: int(x) + 1, numbers_file.readline().strip().split(',')))

    # Extract selected columns from data.txt
    with open(data_file, 'r') as data_file:
        for line in data_file:
            fields = line.strip().split()
            selected_fields = [fields[i-1] if 1 <= i <= len(fields) else '' for i in selected_columns]
            print('\t'.join(selected_fields))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extract_columns.py numbers.txt data.txt")
        sys.exit(1)

    numbers_file_path = sys.argv[1]
    data_file_path = sys.argv[2]

    extract_columns(numbers_file_path, data_file_path)
