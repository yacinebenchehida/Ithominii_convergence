import sys

def extract_and_compare(input_file, species_name):
	first_line_values = None
	last_line_values = None

	with open(input_file, 'r') as file:
		for line in file:
			if species_name in line:
				values = line.strip().split('\t')
				if len(values) >= 5:
					try:
						start = int(values[3])
						end = int(values[4])
					except ValueError:
						continue  # Skip lines where start or end positions are not valid integers
					
					if not first_line_values:
						first_line_values = [start, end]
					
					last_line_values = [start, end]

	if first_line_values and last_line_values:
		all_values = first_line_values + last_line_values
		min_value = min(all_values)
		max_value = max(all_values)
		if(last_line_values > first_line_values):
			return min_value, max_value
		else:
			return max_value, min_value
	else:
		return None, None

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("Usage: python find_min_max.py <input_file> <species_name>")
		sys.exit(1)

	input_file = sys.argv[1]
	species_name = sys.argv[2]

	min_value, max_value = extract_and_compare(input_file, species_name)

	if min_value is not None and max_value is not None:
		print(str(min_value) + "\t" + str(max_value))
	else:
		print(f"No data found for species: {species_name}")
