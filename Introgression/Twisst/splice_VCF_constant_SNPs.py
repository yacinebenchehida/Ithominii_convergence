import argparse
import gzip

#############################################
# Function to extract the header of the VCF #
#############################################
def get_vcf_header(input_vcf):
	header_lines = []
	with gzip.open(input_vcf, 'rt') as f:
		for line in f:
			if line.startswith('#'):
				header_lines.append(line)
			else:
				break
	return header_lines

######################################################################
# Function to determine if a postion in the VCF variant or invariant #
######################################################################
def is_variant(line):
	fields = line.split('\t')
	ref = fields[3]
	alt = fields[4]
	
	# If the ALT field is not a ".", it's a proper SNP
	return alt != '.'

###############################################################
# Function that logs the genomic range of each phylogeny part #
###############################################################
def log_phylogeny_ranges(log_file, file_count, start_position, end_position):
	with open(log_file, 'a') as log:
		log.write(f"{file_count}\t{start_position}\t{end_position}\n")

######################################################
# Function that splice VCF into slices of n variants #
######################################################
def slice_vcf(input_vcf, output_prefix, snps_per_file, log_file):
	# Extract header of the VCF
	header = get_vcf_header(input_vcf)
	
	# Counters
	variant_count = 0  # To count only variant sites
	file_count = 1	 # To number output files
	start_position = None  # Start position of the current phylogeny range
	
	# Open the input VCF for reading (bgzipped)
	with gzip.open(input_vcf, 'rt') as vcf_in:
		# Open the first output VCF (plain text, not gzipped)
		output_vcf = f"{output_prefix}_part_{file_count}.vcf"
		vcf_out = open(output_vcf, 'w')
		
		# Write the header to the first output file
		vcf_out.writelines(header)
		
		for line in vcf_in:
			if line.startswith('#'):
				# Skip header lines (already written)
				continue
			
			# Extract the position (from column 2 of the VCF line)
			position = int(line.split('\t')[1])
			
			# Set start_position for the first variant in the current file
			if start_position is None:
				start_position = position
			
			# Write the SNP or invariant site (variant or non-variant line) to the current file
			vcf_out.write(line)
			
			# If the line is a variant (SNP or other), count it
			if is_variant(line):
				variant_count += 1
			
			# If we've written 'snps_per_file' variants, close this file and open a new one
			if variant_count >= snps_per_file:
				end_position = position  # Record the end position of this phylogeny
				
				# Log the phylogeny range
				log_phylogeny_ranges(log_file, file_count, start_position, end_position)
				
				# Close the current file
				vcf_out.close()
				
				# Increment the file count and reset variant counter
				file_count += 1
				variant_count = 0
				start_position = None  # Reset for the next file
				
				# Open a new output VCF file for the next batch
				output_vcf = f"{output_prefix}_part_{file_count}.vcf"
				vcf_out = open(output_vcf, 'w')
				
				# Write the header to the new file
				vcf_out.writelines(header)
		
		# Close the last open file
		vcf_out.close()
		
		# If there are still variants in the last file, log the final range
		if start_position is not None:
			end_position = position
			log_phylogeny_ranges(log_file, file_count, start_position, end_position)
		
	print(f"Finished splitting the VCF into {file_count} parts, and logged ranges to {log_file}.")

########
# Main #
########
if __name__ == "__main__":
	# Argument parser
	parser = argparse.ArgumentParser(description="Split a bgzipped VCF file into smaller files with a constant number of variant SNPs while including invariant sites.")
	parser.add_argument("input_vcf", help="Path to the input bgzipped VCF (.vcf.gz) file")
	parser.add_argument("output_prefix", help="Prefix for the output VCF files")
	parser.add_argument("snps_per_file", type=int, help="Number of variant SNPs per output VCF file")
	parser.add_argument("log_file", help="Path to the log file that will track the start and end positions of each part")

	# Parse arguments
	args = parser.parse_args()

	# Call the function with command line arguments
	slice_vcf(args.input_vcf, args.output_prefix, args.snps_per_file, args.log_file)
