import argparse
import subprocess
import os

def extract_variant_positions(input_vcf, variant_positions_file):
    """Extract the positions of variant sites from the VCF file."""
    command = f"bcftools query -f '%CHROM\\t%POS\\n' -v snps {input_vcf} > {variant_positions_file}"
    subprocess.run(command, shell=True, check=True)

def split_vcf(input_vcf, variant_positions_file, chunk_size, output_dir):
    """Split the VCF file into smaller chunks each containing `chunk_size` variants."""
    with open(variant_positions_file, 'r') as file:
        positions = file.readlines()
    
    total_variants = len(positions)
    start = 0

    for i in range(0, total_variants, chunk_size):
        end = i + chunk_size
        if end > total_variants:
            end = total_variants

        start_pos = positions[i].strip().split('\t')[1]
        end_pos = positions[end - 1].strip().split('\t')[1]

        output_vcf = os.path.join(output_dir, f"chunk_{i + 1}_to_{end}.vcf")

        command = f"bcftools view --regions {start_pos}-{end_pos} {input_vcf} -o {output_vcf}"
        subprocess.run(command, shell=True, check=True)
        print(f"Created {output_vcf}")

def main():
    parser = argparse.ArgumentParser(description="Split a VCF file into smaller chunks, each containing a specified number of variant sites.")
    parser.add_argument('-i', '--input', required=True, help="Input VCF file")
    parser.add_argument('-c', '--chunk_size', type=int, required=True, help="Number of variants per chunk")
    parser.add_argument('-o', '--output_dir', required=True, help="Directory to store the output VCF chunks")
    args = parser.parse_args()

    input_vcf = args.input
    chunk_size = args.chunk_size
    output_dir = args.output_dir

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    variant_positions_file = os.path.join(output_dir, "variant_positions.txt")
    
    # Step 1: Extract variant positions
    extract_variant_positions(input_vcf, variant_positions_file)

    # Step 2: Split the VCF into smaller chunks
    split_vcf(input_vcf, variant_positions_file, chunk_size, output_dir)

    # Cleanup
    os.remove(variant_positions_file)
    print("Finished splitting VCF.")

if __name__ == "__main__":
    main()
