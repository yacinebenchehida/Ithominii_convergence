import argparse
import subprocess
import pysam
import math
import os
import matplotlib.pyplot as plt

def read_population_file(pop_file):
    """Reads the population file and returns a dictionary mapping samples to species."""
    pop_dict = {}
    with open(pop_file, 'r') as f:
        for line_num, line in enumerate(f, start=1):
            parts = line.strip().split("\t")  # Expecting tab-separated values
            if len(parts) != 2:
                print(f"Warning: Skipping line {line_num} in {pop_file} (Expected 2 tab-separated values, got {len(parts)}): '{line.strip()}'")
                continue
            sample, species = parts
            pop_dict[sample] = species
    return pop_dict

def extract_chromosome(vcf_file, chromosome, output_chromosome_file):
    """Extracts a specific chromosome from the VCF file using bcftools."""
    command = f"module load BCFtools/1.19-GCC-13.2.0 && bcftools view -r {chromosome} {vcf_file} -Oz -o {output_chromosome_file}"
    subprocess.run(command, shell=True, check=True)

def extract_transpolymorphic_sites(vcf_file, pop_file, chromosome, output_file, working_dir):
    """Extract sites polymorphic in all species from the VCF based on the population file."""
    
    os.makedirs(working_dir, exist_ok=True)
    output_chromosome_file = os.path.join(working_dir, "extracted_chromosome.vcf.gz")
    extract_chromosome(vcf_file, chromosome, output_chromosome_file)
    
    vcf = pysam.VariantFile(output_chromosome_file)
    populations = read_population_file(pop_file)
    species_samples = {species: [] for species in set(populations.values())}

    for sample, species in populations.items():
        species_samples[species].append(sample)

    transpolymorphic_sites = []
    all_polymorphic_sites = set()
    total_sites = set()

    for record in vcf.fetch():
        relevant_samples = [s for s in record.samples if s in populations]
        species_allele_counts = {species: {} for species in species_samples.keys()}

        for sample in relevant_samples:
            species = populations[sample]
            alleles = record.samples[sample].alleles
            if None in alleles:
                continue
            
            for allele in alleles:
                if allele not in species_allele_counts[species]:
                    species_allele_counts[species][allele] = 0
                species_allele_counts[species][allele] += 1

        num_species = len(species_allele_counts)
        polymorphic_in_required_taxa = False
        
        if num_species > 3:
            polymorphic_in_required_taxa = sum(len(allele_counts) >= 2 for allele_counts in species_allele_counts.values()) >= 4
        elif num_species <= 3:
            polymorphic_in_required_taxa = all(len(allele_counts) >= 2 for allele_counts in species_allele_counts.values())

        if any(len(allele_counts) >= 2 for allele_counts in species_allele_counts.values()):
            all_polymorphic_sites.add(record.pos)
        
        total_sites.add(record.pos)

        if polymorphic_in_required_taxa and all(all(count >= 4 for count in allele_counts.values()) for allele_counts in species_allele_counts.values()):
            output_line = f"{record.pos}"
            for species, allele_counts in species_allele_counts.items():
                if allele_counts:
                    alleles_info = " ".join([f"{allele}:{count}" for allele, count in allele_counts.items()])
                    output_line += f" {species} [{alleles_info}]"
            transpolymorphic_sites.append(output_line)

    output_file = os.path.join(working_dir, os.path.basename(output_file))
    with open(output_file, 'w') as f:
        for line in transpolymorphic_sites:
            f.write(line + "\n")

    return total_sites, all_polymorphic_sites

def count_sites_in_windows(output_file, total_sites, all_polymorphic_sites, window_size):
    with open(output_file, 'r') as f:
        positions = [int(line.split()[0]) for line in f.readlines()]

    min_position = min(positions)
    max_position = max(positions)
    total_windows = math.ceil((max_position - min_position) / window_size)

    window_counts = [0] * total_windows
    available_site_ratios = []
    polymorphic_site_ratios = []
    mid_points = []

    for pos in positions:
        window_index = (pos - min_position) // window_size
        window_counts[window_index] += 1

    window_output_file = os.path.join(os.path.dirname(output_file), "window_counts.txt")
    with open(window_output_file, 'w') as f:
        f.write(f"Transpolymorphic site counts per {window_size}-bp window:\n")
        for i, count in enumerate(window_counts):
            start = min_position + i * window_size
            end = start + window_size - 1
            mid_point = (start + end) // 2
            mid_points.append(mid_point)

            available_sites = len([s for s in total_sites if start <= s <= end])
            polymorphic_sites = len([s for s in all_polymorphic_sites if start <= s <= end])

            ratio_available = count / available_sites if available_sites else 0
            ratio_polymorphic = count / polymorphic_sites if polymorphic_sites else 0

            available_site_ratios.append(ratio_available)
            polymorphic_site_ratios.append(ratio_polymorphic)

            f.write(f"Window {i + 1}: {start}-{end} bp, Sites: {count}, Ratio_Available: {ratio_available:.4f}, Ratio_Polymorphic: {ratio_polymorphic:.4f}\n")
    
    print(f"Sliding window results written to {window_output_file}")

    return mid_points, window_counts, available_site_ratios, polymorphic_site_ratios

def plot_results(mid_points, transpolymorphic_counts, available_ratios, polymorphic_ratios, working_dir):
    # Plot Transpolymorphic Sites
    plt.figure(figsize=(10, 6))
    plt.plot(mid_points, transpolymorphic_counts, label="Transpolymorphic Sites", color="blue", marker="o")
    plt.xlabel("Window Mid Position (bp)")
    plt.ylabel("Number of Transpolymorphic Sites")
    plt.title("Transpolymorphic Sites by Window")
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(working_dir, "transpolymorphic_sites_plot.pdf"))
    plt.close()

    # Plot Ratio to Available Sites
    plt.figure(figsize=(10, 6))
    plt.plot(mid_points, available_ratios, label="Ratio to Available Sites", color="orange", marker="x")
    plt.xlabel("Window Mid Position (bp)")
    plt.ylabel("Ratio to Available Sites")
    plt.title("Ratio of Transpolymorphic Sites to Available Sites by Window")
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(working_dir, "available_sites_ratio_plot.pdf"))
    plt.close()

    # Plot Ratio to Polymorphic Sites
    plt.figure(figsize=(10, 6))
    plt.plot(mid_points, polymorphic_ratios, label="Ratio to Polymorphic Sites", color="green", marker="s")
    plt.xlabel("Window Mid Position (bp)")
    plt.ylabel("Ratio to Polymorphic Sites")
    plt.title("Ratio of Transpolymorphic Sites to Polymorphic Sites by Window")
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(working_dir, "polymorphic_sites_ratio_plot.pdf"))
    plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract transpolymorphic sites and analyze in windows.")
    parser.add_argument("-v", "--vcf_file", required=True, help="Path to the input VCF file.")
    parser.add_argument("-p", "--pop_file", required=True, help="Path to the population file.")
    parser.add_argument("-c", "--chromosome", required=True, help="Chromosome to analyze.")
    parser.add_argument("-o", "--output_file", required=True, help="Path to the output file for transpolymorphic sites.")
    parser.add_argument("-w", "--window_size", type=int, required=True, help="Window size in base pairs for counting transpolymorphic sites.")
    parser.add_argument("-d", "--working_dir", required=True, help="Directory to save intermediate VCF files.")

    args = parser.parse_args()

    total_sites, all_polymorphic_sites = extract_transpolymorphic_sites(args.vcf_file, args.pop_file, args.chromosome, args.output_file, args.working_dir)
    mid_points, window_counts, available_ratios, polymorphic_ratios = count_sites_in_windows(
        os.path.join(args.working_dir, os.path.basename(args.output_file)), total_sites, all_polymorphic_sites, args.window_size
    )

    plot_results(mid_points, window_counts, available_ratios, polymorphic_ratios, args.working_dir)
