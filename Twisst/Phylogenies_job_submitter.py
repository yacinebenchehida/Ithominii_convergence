import os
import glob
import argparse

def create_slurm_script(group_number, vcf_files, phylo_method, job_dir, output_path, job_name):
    # Define the job script name
    job_script_name = f"job_{group_number}.sh"
    job_script_path = os.path.join(job_dir, job_script_name)
    
    # Create the SLURM job script
    with open(job_script_path, 'w') as job_file:
        job_file.write("#!/bin/bash\n")
        job_file.write(f"#SBATCH --job-name={job_name}_{group_number}\n")  # Use the provided job name
        job_file.write("#SBATCH --time=02:00:00\n")
        job_file.write("#SBATCH --ntasks=1\n") 
        job_file.write("#SBATCH --cpus-per-task=8\n")
        job_file.write("#SBATCH --mem=3G\n") 
        job_file.write("#SBATCH --account=BIOL-SPECGEN-2018\n")
        
        # Load necessary modules if needed (uncomment and modify)
        job_file.write("module load BCFtools/1.19-GCC-13.2.0\n")
        job_file.write("module load R/4.2.1-foss-2022a\n")
        job_file.write("module load RAxML-NG/1.0.2-gompi-2020b\n")
        
        job_file.write("cd $SLURM_SUBMIT_DIR\n") 
        
        # bgzip, tabix, convert to PHYLIP, and run phylogenetic analysis for each VCF file
        for vcf in vcf_files:
            prefix = os.path.splitext(os.path.basename(vcf))[0]  # Get the base name of the VCF file without extension
            result_dir = os.path.join(output_path, prefix)  # Create a specific output directory for this VCF
            
            # Create the directory for this VCF's results
            job_file.write(f"mkdir -p {result_dir}\n")
            
            job_file.write(f"bgzip -c {vcf} > {vcf}.gz\n")
            job_file.write(f"tabix -p vcf {vcf}.gz\n")

            # Convert VCF to PHYLIP
            job_file.write(f"python3 vcf2phylip.py -i {vcf} -r -m 0 --output-folder {result_dir} --output-prefix {prefix}\n")

            # Run the phylogenetic analysis based on the specified method
            if phylo_method == "NJ":
                job_file.write(f"Rscript NJ_tree.R {result_dir}/{prefix}*.phy {result_dir}\n")
                job_file.write(f"find {result_dir} ! -name '*.newick' -type f -delete\n")
            elif phylo_method == "ML":
                job_file.write(f"raxml-ng-mpi --search --msa {result_dir}/{prefix}*.phy --model GTR+G+I --threads 8 --force perf_threads --prefix {result_dir}/{prefix}\n")
                job_file.write(f"find {result_dir} ! -name '*.bestTree' -type f -delete\n")
            
            # Remove the VCF and index files after processing
            job_file.write(f"rm {vcf} {vcf}.gz {vcf}.gz.tbi\n")
        
    return job_script_path


def submit_jobs(vcf_directory, phylo_method, job_dir, output_path, job_name):
    # Create output directory for job scripts if it doesn't exist
    os.makedirs(job_dir, exist_ok=True)

    # List all VCF files
    vcf_files = sorted(glob.glob(os.path.join(vcf_directory, '*.vcf')))
    total_files = len(vcf_files)
    
    # Group files into batches of 500
    for i in range(0, total_files, 500):
        group_number = (i // 500) + 1  # Batch number
        batch_files = vcf_files[i:i + 500]  # Get the next 100 files (or less for the last batch)

        # Create SLURM job script for this batch
        job_script_path = create_slurm_script(group_number, batch_files, phylo_method, job_dir, output_path, job_name)

        # Submit the job to SLURM
        os.system(f"sbatch {job_script_path}")  # Submit job script


# Main function
if __name__ == "__main__":
    # Argument parser
    parser = argparse.ArgumentParser(description="Submit batch jobs for VCF processing and phylogenetic analysis.")
    parser.add_argument("vcf_directory", help="Directory containing VCF files.")
    parser.add_argument("phylo_method", choices=["NJ", "ML"], help="Phylogenetic method to use ('NJ' or 'ML').")
    parser.add_argument("job_dir", help="Directory to save SLURM job scripts.")
    parser.add_argument("output_path", help="Directory to save output PHYLIP files.")
    parser.add_argument("job_name", help="Name for the SLURM job.")

    # Parse arguments
    args = parser.parse_args()

    # Call the function with command line arguments
    submit_jobs(args.vcf_directory, args.phylo_method, args.job_dir, args.output_path, args.job_name)
