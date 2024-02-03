#!/bin/bash
# Author: Yacine Ben Chehida

#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=15G 
#SBATCH --account=BIOL-SPECGEN-2018 
#SBATCH --job-name=CALL_SV

module load Python/2.7.18-GCCcore-11.3.0-bare

MANTA_path="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/manta-1.6.0.centos6_x86_64/bin"
INPUT_path="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/8_SV/Inputs"
RESULTS_path="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/8_SV/Results"
SCRIPT_path="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/8_SV/Script"

cd $INPUT_path
BAM=$(for i in *bam; do echo -e "--bam $INPUT_path/$i"; done|perl -pe 's/\n/ /g')
echo $BAM

cd $SCRIPT_path

$MANTA_path/configManta.py --config $INPUT_path/configManta.py.ini --runDir $RESULTS_path/Combined

/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/8_SV/Results/Combined/runWorkflow.py -j 4
