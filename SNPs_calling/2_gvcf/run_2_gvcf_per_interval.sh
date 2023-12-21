PATH_SCRIPT="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/2_gvcf/Script"
PATH_REF="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/0_Data/reference_genomes"
PATH_DATA=" /mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/0_Data/Fastq_files"

cd $PATH_DATA

for batch in NEOF_November_2023
do cd $batch
	for i in Sample* 
	do 
		ref=$(grep $i $PATH_SCRIPT/ref_subset.txt| awk '{print $2}' )
		if [[ $ref == "Mechanitis_messenoides" ]]; then
			echo $i $(ls $PATH_REF/$ref/Mech.mess.mm15-0374.mother.fa)
			sbatch $PATH_SCRIPT/2_gvcf_per_interval.sh $i $ref $batch	 
		fi
		if [[ $ref == "Melinaea_marsaeus" ]]; then
			cd $PATH_SCRIPT
			echo $i $(ls $PATH_REF/$ref/CAM015037_Melinaea_marsaeus_rileyi.fa)
			sbatch $PATH_SCRIPT/2_gvcf_per_interval.sh $i $ref $batch
		fi
		if [[ $ref == "Melinaea_menophilus" ]]; then
			echo $i $(ls $PATH_REF/$ref/CAM015033_Melinaea_menophilus_ssp_nov.fa)
			sbatch $PATH_SCRIPT/2_gvcf_per_interval.sh $i $ref $batch
		fi
		if [[ $ref == "Heliconius_demeter" ]]; then
			echo $i $(ls $PATH_REF/$ref/Hdem.assembly.v1.5.fasta)
			sbatch $PATH_SCRIPT/2_gvcf_per_interval.sh $i $ref $batch
		fi
		if [[ $ref == "Heliconius_neruda" ]]; then
			echo $i $(ls $PATH_REF/$ref/Haoe.assembly.v1.2.fasta)
			sbatch $PATH_SCRIPT/2_gvcf_per_interval.sh $i $ref $batch
		fi
		if [[ $ref == "Chetone_histrio" ]]; then
			cd $PATH_SCRIPT
      echo $i $(ls $PATH_REF/$ref/chetone_histrio_mtDNA_05_02_23.fasta)
      sbatch $PATH_SCRIPT/2_gvcf_per_interval.sh $i $ref $batch	
			cd $PATH_DATA
    fi
	done
	cd $PATH_DATA
done
