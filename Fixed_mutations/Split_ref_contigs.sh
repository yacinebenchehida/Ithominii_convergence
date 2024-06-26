RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/Fixed_mutations/Split_reference_genomes/Results"
WD="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/0_Data/reference_genomes"

cd $WD

for i in $1
	do cd $i
                FASTA_REF=$(ls $WD/"$i"|grep -E "*.fa$|*.fasta$")
                mkdir -p $RESULTS/$i
                cat $FASTA_REF|grep ">"|sed 's/>//g' |awk '{print $1}' > $RESULTS/$i/contigs.txt
                cd $RESULTS/$i
                Num_Contigs=$(cat contigs.txt| wc -l)
                echo $Num_Contigs
                if [[ $Num_Contigs -lt 99 ]]
                then
                    	file=contigs.txt; nsplit=$Num_Contigs; len=$(wc -l < $file); split -l$(($len/$nsplit)) "$file" contigs_"$i"_part_0 -a 3 -d
                else
                    	file=contigs.txt; nsplit=99; len=$(wc -l < $file); split -l$(($len/$nsplit)) "$file" contigs_"$i"_part_0 -a 3 -d
                fi
        cd $WD
done
