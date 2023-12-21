#!/bin/bash
#Author: Yacine Ben Chehida

#SBATCH --time=0-00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --account=BIOL-SPECGEN-2018

cd $1
COUNTER=1
if compgen -G "*.txt" > /dev/null; then rm *.txt; fi
for i in $(ls|grep -v txt)
	do 
	if [[ $COUNTER == 1 ]]; then
		cat  $i/Blast_hits.txt|awk 'NR==1' > Final_Blast_results.txt
		cat $i/Blast_hits.txt|awk 'NR>1'|sort -V -k15 >> Final_Blast_results.txt
#		cat $i/Blast_hits.txt|awk 'NR>1'|sort -V -k15|grep "Bacteria" >> Final_Contaminated_windows.txt 
		cat $i/Blast_hits.txt|awk '$14 ~ /^[0-9]+$/' >> Edirect_failed_windows.txt
		let COUNTER++
	else
		cat $i/Blast_hits.txt|awk 'NR>1'|sort -V -k15 >> Final_Blast_results.txt
#		cat $i/Blast_hits.txt|awk 'NR>1'|sort -V -k15| grep "Bacteria" >> Final_Contaminated_windows.txt
		cat $i/Blast_hits.txt|awk '$14 ~ /^[0-9]+$/' >> Edirect_failed_windows.txt
                let COUNTER++
	fi
done


../../Script/Taxonomy_corrector.sh > corrected_COMPLETE_Blast_results.txt
cat corrected_COMPLETE_Blast_results.txt |grep -v -E "Lepi|Insecta"|awk  'NF==16 {print}' > Final_list_to_check.txt
cat corrected_COMPLETE_Blast_results.txt |grep Virus >> Final_list_to_check.txt
