#! /bin/bash

####################
# Set useful paths #
####################
SCRIPT="/shared/biology/bioldata1/bl-kd684/yacine/Conv_Evol/Analyses_bio/BUSCO/Scripts"
PAL2NAL="/shared/biology/bioldata1/bl-kd684/yacine/Conv_Evol/Tools/pal2nal.v14/pal2nal.pl"

#########################
# Set working directory #
#########################
cd /shared/biology/bioldata1/bl-kd684/yacine/Conv_Evol/Analyses_bio/BUSCO/Results

########################################
# Identify common genes to all species #
########################################
for i in Mechanitis_mazaeus Melinaea_isocomma Melinaea_mothone bicyclus Biston Bombyx Monarch Pieris Plutella Chetone_histrio  Mechanitis_messenoides  Melinaea_marsaeus  Melinaea_menophilus Heliconius_erato Heliconius_melpomene Heliconius_numata Heliconius_pardalinus Hypothyris_anastasia Ithomia_salapia; do cd $i/$i/run_lepidoptera_odb10/busco_sequences/single*; ls *faa > ../../../../../"$i"_busco_gene.txt ; cd ../../../../..; done
paste *txt > common_busco_genes.txt
Rscript $SCRIPT/common_genes.R common_busco_genes.txt
rm common_busco_genes.txt *_busco_gene.txt

#######################################################
# Extract common sequences to all species in a folder #
#######################################################
[ -e concatanated_genes_prot ] && rm -r concatanated_genes_prot
[ -e concatanated_genes_nucl ] && rm -r concatanated_genes_nucl
mkdir -p concatanated_genes_prot concatanated_genes_nucl

for i in  Mechanitis_mazaeus Melinaea_isocomma Melinaea_mothone bicyclus Biston Bombyx Monarch Pieris Plutella Chetone_histrio Mechanitis_messenoides  Melinaea_marsaeus  Melinaea_menophilus Heliconius_erato Heliconius_melpomene Heliconius_numata Heliconius_pardalinus Hypothyris_anastasia Ithomia_salapia
do 
	cat common_genes.txt|perl -pe 's/\.faa//g'| while read line; do cat $i/$i/run_lepidoptera_odb10/busco_sequences/single*/$line*faa >> concatanated_genes_prot/$line.faa; cat $i/$i/run_lepidoptera_odb10/busco_sequences/single*/$line*fna >> concatanated_genes_nucl/$line.fna
	done
done

#######################################################################################
# Use muscle to align proteins and then pal2nal to reverse translate it to nucleotide #
#######################################################################################
mkdir -p aligned_concatanated_genes_prot aligned_concatanated_genes_nucl
cd concatanated_genes_prot
#module purge
#module load bio/libMuscle/3.7-foss-2016b
MUSCLE="/shared/biology/bioldata1/bl-kd684/yacine/Conv_Evol/Tools/muscle/muscle3.8.31_i86linux32"

for i in *faa
	do 
	NAME=$(echo $i| perl -pe 's/faa/aln/g')
	NAME2=$(echo $i| perl -pe 's/faa/fna/g')
	cat $i |grep ">"|perl -pe 's/>//g' > order_"$NAME"
	$MUSCLE  -in $i -out ../aligned_concatanated_genes_prot/aligned_"$NAME"
	python2 $SCRIPT/reorder_fasta.py ../aligned_concatanated_genes_prot/aligned_"$NAME" order_"$NAME"  > ../aligned_concatanated_genes_prot/tmp_"$i"
	rm order_"$NAME" ../aligned_concatanated_genes_prot/aligned_"$NAME"
	mv ../aligned_concatanated_genes_prot/tmp_"$i" ../aligned_concatanated_genes_prot/aligned_"$NAME"
	$PAL2NAL ../aligned_concatanated_genes_prot/aligned_"$NAME" ../concatanated_genes_nucl/$NAME2 -output fasta -codontable 1 > ../aligned_concatanated_genes_nucl/aligned_"$NAME2"
	done

####################################################
# remove sequences with an excessive number of "-" #
####################################################
module purge
module load bio/Biopython/1.79-foss-2022a
cd ../aligned_concatanated_genes_nucl
for i in *fna; do INDELS=$(python -W ignore $SCRIPT/indels_count.py ../aligned_concatanated_genes_nucl/$i); if [[ "$INDELS" -gt 2 ]]; then echo $i; echo $INDELS; rm ../aligned_concatanated_genes_nucl/$i; fi; done

#########################################################################################################################
# replace fasta headers so they are all the same (essential later the sequences can be merged into a single fasta file) #
#########################################################################################################################
cd ../aligned_concatanated_genes_nucl
for j in *.fna; do for i in Mechanitis_mazaeus Melinaea_isocomma Melinaea_mothone bicyclus Biston Bombyx Monarch Pieris Plutella Chetone_histrio Mechanitis_messenoides  Melinaea_marsaeus  Melinaea_menophilus Heliconius_erato Heliconius_melpomene Heliconius_numata Heliconius_pardalinus Hypothyris_anastasia Ithomia_salapia; do sed -i -E ':a;N;$!ba; s/(>)([a-zA-Z0-9\_\.]+):([0-9]+)-([0-9]+)/\1'"${i}"'/1' $j; done; done

####################################################
# Concatenate all the genes in a single fasta file #
####################################################
#module load bio/SeqKit/0.12.0
SEQKIT="/shared/biology/bioldata1/bl-kd684/yacine/Conv_Evol/Tools/SeqKit/seqkit"
$SEQKIT concat *.fna  > ../final.fna
