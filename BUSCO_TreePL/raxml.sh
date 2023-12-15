module load  bio/RAxML/8.2.12-intel-2019a-pthreads-avx2
module load numlib/Armadillo/10.7.5-foss-2022a

WD="/shared/biology/bioldata1/bl-kd684/yacine/Conv_Evol/Analyses_bio/BUSCO/Results"
SCRIPT="/shared/biology/bioldata1/bl-kd684/yacine/Conv_Evol/Analyses_bio/BUSCO/Scripts"
RAXML="/shared/biology/bioldata1/bl-kd684/yacine/Conv_Evol/Tools/standard-RAxML/raxmlHPC-PTHREADS"

cd $WD
raxmlHPC -s final.fna -m GTRGAMMAI -n BUSCO_PHYLO -p 1234 -T 4
sed -r 's/[0-9:.]//g' RAxML_bestTree.BUSCO_PHYLO > output_topo.tre
raxmlHPC -f j -m GTRGAMMAI -s final.fna -n BS -# 100 -b $RANDOM -T 4
for i in {0..99}; do raxmlHPC -f e -t output_topo.tre -m GTRGAMMAI -s final.fna.BS"$i" -n BS_"$i" -T 4; done
cat RAxML_result.BS* > RAxML_bootstrap.bootstrap_all.tre

/users/ybc502/scratch/Conv_Evol/Tools/phyx/src/pxrr -t RAxML_bootstrap.bootstrap_all.tre -g Plutella -o RAxML_bootstrap.bootstrap_all_rooted.tre

mkdir -p $WD/100_trees
split -l 1 RAxML_bootstrap.bootstrap_all_rooted.tre 100_trees/RAxML_bestTree.BSone_tree
for i in 100_trees/*; do mv $i $i.tre; done
