#Load modules
module load  bio/RAxML/8.2.12-intel-2019a-pthreads-avx2
module load numlib/Armadillo/11.4.3-foss-2022b

#Set paths
WD="/shared/biology/bioldata1/bl-kd684/yacine/Conv_Evol/Analyses_bio/BUSCO/Results"
SCRIPT="/shared/biology/bioldata1/bl-kd684/yacine/Conv_Evol/Analyses_bio/BUSCO/Scripts"
RAXML="/shared/biology/bioldata1/bl-kd684/yacine/Conv_Evol/Tools/standard-RAxML/raxmlHPC-PTHREADS"

cd $WD

# Run raxml
$RAXML -s final.fna -m GTRGAMMAI -n BUSCO_PHYLO -p 1234 -T 4
sed -r 's/[0-9:.]//g' RAxML_bestTree.BUSCO_PHYLO > output_topo.tre
$RAXML -f j -m GTRGAMMAI -s final.fna -n BS -# 100 -b $RANDOM -T 4
for i in {0..99}; do $RAXML -f e -t output_topo.tre -m GTRGAMMAI -s final.fna.BS"$i" -n BS_"$i" -T 4; done
cat RAxML_result.BS* > RAxML_bootstrap.bootstrap_all.tre

/shared/biology/bioldata1/bl-kd684/yacine/Conv_Evol/Tools/phyx/src -t RAxML_bootstrap.bootstrap_all.tre -g Plutella_xylostella -o RAxML_bootstrap.bootstrap_all_rooted.tre

mkdir -p $WD/100_trees
split -l 1 RAxML_bootstrap.bootstrap_all_rooted.tre 100_trees/RAxML_bestTree.BSone_tree
for i in 100_trees/*; do mv $i $i.tre; done
