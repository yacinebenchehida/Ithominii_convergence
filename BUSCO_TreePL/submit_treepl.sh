#!/bin/bash


TREEPL="/shared/biology/bioldata1/bl-kd684/yacine/Conv_Evol/Tools/treePL/src/treePL"

module purge
module load toolchain/foss/2020a
module load system/ADOL-C/2.7.2-gompi-2020a
module load numlib/NLopt/2.6.1-GCCcore-9.3.0
module load  bio/Beast/2.5.2-foss-2018b
module load numlib/NLopt/2.7.1-GCCcore-12.2.0

mkdir -p ../Results/step1_output
for i in {1..100}; do python3 ./run_treepl.py ../Results/100_trees $i ../Results/step1_output; done
for i in {1..100}; do ./run_treepl2.py ../Results/100_trees $i ../Results/step2.1_output; done
for i in {1..100}; do python3 ./run_treepl3.py ../Results/100_trees $i ../Results/step3_output; done


cd ../Results/step3_output
cat *.tre > treeall_dated.tre
treeannotator -burnin 10 treeall_dated.tre out.txt
