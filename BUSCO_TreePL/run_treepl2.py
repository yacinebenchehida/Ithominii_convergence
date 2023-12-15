#! /usr/bin/env python
# Modified from https://github.com/sunray1/treepl/blob/master/run_treepl2.py
import sys, os, re

getfirstline = "grep -n 'PLACE THE LINES' %s"
getlastline = "wc -l %s"
getlinestoadd = "sed -n '%s,%sp' %s"

mkdircommd = "mkdir -p %s/%s"
cpconfig = "cp ../Inputs/treepl_step2.config %s/%s"
sandr = "sed -i 's/%s/%s/g' %s"
sandr_first = "sed -i '0,/%s/{s/%s/%s/}' %s"
runtreepl = "/users/ybc502/scratch/Conv_Evol/Tools/treePL/src/treePL %s"
mvcmnd = "mv %s %s"

#./prog.py tree_directory num out_dir
#num starts with 1
args = sys.argv
tree_directory = args[1]
num = args[2]
outdir = args[3]

#get treesuffix
tree_files = os.listdir(tree_directory)
ind_tree_file = tree_files[int(num)-1]
ind_tree_suff = ind_tree_file.split("tree")[1].split(".")[0]
tree_directory = re.sub("/","\/",tree_directory)
outdir = re.sub("/","\/",outdir)

#make tree directory and cp config file over
if os.path.exists(outdir + "/tree" + ind_tree_suff + ".out") == False:
	os.system(mkdircommd % (outdir, ind_tree_suff))
	os.system(cpconfig % (outdir, ind_tree_suff))

	#edit config file
	os.system(sandr % ("treefile =", "treefile = ..\/..\/" + tree_directory + "\/" + ind_tree_file, outdir + "/" + ind_tree_suff + "/" + "treepl_step2.config"))
	os.system(sandr % ("outfile =", "outfile = doesntmatter.tre", outdir + "/" + ind_tree_suff + "/" + "treepl_step2.config"))
	#get prime lines
	firstline = os.popen(getfirstline % ("../Results/step1_output/tree" + ind_tree_suff + "_step1.out")).read()
	startnum = int(firstline.split(":")[0]) + 1
	lastline = os.popen(getlastline % ("../Results/step1_output/tree" + ind_tree_suff + "_step1.out")).read()
	endnum = lastline.split()[0]

	primeadd = os.popen(getlinestoadd % (str(startnum), endnum, "../Results/step1_output/tree" + ind_tree_suff + "_step1.out"))
	for i in primeadd:
		os.system(sandr_first % ("#|", "#|", i.strip(), outdir + "/" + ind_tree_suff + "/" + "treepl_step2.config"))

	#runtreepl
	outdir = re.sub("\\\\","",outdir)
	os.chdir(outdir + "/" + ind_tree_suff)
	
	os.system(runtreepl % ("treepl_step2.config"))
	os.system(mvcmnd % ("cv.out", "../tree" + ind_tree_suff + ".out"))
