There are two ASTRAL pipelines. 
1) The first uses ANGSD to create a fasta file from each bam file (using dofasta), then fasta files are sliced into windows of 100000kb using trimal. A phylogeny is inferred in each window using raxml. All the generated trees are then fed to ASTRAL. 
2) The second uses a VCF file. The vcf are sliced using bcftools. A phylogeny is inferred in each window using raxml-ng. All th egenerated trees are then fed to ASTRAL. 
