#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Script taken from: https://github.com/edgardomortiz/vcf2phylip/blob/master/vcf2phylip.py

"""
The script converts a collection of SNPs in VCF format into a PHYLIP, FASTA,
NEXUS, or binary NEXUS file for phylogenetic analysis. The code is optimized
to process VCF files with sizes >1GB. For small VCF files the algorithm slows
down as the number of taxa increases (but is still fast).

Any ploidy is allowed, but binary NEXUS is produced only for diploid VCFs.
"""

__author__      = "Edgardo M. Ortiz"
__credits__     = "Juan D. Palacio-Mejía"
__version__     = "2.8"
__email__       = "e.ortiz.v@gmail.com"
__date__        = "2021-08-10"

import argparse
import gzip
import random
import sys
from pathlib import Path

# Dictionary of IUPAC ambiguities for nucleotides
# '*' is a deletion in GATK, deletions are ignored in consensus, lowercase consensus is used when an
# 'N' or '*' is part of the genotype. Capitalization is used by some software but ignored by Geneious
# for example
AMBIG = {
    "A"    :"A", "C"    :"C", "G"    :"G", "N"    :"N", "T"    :"T",
    "*A"   :"a", "*C"   :"c", "*G"   :"g", "*N"   :"n", "*T"   :"t",
    "AC"   :"M", "AG"   :"R", "AN"   :"a", "AT"   :"W", "CG"   :"S",
    "CN"   :"c", "CT"   :"Y", "GN"   :"g", "GT"   :"K", "NT"   :"t",
    "*AC"  :"m", "*AG"  :"r", "*AN"  :"a", "*AT"  :"w", "*CG"  :"s",
    "*CN"  :"c", "*CT"  :"y", "*GN"  :"g", "*GT"  :"k", "*NT"  :"t",
    "ACG"  :"V", "ACN"  :"m", "ACT"  :"H", "AGN"  :"r", "AGT"  :"D",
    "ANT"  :"w", "CGN"  :"s", "CGT"  :"B", "CNT"  :"y", "GNT"  :"k",
    "*ACG" :"v", "*ACN" :"m", "*ACT" :"h", "*AGN" :"r", "*AGT" :"d",
    "*ANT" :"w", "*CGN" :"s", "*CGT" :"b", "*CNT" :"y", "*GNT" :"k",
    "ACGN" :"v", "ACGT" :"N", "ACNT" :"h", "AGNT" :"d", "CGNT" :"b",
    "*ACGN":"v", "*ACGT":"N", "*ACNT":"h", "*AGNT":"d", "*CGNT":"b",
    "*"    :"-", "*ACGNT":"N",
}

# Dictionary for translating biallelic SNPs into SNAPP, only for diploid VCF
# 0 is homozygous reference
# 1 is heterozygous
# 2 is homozygous alternative
GEN_BIN = {
    "./.":"?",
    ".|.":"?",
    "0/0":"0",
    "0|0":"0",
    "0/1":"1",
    "0|1":"1",
    "1/0":"1",
    "1|0":"1",
    "1/1":"2",
    "1|1":"2",
}


def extract_sample_names(vcf_file):
    """
    Extract sample names from VCF file
    """
    if vcf_file.lower().endswith(".gz"):
        opener = gzip.open
    else:
        opener = open
    sample_names = []
    with opener(vcf_file, "rt") as vcf:
        for line in vcf:
            line = line.strip("\n")
            if line.startswith("#CHROM"):
                record = line.split("\t")
                sample_names = [record[i].replace("./", "") for i in range(9, len(record))]
                break
    return sample_names


def is_anomalous(record, num_samples):
    """
    Determine if the number of samples in current record corresponds to number of samples described
    in the line '#CHROM'
    """
    return bool(len(record) != num_samples + 9)


def is_snp(record):
    """
    Determine if current VCF record is a SNP (single nucleotide polymorphism) as opposed to MNP
    (multinucleotide polymorphism)
    """
    # <NON_REF> must be replaced by the REF in the ALT field for GVCFs from GATK
    alt = record[4].replace("<NON_REF>", record[3])
    return bool(len(record[3]) == 1 and len(alt) - alt.count(",") == alt.count(",") + 1)


def num_genotypes(record, num_samples):
    """
    Get number of genotypes in VCF record, total number of samples - missing genotypes
    """
    missing = 0
    for i in range(9, num_samples + 9):
        if record[i].startswith("."):
            missing += 1
    return num_samples - missing


def get_matrix_column(record, num_samples, resolve_IUPAC):
    """
    Transform a VCF record into a phylogenetic matrix column with nucleotides instead of numbers
    """
    nt_dict = {str(0): record[3].replace("-","*").upper(), ".": "N"}
    # <NON_REF> must be replaced by the REF in the ALT field for GVCFs from GATK
    alt = record[4].replace("-", "*").replace("<NON_REF>", nt_dict["0"])
    alt = alt.split(",")
    for n in range(len(alt)):
        nt_dict[str(n+1)] = alt[n]
    column = ""
    for i in range(9, num_samples + 9):
        geno_num = record[i].split(":")[0].replace("/", "").replace("|", "")
        try:
            geno_nuc = "".join(sorted(set([nt_dict[j] for j in geno_num])))
        except KeyError:
            return "malformed"
        if resolve_IUPAC is False:
            column += AMBIG[geno_nuc]
        else:
            column += AMBIG[nt_dict[random.choice(geno_num)]]
    return column


def get_matrix_column_bin(record, num_samples):
    """
    Return an alignment column in NEXUS binary from a VCF record, if genotype is not diploid with at
    most two alleles it will return '?' as state
    """
