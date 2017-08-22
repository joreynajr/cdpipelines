#!/home/joreyna/anaconda2/envs/hla/bin/python

import os 
import argparse
import hlapy
import subprocess

parser = argparse.ArgumentParser(description='Make a diff comparing allele1 and allele2 alignments from ClustalW results.')
parser.add_argument('allele1', help='The name of the first allele.', type=str)
parser.add_argument('allele2', help='The name of the second allele.', type=str)
parser.add_argument('outdir', help='Path for the output dir.', type=str)
args = parser.parse_args()

# EXTRACTING command line data 
allele1 = args.allele1
allele2 = args.allele2 
outdir = args.outdir
gene = allele1.split('*')[0]

tmp_dir = '/home/joreyna/trash/'

gene_clustalw = '/repos/cardips-pipelines/HLA_Typing/sources/IPD_IMGT_HLA_Release_3.29.0_2017_07_27/alignments/{}_gen.txt'.format(gene)

allele1_name = allele1.replace('*', '_').replace(':', '_')
allele1_grep = allele1.replace('*', '\\\*').replace(':', '\:')
allele1_alignment = os.path.join(tmp_dir, '{}_alignment.txt'.format(allele1_name))
cmd = "grep {} {} | cut -f 3- -d ' '  | sed 's/ //' > {}".format(allele1_grep, gene_clustalw, allele1_alignment)
subprocess.call(cmd, shell=True)

allele2_name = allele2.replace('*', '_').replace(':', '_')
allele2_grep = allele2.replace('*', '\\\*').replace(':', '\:')
allele2_alignment = os.path.join(tmp_dir, '{}_alignment.txt'.format(allele2_name))
cmd = "grep {} {} | cut -f 3- -d ' '  | sed 's/ //' > {}".format(allele2_grep, gene_clustalw, allele2_alignment)
subprocess.call(cmd, shell=True)

#diff_fn = os.path.join(tmp_dir, 'diff_{}_v_{}.txt')
out_fn = os.path.join(outdir, '{}_v_{}.diff'.format(allele1_name, allele2_name))
cmd = 'diff -s {} {} > {}'.format(allele1_alignment, allele2_alignment, out_fn)
subprocess.call(cmd, shell=True)


#os.remove(allele1_alignment)
#os.remove(allele2_alignment)

print out_fn

