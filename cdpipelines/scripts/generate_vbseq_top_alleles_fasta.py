#!/home/joreyna/anaconda2/envs/hla/bin/python

import os 
import argparse
import pandas as pd 
import hlapy
import subprocess
import warnings
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser(description='Generate a lower resolution HLA type fasta file')

parser.add_argument('sample_name', help='The name of the sample.', type=str)
parser.add_argument('resolution', help='HLA type resolution to consider', type=int)
parser.add_argument('num_top_alleles', help='Number of top alleles to consider.', type=int)
parser.add_argument('hla_bed', help='Bed file derived from hla_gen.fasta', type=str)
parser.add_argument('hla_fasta', help='HLA reference file from IPD/IMGT HLA database, (e.g. hla_gen.fasta)', type=str)
parser.add_argument('output_dir', type=str)
parser.add_argument('raw_vbseq', help='Raw output from HLA-VBSeq perl parser.', type=str)

args = parser.parse_args()

sample_name = args.sample_name
top = args.num_top_alleles 
resolution = args.resolution  
hla_bed = args.hla_bed
hla_fasta = args.hla_fasta
output_dir = args.output_dir
raw_vbseq = args.raw_vbseq 

# If the resolution value has been set to something other than 8 then we will 
# filter the allele sequences for the lower resolution types. Else we will filter for 
# just the top 5 resolutions types which correspond to the 8-digit resolution.
output_fasta = os.path.join(output_dir, '{}_res_{}_top_{}.fasta'.format(sample_name, resolution, top))
output_bed = os.path.join(output_dir, '{}_res_{}_top_{}.bed'.format(sample_name, resolution, top))

# DETERMINE top alleles for each HLA gene 
print('Determining the top allele for each HLA gene.')
raw_vbseq = pd.read_table(raw_vbseq, header=None, names=['allele', 'mean_read_depth'])
raw_vbseq['gene'] = [x[0] for x in raw_vbseq.allele.str.split('*')]    
gene_grps = raw_vbseq.groupby('gene')
top_alleles = gene_grps.apply(\
      lambda x: x.sort_values('mean_read_depth', ascending=False).iloc[0:top])
top_alleles.reset_index(drop=True, inplace=True)
top_alleles.set_index('allele', inplace=True)

if resolution != 8: 
    # IDENTIFY the lower digit resolution of the top alleles
    print('Identifying the lower resolution, {}, of the top alleles.'.format(resolution))
    resolution_name = '{}_digit'.format(resolution)
    top_alleles[resolution_name] = \
        [hlapy.get_lower_resolution_by_trimming(x, resolution) for x in 
             top_alleles.index.tolist()]
    query_alleles = set(top_alleles[resolution_name].tolist())
    query_genes = set([x.split('*')[0] for x in query_alleles])

    # EXTRACT top lower resolution alleles
    print('Extracting the lower resolution alleles.') 
    filter_allele_list = []
    for x in open(hla_bed):
        accession_id, start, end, allele = x.split()
        gene = allele.split('*')[0]

        lower_resolution_type = \
            hlapy.get_lower_resolution_by_trimming(allele, resolution)
            
        if lower_resolution_type in query_alleles or gene not in query_genes:
            filter_allele_list.append([accession_id, start, end, allele ])
else:
    filter_allele_list = []
    top_alleles_list = top_alleles.index.tolist()
    query_genes = set([x.split('*')[0] for x in top_alleles_list])
    for x in open(hla_bed):
        accession_id, start, end, allele = x.split()
        gene = allele.split('*')[0]
        if allele in top_alleles_list or gene not in query_genes:
            filter_allele_list.append([accession_id, start, end, allele ])
        
# GENERATE a bed file of alleles
print('Generating a bed file of the lower resolutions')
with open(output_bed, 'w') as bed:
    for accession_id, start, end, allele in filter_allele_list:
        name = ' '.join([accession_id, allele, end, 'bp'])
        allele = '\t'.join([accession_id, start, end, name])
        bed.write( allele + '\n')

# GENERATE a fasta file of the alleles
cmd = 'bedtools getfasta -fi {} -bed {} -fo {} -name'.\
    format(hla_fasta, output_bed, output_fasta)
subprocess.check_output(cmd, shell=True)
