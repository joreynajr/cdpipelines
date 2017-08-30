#!/home/joreyna/anaconda2/envs/hla/bin/python

import os 
import argparse
import hlapy
import subprocess

testing = False 
def determine_unknown_indel_length(seq1, seq2):  
    """
    Determine how much unknown sequence is found in one allele versus another. 
    
    Parameters
    ----------
    seq1: str
        First sequence
    
    seq2: str 
        Second sequence 
        
    Returns
    -------
    
    """
    
    bp1 = seq1[0]
    bp2 = seq2[0]
    skip = 0 

    if len(seq1) != len(seq2):
        return 'error - cannot compare', 0, 0

    if bp1 != '*' and bp2 != '*':
        return 'no diff', 0, 0

    if bp1 != '*' and bp2 == '*':
        five_prime_diff = 'deletion'

    elif bp1 == '*' and bp2 != '*':
        five_prime_diff = 'insertion'

    else:
        five_prime_diff = 'tbd'

    variant_len = 0
    for i in range(len(seq1)):

        bp1 = seq1[i]
        bp2 = seq2[i]

        if five_prime_diff == 'deletion':
            if bp1 != '*' and bp2 == '*':
                variant_len +=1 

            elif bp2 != '*':
                break 

        elif five_prime_diff == 'insertion':
            if bp1 == '*' and bp2 != '*':
                variant_len +=1 

            elif bp1 != '*':
                break

        elif five_prime_diff == 'tbd':

            if bp1 == '*' and bp2 != '*':
                five_prime_diff = 'insertion'
                variant_len += 1 

            elif bp1 != '*' and bp2 == '*':
                five_prime_diff = 'deletion'
                variant_len += 1    

            elif bp1 != '*' and bp2 != '*':
                five_prime_diff = 'no diff'

            else:
                skip += 1
                
    return five_prime_diff, variant_len, skip

def report_unknown_indel_results(allele1, allele2, \
             five_prime_diff, five_prime_variant_len, three_prime_diff, three_prime_variant_len):
    
    """
    Reports the sequence alignment between allele1 and allele2. 
    
    Parameters
    ----------
    allele1: str
        Full HLA allele name of the first allele. 
    allele2: str
        Full HLA allele name of the second allele. 
    five_prime_diff: str
        Flag that denotes the differences between the first and second allele at the 5' end end.
    five_prime_variant_len: int
        The length of the variant found at the 5' end prime end. 
    three_prime_diff: str
        Flag that denotes the differences between the first and second allele at the 5' end end.
    three_prime_variant_len: int
        The length of the variant found at the 3' end prime end. 
    
    Returns
    -------
    relative_sequence_alignment: str
        Report describing the differences between the alleles. 
    
    """
    
    if five_prime_diff == 'error - cannot compare' or three_prime_diff == 'error - cannot compare':
        relative_sequence_alignment = '{} and {} cannot be compared using the MSF file'.format(allele1, allele2)
        return relative_sequence_alignment
    
    if five_prime_diff == 'deletion':
        if three_prime_diff == 'deletion':
            relative_sequence_alignment = "{} is a superset of {}; 5' end {}bp overhang; 3' end {}bp overhang".\
                format(allele1, allele2, five_prime_variant_len, three_prime_variant_len)

        elif three_prime_diff == 'insertion':
            relative_sequence_alignment = "{} has a staggered 3' end overlap with {}; 5' end {}bp overhang; 3' {}bp overhang".\
                format(allele1, allele2, five_prime_variant_len, three_prime_variant_len)  
            
        elif three_prime_diff == 'no diff':
            relative_sequence_alignment = "{} is a superset of {}; 5' end {}bp overhang; 3' end {}bp overhang".\
                format(allele1, allele2, five_prime_variant_len, three_prime_variant_len) 
            
    elif five_prime_diff == 'insertion':
        if three_prime_diff == 'insertion':
            relative_sequence_alignment = "{} is a subset of {}; 5' end {}bp overhang; 3' end {}bp overhang".\
                format(allele1, allele2, five_prime_variant_len, three_prime_variant_len)  
                
        elif three_prime_diff == 'deletion':
            relative_sequence_alignment = "{} has a staggered 5' end overlap with {}; 5' {}bp overhang; 3' end {}bp overhang".\
                format(allele1, allele2, five_prime_variant_len, three_prime_variant_len) 
                
        elif three_prime_diff == 'no diff':
            relative_sequence_alignment = "{} is a subset of {}; 5' end {}bp overhang; 3' end {}bp overhang".\
                format(allele1, allele2, five_prime_variant_len, three_prime_variant_len) 
            
            
    elif five_prime_diff == 'no diff':
        if three_prime_diff == 'deletion':
            relative_sequence_alignment = "{} is a superset of {}; 5' end {}bp overhang; 3' end {}bp overhang".\
                format(allele1, allele2, five_prime_variant_len, three_prime_variant_len) 
        
        elif three_prime_diff == 'insertion':
            relative_sequence_alignment = "{} is a subset of {}; 5' end {}bp overhang; 3' end {}bp overhang".\
                format(allele1, allele2, five_prime_variant_len, three_prime_variant_len) 
                
        elif three_prime_diff == 'no diff':
            relative_sequence_alignment = "{} is aligned with {}; 5' end {}bp overhang; 3' end {}bp overhang".\
                format(allele1, allele2, five_prime_variant_len, three_prime_variant_len) 


            
    return relative_sequence_alignment

def generate_unknown_indel_report(allele1_fn, allele2_fn):
    # LOADING the sequence data 
    linearize_seq = lambda x: x.replace('\n', '').replace(' ', '').replace('\t', '').replace('|', '')
    seq1 = linearize_seq(open(allele1_fn).read())
    seq2 = linearize_seq(open(allele2_fn).read())
    
    

    # COMPARING the two sequences for any unknown sequences between them.  
    five_prime_diff, five_prime_variant_len, five_prime_skip = determine_unknown_indel_length(seq1, seq2)
    three_prime_diff, three_prime_variant_len, three_prime_skip = determine_unknown_indel_length(list(reversed(seq1)), list(reversed(seq2)))

    return report_unknown_indel_results(allele1, allele2, \
                 five_prime_diff, five_prime_variant_len, three_prime_diff, three_prime_variant_len)



# Determining the number of variants
def determine_num_snps_ins_dels(seq1, seq2, five_prime_variant_len, five_prime_skip, three_prime_variant_len, three_prime_skip):  
    core_seq1 = seq1[(five_prime_variant_len + five_prime_skip) : len(seq1) - (three_prime_variant_len + three_prime_skip)]
    core_seq2 = seq2[(five_prime_variant_len + five_prime_skip) : len(seq2) - (three_prime_variant_len + three_prime_skip)]

    #print 'core_seq1:', core_seq1
    #print 'core_seq2:', core_seq2

    ins = 0
    dels = 0
    snps = 0
    for i in range(len(core_seq1)):

        core1 = core_seq1[i]
        core2 = core_seq2[i]

        if core1 == core2:
            continue

        else:

            if core1 == '.':
                ins += 1 
            elif core2 == '.':
                dels += 1
            else:
                snps += 1 

    return snps, ins, dels

if testing == True:

    print 'Testing the function to generate unknown indels in UTR sequences.'
    
    # Testing subset 
    allele1 = 'DRB1_11_01_01_03'
    allele2 = 'DRB1_11_01_01_01'

    allele1_fn = 'test/DRB1_11_01_01_03_alignment.txt'.format(allele1)
    allele2_fn = 'test/DRB1_11_01_01_01_alignment.txt'.format(allele2)

    print generate_unknown_indel_report(allele1_fn, allele2_fn)

    # Testing superset  
    allele1 = 'DRB1_11_01_01_01'
    allele2 = 'DRB1_11_01_01_03'

    allele1_fn = 'test/DRB1_11_01_01_01_alignment.txt'.format(allele1)
    allele2_fn = 'test/DRB1_11_01_01_03_alignment.txt'.format(allele2)

    print generate_unknown_indel_report(allele1_fn, allele2_fn)

    # Testing exact match
    # Testing superset  
    allele1 = 'DRB1_11_01_01_01'
    allele2 = 'DRB1_11_01_01_01'

    allele1_fn = 'test/DRB1_11_01_01_01_alignment.txt'.format(allele1)
    allele2_fn = 'test/DRB1_11_01_01_01_alignment.txt'.format(allele2)

    print generate_unknown_indel_report(allele1_fn, allele2_fn)

if __name__ == '__main__': 
    parser = argparse.ArgumentParser(description='Make a diff comparing allele1 and allele2 alignments from ClustalW results.')
    parser.add_argument('allele1', help='The name of the first allele.', type=str)
    parser.add_argument('allele2', help='The name of the second allele.', type=str)
    parser.add_argument('--outdir', help='Path for the output dir.', type=str)
    parser.add_argument('--report', help='If set prints the data in report style.', action='store_true')
    parser.add_argument('--test_dir', help='If set serves as the directory for intermediate data outut.', type=str)
    args = parser.parse_args()

    # EXTRACTING command line data 
    allele1 = args.allele1
    allele2 = args.allele2 
    outdir = args.outdir
    test_dir = args.test_dir
    gene = allele1.split('*')[0]
    report = args.report

    tmp_dir = '/home/joreyna/trash/'

    gene_clustalw = '/repos/cardips-pipelines/HLA_Typing/sources/IPD_IMGT_HLA_Release_3.29.0_2017_07_27/alignments/{}_gen.txt'.format(gene)

    # EXTRACTING allele1 sequence data 
    allele1_name = allele1.replace('*', '_').replace(':', '_')
    allele1_grep = allele1.replace('*', '\\\*').replace(':', '\:')
    allele1_alignment = os.path.join(tmp_dir, '{}_alignment.txt'.format(allele1_name))
    cmd = "grep -w {} {} | cut -f 3- -d ' '  | sed 's/^\s*//'".format(allele1_grep, gene_clustalw, allele1_alignment)
    allele1_alignment = subprocess.check_output(cmd, shell=True)

    # EXTRACTING allele2 sequence data 
    allele2_name = allele2.replace('*', '_').replace(':', '_')
    allele2_grep = allele2.replace('*', '\\\*').replace(':', '\:')
    allele2_alignment = os.path.join(tmp_dir, '{}_alignment.txt'.format(allele2_name))
    cmd = "grep -w {} {} | cut -f 3- -d ' '  | sed 's/^\s*//'".format(allele2_grep, gene_clustalw, allele2_alignment)

    allele2_alignment = subprocess.check_output(cmd, shell=True)

    if test_dir:
        fn = os.path.join(test_dir, '{}_alignment.txt'.format(allele1_name))
        with open(fn, 'w') as f: 
            f.write(allele1_alignment)

        fn = os.path.join(test_dir, '{}_alignment.txt'.format(allele2_name))
        with open(fn, 'w') as f: 
            f.write(allele2_alignment)

    linearize_seq = lambda x: x.replace('\n', '').replace(' ', '').replace('\t', '')
    seq1 = linearize_seq(allele1_alignment)
    seq2 = linearize_seq(allele2_alignment)

    # COMPARING the two sequences for any "unknown" UTR sequences between them.  
    five_prime_diff, five_prime_variant_len, five_prime_skip = determine_unknown_indel_length(seq1, seq2)
    three_prime_diff, three_prime_variant_len, three_prime_skip = determine_unknown_indel_length(list(reversed(seq1)), list(reversed(seq2)))

    utr_sequence_report = report_unknown_indel_results(allele1, allele2, \
                 five_prime_diff, five_prime_variant_len, three_prime_diff, three_prime_variant_len)

    # GENERATING a final report 
    if five_prime_diff == 'error - cannot compare' or three_prime_diff == 'error - cannot compare':
    
        if report:
            final_report = ' Error - cannot accurately determine variants.'
        else:
            header = ['snps', 'ins', 'dels', 'description']
            data = ['na', 'na', 'na', 'Error - cannot accurately determine variants.']
            final_report = '\t'.join(header) + '\n' + '\t'.join([str(x) for x in data])

    else:
        description, five_prime_rep, three_prime_rep = utr_sequence_report.split(';')
        # COUNTING the number of snps, ins, and dels relative to the two sequences.  
        snps, ins, dels = determine_num_snps_ins_dels(seq1, seq2, five_prime_variant_len, five_prime_skip, three_prime_variant_len, three_prime_skip)

        if report:

            variant_report = ' Body {} snps; {} ins; {} dels'.format(snps, ins, dels)
            final_report = description + '\n\t-'
            final_report += '\n\t-'.join([five_prime_rep, variant_report, three_prime_rep])

        else:
            header = ['snps', 'ins', 'dels', 'description']
            data = [snps, ins, dels, description]
            final_report = '\t'.join(header) + '\n' + '\t'.join([str(x) for x in data])

    # GENERATE output 
    if not outdir: 
        print final_report 

    else: 
        if report:
            fn = os.path.join(outdir, '{}_v_{}.txt'.format(allele1_name, allele2_name))
        else:
            fn = os.path.join(outdir, '{}_v_{}.tsv'.format(allele1_name, allele2_name))
        with open(fn, 'w') as f:
            f.write(final_report)

        








