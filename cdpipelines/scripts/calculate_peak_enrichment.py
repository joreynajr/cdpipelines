"""
Calculate peak information which is used for quality control of the
ATAC-seq pipeline. This code is celltype specific. Cell types 
include induced pluripotent stem cells (iPSC), cardiomyocytes (CM), 
pancreas progenitors (PP3), and retina cells (RPE). More celltypes 
may be added over time. 

"""

# Imports
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
import subprocess
import logging
import warnings
functional_celltypes = ['iPSC', 'CM', 'PP3', 'RPE']

def findCellTypeAndReadCount(dataID):
	#1.  Query the database for cell type and number of input reads
	subp1=subprocess.Popen(['python','queryDatabase.py','-data','data_id,cell,Number_of_input_reads','-dataType','atacs','--pickle','atac_QC_database_query.pkl'])
	subp1.wait()
	query=pd.read_pickle('atac_QC_database_query.pkl')
	query.index=[i for i in range(len(query))]
	query['data_id']=[ID[:8]+'-'+ID[8:12]+'-'+ID[12:16]+'-'+ID[16:20]+'-'+ID[20:] for ID in query['data_id']]

	dataIDtoCell={query.loc[ind,'data_id']:query.loc[ind,'cell'] for ind in query.index}
	dataIDtoReadNum={query.loc[ind,'data_id']:query.loc[ind,'Number_of_input_reads'] for ind in query.index}
	
	ID=dataID
	
	if ID not in dataIDtoReadNum:
		return 'readCountErr','readCountErr'
	
	
	if dataIDtoCell[ID]!='CM' and dataIDtoCell[ID]!='iPSC' and dataIDtoCell[ID]!='PP3':
		return 'cellTypeErr','cellTypeErr'
	
	return dataIDtoCell[ID],dataIDtoReadNum[ID]

def readPeaks(peakFile):
	Peaks=[]
	first=True
	with open(peakFile,'r') as f:
		for line in f:
			if first:
				first=False
			else:
				line=line.strip().split()
				Peaks.append((line[0],int(line[1]),int(line[2])-int(line[1])))
	return Peaks

def readRoadmapStates(cellType='iPSC'):
	if cellType in ['CM']:
		filePath='/publicdata/roadmap_15_state_20151104/E083_15_coreMarks_mnemonics_sorted.bed'
	elif cellType in ['PP3']:
		filePath='/publicdata/roadmap_15_state_20151104/E087_15_coreMarks_mnemonics_sorted.bed'
	elif cellType in ['iPSC']:
		filePath='/publicdata/roadmap_15_state_20151104/E020_15_coreMarks_mnemonics_sorted.bed'
	else: # Default to iPSC if code for cell type has not been implemented.
		filePath='/publicdata/roadmap_15_state_20151104/E020_15_coreMarks_mnemonics_sorted.bed'
	chromosomeToRegions={}
	activeStates=['7_Enh','12_EnhBiv','6_EnhG','1_TssA','10_TssBiv','2_TssAFlnk']
	with open(filePath,'r') as f:
		for line in f:
			line=line.strip().split()
			if line[3] in activeStates:
				if line[0] in chromosomeToRegions:
					chromosomeToRegions[line[0]].append((int(line[1]),int(line[2])-int(line[1])))
				else:
					chromosomeToRegions[line[0]]=[(int(line[1]),int(line[2])-int(line[1]))]
	return chromosomeToRegions

def countPeakStates(peaks,chromosomeToRegions):
	"""
	Determines the number of peaks which are fully within promoter or enhancer, partially within promoter or enhancer,
	not within promoter or enhancer, wrong chr. 


	Parameters
	----------
	peaks: peak data (generated using the readPeaks function).
	chromosomeToRegions: celltype specific roadmap data (generated using readRoadmapState the countPeakStates function).

	Returns
	-------
	PeakDistribution: tuple
		counts the number of peaks which are fully within promoter or enhancer, partially within promoter or enhancer, 
		not within promoter or enhancer, or on an incorrect chromosome.
		
	"""
	#tuple of [
	peakDistribution=[0,0,0,0]
	for i in range(len(peaks)):
		peak=peaks[i]
		if i==0 or peaks[i][0]!=peaks[i-1][0]:
			roadmapIndex=0
		if peak[0] not in chromosomeToRegions:
			peakDistribution[3]+=1
		else:
			#to make this faster, keep track of the starting index of each search since the peaks and the regions are both ordered.
			#start each next search at the stop (including) of the previous.  Only makes you need to go through each list one time
			possibleAreas=[]
			while roadmapIndex<len(chromosomeToRegions[peak[0]]):
				region=chromosomeToRegions[peak[0]][roadmapIndex]
					 
				if region[0] <= peak[1] and region[0]+region[1] >=peak[1]:
					possibleAreas.append(region)
					
				if peak[1] >= region[0]+region[1]:
					roadmapIndex+=1
				else:
					break

			if len(possibleAreas)==0:
				peakDistribution[2]+=1
			else:
				possibleAreas=[region for region in possibleAreas if region[0]+region[1] >= peak[1]+peak[2]]
				#possibleAreas now only contains the region that houses the element
				if len(possibleAreas)==0:
					peakDistribution[1]+=1
				else:
					peakDistribution[0]+=1
	return peakDistribution


def exportActiveStates(roadmapStates,filePath):
	chromosomes=sorted(roadmapStates.keys())
	with open(filePath,'w+') as f:
		for chrom in chromosomes:
			for startLen in roadmapStates[chrom]:
				f.write(chrom+'\t'+str(startLen[0])+'\t'+str(startLen[0]+startLen[1])+'\n')
				
	return


def QC_function(dataID, region_coverage_path, var_coverage_path, cellType, NumInputReads, peakFile, bamFile):
	"""
	QC-ing the ATAC-seq results. 
	Determined the fraction of (full) peaks in active elements, the fraction of reads in peaks, and the 
	median depth of common variants in cell type/ATAC-seq specific positions. 
	#1.  Read each roadmap type
	#2.  Read the peak file
	#3.  Find the state counts for the peaks
	#4.  Find the coverage of reads within the roadmap sites
	#5.  Find the median depth of common variants within the reference sites


	Parameters
	----------

	Returns
	-------
	results: pd.DataFrame
		A dataframe of the gathered data. 

	"""
	   
	results_cols = ['Data ID', 'Fraction Peaks in Active Elements', 'Fraction Reads in Active Elements', \
					'Mean Coverage of Common Vars', 'Error Code', 'cellType']

	#1.  Find the state counts for the peaks
	logging.info('Finding Peak States')
	if cellType in ['iPSC', 'CM']:
		cellTypeRoadmap = readRoadmapStates(cellType=cellType)
		active_elements_bed = '/projects/CARDIPS/pipeline/ATACseq/reference_files/{}_active_elements.bed'.format(cellType)
		common_variants_position = '/projects/CARDIPS/pipeline/ATACseq/reference_files/{}_ATAC_Common_Vars.positions'.format(cellType)

		#2.  Read the peaks
		logging.info('Reading Peak File')
		peaks = readPeaks(peakFile)

		#3.  Find the state counts for the peaks
		peakDist = countPeakStates(peaks, cellTypeRoadmap)
		fractionActivePeaks = float(peakDist[0])/sum(peakDist)

		#4.  Find the coverage of reads within the roadmap sites
		logging.info('Finding Peak Coverage')
		subp2 = subprocess.Popen('module load bedtools; bedtools coverage -sorted -b {} -a {} > {}'.format(bamFile, active_elements_bed, region_coverage_path), shell=True)
		subp2.wait()
		activeReads = sum([int(line.strip().split()[3]) for line in open(region_coverage_path,'r')])
		fractionActiveReads = float(activeReads)/NumInputReads

		#5.  Find the median depth of common variants within the reference sites
		logging.info('Finding Variant Coverage')
		subp3 = subprocess.Popen('module load samtools; samtools mpileup -A --positions {} -o {} {}'.format(common_variants_position, var_coverage_path, bamFile), shell=True)
		subp3.wait()
		meanVarCov = np.mean([int(line.strip().split()[3]) for line in open(var_coverage_path,'r')])

		results = [dataID, fractionActivePeaks, fractionActiveReads, meanVarCov, 'No Error During QC', cellType]
		results = pd.DataFrame([results], columns=results_cols)

	else:

		results = [dataID, np.nan, np.nan, np.nan, 'Celltype {} not supported'.format(cellType), cellType]
		results = pd.DataFrame([results], columns=results_cols)
		warnings.warn('Error: Celltype {} is not supported. Empty dataframe was output.'.format(cellType))

	return results  


if __name__ == '__main__':

	# Command Line Arguments
	parser=argparse.ArgumentParser()
	parser.add_argument('-dataID',help='Sample Data ID', required=True)
	parser.add_argument('-bam',help='Sorted rm_dup bam file', required=True)
	parser.add_argument('-peak',help='Narrow Peak file', required=True)
	parser.add_argument('-reg_cov',help='File path to create read coverage pileup.', required=True)
	parser.add_argument('-var_cov',help='File path to create variant coverage pileup.', required=False)
	parser.add_argument('-output',help='Output tsv with qc metric file path.', required=True)
	parser.add_argument('-logPath',help='Output log file path',required=True)
	parser.add_argument('--cellType',help='If given, database will not be queried and supplied cell type will be used.  Requires NumInputReads.  Optional Argument', required=False,default='NotGiven')
	parser.add_argument('--NumInputReads',help='If given, database will not be queried and supplied input read count will be used.  Requires cellType.  Optional Argument', required=False,default=-1,type=int)
	args = vars(parser.parse_args())

	# Set Up Logging
	logging.basicConfig(level=logging.DEBUG,
						format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
						datefmt='%m-%d %H:%M',
						filename=args['logPath'],
						filemode='w+')
	# Handler which writes INFO messages or higher to the sys.stderr so that console will display output too
	console = logging.StreamHandler()
	console.setLevel(logging.INFO)
	logging.getLogger('').addHandler(console)

	logging.info('''
	=====================
	ATACseq Peak QC
	Author:  Bill Greenwald
	Version: 1.3
	=====================
	''')

	ID = args['dataID']
	# Parsing the inputs and extracting data for the qc functions 
	results_cols = ['Data ID','Fraction Peaks in Active Elements','Fraction Reads in Active Elements','Mean Coverage of Common Vars','Error Code']

	if args['cellType']=='Not Given' or args['NumInputReads']==-1:
		if args['cellType']=='Not Given' and args['NumInputReads']==-1:
			cellType,NumInputReads=findCellTypeAndReadCount(args['input'])

			if NumInputReads=='readCountErr':
				results=[ID, np.nan, np.nan, np.nan,'Error: Number of input reads was not passed.']
				results=pd.DataFrame([results], columns=results_cols)
				results.to_csv(args['output'], sep='\t', index=False)
				warnings.warn('--NumInputReads must be passed. Empty dataframe was output.')
				exit()

			elif cellType=='cellTypeErr':
				results=[ID, np.nan, np.nan, np.nan,'Error: Celltype information was not passed.']
				results=pd.DataFrame([results], columns=results_cols)
				results.to_csv(args['output'], sep='\t', index=False)
				warnings.warn('--cellType must be used. Empty dataframe was output.')
				exit()
		else:
			results=[ID, np.nan, np.nan, np.nan,'Error: Read count and celltype information was not passed.']
			results=pd.DataFrame([results], columns=results_cols)
			results.to_csv(args['output'], sep='\t', index=False)
			warnings.warn('--cellType and --NumInputReads must be used in tandem. Empty dataframe was output.')
			exit()
	else:
		cellType = args['cellType']
		NumInputReads = args['NumInputReads']
	# Running the qc function 
	results = QC_function(ID, args['reg_cov'], args['var_cov'], cellType, NumInputReads, args['peak'], args['bam'])
	results.to_csv(args['output'], sep='\t', index=False)


