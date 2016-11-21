
# coding: utf-8

# # Command Line Arguments

# In[ ]:

import argparse

parser=argparse.ArgumentParser()
# parser._optionals.title = "Flag Arguments"
parser.add_argument('-dataID',help="Sample Data ID", required=True)
parser.add_argument('-bam',help="Sorted rm_dup bam file", required=True)
parser.add_argument('-peak',help="Narrow Peak file", required=True)
parser.add_argument('-reg_cov',help="File path to create read coverage pileup.", required=True)
parser.add_argument('-var_cov',help="File path to create variant coverage pileup.", required=False)
parser.add_argument('-output',help="Output tsv with qc metric file path.", required=True)
parser.add_argument('-logPath',help="Output log file path",required=True)
parser.add_argument('--cellType',help="If given, database will not be queried and supplied cell type will be used.  Requires NumInputReads.  Optional Argument", required=False,default="NotGiven")
parser.add_argument('--NumInputReads',help="If given, database will not be queried and supplied input read count will be used.  Requires cellType.  Optional Argument", required=False,default=-1,type=int)
args = vars(parser.parse_args())


# # Imports

# In[1]:

import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
import subprocess
import logging


# # Set Up Logging

# In[ ]:

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename=args['logPath'],
                    filemode='w+')

# define a Handler which writes INFO messages or higher to the sys.stderr so that console will display output too
console = logging.StreamHandler()
console.setLevel(logging.INFO)
logging.getLogger('').addHandler(console)


# # Functions

# In[1]:

def readPeaks(peakFile):
    Peaks=[]
    first=True
    with open(peakFile,"r") as f:
        for line in f:
            if first:
                first=False
            else:
                line=line.strip().split()
                Peaks.append((line[0],int(line[1]),int(line[2])-int(line[1])))
    return Peaks

def readRoadmapStates(cellType="CM"):
    if cellType=="CM":
        filePath="/publicdata/roadmap_15_state_20151104/E083_15_coreMarks_mnemonics_sorted.bed"
    elif cellType=="islet":
        filePath="/publicdata/roadmap_15_state_20151104/E087_15_coreMarks_mnemonics_sorted.bed"
    elif cellType=="iPSC":
        filePath="/publicdata/roadmap_15_state_20151104/E020_15_coreMarks_mnemonics_sorted.bed"
    chromosomeToRegions={}
    activeStates=["7_Enh","12_EnhBiv","6_EnhG","1_TssA","10_TssBiv","2_TssAFlnk"]
    with open(filePath,"r") as f:
        for line in f:
            line=line.strip().split()
            if line[3] in activeStates:
                if line[0] in chromosomeToRegions:
                    chromosomeToRegions[line[0]].append((int(line[1]),int(line[2])-int(line[1])))
                else:
                    chromosomeToRegions[line[0]]=[(int(line[1]),int(line[2])-int(line[1]))]
    return chromosomeToRegions

def countPeakStates(peaks,chromosomeToRegions):
    #tuple of [fully within promoter or enhancer, partially within promoter or enhancer, not within promoter or enhancer,wrong chr]
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
            
#             possibleAreas=[region for region in chromosomeToRegions[peak[0]] if region[0] <= peak[1] and region[0]+region[1] >=peak[1]]
            #possibleAreas now contains all peaks that start before or on the site and extend past or to the site
            
#             Tracer()()
            
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
    with open(filePath,"w+") as f:
        for chrom in chromosomes:
            for startLen in roadmapStates[chrom]:
                f.write(chrom+"\t"+str(startLen[0])+"\t"+str(startLen[0]+startLen[1])+"\n")
                
    return

def QC_function(dataID,region_coverage_path,var_coverage_path,cellType,NumInputReads,peakFile,bamFile):
    #6 steps:
    #1.  Query the database for cell type and number of input reads (done in previous function call or passed as cmd line arg)
    #2.  Read each roadmap type
    #3.  Read the peak file
    #4.  Find the state counts for the peaks
    #5.  Find the coverage of reads within the roadmap sites
    #6.  Find the median depth of common variants within the reference sites
       
    ID=dataID
    
    results=[]
    
    #Start sample by sample analysis

    #3.  Read the peaks
    logging.info("Reading Peak File")
    peaks=readPeaks(peakFile)

    #4.  Find the state counts for the peaks
    logging.info("Finding Peak States")
    if cellType=="iPSC":
        cellTypeRoadmap=readRoadmapStates(cellType="iPSC")
    elif cellType=="CM":
        cellTypeRoadmap=readRoadmapStates()
    peakDist=countPeakStates(peaks,cellTypeRoadmap)

    fractionActivePeaks=float(peakDist[0])/sum(peakDist)

    logging.info("Finding Peak Coverage")
    #5.  Call bedtools coverage for active regions
    if cellType=="iPSC":
        subp2=subprocess.Popen(["module load bedtools; bedtools coverage -sorted -b "+bamFile+" -a /projects/CARDIPS/pipeline/ATACseq/reference_files/iPSC_active_elements.bed > "+region_coverage_path],shell=True)
    elif cellType=="CM":
        subp2=subprocess.Popen(["module load bedtools; bedtools coverage -sorted -b "+bamFile+" -a /projects/CARDIPS/pipeline/ATACseq/reference_files/CM_active_elements.bed > "+region_coverage_path],shell=True)
    subp2.wait()
    activeReads=sum([int(line.strip().split()[3]) for line in open(region_coverage_path,"r")])
    fractionActiveReads=float(activeReads)/NumInputReads

    logging.info("Finding Variant Coverage")
    #6. Call bedtools coverage for common variants
    if cellType=="iPSC":
        subp3=subprocess.Popen(["module load samtools; samtools mpileup -A --positions /projects/CARDIPS/pipeline/ATACseq/reference_files/iPSC_ATAC_Common_Vars.positions -o "+var_coverage_path+" "+bamFile],shell=True)
    elif cellType=="CM":
        subp3=subprocess.Popen(["module load samtools; samtools mpileup -A --positions /projects/CARDIPS/pipeline/ATACseq/reference_files/CM_ATAC_Common_Vars.positions -o "+var_coverage_path+" "+bamFile],shell=True)
    subp3.wait()
    meanVarCov=np.mean([int(line.strip().split()[3]) for line in open(var_coverage_path,"r")])

    res=[ID,fractionActivePeaks,fractionActiveReads,meanVarCov,"No Error During QC"]
    results.append(res)

    return pd.DataFrame(results,columns=["Data ID","Fraction Peaks in Active Elements","Fraction Reads in Active Elements","Mean Coverage of Common Vars","Error Code"])

def findCellTypeAndReadCount(dataID):
    #1.  Query the database for cell type and number of input reads
    subp1=subprocess.Popen(["python","queryDatabase.py","-data","data_id,cell,Number_of_input_reads","-dataType","atacs","--pickle","atac_QC_database_query.pkl"])
    subp1.wait()
    query=pd.read_pickle("atac_QC_database_query.pkl")
    query.index=[i for i in range(len(query))]
    query["data_id"]=[ID[:8]+"-"+ID[8:12]+"-"+ID[12:16]+"-"+ID[16:20]+"-"+ID[20:] for ID in query["data_id"]]

    dataIDtoCell={query.loc[ind,"data_id"]:query.loc[ind,"cell"] for ind in query.index}
    dataIDtoReadNum={query.loc[ind,"data_id"]:query.loc[ind,"Number_of_input_reads"] for ind in query.index}
    
    ID=dataID
    
    if ID not in dataIDtoReadNum:
        return "readCountErr","readCountErr"
    
    
    if dataIDtoCell[ID]!="CM" and dataIDtoCell[ID]!="iPSC":
        return "cellTypeErr","cellTypeErr"
    
    return dataIDtoCell[ID],dataIDtoReadNum[ID]


# In[ ]:

logging.info("""
=====================
ATACseq Peak QC
Author:  Bill Greenwald
Version: 1.3
=====================
""")


# In[2]:
results = []
if args['cellType']=="Not Given" or args['NumInputReads']==-1:
    if args['cellType']=="Not Given" and args['NumInputReads']==-1:
        cellType,NumInputReads=findCellTypeAndReadCount(args['input'])
        if NumInputReads=="readCountErr":
            res=[ID,-1,-1,-1,"Data ID not linked to sample"]
            results.append(res)
            temp=pd.DataFrame(results,columns=["Data ID","Fraction Peaks in Active Elements","Fraction Reads in Active Elements","Mean Coverage of Common Vars","Error Code"])
            temp.to_csv(args['output'],sep="\t",index=False)
            exit()
        elif cellType=="cellTypeErr":
            res=[ID,-1,-1,-1, "Cell type is not iPSC or CM"]
            results.append(res)
            temp=pd.DataFrame(results,columns=["Data ID","Fraction Peaks in Active Elements","Fraction Reads in Active Elements","Mean Coverage of Common Vars","Error Code"])
            temp.to_csv(args['output'],sep="\t",index=False)
            exit()
    else:
        print("--cellType and --NumInputReads must be used in tandem.")
        exit()
else:
    cellType=args['cellType']
    if cellType!="iPSC" and cellType!="CM":
        res = [args['dataID'], '', '','','Celltype {} not supported'.format(cellType)]
        results.append(res)
        temp=pd.DataFrame(results,columns=["Data ID","Fraction Peaks in Active Elements","Fraction Reads in Active Elements","Mean Coverage of Common Vars","Error Code"])
        temp.to_csv(args['output'],sep="\t",index=False)
        exit()
    NumInputReads=args['NumInputReads']
QC_function(args['dataID'],args['reg_cov'],args['var_cov'],cellType,NumInputReads,args['peak'],args['bam']).to_csv(args['output'],sep="\t",index=False)


