
# coding: utf-8

# # Command Line Arguments

# In[132]:

import argparse


# In[ ]:

parser=argparse.ArgumentParser()
parser._optionals.title = "Flag Arguments"
parser.add_argument('-data',help="Comma separated list of data information desired.  Required Argument", required=True)
parser.add_argument('-dataType',help="Type of desired data.  Required Argument", required=True)
parser.add_argument('--filters',help="comma separated list of projects to search in.  Defaults to all.  Optional Argument.", required=False,default="All")
parser.add_argument('--pickle',help="Write output to specified pickle file instead of a tsv to std out.  Optional Argument", required=False,default="%#$")
args = vars(parser.parse_args())


# # Imports

# In[ ]:

import pandas as pd
import MySQLdb
from IPython.core.debugger import Tracer


# # Define Database Architecture
# 
# ### We will define the following tree structure:
# 
# The subject_subject table is the root
# 
# Each child of the root, no matter what table it actually is, is the project table.  These can have separate children metadata tables themselves, but they also are the top level parent node

# In[37]:

class Node(object):
    def __init__(self, name, tableType, table_members):
        self.name=name
        self.tableType = tableType
        self.tableMembers=table_members
        self.children = []
        if tableType=="Root":
            self.parents=None
        else:
            self.parents=[]
        
    def add_child(self, obj):
        self.children.append(obj)
        obj.parents.append(self)
        
    def describe(self):
        print("Table Name:"+str(self.name))
        print("Table Type:"+str(self.tableType))
        print("Table Members:"+str(self.tableMembers))
        print("Children:"+str([x.name for x in self.children]))
        if self.parents==None:
            print("Parents: Root of database")
        else:
            print("Parents:"+str([x.name for x in self.parents]))


# In[261]:

#"\""+"\",\"".join(list(queryDatabase("describe subject_subject")["Field"]))+"\""


# In[262]:

subjSubj=Node("Subject Subject","Root",["family_id","subj.name","ext_name","subj.id","sex","disease","relation","father_id","mother_id","twin_id","self_reported_race","ethnicity_group","expected_super_pop","ancestry_comment","age","height","weight","diag_physician","diag_self","comment","cardips_consent","collection_mechanism","consent_status","hipaa_consent"])

#define next layer of tree

fam1070_tissue=Node("Family1070_tissue","Project",["subject_id","tis.id","clone","cell","passage","day","differentiation_success","differentiation_protocol","pct_ctnt_positive_cells","pct_ctnt_highly_positive_cells","pct_ctni_positive_cells","harvest_day","operator1","operator2","harvest_id","culture_vessels","tissue_unique_identifier"])
timecourse_tissue=Node("Timecourse_tissue","Project",["subject_id","tis.id","name","cell","clone","passage","day","diff_success","diff_protocol","harvest_day","operator1","operator2","culture_vessel","harvest_id"])

subjSubj.add_child(fam1070_tissue)
subjSubj.add_child(timecourse_tissue)

#define next layer of tree:

fam1070_rnas=Node("Family1070_rnas","Sequencing",["id","tissue_id","fastq_or_file_name","igm_manifest_sample_name","igm_manifest_sample_code","library_sequencing_manifest_date","data_type","assay_name","assay_done_by","test_number","assay_replicate","library_replicate","sequencing_replicate","standardized_seq_id","comment","rna_extraction_date","rin","rna_concentration"])
fam1070_chips=Node("Family1070_chips","Sequencing",["id","tissue_id","fastq_or_file_name","igm_manifest_sample_name","igm_manifest_sample_code","library_sequencing_manifest_date","data_type","assay_name","assay_done_by","test_number","assay_replicate","library_replicate","sequencing_replicate","standardized_seq_id","comment","antibody_cat_no","antibody_lot","antibody_micrograms","sonication_cycles","date_of_ip_assay","micrograms_chromatin_input","nanograms_dna_ip_yield"])
fam1070_atacs=Node("Family1070_atacs","Sequencing",["id","tissue_id","fastq_or_file_name","igm_manifest_sample_name","igm_manifest_sample_code","library_sequencing_manifest_date","data_type","assay_name","assay_done_by","test_number","assay_replicate","library_replicate","sequencing_replicate","standardized_seq_id","comment","tagmentation_date","avg_mol_weight"])
fam1070_hic=Node("Family1070_hic","Sequencing",["id","sequencing_run_flowcell","original_sequence_name","manifest_sample_name","sample_code","manifest_date","assay_date","assay_done_by","test_number","assay_replicate","library_replicate","sequencing_replicate","status","comment","tissue_id"])
# fam1070_gros=Node("Family1070_gros","Sequencing",[])
fam1070_tissue.add_child(fam1070_rnas)
fam1070_tissue.add_child(fam1070_chips)
fam1070_tissue.add_child(fam1070_atacs)
fam1070_tissue.add_child(fam1070_hic)
# fam1070_tissue.add_child(fam1070_gros)

timecourse_rnas=Node("Timecourse_rnas","Sequencing",["seq.id","tissue_id","igm_manifest_sample_name","igm_manifest_sample_code","lib_seq_manifest_date","data_type","assay_name","assay_done_by","assay_replicate","library_replicate","sequencing_replicate","standardized_seq_id","comment","rna_extraction_date","rin","rna_concentration"])
timecourse_chips=Node("Timecourse_chips","Sequencing",["seq.id","tissue_id","igm_manifest_sample_name","igm_manifest_sample_code","lib_seq_manifest_date","data_type","assay_name","assay_done_by","assay_replicate","library_replicate","sequencing_replicate","standardized_seq_id","comment","antibody","sonication_cycles","date_of_ip_assay","micrograms_chromatin_input","nanograms_dna_ip_yield"])
timecourse_atacs=Node("Timecourse_atacs","Sequencing",["seq.id","tissue_id","igm_manifest_sample_name","igm_manifest_sample_code","lib_seq_manifest_date","data_type","assay_name","assay_done_by","assay_replicate","library_replicate","sequencing_replicate","standardized_seq_id","comment","date_of_atac_assay","pcr_cycles","barcode1_name","barcode1_sequence","barcode2_name","barcode2_sequence","date_of_size_selection","average_insert_size","dna_concentration"])
timecourse_tissue.add_child(timecourse_rnas)
timecourse_tissue.add_child(timecourse_chips)
timecourse_tissue.add_child(timecourse_atacs)

#define next layer of tree

data_atacs=Node("Data Atacs","Data",["sequence_id","data.id","data.name","status","comment","sample_id","sample_type_id"])
data_chips=Node("Data Chips","Data",["sequence_id","data.id","data.name","status","comment","sample_id","sample_type_id"])
data_rnas=Node("Data Rnas","Data",["sequence_id","data.id","data.name","status","comment","sample_id","sample_type_id"])
data_hic=Node("Data HIC","Data",["sequence_id","data.id","data.name","status","comment","sample_id","sample_type_id"])

fam1070_tissue.add_child(data_atacs)
fam1070_tissue.add_child(data_chips)
fam1070_tissue.add_child(data_rnas)
fam1070_tissue.add_child(data_hic)

timecourse_tissue.add_child(data_atacs)
timecourse_tissue.add_child(data_chips)
timecourse_tissue.add_child(data_rnas)
timecourse_tissue.add_child(data_hic)

#define next layer of tree

qc_atacs=Node("QC Atacs","QC",["qc.id","total_read","pct_npropealigned","pct_duplicate","pct_mapqlt30","pct_minorchr","pct_blacklist","filtered_read","nrf","pbc","est_flag_len","nsc","rsc","quality_tag","total_peak","peak_with_q001","pct_collapsed_region","pct_multi_peak_region","pct_region_uniq_peak","pct_frip","comment","data_id"])
qc_chips=Node("QC Chips","QC",["qc.id","number_of_input_reads","pct_uniquely_mapped_reads","median_insert_size","pct_duplication","pct_chrx","pct_chry","pct_chrm","pct_promoter","pct_passed_filter","pct_lt140","efficiency","data_id","number_of_peaks"])
qc_rnas=Node("QC Atacs","QC",["qc.id","number_of_input_reads","pct_uniquely_mapped_reads","median_insert_size","median_5prime_to_3prime_bias","pct_duplication","pct_ribosomal_bases","pct_coding_bases","pct_utr_bases","pct_intronic_bases","pct_intergenic_bases","pct_mrna_bases","pct_correct_strand_reads","pct_chrx","pct_chry","identity_conc","data_id"])


data_atacs.add_child(qc_atacs)
data_chips.add_child(qc_chips)
data_rnas.add_child(qc_rnas)


#create dict of data type to children nodes

qcNodes={"atacs":qc_atacs,"chips":qc_chips,"rnas":qc_rnas,"hic":data_hic}


# # Functions

# In[248]:

def queryDatabase(Query):
    mysql_cn= MySQLdb.connect(host='fl-hn1', 
                port=3306,user='cardips', passwd='G00dC@kes', 
                db='cardips')
    return pd.read_sql(Query,con=mysql_cn)

def findMetadata(metaDatum,currentNode):
    if currentNode.parents==None:
        if metaDatum in currentNode.tableMembers:
            return [currentNode]
        else:
            return []
    ParentMetadataResults=[]
    for parent in currentNode.parents:
        ParentMetadataResults.extend(findMetadata(metaDatum,parent))
    if metaDatum in currentNode.tableMembers:
        results=[currentNode]
        results.extend(ParentMetadataResults)
        return results
    else:
        return ParentMetadataResults

def queryAllMetadata(dataNodes,dataType,metaData):
    startNode=dataNodes[dataType]
    potentialNodes=[]
    for metaDatum in metaData:
        res=findMetadata(metaDatum,startNode)
        potentialNodes.extend(res)
        
        
    return list(set(potentialNodes))
        

def getHighestLevelNodes(potentialNodes):
    temp=[x for x in potentialNodes if x.tableType=="Root"]
    if len(temp)>0:
        return temp
    
    temp=[x for x in potentialNodes if x.tableType=="Project"]
    if len(temp)>0:
        return temp
    
    temp=[x for x in potentialNodes if x.tableType=="Sequencing"]
    if len(temp)>0:
        return temp
    
    temp=[x for x in potentialNodes if x.tableType=="Data"]
    if len(temp)>0:
        return temp
    
    temp=[x for x in potentialNodes if x.tableType=="QC"]
    if len(temp)>0:
        return temp

    
def filterNodes(nodes,filters):
    keepNodes=[]
    for node in nodes:
        flag=False
        for filt in filters:
            if filt in node.name:
                flag=True
        if flag==False:
            keepNodes.append(node)
    return keepNodes
    
    
def makeQuerries(nodes,data,dataType):
    results=pd.DataFrame()
    for node in nodes:
        if node.tableType=="QC":
            results=results.append(queryDatabase("SELECT "+data+" FROM qc_"+dataType+""))
        elif node.tableType=="Data":
            results=results.append(queryDatabase("SELECT "+data+" FROM qc_"+dataType+" qc LEFT JOIN data_"+dataType+" data on qc.data_id=data.id"))
        elif node.tableType=="Sequencing":
            if "Family1070" in node.name:
                results=results.append(queryDatabase("SELECT "+data+" FROM qc_"+dataType+" qc LEFT JOIN data_"+dataType+" data on qc.data_id=data.id LEFT JOIN family1070_"+dataType+" fam on data.sample_id=fam.id"))
            elif "Timecourse" in node.name:
                results=results.append(queryDatabase("SELECT "+data+" FROM qc_"+dataType+" qc LEFT JOIN data_"+dataType+" data on qc.data_id=data.id LEFT JOIN timecourse_"+dataType+" time on data.sample_id=time.id"))                
        elif node.tableType=="Project":
            if "Family1070" in node.name:
                results=results.append(queryDatabase("SELECT "+data+" FROM qc_"+dataType+" qc LEFT JOIN data_"+dataType+" data on qc.data_id=data.id LEFT JOIN family1070_"+dataType+" fam on data.sample_id=fam.id LEFT JOIN family1070_tissue tis on fam.tissue_id=tis.id WHERE data.sample_type_id=(SELECT id FROM django_content_type WHERE app_label='family1070' AND model='"+dataType+"')"))
            elif "Timecourse" in node.name:
                results=results.append(queryDatabase("SELECT "+data+" FROM qc_"+dataType+" qc LEFT JOIN data_"+dataType+" data on qc.data_id=data.id LEFT JOIN timecourse_"+dataType+" time on data.sample_id=time.id LEFT JOIN timecourse_tissue tis on time.tissue_id=tis.id WHERE data.sample_type_id=(SELECT id FROM django_content_type WHERE app_label='timecourse' AND model='"+dataType+"')"))
        elif node.tableType=="Root":
                results=results.append(queryDatabase("SELECT "+data+" FROM qc_"+dataType+" qc LEFT JOIN data_"+dataType+" data on qc.data_id=data.id LEFT JOIN family1070_"+dataType+" fam on data.sample_id=fam.id LEFT JOIN family1070_tissue tis on fam.tissue_id=tis.id LEFT JOIN subject_subject subj on tis.subject_id=subj.id WHERE data.sample_type_id=(SELECT id FROM django_content_type WHERE app_label='family1070' AND model='"+dataType+"')"))
                results=results.append(queryDatabase("SELECT "+data+" FROM qc_"+dataType+" qc LEFT JOIN data_"+dataType+" data on qc.data_id=data.id LEFT JOIN timecourse_"+dataType+" time on data.sample_id=time.id LEFT JOIN timecourse_tissue tis on time.tissue_id=tis.id LEFT JOIN subject_subject subj on tis.subject_id=subj.id WHERE data.sample_type_id=(SELECT id FROM django_content_type WHERE app_label='timecourse' AND model='"+dataType+"')"))
        
            
    return results


# # Testing

# In[257]:

# args={"data":"data_id,pct_uniquely_mapped_reads,cell,subj.name","dataType":"rnas","filters":"Family1070","pickle":"%#$"}


# # Run the program

# In[263]:

res=queryAllMetadata(qcNodes,args['dataType'],args['data'].split(","))
res=getHighestLevelNodes(res)
res=filterNodes(res,args['filters'].split(","))
if args['pickle']=="%#$":
    print(makeQuerries(res,args['data'],args['dataType']).to_csv(sep="\t",index=False))
else:
    pick=makeQuerries(res,args['data'],args['dataType'])
    pick.to_pickle(args['pickle'])

