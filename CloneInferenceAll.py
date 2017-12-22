##############################################################################################
### PROJECT: Immune Repertoire. Analysis B cells antibodies
###
### CITATION: 
###
### PROCESS: To obtain the number of clones per sample
###           
### Author: Silvia Pineda
### Date: January, 2017
############################################################################################
###import pandas library to work with data frames
import pandas as pd
import numpy as np
from functools import partial
from multiprocessing import Pool
################################
######### Functions ############
################################

def MatchNucleotides( a,b ):
    count = len(a)
    threshold = len(a)*0.95
    for i in range(0, len(a)):
        if(a[i]!=b[i]):
            count = count-1
            if count < threshold:
                return False
    return True 

def NucleotideBelongsClone(nucleotide,clone):
    for i in range(0,len(clone)):
        matches=MatchNucleotides(nucleotide,clone[i])
        if(matches): 
            return True
    return False      

def AddNucleotide(nucleotide,clones):
    clones_new=[]
    nucleotide_clone = [nucleotide]
    for i in range(0,len(clones)):
        if(NucleotideBelongsClone(nucleotide,clones[i])):
            nucleotide_clone = nucleotide_clone + clones[i]
        else:
            clones_new.append(clones[i])
    clones_new.append(nucleotide_clone)
    return clones_new

def ProcessNucleotides(nucleotides):
    Listclones=[]
    for i in range(0,len(nucleotides)):
        Listclones = AddNucleotide(nucleotides[i], Listclones)
    return Listclones

def ObtainNumberClones(ListClones):
    count = 1
    df_clones = pd.DataFrame([])
    for i in range(0,len(ListClones)):
        for j in range(0,len(ListClones[i])):
            df_clones = df_clones.append({'clone':ListClones[i][j],
                                    'number':count},ignore_index=True)
        count=count+1
    return df_clones


def ProcessGroup(data_clonesInference, V_J_CDR3_unique, j):
    print (j)
    data_clonesInference_V_J_CDR3_unique = data_clonesInference[data_clonesInference['V_J_lenghCDR3'] == V_J_CDR3_unique[j]]
    nucleotides = list(data_clonesInference_V_J_CDR3_unique['cdr3_seq'])
    ##Obtain the number of clones inferred
    ClonesInfered = ObtainNumberClones(ProcessNucleotides(nucleotides))
    data_clonesInference_V_J_CDR3_unique.loc[data_clonesInference_V_J_CDR3_unique.index,'CloneId'] = pd.Series(ClonesInfered['number']).values
    data_clonesInference_V_J_CDR3_unique['CloneId'] = data_clonesInference_V_J_CDR3_unique['V_J_lenghCDR3'] + "_" + data_clonesInference_V_J_CDR3_unique['CloneId'].map(str)
    return data_clonesInference_V_J_CDR3_unique

def ProcessSample(data_clonesInference):
    result = pd.DataFrame([])
    V_J_CDR3_unique = data_clonesInference['V_J_lenghCDR3'].unique()
    ClonesInfered = 0
    ##Loop for each V_J_CDR3 unique 
    pool = Pool(30)
    partialFunction = partial(ProcessGroup, data_clonesInference, V_J_CDR3_unique)
    result = pool.map(partialFunction, range(0,len(V_J_CDR3_unique)))
    return result
    
########################
###### Main program ####
########################
print("Start")
data_clonesInference = pd.read_csv("/home/pinedasans/VDJ/clonsInfered/data_clonesInference.txt",sep="\t")

###Obtain the unique subjects
#sample_unique = data_clonesInference['sample_id'].unique()

##Declare result: total number of clones infered per subject
result_ClonesInfered = pd.DataFrame([])

###Loop for each subject
final_result = ProcessSample(data_clonesInference)
  
result_ClonesInfered = result_ClonesInfered.append(final_result)  

###Result    
result_ClonesInfered.to_csv('/home/pinedasans/VDJ/clonsInfered/ClonesInferedAll.csv')
print("End")