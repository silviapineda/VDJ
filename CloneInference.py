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

################################
######### Functions ############
################################

def MatchNucleotides( a,b ):
    count = len(a)
    threshold = len(a)*0.9
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


    
    
    
def MinReadDepth(data_clonesInference,amplification,min_to_down,sample_to_down):
    if (amplification is not None):
        data_clonesInference_amp = data_clonesInference[data_clonesInference['amplification_template'] == amplification]
        reads = pd.Series(data_clonesInference_amp[sample_to_down]).value_counts()
    else:
        reads = pd.Series(data_clonesInference[sample_to_down]).value_counts()
    
    min_read_depth = min(reads[reads>min_to_down])
    return min_read_depth

def DownSampling(data,min_read_depth):
    data_sample_down = data.sample(min_read_depth) 
    return data_sample_down

    
def DataToDownSampled(data_clonesInference,sample_by):
    data_clonesInference_amp = data_clonesInference[data_clonesInference['amplification_template'] == amplification]
    reads = pd.Series(data_clonesInference_amp[sample_to_down]).value_counts()
        return data_clonesInference[data_clonesInference['amplification_template'] == amplification]
    else:
        return data_clonesInference

        
        
     
####To obtain the clones by downsampling by individuals

min_read_depth = MinReadDepth(data_clonesInference,None,100,'sample_id')
data_clones_to_downSampled = DataToDownSampled(data_clonesInference,None)

min_read_depth = MinReadDepth(data_clonesInference,"cDNA",100,'specimen_label')
data_clones_to_downSampled_cDNA = DataToDownSampled(data_clonesInference,"cDNA")

min_read_depth = MinReadDepth(data_clonesInference,"gDNA",1000,'specimen_label')
data_clones_to_downSampled_gDNA = DataToDownSampled(data_clonesInference,"gDNA")

data_clonesInference_sample_unique = data_clones_to_downSampled[data_clones_to_downSampled['sample_id'] == sample_unique[i]]
data_clonesInference_sample_unique_down = DownSampling(data_clonesInference_sample_unique,min_read_depth)
    
        
def ProcessSample(data_clonesInference_sample_unique):
     ##To calculate the minimum read depth per individual to calculate the clones
   
    V_J_CDR3_unique = data_clonesInference_sample_unique_down['V_J_lenghCDR3'].unique()
    ClonesInfered = 0
       
    ##Loop for each V_J_CDR3 unique 
    result = data_clonesInference_sample_unique_down
    for j in range(0,len(V_J_CDR3_unique)):
        #print (j)
        nucleotides = list(data_clonesInference_V_J_CDR3_unique['cdr3_seq'])
        ##Obtain the number of clones infered per sample
        ClonesInfered = ObtainNumberClones(ProcessNucleotides(nucleotides))
        result.loc[data_clonesInference_V_J_CDR3_unique.index,'numberClone'] = pd.Series(ClonesInfered['number']).values
    return result
  

    
########################
###### Main program ####
########################

data_clonesInference = pd.read_csv("/Users/Pinedasans/Data/VDJ/data_clonesInference.txt",sep="\t")

###Obtain the unique subjects
sample_unique = data_clonesInference['sample_id'].unique()

##Declare result: total number of clones infered per subject
result_ClonesInfered = pd.DataFrame([])
  
###Loop for each subject
for i in range(0,len(sample_unique)):
    print (i)
    data_clonesInference_sample_unique = data_clones_to_downSampled[data_clones_to_downSampled['sample_id'] == sample_unique[i]]
    data_clonesInference_sample_unique_down = DownSampling(data_clonesInference_sample_unique,min_read_depth)
    final_result = ProcessSample(data_clonesInference_sample_unique)
    result_ClonesInfered = result_ClonesInfered.append(final_result)  

###Result    
result_ClonesInfered.to_csv('/Users/Pinedasans/Data/VDJ/ClonesInfered_downsampled.csv')

