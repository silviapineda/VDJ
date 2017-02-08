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
    clones=[]
    for i in range(0,len(nucleotides)):
        clones = AddNucleotide(nucleotides[i], clones)
    return clones



########################
###### Main program ####
########################

###import pandas library to work with data frames
import pandas
data_clonesInference = pandas.read_csv("/Users/Pinedasans/Data/VDJ/data_clonesInference.txt",sep="\t")

###Obtain the unique samples
sample_unique = data_clonesInference['specimen_label'].unique()

##Declare result: total number of clones infered per sample
ClonesInfered_sample = []
ReadsPerClone = []

###Loop for each sample
for i in range(0,len(sample_unique)):
    print (i)
    data_clonesInference_sample_unique = data_clonesInference[data_clonesInference['specimen_label'] == sample_unique[i]]
    V_J_CDR3_unique = data_clonesInference_sample_unique['V_J_lenghCDR3'].unique()
    ClonesInfered = 0
    
    ##Loop for each V_J_CDR3 unique 
    for j in range(0,len(V_J_CDR3_unique)):
        data_clonesInference_V_J_CDR3_unique = data_clonesInference_sample_unique[data_clonesInference_sample_unique['V_J_lenghCDR3'] == V_J_CDR3_unique[j]]
        nucleotides = list(data_clonesInference_V_J_CDR3_unique['cdr3_seq'])
        ##Obtain the number of clones infered per sample
        ClonesInfered = ClonesInfered + len(ProcessNucleotides(nucleotides))
    ClonesInfered_sample.append(ClonesInfered)  
    ReadsPerClone.append(len(data_clonesInference_sample_unique)/ClonesInfered_sample[0])

###Result    
clones_count_sample = pandas.DataFrame(
{'specimen_label':sample_unique,
'num_clones':ClonesInfered_sample,
'reads_per_clone':ReadsPerClone
}) 
clones_count_sample.to_csv('/Users/Pinedasans/Data/VDJ/clones_count_sample.csv')