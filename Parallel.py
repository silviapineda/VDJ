def ProcessSample(data_clonesInference_sample_unique):
    min_read_depth = min(pd.Series(data_clonesInference['sample_id']).value_counts())
    ##To calculate the minimum read depth per individual to calculate the clones
    data_clonesInference_sample_unique_down = data_clonesInference_sample_unique.sample(min_read_depth) 
    V_J_CDR3_unique = data_clonesInference_sample_unique_down['V_J_lenghCDR3'].unique()
    ClonesInfered = 0
       
    ##Loop for each V_J_CDR3 unique 
    result = data_clonesInference_sample_unique_down
    for j in range(0,len(V_J_CDR3_unique)):
        #print (j)
        data_clonesInference_V_J_CDR3_unique = data_clonesInference_sample_unique_down[data_clonesInference_sample_unique_down['V_J_lenghCDR3'] == V_J_CDR3_unique[j]]
        nucleotides = list(data_clonesInference_V_J_CDR3_unique['cdr3_seq'])
        ##Obtain the number of clones infered per sample
        ClonesInfered = ObtainNumberClones(ProcessNucleotides(nucleotides))
        result.loc[data_clonesInference_V_J_CDR3_unique.index,'numberClone'] = pd.Series(ClonesInfered['number']).values
    return result
########################
###### Main program ####
########################

data_clonesInference = pandas.read_csv("/Users/Pinedasans/Data/VDJ/data_clonesInference.txt",sep="\t")

###Obtain the unique subjects
sample_unique = data_clonesInference['sample_id'].unique()

##Declare result: total number of clones infered per subject
result_ClonesInfered = []

rangeSamples = range(0,len(sample_unique))
###Loop for each subject
for i in range(0,len(sample_unique)):
    print (i)
    data_clonesInference_sample_unique = data_clonesInference[data_clonesInference['specimen_label'] == sample_unique[i]]
    
    final_result = []
    with Pool(10) as p:
        final_result = p.map(ProcessSample, data_clonesInference_sample_unique)
      
for i in range(0,len(final_result)):
    result_ClonesInfered = result_ClonesInfered.append(final_result[i])  

###Result    
result_ClonesInfered.to_csv('/Users/Pinedasans/Data/VDJ/ClonesInfered_downsampled-test.csv')