import os
import sys
import numpy as np
import math
from copy import deepcopy
import matplotlib.pyplot as plt
import pandas as pd
import pickle


from scipy import stats
from random import sample
from Bio import SeqIO



if __name__ == '__main__':
    
    sys.path.append(os.path.abspath('directory'))
    import funcFile
    import funcPeakAlign
    import funcSeqAll
    import funcToolsAll
    import funcByRef
    import funcGeneral
    import funcTimeWarp
    import BoXFP
    
    np.set_printoptions(threshold=sys.maxsize)
    
    #read in data list
    with open('data_files_TCV_4.fl' ,'r') as f:
        file_list =f.readlines()
    
    #open data file
    file_1= open('190902_4_0.obj','rb')
    
    #load data
    data_arrB=pickle.load(file_1)
    
    #read in genome sequence
    for record in SeqIO.parse('TCV_ref_genome.fasta','fasta'):
        seq_arr = record.seq
        
    #print section of sequence (mainly for use in annotating plots)    
    print seq_arr[23:373]
                
    #specify reference sequence
    ref_seq='TCV_ref_genome.fasta'    

    #partition sequence data
    partition_SB=BoXFP.S1_partitioning(data_arrB,3)

    
    
    #plot partitioned data
    col=iter(plt.cm.rainbow(np.linspace(0,1,50)))
    
    for i in range(len(partition_SB)):
        c=next(col)
        plt.plot(np.arange(1,351),partition_SB[i][:,1],c=c)
    plt.show()

    #calculate correlations between partitioned sequences
    SB_cov_list,SB_cov_matrix=BoXFP.correl_assessor(partition_SB,1)

    #calculate mean correlations for each sequence
    sum_arrayB=SB_cov_matrix.mean(axis=1)
    
    #find those sequences with a high mean correlation
    keep_seqs=np.where(sum_arrayB>0.55)[0]
    
    #select only high cor
    partition_SB=[partition_SB[i] for i in keep_seqs]
    
    
    plt.hist(sum_arrayB)
    plt.show()
    
   
    #create array on values to iterate over
    x = np.linspace(0,1,20)
    
    #create empty arrays to store information
    pairwise_test_mat = np.empty([20])
    pos_test_mat = np.empty([20,20])
    corr_test_mat = np.empty([20,20])
    nuc_count_mat = np.empty([20,20])
    signif_mat = np.empty([20,20])

        
    #iterate over cutoff point one.     
    for i,x1 in enumerate(x):
        
        #perform postion vote
        correl_val,vote_arr = BoXFP.position_vote(partition_SB,x1,0.5)
        
        
        if correl_val<1:
            pairwise_test_mat[i] = correl_val
        else:
            pairwise_test_mat[i] = np.nan
    #find max correlation cutoff1 value
    max_corr = np.where(pairwise_test_mat == np.nanmax(pairwise_test_mat))
    
    print pairwise_test_mat
    max_corr=np.array(max_corr)
    max_x1=x[max_corr][0]
    
    print max_x1
    
    #print max pairwise correlation
    print np.nanmax(pairwise_test_mat) 
    
    #iterate over all values for cutoff1 and cutoff2
    for i,x1 in enumerate(x):
        for j,x2 in enumerate(x):
            
            #generate consensus sequence
            corel_val,vote_arr = BoXFP.position_vote(partition_SB,x1,x2)
            
            #get nucleotide count
            nuc_count_mat[j] = np.count_nonzero(vote_arr)/float(len(vote_arr))*100
            
            #perform sequence search
            signif_mat[i,j],pos_test_mat[i,j],corr_test_mat[i,j]= BoXFP.sequence_search_area(ref_seq,vote_arr)
        

    #print out metric matrices
    print corr_test_mat
    print signif_mat
    print pos_test_mat.astype(int)
    
    #get max correlation values of cutoff1 and cutoff2
    max_corr = np.where(corr_test_mat == np.nanmax(corr_test_mat))
    
    #print the position of this highest correlation 
    print pos_test_mat[max_corr[0],max_corr[1]]
 
 
   


        
  
