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
from scipy.stats.stats import pearsonr




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
    
    #read data file list  
    with open('data_files_TCV_4.fl' ,'r') as f:
        file_list =f.readlines()
    
    
    #initialise cor and aver arrrays
    corAB_R_25=[]
    corAC_R_25=[]
    corBC_R_25=[]  
    
    corAB_R_50=[]
    corAC_R_50=[]
    corBC_R_50=[] 
    
    corAB_R_100=[]
    corAC_R_100=[]
    corBC_R_100=[]    
    

    corAB_V_25=[]
    corAC_V_25=[]
    corBC_V_25=[]  
    
    corAB_V_50=[]
    corAC_V_50=[]
    corBC_V_50=[]        
    
    corAB_V_100=[]
    corAC_V_100=[]
    corBC_V_100=[]    
    
    avers_V_100=[]
    avers_V_25=[]
    avers_V_50=[]
  
    avers_R_100=[]
    avers_R_25=[]
    avers_R_50=[]
    
    #carry out analysis over windows
    for i in range(10):
        #print the window number
        print 'window '+str(i)
        #get file path for .obj file
        file_path="190902_4_"+str(i)+".obj"
        
        #open data file
        file_1= open(file_path,'rb')
        
        #load data
        data_arr2=pickle.load(file_1)
       
       #find peaks in SM 
        peaksTM=BoXFP.peak_finder(data_arr2,4,.25,TM=1,lower_limit=1400)

        #plot peaks in SM with position outlined
        BoXFP.sm_plotter(data_arr2,peaksTM,file_list)

        #run the partioning algorithm 
        bin_alloc,partition_RX,peak_info= BoXFP.RX_partitioning_replicates(data_arr2,1,.25,ll=1400)
        #partition_S1=BoXFP.TM_partitioning_v2(data_arr2,2)


        #list indices in data file that correspond to each dataset
        R_0=[0,1,2]
        R_25=[6,7,8]
        R_50=[9,10,11]
        R_100=[3,4,5]
        
        V_0=[12,13,14]
        V_25=[18,19,20]
        V_50=[21,22,23]
        V_100=[15,16,17]
 
       
       #Run relaignment on partitioned data
        part_R_0=BoXFP.RX_partition_realignment(partition_RX,bin_alloc,peak_info,R_0,data_arr2)

        part_R_25=BoXFP.RX_partition_realignment(partition_RX,bin_alloc,peak_info,R_25,data_arr2)

        part_R_50=BoXFP.RX_partition_realignment(partition_RX,bin_alloc,peak_info,R_50,data_arr2)

        part_R_100=BoXFP.RX_partition_realignment(partition_RX,bin_alloc,peak_info,R_100,data_arr2)

        part_V_0=BoXFP.RX_partition_realignment(partition_RX,bin_alloc,peak_info,V_0,data_arr2)

        part_V_25=BoXFP.RX_partition_realignment(partition_RX,bin_alloc,peak_info,V_25,data_arr2)

        part_V_50=BoXFP.RX_partition_realignment(partition_RX,bin_alloc,peak_info,V_50,data_arr2)

        part_V_100=BoXFP.RX_partition_realignment(partition_RX,bin_alloc,peak_info,V_100,data_arr2)


        #produce correlation matrices for replicates
        vbv,cor_mat_R_100=BoXFP.correl_assessor(part_R_100,'amp')
        vbv,cor_mat_R_25=BoXFP.correl_assessor(part_R_25,'amp')
        vbv,cor_mat_R_50=BoXFP.correl_assessor(part_R_50,'amp')
        
        vbv,cor_mat_V_100=BoXFP.correl_assessor(part_V_100,'amp')
        vbv,cor_mat_V_25=BoXFP.correl_assessor(part_V_25,'amp')
        vbv,cor_mat_V_50=BoXFP.correl_assessor(part_V_50,'amp')
                
        
        #store correlations
        corAB_R_100=np.append(corAB_R_100,cor_mat_R_100[0,1])
        corAC_R_100=np.append(corAC_R_100,cor_mat_R_100[0,2])
        corBC_R_100=np.append(corBC_R_100,cor_mat_R_100[2,1])
        
        corAB_R_25=np.append(corAB_R_25,cor_mat_R_25[0,1])
        corAC_R_25=np.append(corAC_R_25,cor_mat_R_25[0,2])
        corBC_R_25=np.append(corBC_R_25,cor_mat_R_25[2,1])
        
        corAB_R_50=np.append(corAB_R_50,cor_mat_R_50[0,1])
        corAC_R_50=np.append(corAC_R_50,cor_mat_R_50[0,2])
        corBC_R_50=np.append(corBC_R_50,cor_mat_R_50[2,1])
        
        corAB_V_100=np.append(corAB_V_100,cor_mat_V_100[0,1])
        corAC_V_100=np.append(corAC_V_100,cor_mat_V_100[0,2])
        corBC_V_100=np.append(corBC_V_100,cor_mat_V_100[2,1])
        
        corAB_V_25=np.append(corAB_V_25,cor_mat_V_25[0,1])
        corAC_V_25=np.append(corAC_V_25,cor_mat_V_25[0,2])
        corBC_V_25=np.append(corBC_V_25,cor_mat_V_25[2,1])
        
        corAB_V_50=np.append(corAB_V_50,cor_mat_V_50[0,1])
        corAC_V_50=np.append(corAC_V_50,cor_mat_V_50[0,2])
        corBC_V_50=np.append(corBC_V_50,cor_mat_V_50[2,1])
        
 


        #Calculate areas under peaks in partitioned data
        amp_av_R_0,area_av_R_0,area_sd_R_0=BoXFP.RX_calculator_replicates(part_R_0,data_arr2,R_0)
        amp_av_R_25,area_av_R_25,area_sd_R_25=BoXFP.RX_calculator_replicates(part_R_25,data_arr2,R_25)
        amp_av_R_50,area_av_R_50,area_sd_R_50=BoXFP.RX_calculator_replicates(part_R_50,data_arr2,R_50)
        amp_av_R_100,area_av_R_100,area_sd_R_100=BoXFP.RX_calculator_replicates(part_R_100,data_arr2,R_100)


        amp_av_V_0,area_av_V_0,area_sd_V_0=BoXFP.RX_calculator_replicates(part_V_0,data_arr2,V_0)
        amp_av_V_25,area_av_V_25,area_sd_V_25=BoXFP.RX_calculator_replicates(part_V_25,data_arr2,V_25)
        amp_av_V_50,area_av_V_50,area_sd_V_50=BoXFP.RX_calculator_replicates(part_V_50,data_arr2,V_50)
        amp_av_V_100,area_av_V_100,area_sd_V_100=BoXFP.RX_calculator_replicates(part_V_100,data_arr2,V_100)






        #scale background based on peak amplitudes
        sf_V_25=funcSeqAll.scaleShapeDataWindow(amp_av_V_25,amp_av_V_0,deg=40,rate=0.25,step=10,fit='linear')

        sf_V_50=funcSeqAll.scaleShapeDataWindow(amp_av_V_50,amp_av_V_0,deg=40,rate=0.25,step=10,fit='linear')

        sf_V_100=funcSeqAll.scaleShapeDataWindow(amp_av_V_100,amp_av_V_0,deg=40,rate=0.25,step=10,fit='linear')
        sf_R_25=funcSeqAll.scaleShapeDataWindow(amp_av_R_25,amp_av_R_0,deg=40,rate=0.25,step=10,fit='linear')

        sf_R_50=funcSeqAll.scaleShapeDataWindow(amp_av_R_50,amp_av_R_0,deg=40,rate=0.25,step=10,fit='linear')

        sf_R_100=funcSeqAll.scaleShapeDataWindow(amp_av_R_100,amp_av_R_0,deg=40,rate=0.25,step=10,fit='linear')




        #calculate the background corrected areas
        ad_R_25,nad_R_25,aver_R_25=BoXFP.RX_correction(area_av_R_25,area_av_R_0,sf_R_25)
        ad_R_50,nad_R_50,aver_R_50=BoXFP.RX_correction(area_av_R_50,area_av_R_0,sf_R_50)
        ad_R_100,nad_R_100,aver_R_100=BoXFP.RX_correction(area_av_R_100,area_av_R_0,sf_R_100)

        ad_V_25,nad_V_25,aver_V_25=BoXFP.RX_correction(area_av_V_25,area_av_V_0,sf_V_25)
        ad_V_50,nad_V_50,aver_V_50=BoXFP.RX_correction(area_av_V_50,area_av_V_0,sf_V_50)
        ad_V_100,nad_V_100,aver_V_100=BoXFP.RX_correction(area_av_V_100,area_av_V_0,sf_V_100)
        
        #bind the normalisation factors (avers)
        
        avers_V_100=np.append(avers_V_100,aver_V_100)
        avers_V_25=np.append(avers_V_25,aver_V_25)
        avers_V_50=np.append(avers_V_50,aver_V_50)
        
        
        avers_R_100=np.append(avers_R_100,aver_R_100)
        avers_R_25=np.append(avers_R_25,aver_R_25)
        avers_R_50=np.append(avers_R_50,aver_R_50)



        #calculate the error propogation 
        ad_se_V_25=BoXFP.error_propagation((area_sd_V_25/len(V_25)),sf_V_25*(area_sd_V_0/len(V_0)))

        ad_se_V_50=BoXFP.error_propagation((area_sd_V_50/len(V_50)),sf_V_50*(area_sd_V_0/len(V_0)))

        ad_se_V_100=BoXFP.error_propagation((area_sd_V_100/len(V_100)),sf_V_100*(area_sd_V_0/len(V_0)))

        ad_se_R_25=BoXFP.error_propagation((area_sd_R_25/len(R_25)),sf_R_25*(area_sd_R_0/len(R_0)))

        ad_se_R_50=BoXFP.error_propagation((area_sd_R_50/len(R_50)),sf_R_50*(area_sd_R_0/len(R_0)))

        ad_se_R_100=BoXFP.error_propagation((area_sd_R_100/len(R_100)),sf_R_100*(area_sd_R_0/len(R_0)))


        #calculate normalised standard errors
        nad_se_V_25=ad_se_V_25/(aver_V_25)
        nad_se_V_50=ad_se_V_50/(aver_V_50)
        nad_se_V_100=ad_se_V_100/(aver_V_100)


        nad_se_R_25=ad_se_R_25/(aver_R_25)
        nad_se_R_50=ad_se_R_50/(aver_R_50)
        nad_se_R_100=ad_se_R_100/(aver_R_100)
        
        ####Difference mapping####

        #merge V and T corrected areas
        diffs_25=np.append(ad_V_25,ad_R_25)

        #calculate opercentages average and outliers 
        PO_25,PA_25=funcSeqAll.findPOutlierBox(diffs_25)

        
        diffs_50=np.append(ad_V_50,ad_R_50)

        PO_50,PA_50=funcSeqAll.findPOutlierBox(diffs_50)


        diffs_100=np.append(ad_V_100,ad_R_100)

        PO_100,PA_100=funcSeqAll.findPOutlierBox(diffs_100)


        #renormalise V and T data 
        nad_V_25_2,aver_V_25_2=funcSeqAll.normSimple(ad_V_25,PO_25,PA_25)
        nad_V_50_2,aver_V_50_2=funcSeqAll.normSimple(ad_V_50,PO_50,PA_50)
        nad_V_100_2,aver_V_100_2=funcSeqAll.normSimple(ad_V_100,PO_100,PA_100)

        nad_R_25_2,aver_R_25_2=funcSeqAll.normSimple(ad_R_25,PO_25,PA_25)
        nad_R_50_2,aver_R_50_2=funcSeqAll.normSimple(ad_R_50,PO_50,PA_50)
        nad_R_100_2,aver_R_100_2=funcSeqAll.normSimple(ad_R_100,PO_100,PA_100)


        #calcute difference mapping
        nad_D_25=nad_V_25_2-nad_R_25_2
        nad_D_50=nad_V_50_2-nad_R_50_2
        nad_D_100=nad_V_100_2-nad_R_100_2

        #calculate new errors
        nad_se_R_25=ad_se_R_25/aver_R_25_2
        nad_se_R_50=ad_se_R_50/aver_R_50_2
        nad_se_R_100=ad_se_R_100/aver_R_100_2

        nad_se_V_25=ad_se_V_25/aver_V_25_2
        nad_se_V_50=ad_se_V_50/aver_V_50_2
        nad_se_V_100=ad_se_V_100/aver_V_100_2

        #propagate errors
        nad_se_D_25=BoXFP.error_propagation(nad_se_V_25,nad_se_R_25)
        nad_se_D_50=BoXFP.error_propagation(nad_se_V_50,nad_se_R_50)
        nad_se_D_100=BoXFP.error_propagation(nad_se_V_100,nad_se_R_100)


        #print replicate correlations
        print 'R_0'
        print BoXFP.correl_assessor(part_R_0,'amp')

        print 'R_25'
        print BoXFP.correl_assessor(part_R_25,'amp')

        print 'R_50'
        print BoXFP.correl_assessor(part_R_50,'amp')

        print 'R_100'
        print BoXFP.correl_assessor(part_R_100,'amp')


        print 'V_0'
        print BoXFP.correl_assessor(part_V_0,'amp')

        print 'V_25'
        print BoXFP.correl_assessor(part_V_25,'amp')

        print 'V_50'
        print BoXFP.correl_assessor(part_V_50,'amp')

        print 'V_100'
        print BoXFP.correl_assessor(part_V_100,'amp')



        #print normalisation factors
        print 'aver_V_25 '+str(aver_V_25)
        print 'aver_V_50 '+str(aver_V_50)
        print 'aver_V_10 '+str(aver_V_100)

        print 'aver_R_25 '+str(aver_R_25)
        print 'aver_R_50 '+str(aver_R_50)
        print 'aver_R_10 '+str(aver_R_100)


        #bin the data output
        if i==0:

                nad_R_100_arr=deepcopy(nad_R_100)
                nad_se_R_100_arr=deepcopy(nad_se_R_100)
                nad_R_25_arr=deepcopy(nad_R_25)
                nad_se_R_25_arr=deepcopy(nad_se_R_25)
                nad_R_50_arr=deepcopy(nad_R_50)
                nad_se_R_50_arr=deepcopy(nad_se_R_50)

                nad_V_100_arr=deepcopy(nad_V_100)
                nad_se_V_100_arr=deepcopy(nad_se_V_100)
                nad_V_25_arr=deepcopy(nad_V_25)
                nad_se_V_25_arr=deepcopy(nad_se_V_25)
                nad_V_50_arr=deepcopy(nad_V_50)
                nad_se_V_50_arr=deepcopy(nad_se_V_50)
                
                
                ad_R_100_arr=deepcopy(ad_R_100)
                ad_se_R_100_arr=deepcopy(ad_se_R_100)
                ad_R_25_arr=deepcopy(ad_R_25)
                ad_se_R_25_arr=deepcopy(ad_se_R_25)
                ad_R_50_arr=deepcopy(ad_R_50)
                ad_se_R_50_arr=deepcopy(ad_se_R_50)

                ad_V_100_arr=deepcopy(ad_V_100)
                ad_se_V_100_arr=deepcopy(ad_se_V_100)
                ad_V_25_arr=deepcopy(ad_V_25)
                ad_se_V_25_arr=deepcopy(ad_se_V_25)
                ad_V_50_arr=deepcopy(ad_V_50)
                ad_se_V_50_arr=deepcopy(ad_se_V_50)

                
                nad_D_100_arr=deepcopy(nad_D_100)
                nad_se_D_100_arr=deepcopy(nad_se_D_100)
                nad_D_25_arr=deepcopy(nad_D_25)
                nad_se_D_25_arr=deepcopy(nad_se_D_25)
                nad_D_50_arr=deepcopy(nad_D_50)
                nad_se_D_50_arr=deepcopy(nad_se_D_50)
        else:
            nad_R_100_arr=np.vstack((nad_R_100_arr,nad_R_100))
            nad_se_R_100_arr=np.vstack((nad_se_R_100_arr,nad_se_R_100))

            nad_R_25_arr=np.vstack((nad_R_25_arr,nad_R_25))
            nad_se_R_25_arr=np.vstack((nad_se_R_25_arr,nad_se_R_25))

            nad_R_50_arr=np.vstack((nad_R_50_arr,nad_R_50))
            nad_se_R_50_arr=np.vstack((nad_se_R_50_arr,nad_se_R_50))

            nad_V_100_arr=np.vstack((nad_V_100_arr,nad_V_100))
            nad_se_V_100_arr=np.vstack((nad_se_V_100_arr,nad_se_V_100))

            nad_V_25_arr=np.vstack((nad_V_25_arr,nad_V_25))
            nad_se_V_25_arr=np.vstack((nad_se_V_25_arr,nad_se_V_25))

            nad_V_50_arr=np.vstack((nad_V_50_arr,nad_V_50))
            nad_se_V_50_arr=np.vstack((nad_se_V_50_arr,nad_se_V_50))
            
            
            
            ad_R_100_arr=np.vstack((ad_R_100_arr,ad_R_100))
            ad_se_R_100_arr=np.vstack((ad_se_R_100_arr,ad_se_R_100))

            ad_R_25_arr=np.vstack((ad_R_25_arr,ad_R_25))
            ad_se_R_25_arr=np.vstack((ad_se_R_25_arr,ad_se_R_25))

            ad_R_50_arr=np.vstack((ad_R_50_arr,ad_R_50))
            ad_se_R_50_arr=np.vstack((ad_se_R_50_arr,ad_se_R_50))

            ad_V_100_arr=np.vstack((ad_V_100_arr,ad_V_100))
            ad_se_V_100_arr=np.vstack((ad_se_V_100_arr,ad_se_V_100))

            ad_V_25_arr=np.vstack((ad_V_25_arr,ad_V_25))
            ad_se_V_25_arr=np.vstack((ad_se_V_25_arr,ad_se_V_25))

            ad_V_50_arr=np.vstack((ad_V_50_arr,ad_V_50))
            ad_se_V_50_arr=np.vstack((ad_se_V_50_arr,ad_se_V_50))
            
            

            nad_D_100_arr=np.vstack((nad_D_100_arr,nad_D_100))
            nad_se_D_100_arr=np.vstack((nad_se_D_100_arr,nad_se_D_100))

            nad_D_25_arr=np.vstack((nad_D_25_arr,nad_D_25))
            nad_se_D_25_arr=np.vstack((nad_se_D_25_arr,nad_se_D_25))

            nad_D_50_arr=np.vstack((nad_D_50_arr,nad_D_50))
            nad_se_D_50_arr=np.vstack((nad_se_D_50_arr,nad_se_D_50))
    
    #set up empty correlation matrices for comparing data from different windows
    ca_D_25=np.zeros([len(nad_D_50_arr),len(nad_D_50_arr)])
    ca_D_50=np.zeros([len(nad_D_50_arr),len(nad_D_50_arr)])
    ca_D_100=np.zeros([len(nad_D_50_arr),len(nad_D_50_arr)])
    ca_R_25=np.zeros([len(nad_D_50_arr),len(nad_D_50_arr)])
    ca_R_50=np.zeros([len(nad_D_50_arr),len(nad_D_50_arr)])
    ca_R_100=np.zeros([len(nad_D_50_arr),len(nad_D_50_arr)])
    ca_V_25=np.zeros([len(nad_D_50_arr),len(nad_D_50_arr)])
    ca_V_50=np.zeros([len(nad_D_50_arr),len(nad_D_50_arr)])
    ca_V_100=np.zeros([len(nad_D_50_arr),len(nad_D_50_arr)])
    
    
    #calculate different correlations
    for i in range(len(nad_D_50_arr)):
        for j in range(len(nad_D_50_arr)):
            ca_R_25[i,j]=pearsonr(nad_R_25_arr[i,:],nad_R_25_arr[j,:])[0]
            ca_R_50[i,j]=pearsonr(nad_R_50_arr[i,:],nad_R_50_arr[j,:])[0]
            ca_R_100[i,j]=pearsonr(nad_R_100_arr[i,:],nad_R_100_arr[j,:])[0]
            
            ca_D_25[i,j]=pearsonr(nad_D_25_arr[i,:],nad_D_25_arr[j,:])[0]
            ca_D_50[i,j]=pearsonr(nad_D_50_arr[i,:],nad_D_50_arr[j,:])[0]
            ca_D_100[i,j]=pearsonr(nad_D_100_arr[i,:],nad_D_100_arr[j,:])[0]
            
            
            ca_V_25[i,j]=pearsonr(nad_V_25_arr[i,:],nad_V_25_arr[j,:])[0]
            ca_V_50[i,j]=pearsonr(nad_V_50_arr[i,:],nad_V_50_arr[j,:])[0]
            ca_V_100[i,j]=pearsonr(nad_V_100_arr[i,:],nad_V_100_arr[j,:])[0]
    
    #calculate mean correlation for each window
    mean_ca_R_25=np.mean(ca_R_25,axis=1)
    mean_ca_R_50=np.mean(ca_R_50,axis=1)
    mean_ca_R_100=np.mean(ca_R_100,axis=1)
    
    mean_ca_V_25=np.mean(ca_V_25,axis=1)
    mean_ca_V_50=np.mean(ca_V_50,axis=1)
    mean_ca_V_100=np.mean(ca_V_100,axis=1)
    
    mean_ca_D_25=np.mean(ca_D_25,axis=1)
    mean_ca_D_50=np.mean(ca_D_50,axis=1)
    mean_ca_D_100=np.mean(ca_D_100,axis=1)
    
    print mean_ca_D_50
    print mean_ca_D_25
    print mean_ca_D_100
    
    

    #discard any windows that have poor weak correlation
    keep_ca_D_25=np.where(mean_ca_D_25>0.75)[0]
    keep_ca_D_50=np.where(mean_ca_D_50>0.75)[0]
    keep_ca_D_100=np.where(mean_ca_D_100>0.75)[0]
    
    
    #calculate mean of data over kept windows
    nad_R_25=np.mean(nad_R_25_arr[keep_ca_D_25,:],axis=0)
    nad_R_50=np.mean(nad_R_50_arr[keep_ca_D_50,:],axis=0)
    nad_R_100=np.mean(nad_R_100_arr[keep_ca_D_100,:],axis=0)
    
    nad_se_R_25=np.mean(nad_se_R_25_arr[keep_ca_D_25,:],axis=0)
    nad_se_R_50=np.mean(nad_se_R_50_arr[keep_ca_D_50,:],axis=0)
    nad_se_R_100=np.mean(nad_se_R_100_arr[keep_ca_D_100,:],axis=0)
    
    nad_V_25=np.mean(nad_V_25_arr[keep_ca_D_25,:],axis=0)
    nad_V_50=np.mean(nad_V_50_arr[keep_ca_D_50,:],axis=0)
    nad_V_100=np.mean(nad_V_100_arr[keep_ca_D_100,:],axis=0)
    
    nad_se_V_25=np.mean(nad_se_V_25_arr[keep_ca_D_25,:],axis=0)
    nad_se_V_50=np.mean(nad_se_V_50_arr[keep_ca_D_50,:],axis=0)
    nad_se_V_100=np.mean(nad_se_V_50_arr[keep_ca_D_100,:],axis=0)
    
    
    ad_R_25=np.mean(ad_R_25_arr[keep_ca_D_25,:],axis=0)
    ad_R_50=np.mean(ad_R_50_arr[keep_ca_D_50,:],axis=0)
    ad_R_100=np.mean(ad_R_100_arr[keep_ca_D_100,:],axis=0)
    
    ad_se_R_25=np.mean(ad_se_R_25_arr[keep_ca_D_25,:],axis=0)
    ad_se_R_50=np.mean(ad_se_R_50_arr[keep_ca_D_50,:],axis=0)
    ad_se_R_100=np.mean(ad_se_R_100_arr[keep_ca_D_100,:],axis=0)
    
    ad_V_25=np.mean(ad_V_25_arr[keep_ca_D_25,:],axis=0)
    ad_V_50=np.mean(ad_V_50_arr[keep_ca_D_50,:],axis=0)
    ad_V_100=np.mean(ad_V_100_arr[keep_ca_D_100,:],axis=0)
    
    ad_se_V_25=np.mean(ad_se_V_25_arr[keep_ca_D_25,:],axis=0)
    ad_se_V_50=np.mean(ad_se_V_50_arr[keep_ca_D_50,:],axis=0)
    ad_se_V_100=np.mean(ad_se_V_50_arr[keep_ca_D_100,:],axis=0)
    
    
    nad_D_25=np.mean(nad_D_25_arr[keep_ca_D_25,:],axis=0)
    nad_D_50=np.mean(nad_D_50_arr[keep_ca_D_50,:],axis=0)
    nad_D_100=np.mean(nad_D_100_arr[keep_ca_D_100,:],axis=0)
    
    nad_se_D_25=np.mean(nad_se_D_25_arr[keep_ca_D_25,:],axis=0)
    nad_se_D_50=np.mean(nad_se_D_50_arr[keep_ca_D_50,:],axis=0)
    nad_se_D_100=np.mean(nad_se_D_100_arr[keep_ca_D_100,:],axis=0)
    
    
    #### write normalisation factors to a file###
    avers_arr_V=np.transpose([avers_V_25,avers_V_50,avers_V_100])
    
    avers_labs_V=['aver_V_25','aver_V_50','aver_V_100']
    
    avers_df_V=pd.DataFrame(avers_arr_V,columns=avers_labs_V)
    
    avers_df_V.to_csv('TCV_pri1_avers_V.csv',sep=',',index=False)
    
    avers_arr_R=np.transpose([avers_R_25,avers_R_50,avers_R_100])
    
    avers_labs_R=['aver_R_25','aver_R_50','aver_R_100']
    
    avers_df_R=pd.DataFrame(avers_arr_R,columns=avers_labs_R)
    
    avers_df_R.to_csv('TCV_pri1_avers_R.csv',sep=',',index=False)
    
    
    ##### write window correlation matrices to file
    np.savetxt('cor_mat_D_25.csv',ca_D_25)
    np.savetxt('cor_mat_D_50.csv',ca_D_50)
    np.savetxt('cor_mat_D_100.csv',ca_D_100)
    
    
    #### write mean window correlations to file
    mean_cor_arr=np.transpose([mean_ca_D_25,mean_ca_D_50,mean_ca_D_100])
    
    mean_cor_labs=['cor_D_25','cor_D_50','cov_D_100']
   
    mean_cor_df=pd.DataFrame(mean_cor_arr,columns=mean_cor_labs)
    
    mean_cor_df.to_csv('TCV_pri1_mean_cors.csv',sep=',',index=False)
    
    ##### write replicate correlations to a file
    cor_labs=['corAB_100','corBC_100','corAC_100','corAB_25','corBC_25','corAC_25','corAB_50','corBC_50','corAC_50'] 
    cor_arr_R=np.transpose([corAB_R_100,corBC_R_100,corAC_R_100,corAB_R_25,corBC_R_25,corAC_R_25,corAB_R_50,corBC_R_50,corAC_R_50])
    
    cor_arr_V=np.transpose([corAB_V_100,corBC_V_100,corAC_V_100,corAB_V_25,corBC_V_25,corAC_V_25,corAB_V_50,corBC_V_50,corAC_V_50])
    
    
    
    cor_df_R=pd.DataFrame(cor_arr_R,columns=cor_labs)
    
    cor_df_V=pd.DataFrame(cor_arr_V,columns=cor_labs)
    
    
    
    cor_df_R.to_csv('TCV_pri1_transcript_cor.csv',sep=',',index=False)
    
    cor_df_V.to_csv('TCV_pri1_virion_cor.csv',sep=',',index=False)
    
 
    #specify nucleotide position
    nuc_pos=np.arange(350)+24
    
   ####write normalised data to files
    full_arr_R=np.transpose([nuc_pos,nad_R_25[::-1],nad_se_R_25[::-1],nad_R_50[::-1],nad_se_R_50[::-1],nad_se_R_100[::-1],nad_se_R_100[::-1]])
    
    col_names_R=['genome_position','norm_react_R_25','norm_react_se_R_25','norm_react_R_50','norm_react_se_R_50','norm_react_R_100','norm_react_se_R_100']
    
    full_df_R=pd.DataFrame(full_arr_R,columns=col_names_R)
    
    full_df_R.to_csv('TCV_pri1_transcript_data.csv',sep=',',index=False)
    
   
    full_arr_V=np.transpose([nuc_pos,nad_V_25[::-1],nad_se_V_25[::-1],nad_V_50[::-1],nad_se_V_50[::-1],nad_V_100[::-1],nad_se_V_100[::-1]])
    
    col_names_V=['genome_position','norm_react_V_25','norm_react_se_V_25','norm_react_V_50','norm_react_se_V_50','norm_react_V_100','norm_react_se_V_100']
    
    full_df_V=pd.DataFrame(full_arr_V,columns=col_names_V)
    
    full_df_V.to_csv('TCV_pri1_virion_data.csv',sep=',',index=False)
    
   #####write unnormalised data to files
    full_arr_R=np.transpose([nuc_pos,ad_R_25[::-1],ad_se_R_25[::-1],ad_R_50[::-1],ad_se_R_50[::-1],ad_se_R_100[::-1],ad_se_R_100[::-1]])
    
    col_names_R=['genome_position','ad_react_R_25','ad_react_se_R_25','ad_react_R_50','ad_react_se_R_50','ad_react_R_100','ad_react_se_R_100']
    
    full_df_R=pd.DataFrame(full_arr_R,columns=col_names_R)
    
    full_df_R.to_csv('TCV_pri1_transcript_data_ad.csv',sep=',',index=False)
    
    full_arr_V=np.transpose([nuc_pos,ad_V_25[::-1],ad_se_V_25[::-1],ad_V_50[::-1],ad_se_V_50[::-1],ad_V_100[::-1],ad_se_V_100[::-1]])
    
    col_names_V=['genome_position','ad_react_V_25','ad_react_se_V_25','ad_react_V_50','ad_react_se_V_50','ad_react_V_100','ad_react_se_V_100']
    
    full_df_V=pd.DataFrame(full_arr_V,columns=col_names_V)
    
    full_df_V.to_csv('TCV_pri1_virion_data_ad.csv',sep=',',index=False)
    
    ####write difference mapping to file
    full_arr_D=np.transpose([nuc_pos,nad_D_25[::-1],nad_se_D_25[::-1],nad_D_50[::-1],nad_se_D_50[::-1],nad_D_100[::-1],nad_se_D_100[::-1]])
    
    col_names_D=['genome_position','norm_react_D_25','norm_react_se_D_25','norm_react_D_50','norm_react_se_D_50','norm_react_D_100','norm_react_se_D_100']
    
    full_df_D=pd.DataFrame(full_arr_D,columns=col_names_D)
    
    full_df_D.to_csv('TCV_pri1_difference_data.csv',sep=',',index=False)
    
   
    ####produce snapshots of images of reactivity profiles
    BoXFP.sequence_snapshots(nuc_pos,nad_R_25,nad_se_R_25,col='g',virus='TCV',primer='pri1',condition='R',treatment='25')
    
    BoXFP.sequence_snapshots(nuc_pos,nad_R_50,nad_se_R_50,col='g',virus='TCV',primer='pri1',condition='R',treatment='50')
    
    BoXFP.sequence_snapshots(nuc_pos,nad_R_100,nad_se_R_100,col='g',virus='TCV',primer='pri1',condition='R',treatment='100')
    
    BoXFP.sequence_snapshots(nuc_pos,nad_V_25,nad_se_V_25,col='b',virus='TCV',primer='pri1',condition='V',treatment='25')
    
    BoXFP.sequence_snapshots(nuc_pos,nad_V_50,nad_se_V_50,col='b',virus='TCV',primer='pri1',condition='V',treatment='50')
    BoXFP.sequence_snapshots(nuc_pos,nad_V_100,nad_se_V_100,col='b',virus='TCV',primer='pri1',condition='V',treatment='100')
   
    BoXFP.sequence_snapshots(nuc_pos,nad_D_25,nad_se_D_25,col='m',virus='TCV',primer='pri1',condition='V-T',treatment='25',diff=True)
    
    BoXFP.sequence_snapshots(nuc_pos,nad_D_50,nad_se_D_50,col='m',virus='TCV',primer='pri1',condition='V-T',treatment='50',diff=True)

    BoXFP.sequence_snapshots(nuc_pos,nad_D_100,nad_se_D_100,col='m',virus='TCV',primer='pri1',condition='V-T',treatment='100',diff=True)

   
    
