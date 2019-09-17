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
    
    #read in data file list    
    with open('data_files.fl' ,'r') as f:
        file_list =f.readlines()
  
    #initialise replicate correlation storage arrays
    corAB_R_50=[]
    corAC_R_50=[]
    corBC_R_50=[]    
    
    corAB_R_100=[]
    corAC_R_100=[]
    corBC_R_100=[]
    
    corAB_V_25=[]
    corAC_V_25=[]
    corBC_V_25=[]
    
    corAB_R_25=[]
    corAC_R_25=[]
    corBC_R_25=[]
    
    corAB_V_50=[]
    corAC_V_50=[]
    corBC_V_50=[]    
    
    corAB_V_100=[]
    corAC_V_100=[]
    corBC_V_100=[]
    
    corAB_G_25=[]
    corAC_G_25=[]
    corBC_G_25=[]
    
    #initialise normalisation factor storage arrays
    avers_V_100=[]
    avers_V_25=[]
    avers_V_50=[]
  
    avers_R_100=[]
    avers_R_25=[]
    avers_R_50=[]
    
    #perform calculations over windows
    for i in range(10):
        #output window number
        print 'window '+str(i)
        
        #open data for window
        file_path="190822_"+str(i)+".obj"

        file_1= open(file_path,'rb')
        
        data_arr2=pickle.load(file_1)

        #calculate size marker peaks positions 
        peka,peaksTM=sam_funcs.peak_finder(data_arr2,4,.25,lower_limit=1400)
        
        #plot size marker positions on size marker traces. 
        sam_funcs.sm_plotter(data_arr2,peaksTM,file_list)

        
        #run partition function on the data
        
        ####NOTE THE CUTOFF VALUE TO NEGATE THE EXIT PEAK.####
        bin_alloc,partition_RX,peak_info=sam_funcs.RX_partitioning_replicates(data_arr2,1,.25,Cap=5000,ll=1400,tm=0,tm_cutoff=19)
        

        #index the data according to treatment conditions
        pr1_V_0=[0,1,2]
        pr1_V_25=[6,8]
        pr1_V_50=[9,10,11]
        pr1_V_100=[3,4,5]

        pr1_R_0=[36,37,38]
        pr1_R_25=[42,43,44]
        pr1_R_50=[45,46,47]
        pr1_R_100=[39,40,41]

        pr19_V_0=[12,14]
        pr19_V_25=[18,19,20]
        pr19_V_50=[21,22,23]
        pr19_V_100=[15,16,17]

        pr19_R_0=[48,49,50]
        pr19_R_25=[54,55,56]
        pr19_R_50=[57,58,59]
        pr19_R_100=[51,52,53]

        pr21_V_0=[24,25,26]
        pr21_V_25=[30,31,32]
        pr21_V_50=[33,34,35]
        pr21_V_100=[27,28,29]

        pr21_R_0=[60,61,62]
        pr21_R_25=[66,67,68]
        pr21_R_50=[69,70,71]
        pr21_R_100=[63,64,65]

        part_pr1_V_100=sam_funcs.RX_partition_realignment(partition_RX,bin_alloc,peak_info,pr1_V_100,data_arr2,tm=0,tm_cutoff=19)
        part_pr1_R_0=sam_funcs.RX_partition_realignment(partition_RX,bin_alloc,peak_info,pr1_R_0,data_arr2,tm=0,tm_cutoff=19)

        part_pr1_R_25=sam_funcs.RX_partition_realignment(partition_RX,bin_alloc,peak_info,pr1_R_25,data_arr2,tm=0,tm_cutoff=19)
 
        part_pr1_R_50=sam_funcs.RX_partition_realignment(partition_RX,bin_alloc,peak_info,pr1_R_50,data_arr2,tm=0,tm_cutoff=19)
        part_pr1_R_100=sam_funcs.RX_partition_realignment(partition_RX,bin_alloc,peak_info,pr1_R_100,data_arr2,tm=0,tm_cutoff=19)

        part_pr1_V_0=sam_funcs.RX_partition_realignment(partition_RX,bin_alloc,peak_info,pr1_V_0,data_arr2,tm=0,tm_cutoff=19)

        part_pr1_V_25=sam_funcs.RX_partition_realignment(partition_RX,bin_alloc,peak_info,pr1_V_25,data_arr2,tm=0,tm_cutoff=19)

        part_pr1_V_50=sam_funcs.RX_partition_realignment(partition_RX,bin_alloc,peak_info,pr1_V_50,data_arr2,tm=0,tm_cutoff=19)

                
        #print replicate correlation values
        print 'pr1_R_0'
        print sam_funcs.correl_assessor(part_pr1_R_0,'amp')

        print 'pr1_R_25'
        print sam_funcs.correl_assessor(part_pr1_R_25,'amp')

        print 'pr1_R_50'
        print sam_funcs.correl_assessor(part_pr1_R_50,'amp')

        print 'pr1_R_100'
        print sam_funcs.correl_assessor(part_pr1_R_100,'amp')


        print 'pr1_V_0'
        print sam_funcs.correl_assessor(part_pr1_V_0,'amp')

        print 'pr1_V_25'
        print sam_funcs.correl_assessor(part_pr1_V_25,'amp')

        print 'pr1_V_50'
        print sam_funcs.correl_assessor(part_pr1_V_50,'amp')

        print 'pr1_V_100'
        print sam_funcs.correl_assessor(part_pr1_V_100,'amp')


        #calculate average amppitudes and areas of peaks and the standard deviation of the areas
        amp_av_pr1_R_0,area_av_pr1_R_0,area_sd_pr1_R_0=sam_funcs.RX_calculator_replicates(part_pr1_R_0,data_arr2,pr1_R_0)
        amp_av_pr1_R_25,area_av_pr1_R_25,area_sd_pr1_R_25=sam_funcs.RX_calculator_replicates(part_pr1_R_25,data_arr2,pr1_R_25)
        amp_av_pr1_R_50,area_av_pr1_R_50,area_sd_pr1_R_50=sam_funcs.RX_calculator_replicates(part_pr1_R_50,data_arr2,pr1_R_50)
        amp_av_pr1_R_100,area_av_pr1_R_100,area_sd_pr1_R_100=sam_funcs.RX_calculator_replicates(part_pr1_R_100,data_arr2,pr1_R_100)


        amp_av_pr1_V_0,area_av_pr1_V_0,area_sd_pr1_V_0=sam_funcs.RX_calculator_replicates(part_pr1_V_0,data_arr2,pr1_V_0)
        amp_av_pr1_V_25,area_av_pr1_V_25,area_sd_pr1_V_25=sam_funcs.RX_calculator_replicates(part_pr1_V_25,data_arr2,pr1_V_25)
        amp_av_pr1_V_50,area_av_pr1_V_50,area_sd_pr1_V_50=sam_funcs.RX_calculator_replicates(part_pr1_V_50,data_arr2,pr1_V_50)
        amp_av_pr1_V_100,area_av_pr1_V_100,area_sd_pr1_V_100=sam_funcs.RX_calculator_replicates(part_pr1_V_100,data_arr2,pr1_V_100)


        #calculate replicate correlation values
        vbv,cor_mat_R_100=sam_funcs.correl_assessor(part_pr1_R_100,'amp')
        vbv,cor_mat_R_25=sam_funcs.correl_assessor(part_pr1_R_25,'amp')
        vbv,cor_mat_R_50=sam_funcs.correl_assessor(part_pr1_R_50,'amp')
        vbv,cor_mat_V_100=sam_funcs.correl_assessor(part_pr1_V_100,'amp')
        vbv,cor_mat_V_25=sam_funcs.correl_assessor(part_pr1_V_25,'amp')
        vbv,cor_mat_V_50=sam_funcs.correl_assessor(part_pr1_V_50,'amp')
        
        #store replicate values
        corAB_R_25=np.append(corAB_R_25,cor_mat_R_25[0,1])
        corAC_R_25=np.append(corAC_R_25,cor_mat_R_25[0,2])
        corBC_R_25=np.append(corBC_R_25,cor_mat_R_25[2,1])
        
        corAB_R_50=np.append(corAB_R_50,cor_mat_R_50[0,1])
        corAC_R_50=np.append(corAC_R_50,cor_mat_R_50[0,2])
        corBC_R_50=np.append(corBC_R_50,cor_mat_R_50[2,1])
        
        corAB_R_100=np.append(corAB_R_100,cor_mat_R_100[0,1])
        corAC_R_100=np.append(corAC_R_100,cor_mat_R_100[0,2])
        corBC_R_100=np.append(corBC_R_100,cor_mat_R_100[2,1])
        
        corAB_V_25=np.append(corAB_V_25,cor_mat_V_25[0,1])
        #corAC_V_25=np.append(corAC_V_25,cor_mat_V_25[0,2])
        #corBC_V_25=np.append(corBC_V_25,cor_mat_V_25[2,1])
        
        corAB_V_50=np.append(corAB_V_50,cor_mat_V_50[0,1])
        corAC_V_50=np.append(corAC_V_50,cor_mat_V_50[0,2])
        corBC_V_50=np.append(corBC_V_50,cor_mat_V_50[2,1])
        
        corAB_V_100=np.append(corAB_V_100,cor_mat_V_100[0,1])
        corAC_V_100=np.append(corAC_V_100,cor_mat_V_100[0,2])
        corBC_V_100=np.append(corBC_V_100,cor_mat_V_100[2,1])

        
        #calculate the scaling factors of the background 
        sf_pr1_V_25=funcSeqAll.scaleShapeDataWindow(amp_av_pr1_V_25,amp_av_pr1_V_0,deg=40,rate=0.25,step=10,fit='linear')

        sf_pr1_V_50=funcSeqAll.scaleShapeDataWindow(amp_av_pr1_V_50,amp_av_pr1_V_0,deg=40,rate=0.25,step=10,fit='linear')

        sf_pr1_V_100=funcSeqAll.scaleShapeDataWindow(amp_av_pr1_V_100,amp_av_pr1_V_0,deg=40,rate=0.25,step=10,fit='linear')
        sf_pr1_R_25=funcSeqAll.scaleShapeDataWindow(amp_av_pr1_R_25,amp_av_pr1_R_0,deg=40,rate=0.25,step=10,fit='linear')

        sf_pr1_R_50=funcSeqAll.scaleShapeDataWindow(amp_av_pr1_R_50,amp_av_pr1_R_0,deg=40,rate=0.25,step=10,fit='linear')

        sf_pr1_R_100=funcSeqAll.scaleShapeDataWindow(amp_av_pr1_R_100,amp_av_pr1_R_0,deg=40,rate=0.25,step=10,fit='linear')


        #calculate the background corrected areas and normalised reactiivities 
        ad_pr1_R_25,nad_pr1_R_25,aver_pr1_R_25=sam_funcs.RX_correction(area_av_pr1_R_25,area_av_pr1_R_0,sf_pr1_R_25)
        ad_pr1_R_50,nad_pr1_R_50,aver_pr1_R_50=sam_funcs.RX_correction(area_av_pr1_R_50,area_av_pr1_R_0,sf_pr1_R_50)
        ad_pr1_R_100,nad_pr1_R_100,aver_pr1_R_100=sam_funcs.RX_correction(area_av_pr1_R_100,area_av_pr1_R_0,sf_pr1_R_100)

        ad_pr1_V_25,nad_pr1_V_25,aver_pr1_V_25=sam_funcs.RX_correction(area_av_pr1_V_25,area_av_pr1_V_0,sf_pr1_V_25)
        ad_pr1_V_50,nad_pr1_V_50,aver_pr1_V_50=sam_funcs.RX_correction(area_av_pr1_V_50,area_av_pr1_V_0,sf_pr1_V_50)
        ad_pr1_V_100,nad_pr1_V_100,aver_pr1_V_100=sam_funcs.RX_correction(area_av_pr1_V_100,area_av_pr1_V_0,sf_pr1_V_100)


        #print the normalisation factors

        print 'aver_pr1_V_25 '+str(aver_pr1_V_25)
        print 'aver_pr1_V_50 '+str(aver_pr1_V_50)
        print 'aver_pr1_V_100 '+str(aver_pr1_V_100)

        print 'aver_pr1_R_25 '+str(aver_pr1_R_25)
        print 'aver_pr1_R_50 '+str(aver_pr1_R_50)
        print 'aver_pr1_R_100 '+str(aver_pr1_R_100)
        

        #store normalisation factors
        avers_V_100=np.append(avers_V_100,aver_pr1_V_100)
        avers_V_25=np.append(avers_V_25,aver_pr1_V_25)
        avers_V_50=np.append(avers_V_50,aver_pr1_V_50)
        
        avers_R_100=np.append(avers_R_100,aver_pr1_R_100)
        avers_R_25=np.append(avers_R_25,aver_pr1_R_25)
        avers_R_50=np.append(avers_R_50,aver_pr1_R_50)
       
        #error propagation
        ad_se_pr1_V_25=sam_funcs.error_propagation((area_sd_pr1_V_25/len(pr1_V_25)),sf_pr1_V_25*(area_sd_pr1_V_0/len(pr1_V_0)))

        ad_se_pr1_V_50=sam_funcs.error_propagation((area_sd_pr1_V_50/len(pr1_V_50)),sf_pr1_V_50*(area_sd_pr1_V_0/len(pr1_V_0)))

        ad_se_pr1_V_100=sam_funcs.error_propagation((area_sd_pr1_V_100/len(pr1_V_100)),sf_pr1_V_100*(area_sd_pr1_V_0/len(pr1_V_0)))

        ad_se_pr1_R_25=sam_funcs.error_propagation((area_sd_pr1_R_25/len(pr1_R_25)),sf_pr1_R_25*(area_sd_pr1_R_0/len(pr1_R_0)))

        ad_se_pr1_R_50=sam_funcs.error_propagation((area_sd_pr1_R_50/len(pr1_R_50)),sf_pr1_R_50*(area_sd_pr1_R_0/len(pr1_R_0)))

        ad_se_pr1_R_100=sam_funcs.error_propagation((area_sd_pr1_R_100/len(pr1_R_100)),sf_pr1_R_100*(area_sd_pr1_R_0/len(pr1_R_0)))

        #calculate the normalised standard errors

        nad_se_pr1_V_25=ad_se_pr1_V_25/(aver_pr1_V_25)
        nad_se_pr1_V_50=ad_se_pr1_V_50/(aver_pr1_V_50)
        nad_se_pr1_V_100=ad_se_pr1_V_100/(aver_pr1_V_100)


        nad_se_pr1_R_25=ad_se_pr1_R_25/(aver_pr1_R_25)
        nad_se_pr1_R_50=ad_se_pr1_R_50/(aver_pr1_R_50)
        nad_se_pr1_R_100=ad_se_pr1_R_100/(aver_pr1_R_100)


        ####difference mapping

        diffs_pr1_25=np.append(ad_pr1_V_25,ad_pr1_R_25)

        PO_pr1_25,PA_pr1_25=funcSeqAll.findPOutlierBox(diffs_pr1_25)

        diffs_pr1_50=np.append(ad_pr1_V_50,ad_pr1_R_50)

        PO_pr1_50,PA_pr1_50=funcSeqAll.findPOutlierBox(diffs_pr1_50)

        diffs_pr1_100=np.append(ad_pr1_V_100,ad_pr1_R_100)

        PO_pr1_100,PA_pr1_100=funcSeqAll.findPOutlierBox(diffs_pr1_100)

        nad_pr1_V_25_2,aver_pr1_V_25_2=funcSeqAll.normSimple(ad_pr1_V_25,PO_pr1_25,PA_pr1_25)
        nad_pr1_V_50_2,aver_pr1_V_50_2=funcSeqAll.normSimple(ad_pr1_V_50,PO_pr1_50,PA_pr1_50)
        nad_pr1_V_100_2,aver_pr1_V_100_2=funcSeqAll.normSimple(ad_pr1_V_100,PO_pr1_100,PA_pr1_100)

        nad_pr1_R_25_2,aver_pr1_R_25_2=funcSeqAll.normSimple(ad_pr1_R_25,PO_pr1_25,PA_pr1_25)
        nad_pr1_R_50_2,aver_pr1_R_50_2=funcSeqAll.normSimple(ad_pr1_R_50,PO_pr1_50,PA_pr1_50)
        nad_pr1_R_100_2,aver_pr1_R_100_2=funcSeqAll.normSimple(ad_pr1_R_100,PO_pr1_100,PA_pr1_100)


        
        nad_pr1_D_25=nad_pr1_V_25_2-nad_pr1_R_25_2
        nad_pr1_D_50=nad_pr1_V_50_2-nad_pr1_R_50_2
        nad_pr1_D_100=nad_pr1_V_100_2-nad_pr1_R_100_2

        nad_se_pr1_R_25=ad_se_pr1_R_25/aver_pr1_R_25_2
        nad_se_pr1_R_50=ad_se_pr1_R_50/aver_pr1_R_50_2
        nad_se_pr1_R_100=ad_se_pr1_R_100/aver_pr1_R_100_2

        nad_se_pr1_V_25=ad_se_pr1_V_25/aver_pr1_V_25_2
        nad_se_pr1_V_50=ad_se_pr1_V_50/aver_pr1_V_50_2
        nad_se_pr1_V_100=ad_se_pr1_V_100/aver_pr1_V_100_2





        nad_se_pr1_D_25=sam_funcs.error_propagation(nad_se_pr1_V_25,nad_se_pr1_R_25)
        nad_se_pr1_D_50=sam_funcs.error_propagation(nad_se_pr1_V_50,nad_se_pr1_R_50)
        nad_se_pr1_D_100=sam_funcs.error_propagation(nad_se_pr1_V_100,nad_se_pr1_R_100)

        
        
        #store data
        if i==0:

            nad_R_100_arr=deepcopy(nad_pr1_R_100)
            nad_se_R_100_arr=deepcopy(nad_se_pr1_R_100)
            nad_R_25_arr=deepcopy(nad_pr1_R_25)
            nad_se_R_25_arr=deepcopy(nad_se_pr1_R_25)
            nad_R_50_arr=deepcopy(nad_pr1_R_50)
            nad_se_R_50_arr=deepcopy(nad_se_pr1_R_50)

            nad_V_100_arr=deepcopy(nad_pr1_V_100)
            nad_se_V_100_arr=deepcopy(nad_se_pr1_V_100)
            nad_V_25_arr=deepcopy(nad_pr1_V_25)
            nad_se_V_25_arr=deepcopy(nad_se_pr1_V_25)
            nad_V_50_arr=deepcopy(nad_pr1_V_50)
            nad_se_V_50_arr=deepcopy(nad_se_pr1_V_50)
            
            ad_R_100_arr=deepcopy(ad_pr1_R_100)
            ad_se_R_100_arr=deepcopy(ad_se_pr1_R_100)
            ad_R_25_arr=deepcopy(ad_pr1_R_25)
            ad_se_R_25_arr=deepcopy(ad_se_pr1_R_25)
            ad_R_50_arr=deepcopy(ad_pr1_R_50)
            ad_se_R_50_arr=deepcopy(ad_se_pr1_R_50)

            ad_V_100_arr=deepcopy(ad_pr1_V_100)
            ad_se_V_100_arr=deepcopy(ad_se_pr1_V_100)
            ad_V_25_arr=deepcopy(ad_pr1_V_25)
            ad_se_V_25_arr=deepcopy(ad_se_pr1_V_25)
            ad_V_50_arr=deepcopy(ad_pr1_V_50)
            ad_se_V_50_arr=deepcopy(ad_se_pr1_V_50)

            nad_D_100_arr=deepcopy(nad_pr1_D_100)
            nad_se_D_100_arr=deepcopy(nad_se_pr1_D_100)
            nad_D_25_arr=deepcopy(nad_pr1_D_25)
            nad_se_D_25_arr=deepcopy(nad_se_pr1_D_25)
            nad_D_50_arr=deepcopy(nad_pr1_D_50)
            nad_se_D_50_arr=deepcopy(nad_se_pr1_D_50)
        else:
            nad_R_100_arr=np.vstack((nad_R_100_arr,nad_pr1_R_100))
            nad_se_R_100_arr=np.vstack((nad_se_R_100_arr,nad_se_pr1_R_100))

            nad_R_25_arr=np.vstack((nad_R_25_arr,nad_pr1_R_25))
            nad_se_R_25_arr=np.vstack((nad_se_R_25_arr,nad_se_pr1_R_25))

            nad_R_50_arr=np.vstack((nad_R_50_arr,nad_pr1_R_50))
            nad_se_R_50_arr=np.vstack((nad_se_R_50_arr,nad_se_pr1_R_50))

            nad_V_100_arr=np.vstack((nad_V_100_arr,nad_pr1_V_100))
            nad_se_V_100_arr=np.vstack((nad_se_V_100_arr,nad_se_pr1_V_100))

            nad_V_25_arr=np.vstack((nad_V_25_arr,nad_pr1_V_25))
            nad_se_V_25_arr=np.vstack((nad_se_V_25_arr,nad_se_pr1_V_25))

            nad_V_50_arr=np.vstack((nad_V_50_arr,nad_pr1_V_50))
            nad_se_V_50_arr=np.vstack((nad_se_V_50_arr,nad_se_pr1_V_50))
            
            
            ad_R_100_arr=np.vstack((ad_R_100_arr,ad_pr1_R_100))
            ad_se_R_100_arr=np.vstack((ad_se_R_100_arr,ad_se_pr1_R_100))

            ad_R_25_arr=np.vstack((ad_R_25_arr,ad_pr1_R_25))
            ad_se_R_25_arr=np.vstack((ad_se_R_25_arr,ad_se_pr1_R_25))

            ad_R_50_arr=np.vstack((ad_R_50_arr,ad_pr1_R_50))
            ad_se_R_50_arr=np.vstack((ad_se_R_50_arr,ad_se_pr1_R_50))

            ad_V_100_arr=np.vstack((ad_V_100_arr,ad_pr1_V_100))
            ad_se_V_100_arr=np.vstack((ad_se_V_100_arr,ad_se_pr1_V_100))

            ad_V_25_arr=np.vstack((ad_V_25_arr,ad_pr1_V_25))
            ad_se_V_25_arr=np.vstack((ad_se_V_25_arr,ad_se_pr1_V_25))

            ad_V_50_arr=np.vstack((ad_V_50_arr,ad_pr1_V_50))
            ad_se_V_50_arr=np.vstack((ad_se_V_50_arr,ad_se_pr1_V_50))

            nad_D_100_arr=np.vstack((nad_D_100_arr,nad_pr1_D_100))
            nad_se_D_100_arr=np.vstack((nad_se_D_100_arr,nad_se_pr1_D_100))

            nad_D_25_arr=np.vstack((nad_D_25_arr,nad_pr1_D_25))
            nad_se_D_25_arr=np.vstack((nad_se_D_25_arr,nad_se_pr1_D_25))

            nad_D_50_arr=np.vstack((nad_D_50_arr,nad_pr1_D_50))
            nad_se_D_50_arr=np.vstack((nad_se_D_50_arr,nad_se_pr1_D_50))
    
    #intialise correlation matrices for window comparisons
    ca_D_25=np.zeros([len(nad_D_50_arr),len(nad_D_50_arr)])
    ca_D_50=np.zeros([len(nad_D_50_arr),len(nad_D_50_arr)])
    ca_D_100=np.zeros([len(nad_D_50_arr),len(nad_D_50_arr)])
    ca_R_25=np.zeros([len(nad_D_50_arr),len(nad_D_50_arr)])
    ca_R_50=np.zeros([len(nad_D_50_arr),len(nad_D_50_arr)])
    ca_R_100=np.zeros([len(nad_D_50_arr),len(nad_D_50_arr)])
    ca_V_25=np.zeros([len(nad_D_50_arr),len(nad_D_50_arr)])
    ca_V_50=np.zeros([len(nad_D_50_arr),len(nad_D_50_arr)])
    ca_V_100=np.zeros([len(nad_D_50_arr),len(nad_D_50_arr)])
    
    
    #calculate correlations between windows
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
            
            
    # calculate mean correlation value for each window
    mean_ca_R_25=np.mean(ca_R_25,axis=1)
    mean_ca_R_50=np.mean(ca_R_50,axis=1)
    mean_ca_R_100=np.mean(ca_R_100,axis=1)
    
    mean_ca_V_25=np.mean(ca_V_25,axis=1)
    mean_ca_V_50=np.mean(ca_V_50,axis=1)
    mean_ca_V_100=np.mean(ca_V_100,axis=1)
    
    mean_ca_D_25=np.mean(ca_D_25,axis=1)
    mean_ca_D_50=np.mean(ca_D_50,axis=1)
    mean_ca_D_100=np.mean(ca_D_100,axis=1)
    
    #keep those windows with a mean correlation above a threshold
    keep_ca_R_25=np.where(mean_ca_R_25>0.7)
    keep_ca_R_50=np.where(mean_ca_R_50>0.7)[0]
    keep_ca_R_100=np.where(mean_ca_R_100>0.7)[0]
    
    keep_ca_V_25=np.where(mean_ca_V_25>0.7)[0]
    keep_ca_V_50=np.where(mean_ca_V_50>0.7)[0]
    keep_ca_V_100=np.where(mean_ca_V_100>0.7)[0]
    
    keep_ca_D_25=np.where(mean_ca_D_25>0.7)[0]
    keep_ca_D_50=np.where(mean_ca_D_50>0.7)[0]
    keep_ca_D_100=np.where(mean_ca_D_100>0.7)[0]
    
    print np.count_nonzero(mean_ca_R_25>0.7)
    print np.count_nonzero(mean_ca_R_50>0.7)
    
    print np.count_nonzero(mean_ca_V_25>0.7)
    print np.count_nonzero(mean_ca_V_50>0.7)
    
    print np.count_nonzero(mean_ca_D_25>0.7)
    print np.count_nonzero(mean_ca_D_50>0.7)
    
    ####write replicate correlations to file
    cor_labs_R=['corAB_25','corAC_25','corBC_25','corAB_50','corAC_50','corBC_50','corAB_100','corBC_100','corAC_100']
    
    cor_labs_V=['corAB_25','corAB_50','corBC_50','corAC_50','corAB_100','corBC_100','corAC_100']
    
    cor_arr_R=np.transpose([corAB_R_25,corBC_R_25,corAC_R_25,corAB_R_50,corBC_R_50,corAC_R_50,corAB_V_100,corBC_V_100,corAC_V_100])
    cor_arr_V=np.transpose([corAB_V_25,corAB_V_50,corBC_V_50,corAC_V_50,corAB_V_100,corBC_V_100,corAC_V_100])

    
    cor_df_R=pd.DataFrame(cor_arr_R,columns=cor_labs_R)
    
    cor_df_V=pd.DataFrame(cor_arr_V,columns=cor_labs_V)
        
    cor_df_R.to_csv('MS2_pri1_transcript_cor.csv',sep=',',index=False)
    
    cor_df_V.to_csv('MS2_pri1_virion_cor.csv',sep=',',index=False)
    
    ####write windowing correlations to file
    np.savetxt('cor_mat_D_25.csv',ca_D_25)
    np.savetxt('cor_mat_D_50.csv',ca_D_50)
    np.savetxt('cor_mat_D_100.csv',ca_D_100)
    
    mean_cor_arr=np.transpose([mean_ca_D_100,mean_ca_D_50])
    
    mean_cor_labs=['cor_D_100','cov_D_50']
   
    mean_cor_df=pd.DataFrame(mean_cor_arr,columns=mean_cor_labs)
    
    mean_cor_df.to_csv('MS2_pri1_mean_cors.csv',sep=',',index=False)
    
    ####write normalisation factors to file
    avers_arr_V=np.transpose([avers_V_100,avers_V_25,avers_V_50])
    
    avers_labs_V=['aver_V_100','aver_V_25','aver_V_50']
    
    avers_df_V=pd.DataFrame(avers_arr_V,columns=avers_labs_V)
    
    avers_df_V.to_csv('MS2_pri1_avers_V.csv',sep=',',index=False)
    
    avers_arr_R=np.transpose([avers_R_100,avers_R_50])
    
    avers_labs_R=['aver_R_100','aver_R_50']
    
    avers_df_R=pd.DataFrame(avers_arr_R,columns=avers_labs_R)
    
    avers_df_R.to_csv('MS2_pri1_avers_R.csv',sep=',',index=False)
    
    
    #calculate mean of profiles for windows kept
    nad_pr1_R_25=np.mean(nad_R_25_arr[keep_ca_D_25,:],axis=0)
    nad_pr1_R_50=np.mean(nad_R_50_arr[keep_ca_D_50,:],axis=0)
    nad_pr1_R_100=np.mean(nad_R_100_arr[keep_ca_D_100,:],axis=0)
    
    nad_se_pr1_R_25=np.mean(nad_se_R_25_arr[keep_ca_D_25,:],axis=0)
    nad_se_pr1_R_50=np.mean(nad_se_R_50_arr[keep_ca_D_50,:],axis=0)
    nad_se_pr1_R_100=np.mean(nad_se_R_100_arr[keep_ca_D_100,:],axis=0)
    
    nad_pr1_V_25=np.mean(nad_V_25_arr[keep_ca_V_25,:],axis=0)
    nad_pr1_V_50=np.mean(nad_V_50_arr[keep_ca_V_50,:],axis=0)
    nad_pr1_V_100=np.mean(nad_V_100_arr[keep_ca_V_100,:],axis=0)
    
    nad_se_pr1_V_25=np.mean(nad_se_V_25_arr[keep_ca_V_25,:],axis=0)
    nad_se_pr1_V_50=np.mean(nad_se_V_50_arr[keep_ca_V_50,:],axis=0)
    nad_se_pr1_V_100=np.mean(nad_se_V_50_arr[keep_ca_V_100,:],axis=0)
    
    
    ad_pr1_R_25=np.mean(ad_R_25_arr[keep_ca_D_25,:],axis=0)
    ad_pr1_R_50=np.mean(ad_R_50_arr[keep_ca_D_50,:],axis=0)
    ad_pr1_R_100=np.mean(ad_R_100_arr[keep_ca_D_100,:],axis=0)
    
    ad_se_pr1_R_25=np.mean(ad_se_R_25_arr[keep_ca_D_25,:],axis=0)
    ad_se_pr1_R_50=np.mean(ad_se_R_50_arr[keep_ca_D_50,:],axis=0)
    ad_se_pr1_R_100=np.mean(ad_se_R_100_arr[keep_ca_D_100,:],axis=0)
    
    ad_pr1_V_25=np.mean(ad_V_25_arr[keep_ca_V_25,:],axis=0)
    ad_pr1_V_50=np.mean(ad_V_50_arr[keep_ca_V_50,:],axis=0)
    ad_pr1_V_100=np.mean(ad_V_100_arr[keep_ca_V_100,:],axis=0)
    
    ad_se_pr1_V_25=np.mean(ad_se_V_25_arr[keep_ca_V_25,:],axis=0)
    ad_se_pr1_V_50=np.mean(ad_se_V_50_arr[keep_ca_V_50,:],axis=0)
    ad_se_pr1_V_100=np.mean(ad_se_V_50_arr[keep_ca_V_100,:],axis=0)
    
    nad_pr1_D_25=np.mean(nad_D_25_arr[keep_ca_D_25,:],axis=0)
    nad_pr1_D_50=np.mean(nad_D_50_arr[keep_ca_D_50,:],axis=0)
    nad_pr1_D_100=np.mean(nad_D_100_arr[keep_ca_D_100,:],axis=0)
    
    nad_se_pr1_D_25=np.mean(nad_se_D_25_arr[keep_ca_D_25,:],axis=1)
    nad_se_pr1_D_50=np.mean(nad_se_D_50_arr[keep_ca_D_50,:],axis=0)
    nad_se_pr1_D_100=np.mean(nad_se_D_100_arr[keep_ca_D_100,:],axis=0)
    
    #specify genome position
    nuc_pos1=np.arange(310)+32
    
    ####write normalised data to files####
    full_arr_pr1_R=np.transpose([nuc_pos1,nad_pr1_R_50[::-1],nad_se_pr1_R_50[::-1],nad_pr1_R_100[::-1],nad_se_pr1_R_100[::-1]])
    
    col_names_pr1_R=['genome_position','norm_react_R_50','norm_react_se_R_50','norm_react_R_100','norm_react_se_R_100']
    
    full_df_pr1_R=pd.DataFrame(full_arr_pr1_R,columns=col_names_pr1_R)
    
    full_df_pr1_R.to_csv('MS2_pri1_transcript_data_3.csv',sep=',',index=False)
    
    full_arr_pr1_V=np.transpose([nuc_pos1,nad_pr1_V_25[::-1],nad_se_pr1_V_25[::-1],nad_pr1_V_50[::-1],nad_se_pr1_V_50[::-1],nad_pr1_V_100[::-1],nad_se_pr1_V_100[::-1]])
    
    col_names_pr1_V=['genome_position','norm_react_V_25','norm_react_se_V_25','norm_react_V_50','norm_react_se_V_50','norm_react_V_100','norm_react_se_V_100']
    
    full_df_pr1_V=pd.DataFrame(full_arr_pr1_V,columns=col_names_pr1_V)
    
    full_df_pr1_V.to_csv('MS2_pri1_virion_data_3.csv',sep=',',index=False)
    
    
    #write unnormalised data to files
    full_arr_pr1_R=np.transpose([nuc_pos1,ad_pr1_R_50[::-1],ad_se_pr1_R_50[::-1],ad_pr1_R_100[::-1],ad_se_pr1_R_100[::-1]])
    
    col_names_pr1_R=['genome_position','norm_react_R_50','norm_react_se_R_50','norm_react_R_100','norm_react_se_R_100']
    
    full_df_pr1_R=pd.DataFrame(full_arr_pr1_R,columns=col_names_pr1_R)
    
    full_df_pr1_R.to_csv('MS2_pri1_transcript_data_ad.csv',sep=',',index=False)
    
    full_arr_pr1_V=np.transpose([nuc_pos1,ad_pr1_V_25[::-1],ad_se_pr1_V_25[::-1],ad_pr1_V_50[::-1],ad_se_pr1_V_50[::-1],ad_pr1_V_100[::-1],ad_se_pr1_V_100[::-1]])
    
    col_names_pr1_V=['genome_position','norm_react_V_25','norm_react_se_V_25','norm_react_V_50','norm_react_se_V_50','norm_react_V_100','norm_react_se_V_100']
    
    full_df_pr1_V=pd.DataFrame(full_arr_pr1_V,columns=col_names_pr1_V)
    
    full_df_pr1_V.to_csv('MS2_pri1_virion_data_ad.csv',sep=',',index=False)
    
    ####write difference maps to file
    full_arr_pr1_D=np.transpose([nuc_pos1,nad_pr1_D_50[::-1],nad_se_pr1_D_50[::-1],nad_pr1_D_100[::-1],nad_se_pr1_D_100[::-1]])
    
    col_names_pr1_D=['genome_position','norm_react_D_50','norm_react_se_D_50','norm_react_D_100','norm_react_se_D_100']
    
    full_df_pr1_D=pd.DataFrame(full_arr_pr1_D,columns=col_names_pr1_D)
    
    full_df_pr1_D.to_csv('MS2_pri1_difference_data_3.csv',sep=',',index=False)

    ####produce normalised reactivity snapshots
    sam_funcs.sequence_snapshots(nuc_pos1,nad_pr1_R_25,nad_se_pr1_R_25,col='g',virus='MS2',primer='pri1',condition='R',treatment='25')
    
    sam_funcs.sequence_snapshots(nuc_pos1,nad_pr1_R_50,nad_se_pr1_R_50,col='g',virus='MS2',primer='pri1',condition='R',treatment='50')
    
    sam_funcs.sequence_snapshots(nuc_pos1,nad_pr1_R_100,nad_se_pr1_R_100,col='g',virus='MS2',primer='pri1',condition='R',treatment='100')
    
    sam_funcs.sequence_snapshots(nuc_pos1,nad_pr1_V_25,nad_se_pr1_V_25,col='b',virus='MS2',primer='pri1',condition='V',treatment='25')
    
    sam_funcs.sequence_snapshots(nuc_pos1,nad_pr1_V_50,nad_se_pr1_V_50,col='b',virus='MS2',primer='pri1',condition='V',treatment='50')

    
    sam_funcs.sequence_snapshots(nuc_pos1,nad_pr1_V_100,nad_se_pr1_V_100,col='b',virus='MS2',primer='pri1',condition='V',treatment='100')
   
   
    sam_funcs.sequence_snapshots(nuc_pos1,nad_pr1_D_50,nad_se_pr1_D_50,col='m',virus='MS2',primer='pri1',condition='V-T',treatment='50',diff=True)

    sam_funcs.sequence_snapshots(nuc_pos1,nad_pr1_D_100,nad_se_pr1_D_100,col='m',virus='MS2',primer='pri1',condition='V-T',treatment='100',diff=True)
