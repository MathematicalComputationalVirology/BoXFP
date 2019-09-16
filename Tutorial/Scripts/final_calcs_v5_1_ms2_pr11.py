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
    
    #sys.path.append(os.path.abspath('/Users/mqbppsc3/Desktop/externalLibraries')) 
    sys.path.append(os.path.abspath('/home/samclark/sclark/Cap_elec/software/externalLibraries'))
    import funcFile
    import funcPeakAlign
    import funcSeqAll
    import funcToolsAll
    import funcByRef
    import funcGeneral
    import funcTimeWarp
    import sam_funcs
    
    np.set_printoptions(threshold=sys.maxsize)
    #read in file list
    with open('data_files.fl' ,'r') as f:
        file_list =f.readlines()
    """
    #read in data initially to get the SM peak positions
    data_arr0 = sam_funcs.data_reader(file_list,6000,1400)
        
    n=len(data_arr0)

    #plot the size marker traces to ensure that all peaks have been captured
    color=iter(plt.cm.rainbow(np.linspace(0,1,n)))
    for i in range(n):
        col=next(color)
        plt.plot(data_arr0[i]['Position'],data_arr0[i]['SizeMarker'],c=col)
    plt.show()

    #run preprocessing
    data_arr = sam_funcs.preprocess(data_arr0)

    #run mobility shift
    data_arr1=sam_funcs.mobility_shift(data_arr)
    
    #find peaks in TM traces
    peka,peaksTM=sam_funcs.peak_finder_v2(data_arr1,4,.25,cap=7000)
    data_arr=[]
    
    #list all those peask that have a disproportionate number of SM peaks
    for i in range(len(peaksTM)):
        
        lenTM=len(peaksTM[i])
        
        if lenTM!=21:
            print file_list[i]+' '+str(i)
    
    #plot the SM traces with the peak positions marked to make sure the peak finder function has found the right peaks.
    sam_funcs.sm_plotter(data_arr1,peaksTM,file_list)
    
    #run the data reader version two that carries out the windowing and stores the windows in a pickle .obj file
    sam_funcs.data_reader_v2(file_list,peaksTM,'190902_11')
    """
    #initialise replicate correlations
    corAB_R_25=[]
    corAC_R_25=[]
    corBC_R_25=[]  
    
    corAB_R_50=[]
    corAC_R_50=[]
    corBC_R_50=[]    
    
    corAB_R_10=[]
    corAC_R_10=[]
    corBC_R_10=[]
    
    corAB_Re_25=[]
    corAC_Re_25=[]
    corBC_Re_25=[]  
    
    corAB_Re_50=[]
    corAC_Re_50=[]
    corBC_Re_50=[]    
    
    corAB_Re_10=[]
    corAC_Re_10=[]
    corBC_Re_10=[]
    
    corAB_V_25=[]
    corAC_V_25=[]
    corBC_V_25=[]  
    
    corAB_V_50=[]
    corAC_V_50=[]
    corBC_V_50=[]    
    
    corAB_V_10=[]
    corAC_V_10=[]
    corBC_V_10=[]
    
    corAB_G_25=[]
    corAC_G_25=[]
    corBC_G_25=[]  
    
    corAB_G_50=[]
    corAC_G_50=[]
    corBC_G_50=[]    
    
    corAB_G_10=[]
    corAC_G_10=[]
    corBC_G_10=[]
    
    avers_V_10=[]
    avers_V_25=[]
    avers_V_50=[]

    avers_G_10=[]
    avers_G_25=[]
    avers_G_50=[]
    
    avers_R_10=[]
    avers_R_25=[]
    avers_R_50=[]
    
    
    avers_Re_10=[]
    avers_Re_25=[]
    avers_Re_50=[]
    
    cor_GV_10=[]
    cor_GV_25=[]
    cor_GV_50=[]
    
    cor_GR_10=[]
    cor_GR_25=[]
    cor_GR_50=[]
    
    
    
    #iterate through windows
    for i in range(10):
        
        #print window number
        print 'window '+str(i)
        
        #open window data file
        file_path="190732_pp_"+str(i)+".obj"
        file_1= open(file_path,'rb')
        
        #load data
        data_arr2=pickle.load(file_1)
       
        #run partition for data with no replicates
        partition_RX,partition_RX2=sam_funcs.RX_partitioning_single(data_arr2,1,file_list)

        #run partition for data with replicates
        bin_alloc,partition_RX,peak_info= sam_funcs.RX_partitioning_replicates(data_arr2,1,0.25)


        #perform peak finding and plot peaks
        peka,peaksTM=sam_funcs.peak_finder_v2(data_arr2,4,.25,TM=1)

        sam_funcs.sm_plotter(data_arr2,peaksTM,file_list)

        #index data
        R_0=[27]
        R_10=[28,27,29]
        R_25=[32,33,34]
        R_50=[36,37,38]
        V_0=[13]
        V_10=[14,15]
        V_25=[16,17]
        V_50=[18,19]
        
        Re_0=[20]
        Re_10=[21,22]
        Re_25=[23,24]
        Re_50=[25,26]
        
        G_0=[0]
        G_10=[1,2,3]
        G_25=[5,6,7]
        G_50=[9,10,11]

        #run partition alignment for replicated data
        part_R_10=sam_funcs.RX_partition_realignment(partition_RX,bin_alloc,peak_info,R_10,data_arr2)

        part_R_25=sam_funcs.RX_partition_realignment(partition_RX,bin_alloc,peak_info,R_25,data_arr2)

        part_R_50=sam_funcs.RX_partition_realignment(partition_RX,bin_alloc,peak_info,R_50,data_arr2)
        
        
        part_Re_10=sam_funcs.RX_partition_realignment(partition_RX,bin_alloc,peak_info,Re_10,data_arr2)

        part_Re_25=sam_funcs.RX_partition_realignment(partition_RX,bin_alloc,peak_info,Re_25,data_arr2)

        part_Re_50=sam_funcs.RX_partition_realignment(partition_RX,bin_alloc,peak_info,Re_50,data_arr2)

        part_G_10=sam_funcs.RX_partition_realignment(partition_RX,bin_alloc,peak_info,G_10,data_arr2)

        part_G_25=sam_funcs.RX_partition_realignment(partition_RX,bin_alloc,peak_info,G_25,data_arr2)

        part_G_50=sam_funcs.RX_partition_realignment(partition_RX,bin_alloc,peak_info,G_50,data_arr2)

        part_V_10=sam_funcs.RX_partition_realignment(partition_RX,bin_alloc,peak_info,V_10,data_arr2)

        part_V_25=sam_funcs.RX_partition_realignment(partition_RX,bin_alloc,peak_info,V_25,data_arr2)

        part_V_50=sam_funcs.RX_partition_realignment(partition_RX,bin_alloc,peak_info,V_50,data_arr2)
        

        #get replicate correlations
        vbv,cor_mat_R_10=sam_funcs.correl_assessor(part_R_10,'amp')
        vbv,cor_mat_R_25=sam_funcs.correl_assessor(part_R_25,'amp')
        vbv,cor_mat_R_50=sam_funcs.correl_assessor(part_R_50,'amp')
        
        vbv,cor_mat_Re_10=sam_funcs.correl_assessor(part_Re_10,'amp')
        vbv,cor_mat_Re_25=sam_funcs.correl_assessor(part_Re_25,'amp')
        vbv,cor_mat_Re_50=sam_funcs.correl_assessor(part_Re_50,'amp')
        
        vbv,cor_mat_V_10=sam_funcs.correl_assessor(part_V_10,'amp')
        vbv,cor_mat_V_25=sam_funcs.correl_assessor(part_V_25,'amp')
        vbv,cor_mat_V_50=sam_funcs.correl_assessor(part_V_50,'amp')
        
        vbv,cor_mat_G_10=sam_funcs.correl_assessor(part_G_10,'amp')
        vbv,cor_mat_G_25=sam_funcs.correl_assessor(part_G_25,'amp')
        vbv,cor_mat_G_50=sam_funcs.correl_assessor(part_G_50,'amp')
        
        #bin replicate correlations
        corAB_R_10=np.append(corAB_R_10,cor_mat_R_10[0,1])
        corAC_R_10=np.append(corAC_R_10,cor_mat_R_10[0,2])
        corBC_R_10=np.append(corBC_R_10,cor_mat_R_10[2,1])
        
        corAB_R_25=np.append(corAB_R_25,cor_mat_R_25[0,1])
        corAC_R_25=np.append(corAC_R_25,cor_mat_R_25[0,2])
        corBC_R_25=np.append(corBC_R_25,cor_mat_R_25[2,1])
        
        corAB_R_50=np.append(corAB_R_50,cor_mat_R_50[0,1])
        corAC_R_50=np.append(corAC_R_50,cor_mat_R_50[0,2])
        corBC_R_50=np.append(corBC_R_50,cor_mat_R_50[2,1])
        
        corAB_Re_10=np.append(corAB_Re_10,cor_mat_Re_10[0,1])
        #corAC_Re_10=np.append(corAC_Re_10,cor_mat_Re_10[0,2])
        #corBC_Re_10=np.append(corBC_Re_10,cor_mat_Re_10[2,1])
        
        corAB_Re_25=np.append(corAB_Re_25,cor_mat_Re_25[0,1])
        #corAC_Re_25=np.append(corAC_Re_25,cor_mat_Re_25[0,2])
        #corBC_Re_25=np.append(corBC_Re_25,cor_mat_Re_25[2,1])
        
        corAB_Re_50=np.append(corAB_Re_50,cor_mat_Re_50[0,1])
        #corAC_Re_50=np.append(corAC_Re_50,cor_mat_Re_50[0,2])
        #corBC_Re_50=np.append(corBC_Re_50,cor_mat_Re_50[2,1])
        
        corAB_V_10=np.append(corAB_V_10,cor_mat_V_10[0,1])
        #corAC_V_10=np.append(corAC_V_10,cor_mat_V_10[0,2])
        #corBC_V_10=np.append(corBC_V_10,cor_mat_V_10[2,1])
        
        corAB_V_25=np.append(corAB_V_25,cor_mat_V_25[0,1])
        #corAC_V_25=np.append(corAC_V_25,cor_mat_V_25[0,2])
        #corBC_V_25=np.append(corBC_V_25,cor_mat_V_25[2,1])
        
        corAB_V_50=np.append(corAB_V_50,cor_mat_V_50[0,1])
        #corAC_V_50=np.append(corAC_V_50,cor_mat_V_50[0,2])
        #corBC_V_50=np.append(corBC_V_50,cor_mat_V_50[2,1])
        
        corAB_G_10=np.append(corAB_G_10,cor_mat_G_10[0,1])
        corAC_G_10=np.append(corAC_G_10,cor_mat_G_10[0,2])
        corBC_G_10=np.append(corBC_G_10,cor_mat_G_10[2,1])
        
        corAB_G_25=np.append(corAB_G_25,cor_mat_G_25[0,1])
        corAC_G_25=np.append(corAC_G_25,cor_mat_G_25[0,2])
        corBC_G_25=np.append(corBC_G_25,cor_mat_G_25[2,1])
        
        corAB_G_50=np.append(corAB_G_50,cor_mat_G_50[0,1])
        corAC_G_50=np.append(corAC_G_50,cor_mat_G_50[0,2])
        corBC_G_50=np.append(corBC_G_50,cor_mat_G_50[2,1])
        
        
        print'Re'
        #calculate areas of peaks in replicated data
        amp_av_Re_10,area_av_Re_10,area_sd_Re_10=sam_funcs.area_calcs(part_Re_10,data_arr2,Re_10)
        amp_av_Re_25,area_av_Re_25,area_sd_Re_25=sam_funcs.area_calcs(part_Re_25,data_arr2,Re_25)
        amp_av_Re_50,area_av_Re_50,area_sd_Re_50=sam_funcs.area_calcs(part_Re_50,data_arr2,Re_50)
        
        
        print'R'
        amp_av_R_10,area_av_R_10,area_sd_R_10=sam_funcs.area_calcs(part_R_10,data_arr2,R_10)
        amp_av_R_25,area_av_R_25,area_sd_R_25=sam_funcs.area_calcs(part_R_25,data_arr2,R_25)
        amp_av_R_50,area_av_R_50,area_sd_R_50=sam_funcs.area_calcs(part_R_50,data_arr2,R_50)

        print'G'
        amp_av_G_10,area_av_G_10,area_sd_G_10=sam_funcs.area_calcs(part_G_10,data_arr2,G_10)
        amp_av_G_25,area_av_G_25,area_sd_G_25=sam_funcs.area_calcs(part_G_25,data_arr2,G_25)
        amp_av_G_50,area_av_G_50,area_sd_G_50=sam_funcs.area_calcs(part_G_50,data_arr2,G_50)

        print'V' 
        amp_av_V_10,area_av_V_10,area_sd_V_10=sam_funcs.area_calcs(part_V_10,data_arr2,V_10)
        amp_av_V_25,area_av_V_25,area_sd_V_25=sam_funcs.area_calcs(part_V_25,data_arr2,V_25)
        amp_av_V_50,area_av_V_50,area_sd_V_50=sam_funcs.area_calcs(part_V_50,data_arr2,V_50)


        #calculate areas in non replicated data
        peak_list_R_0=sam_funcs.(partition_RX2[27],data_arr2,27)
        
        peak_list_Re_0=sam_funcs.reac_calculator_2(partition_RX2[20],data_arr2,27)

        peak_list_G_0=sam_funcs.reac_calculator_2(partition_RX2[0],data_arr2,0)

        peak_list_V_0=sam_funcs.reac_calculator_2(partition_RX2[13],data_arr2,13)

        #calculate scaling factor for background correction 
        sf_R_10=funcSeqAll.scaleShapeDataWindow(amp_av_R_10,peak_list_R_0['amp'],deg=40,rate=0.25,step=10,fit='linear')

        sf_R_25=funcSeqAll.scaleShapeDataWindow(amp_av_R_25,peak_list_R_0['amp'],deg=40,rate=0.25,step=10,fit='linear')

        sf_R_50=funcSeqAll.scaleShapeDataWindow(amp_av_R_50,peak_list_R_0['amp'],deg=40,rate=0.25,step=10,fit='linear')
        
        sf_Re_10=funcSeqAll.scaleShapeDataWindow(amp_av_Re_10,peak_list_Re_0['amp'],deg=40,rate=0.25,step=10,fit='linear')

        sf_Re_25=funcSeqAll.scaleShapeDataWindow(amp_av_Re_25,peak_list_Re_0['amp'],deg=40,rate=0.25,step=10,fit='linear')

        sf_Re_50=funcSeqAll.scaleShapeDataWindow(amp_av_Re_50,peak_list_Re_0['amp'],deg=40,rate=0.25,step=10,fit='linear')

        sf_G_10=funcSeqAll.scaleShapeDataWindow(amp_av_G_10,peak_list_G_0['amp'],deg=40,rate=0.25,step=10,fit='linear')

        sf_G_25=funcSeqAll.scaleShapeDataWindow(amp_av_G_25,peak_list_G_0['amp'],deg=40,rate=0.25,step=10,fit='linear')

        sf_G_50=funcSeqAll.scaleShapeDataWindow(amp_av_G_50,peak_list_G_0['amp'],deg=40,rate=0.25,step=10,fit='linear')

        sf_V_10=funcSeqAll.scaleShapeDataWindow(amp_av_V_10,peak_list_V_0['amp'],deg=40,rate=0.25,step=10,fit='linear')

        sf_V_25=funcSeqAll.scaleShapeDataWindow(amp_av_V_25,peak_list_V_0['amp'],deg=40,rate=0.25,step=10,fit='linear')

        sf_V_50=funcSeqAll.scaleShapeDataWindow(amp_av_V_50,peak_list_V_0['amp'],deg=40,rate=0.25,step=10,fit='linear')

        
        #calculate the corrected area differences, normalised area differences
        ad_R_10,nad_R_10,aver_R_10=sam_funcs.area_differences(area_av_R_10,peak_list_R_0['area'],sf_R_10)

        ad_R_25,nad_R_25,aver_R_25=sam_funcs.area_differences(area_av_R_25,peak_list_R_0['area'],sf_R_25)
        ad_R_50,nad_R_50,aver_R_50=sam_funcs.area_differences(area_av_R_50,peak_list_R_0['area'],sf_R_50)
        
        ad_Re_10,nad_Re_10,aver_Re_10=sam_funcs.area_differences(area_av_Re_10,peak_list_Re_0['area'],sf_Re_10)

        ad_Re_25,nad_Re_25,aver_Re_25=sam_funcs.area_differences(area_av_Re_25,peak_list_Re_0['area'],sf_Re_25)
        ad_Re_50,nad_Re_50,aver_Re_50=sam_funcs.area_differences(area_av_Re_50,peak_list_Re_0['area'],sf_Re_50)

        ad_G_10,nad_G_10,aver_G_10=sam_funcs.area_differences(area_av_G_10,peak_list_G_0['area'],sf_G_10)

        ad_G_25,nad_G_25,aver_G_25=sam_funcs.area_differences(area_av_G_25,peak_list_G_0['area'],sf_G_25)
        ad_G_50,nad_G_50,aver_G_50=sam_funcs.area_differences(area_av_G_50,peak_list_G_0['area'],sf_G_50)


        ad_V_10,nad_V_10,aver_V_10=sam_funcs.area_differences(area_av_V_10,peak_list_V_0['area'],sf_V_10)
        ad_V_25,nad_V_25,aver_V_25=sam_funcs.area_differences(area_av_V_25,peak_list_G_0['area'],sf_G_25)
        ad_V_50,nad_V_50,aver_V_50=sam_funcs.area_differences(area_av_V_50,peak_list_V_0['area'],sf_V_50)
        
        
        #bin normalisation factors
        avers_G_10=np.append(avers_G_10,aver_G_10)
        avers_G_25=np.append(avers_G_25,aver_G_25)
        avers_G_50=np.append(avers_G_50,aver_G_50)

        avers_V_10=np.append(avers_V_10,aver_V_10)
        avers_V_25=np.append(avers_V_25,aver_V_25)
        avers_V_50=np.append(avers_V_50,aver_V_50)
        
        avers_R_10=np.append(avers_R_10,aver_R_10)
        avers_R_25=np.append(avers_R_25,aver_R_25)
        avers_R_50=np.append(avers_R_50,aver_R_50)
        
        avers_Re_10=np.append(avers_Re_10,aver_Re_10)
        avers_Re_25=np.append(avers_Re_25,aver_Re_25)
        avers_Re_50=np.append(avers_Re_50,aver_Re_50)

        #print normalisation factors
        print 'norm_R_10 '+str(aver_R_10)
        print 'norm_R_25 '+str(aver_R_25)
        print 'norm_R_50 '+str(aver_R_50)
        
        print 'norm_Re_10 '+str(aver_Re_10)
        print 'norm_Re_25 '+str(aver_Re_25)
        print 'norm_Re_50 '+str(aver_Re_50)

        print 'norm_G_10 '+str(aver_G_10)
        print 'norm_G_25 '+str(aver_G_25)
        print 'norm_G_50 '+str(aver_G_50)

        print 'norm_V_10 '+str(aver_V_10)
        print 'norm_V_25 '+str(aver_V_25)
        print 'norm_V_50 '+str(aver_V_50)

        #calculate normalised standard errors
        nad_se_R_10=area_sd_R_10/(aver_R_10*3)

        nad_se_R_25=area_sd_R_25/(aver_R_25*3)

        nad_se_R_50=area_sd_R_50/(aver_R_50*3)
        
        nad_se_Re_10=area_sd_Re_10/(aver_Re_10*3)

        nad_se_Re_25=area_sd_Re_25/(aver_Re_25*3)

        nad_se_Re_50=area_sd_Re_50/(aver_Re_50*3)

        nad_se_G_10=area_sd_G_10/(aver_G_10*3)

        nad_se_G_25=area_sd_G_25/(aver_G_25*3)

        nad_se_G_50=area_sd_G_50/(aver_G_50*3)

        nad_se_V_10=area_sd_V_10/(aver_V_10*2)

        nad_se_V_25=area_sd_V_25/(aver_V_25*2)

        nad_se_V_50=area_sd_V_50/(aver_V_50*2)
        
        #calculate non normalised standard errors
        ad_se_R_10=area_sd_R_10/(3)

        ad_se_R_25=area_sd_R_25/(3)

        ad_se_R_50=area_sd_R_50/(3)
        
        ad_se_Re_10=area_sd_Re_10/(3)

        ad_se_Re_25=area_sd_Re_25/(3)

        ad_se_Re_50=area_sd_Re_50/(3)

        ad_se_G_10=area_sd_G_10/(3)

        ad_se_G_25=area_sd_G_25/(3)

        ad_se_G_50=area_sd_G_50/(3)

        ad_se_V_10=area_sd_V_10/(2)

        ad_se_V_25=area_sd_V_25/(2)

        ad_se_V_50=area_sd_V_50/(2)

        ####calcualting difference maps
        
        #merge data
        diffs_10_3=np.concatenate((ad_V_10,ad_R_10,ad_G_10),axis=None)

        #find percentage outliers and average
        PO_10_3,PA_10_3=funcSeqAll.findPOutlierBox(diffs_10_3)
        
        diffs_50_3=np.concatenate((ad_V_50,ad_R_50,ad_G_50),axis=None)

        PO_50_3,PA_50_3=funcSeqAll.findPOutlierBox(diffs_50_3)

        diffs_25_3=np.concatenate((ad_V_25,ad_R_25,ad_G_25),axis=None)

        PO_25_3,PA_25_3=funcSeqAll.findPOutlierBox(diffs_25_3)
        
        #renormalise data
        nad_V_50_3,aver_V_50_3_2=funcSeqAll.normSimple(ad_V_50,PO_50_3,PA_50_3)
        nad_R_50_3,aver_R_50_3_2=funcSeqAll.normSimple(ad_R_50,PO_50_3,PA_50_3)
        nad_G_50_3,aver_G_50_3_2=funcSeqAll.normSimple(ad_G_50,PO_50_3,PA_50_3)
       
        nad_V_10_3,aver_V_10_3_2=funcSeqAll.normSimple(ad_V_10,PO_10_3,PA_10_3)
        nad_R_10_3,aver_R_10_3_2=funcSeqAll.normSimple(ad_R_10,PO_10_3,PA_10_3)
        nad_G_10_3,aver_G_10_3_2=funcSeqAll.normSimple(ad_G_10,PO_10_3,PA_10_3)

        nad_V_25_3,aver_V_25_3_2=funcSeqAll.normSimple(ad_V_25,PO_25_3,PA_25_3)
        nad_R_25_3,aver_R_25_3_2=funcSeqAll.normSimple(ad_R_25,PO_25_3,PA_25_3)
        nad_G_25_3,aver_G_25_3_2=funcSeqAll.normSimple(ad_G_25,PO_25_3,PA_25_3)
        
        cor_GV_10=np.append(cor_GV_10,pearsonr(nad_G_10_3,nad_V_10_3)[0])
        cor_GV_25=np.append(cor_GV_25,pearsonr(nad_G_25_3,nad_V_25_3)[0])
        cor_GV_50=np.append(cor_GV_50,pearsonr(nad_G_50_3,nad_V_50_3)[0])
        
        cor_GR_10=np.append(cor_GR_10,pearsonr(nad_G_10_3,nad_R_10_3)[0])
        cor_GR_25=np.append(cor_GR_25,pearsonr(nad_G_25_3,nad_R_25_3)[0])
        cor_GR_50=np.append(cor_GR_50,pearsonr(nad_G_50_3,nad_R_50_3)[0])

        diffs_10=np.append(ad_V_10,ad_R_10)
        diffs_25=np.append(ad_V_25,ad_R_25)
        diffs_50=np.append(ad_V_50,ad_R_50)

        PO_10,PA_10=funcSeqAll.findPOutlierBox(diffs_10)
        PO_25,PA_25=funcSeqAll.findPOutlierBox(diffs_25)
        PO_50,PA_50=funcSeqAll.findPOutlierBox(diffs_50)


        nad_V_10_2,aver_V_10_2=funcSeqAll.normSimple(ad_V_10,PO_10,PA_10)

        nad_V_25_2,aver_V_25_2=funcSeqAll.normSimple(ad_V_25,PO_25,PA_25)
        nad_V_50_2,aver_V_50_2=funcSeqAll.normSimple(ad_V_50,PO_50,PA_50)

        nad_R_10_2,aver_R_10_2=funcSeqAll.normSimple(ad_R_10,PO_10,PA_10)
        nad_R_25_2,aver_R_25_2=funcSeqAll.normSimple(ad_R_25,PO_25,PA_25)
        nad_R_50_2,aver_R_50_2=funcSeqAll.normSimple(ad_R_50,PO_50,PA_50)

        nad_D_10=nad_V_10_2-nad_R_10_2
        nad_D_25=nad_V_25_2-nad_R_25_2
        nad_D_50=nad_V_50_2-nad_R_50_2


        #calculate renormalised standard errors
        nad_se_R_10_2=area_sd_R_10/(aver_R_10_2*3)
        nad_se_V_10_2=area_sd_V_10/(aver_V_10_2*2)

        nad_se_R_25_2=area_sd_R_25/(aver_R_25_2*3)
        nad_se_V_25_2=area_sd_V_25/(aver_V_25_2*2)

        nad_se_R_50_2=area_sd_R_50/(aver_R_50_2*3)
        nad_se_V_50_2=area_sd_V_50/(aver_V_50_2*2)


        nad_se_D_10=sam_funcs.error_propagation(nad_se_R_10_2,nad_se_V_10_2)
        nad_se_D_25=sam_funcs.error_propagation(nad_se_R_25_2,nad_se_V_25_2)
        nad_se_D_50=sam_funcs.error_propagation(nad_se_R_50_2,nad_se_V_50_2)
    
        #store data
        if i==0:
            
            nad_R_10_arr=deepcopy(nad_R_10)
            nad_se_R_10_arr=deepcopy(nad_se_R_10)
            nad_R_25_arr=deepcopy(nad_R_25)
            nad_se_R_25_arr=deepcopy(nad_se_R_25)
            nad_R_50_arr=deepcopy(nad_R_50)
            nad_se_R_50_arr=deepcopy(nad_se_R_50)
            
            nad_Re_10_arr=deepcopy(nad_Re_10)
            nad_se_Re_10_arr=deepcopy(nad_se_Re_10)
            nad_Re_25_arr=deepcopy(nad_Re_25)
            nad_se_Re_25_arr=deepcopy(nad_se_Re_25)
            nad_Re_50_arr=deepcopy(nad_Re_50)
            nad_se_Re_50_arr=deepcopy(nad_se_Re_50)
            
            nad_V_10_arr=deepcopy(nad_V_10)
            nad_se_V_10_arr=deepcopy(nad_se_V_10)
            nad_V_25_arr=deepcopy(nad_V_25)
            nad_se_V_25_arr=deepcopy(nad_se_V_25)
            nad_V_50_arr=deepcopy(nad_V_50)
            nad_se_V_50_arr=deepcopy(nad_se_V_50)
            
            ad_R_10_arr=deepcopy(ad_R_10)
            ad_se_R_10_arr=deepcopy(ad_se_R_10)
            ad_R_25_arr=deepcopy(ad_R_25)
            ad_se_R_25_arr=deepcopy(ad_se_R_25)
            ad_R_50_arr=deepcopy(ad_R_50)
            ad_se_R_50_arr=deepcopy(ad_se_R_50)
            
            ad_Re_10_arr=deepcopy(ad_Re_10)
            ad_se_Re_10_arr=deepcopy(ad_se_Re_10)
            ad_Re_25_arr=deepcopy(ad_Re_25)
            ad_se_Re_25_arr=deepcopy(ad_se_Re_25)
            ad_Re_50_arr=deepcopy(ad_Re_50)
            ad_se_Re_50_arr=deepcopy(ad_se_Re_50)
            
            ad_V_10_arr=deepcopy(ad_V_10)
            ad_se_V_10_arr=deepcopy(ad_se_V_10)
            ad_V_25_arr=deepcopy(ad_V_25)
            ad_se_V_25_arr=deepcopy(ad_se_V_25)
            ad_V_50_arr=deepcopy(ad_V_50)
            ad_se_V_50_arr=deepcopy(ad_se_V_50)
            
            nad_D_10_arr=deepcopy(nad_D_10)
            nad_se_D_10_arr=deepcopy(nad_se_D_10)
            nad_D_25_arr=deepcopy(nad_D_25)
            nad_se_D_25_arr=deepcopy(nad_se_D_25)
            nad_D_50_arr=deepcopy(nad_D_50)
            nad_se_D_50_arr=deepcopy(nad_se_D_50)
            
            nad_G_10_arr=deepcopy(nad_G_10)
            nad_se_G_10_arr=deepcopy(nad_se_G_10)
            nad_G_25_arr=deepcopy(nad_G_25)
            nad_se_G_25_arr=deepcopy(nad_se_G_25)
            nad_G_50_arr=deepcopy(nad_G_50)
            nad_se_G_50_arr=deepcopy(nad_se_G_50)
            
            ad_G_10_arr=deepcopy(ad_G_10)
            ad_se_G_10_arr=deepcopy(ad_se_G_10)
            ad_G_25_arr=deepcopy(ad_G_25)
            ad_se_G_25_arr=deepcopy(ad_se_G_25)
            ad_G_50_arr=deepcopy(ad_G_50)
            ad_se_G_50_arr=deepcopy(ad_se_G_50)
            
            
        else:
            nad_R_10_arr=np.vstack((nad_R_10_arr,nad_R_10))
            nad_se_R_10_arr=np.vstack((nad_se_R_10_arr,nad_se_R_10))
            
            nad_R_25_arr=np.vstack((nad_R_25_arr,nad_R_25))
            nad_se_R_25_arr=np.vstack((nad_se_R_25_arr,nad_se_R_25))
            
            nad_R_50_arr=np.vstack((nad_R_50_arr,nad_R_50))
            nad_se_R_50_arr=np.vstack((nad_se_R_50_arr,nad_se_R_50))
            
            
            nad_Re_10_arr=np.vstack((nad_Re_10_arr,nad_Re_10))
            nad_se_Re_10_arr=np.vstack((nad_se_Re_10_arr,nad_se_Re_10))
            
            nad_Re_25_arr=np.vstack((nad_Re_25_arr,nad_Re_25))
            nad_se_Re_25_arr=np.vstack((nad_se_Re_25_arr,nad_se_Re_25))
            
            nad_Re_50_arr=np.vstack((nad_Re_50_arr,nad_Re_50))
            nad_se_Re_50_arr=np.vstack((nad_se_Re_50_arr,nad_se_Re_50))
            
            nad_V_10_arr=np.vstack((nad_V_10_arr,nad_V_10))
            nad_se_V_10_arr=np.vstack((nad_se_V_10_arr,nad_se_V_10))
            
            nad_V_25_arr=np.vstack((nad_V_25_arr,nad_V_25))
            nad_se_V_25_arr=np.vstack((nad_se_V_25_arr,nad_se_V_25))
            
            nad_V_50_arr=np.vstack((nad_V_50_arr,nad_V_50))
            nad_se_V_50_arr=np.vstack((nad_se_V_50_arr,nad_se_V_50))
            
            ad_R_10_arr=np.vstack((ad_R_10_arr,ad_R_10))
            ad_se_R_10_arr=np.vstack((ad_se_R_10_arr,ad_se_R_10))
            
            ad_R_25_arr=np.vstack((ad_R_25_arr,ad_R_25))
            ad_se_R_25_arr=np.vstack((ad_se_R_25_arr,ad_se_R_25))
            
            ad_R_50_arr=np.vstack((ad_R_50_arr,ad_R_50))
            ad_se_R_50_arr=np.vstack((ad_se_R_50_arr,ad_se_R_50))
            
            ad_Re_10_arr=np.vstack((ad_Re_10_arr,ad_Re_10))
            ad_se_Re_10_arr=np.vstack((ad_se_Re_10_arr,ad_se_Re_10))
            
            ad_Re_25_arr=np.vstack((ad_Re_25_arr,ad_Re_25))
            ad_se_Re_25_arr=np.vstack((ad_se_Re_25_arr,ad_se_Re_25))
            
            ad_Re_50_arr=np.vstack((ad_Re_50_arr,ad_Re_50))
            ad_se_Re_50_arr=np.vstack((ad_se_Re_50_arr,ad_se_Re_50))
            
            ad_V_10_arr=np.vstack((ad_V_10_arr,ad_V_10))
            ad_se_V_10_arr=np.vstack((ad_se_V_10_arr,ad_se_V_10))
            
            ad_V_25_arr=np.vstack((ad_V_25_arr,ad_V_25))
            ad_se_V_25_arr=np.vstack((ad_se_V_25_arr,ad_se_V_25))
            
            ad_V_50_arr=np.vstack((ad_V_50_arr,ad_V_50))
            ad_se_V_50_arr=np.vstack((ad_se_V_50_arr,ad_se_V_50))
            
            nad_D_10_arr=np.vstack((nad_D_10_arr,nad_D_10))
            nad_se_D_10_arr=np.vstack((nad_se_D_10_arr,nad_se_D_10))
            
            nad_D_25_arr=np.vstack((nad_D_25_arr,nad_D_25))
            nad_se_D_25_arr=np.vstack((nad_se_D_25_arr,nad_se_D_25))
            
            nad_D_50_arr=np.vstack((nad_D_50_arr,nad_D_50))
            nad_se_D_50_arr=np.vstack((nad_se_D_50_arr,nad_se_D_50))
            
            
            nad_G_10_arr=np.vstack((nad_G_10_arr,nad_G_10))
            nad_se_G_10_arr=np.vstack((nad_se_G_10_arr,nad_se_G_10))
            
            nad_G_25_arr=np.vstack((nad_G_25_arr,nad_G_25))
            nad_se_G_25_arr=np.vstack((nad_se_G_25_arr,nad_se_G_25))
            
            nad_G_50_arr=np.vstack((nad_G_50_arr,nad_G_50))
            nad_se_G_50_arr=np.vstack((nad_se_G_50_arr,nad_se_G_50))
            
            
            ad_G_10_arr=np.vstack((ad_G_10_arr,ad_G_10))
            ad_se_G_10_arr=np.vstack((ad_se_G_10_arr,ad_se_G_10))
            
            ad_G_25_arr=np.vstack((ad_G_25_arr,ad_G_25))
            ad_se_G_25_arr=np.vstack((ad_se_G_25_arr,ad_se_G_25))
            
            ad_G_50_arr=np.vstack((ad_G_50_arr,ad_G_50))
            ad_se_G_50_arr=np.vstack((ad_se_G_50_arr,ad_se_G_50))
   
    
    

    
    #create empty correlation matrices for windowing correlations
    ca_D_25=np.zeros([len(nad_D_50_arr),len(nad_D_50_arr)])
    ca_D_50=np.zeros([len(nad_D_50_arr),len(nad_D_50_arr)])
    ca_D_10=np.zeros([len(nad_D_50_arr),len(nad_D_50_arr)])
    ca_R_25=np.zeros([len(nad_D_50_arr),len(nad_D_50_arr)])
    ca_R_50=np.zeros([len(nad_D_50_arr),len(nad_D_50_arr)])
    ca_R_10=np.zeros([len(nad_D_50_arr),len(nad_D_50_arr)])
    ca_V_25=np.zeros([len(nad_D_50_arr),len(nad_D_50_arr)])
    ca_V_50=np.zeros([len(nad_D_50_arr),len(nad_D_50_arr)])
    ca_V_10=np.zeros([len(nad_D_50_arr),len(nad_D_50_arr)])
    ca_G_25=np.zeros([len(nad_D_50_arr),len(nad_D_50_arr)])
    ca_G_50=np.zeros([len(nad_D_50_arr),len(nad_D_50_arr)])
    ca_G_10=np.zeros([len(nad_D_10_arr),len(nad_D_10_arr)])

    print 'GR'
    
    print cor_GR_10
    print cor_GR_25
    print cor_GR_50
    print 'GV'
    print cor_GV_10
    print cor_GV_25
    print cor_GV_50
    
    #calcuate correlations between data from windows
    for i in range(len(nad_D_50_arr)):
        for j in range(len(nad_D_50_arr)):
            ca_R_25[i,j]=pearsonr(nad_R_25_arr[i,:],nad_R_25_arr[j,:])[0]
            ca_R_50[i,j]=pearsonr(nad_R_50_arr[i,:],nad_R_50_arr[j,:])[0]
            ca_R_10[i,j]=pearsonr(nad_R_10_arr[i,:],nad_R_10_arr[j,:])[0]

            
            ca_D_25[i,j]=pearsonr(nad_D_25_arr[i,:],nad_D_25_arr[j,:])[0]
            ca_D_50[i,j]=pearsonr(nad_D_50_arr[i,:],nad_D_50_arr[j,:])[0]
            ca_D_10[i,j]=pearsonr(nad_D_10_arr[i,:],nad_D_10_arr[j,:])[0]

            
            ca_V_25[i,j]=pearsonr(nad_V_25_arr[i,:],nad_V_25_arr[j,:])[0]
            ca_V_50[i,j]=pearsonr(nad_V_50_arr[i,:],nad_V_50_arr[j,:])[0]
            ca_V_10[i,j]=pearsonr(nad_V_10_arr[i,:],nad_V_10_arr[j,:])[0]
            
            ca_G_25[i,j]=pearsonr(nad_G_25_arr[i,:],nad_G_25_arr[j,:])[0]
            ca_G_50[i,j]=pearsonr(nad_G_50_arr[i,:],nad_G_50_arr[j,:])[0]
            ca_G_10[i,j]=pearsonr(nad_G_10_arr[i,:],nad_G_10_arr[j,:])[0]
    
    #calculate the mean correlation for each window
    mean_ca_R_25=np.mean(ca_R_25,axis=1)
    mean_ca_R_50=np.mean(ca_R_50,axis=1)
    mean_ca_R_10=np.mean(ca_R_10,axis=1)

    mean_ca_V_25=np.mean(ca_V_25,axis=1)
    mean_ca_V_50=np.mean(ca_V_50,axis=1)
    mean_ca_V_10=np.mean(ca_V_10,axis=1)
    
    mean_ca_G_25=np.mean(ca_G_25,axis=1)
    mean_ca_G_50=np.mean(ca_G_50,axis=1)
    mean_ca_G_10=np.mean(ca_G_10,axis=1)
    
    mean_ca_D_25=np.mean(ca_D_25,axis=1)
    mean_ca_D_50=np.mean(ca_D_50,axis=1)
    mean_ca_D_10=np.mean(ca_D_10,axis=1)
    
    
    #keep windows above a certain average correlation
    keep_ca_D_25=np.where(mean_ca_D_25>0.8)[0]
    keep_ca_D_50=np.where(mean_ca_D_50>0.8)[0]
    keep_ca_D_10=np.where(mean_ca_D_10>0.8)[0]
    
    print 'means'
    print mean_ca_D_10
    print mean_ca_D_25
    print mean_ca_D_50
    
    print np.count_nonzero(mean_ca_R_25>0.8)
    print np.count_nonzero(mean_ca_R_50>0.8)
    print np.count_nonzero(mean_ca_R_10>0.8)
    
    print np.count_nonzero(mean_ca_V_25>0.8)
    print np.count_nonzero(mean_ca_V_50>0.8)
    print np.count_nonzero(mean_ca_V_10>0.8)
    
    print np.count_nonzero(mean_ca_G_25>0.8)
    print np.count_nonzero(mean_ca_G_50>0.8)
    print np.count_nonzero(mean_ca_G_10>0.8)
    
    print np.count_nonzero(mean_ca_D_25>0.8)
    print np.count_nonzero(mean_ca_D_50>0.8)
    print np.count_nonzero(mean_ca_D_10>0.8)
    
    

    ####write normalisation factors
    avers_arr_G=np.transpose([avers_G_10,avers_G_25,avers_G_50])
    
    avers_labs_G=['aver_G_10','aver_G_25','aver_G_50']
    
    avers_df_G=pd.DataFrame(avers_arr_G,columns=avers_labs_G)
    
    avers_df_G.to_csv('MS2_pri11_avers_G.csv',sep=',',index=False)
    
    
    avers_arr_V=np.transpose([avers_V_10,avers_V_25,avers_V_50])
    
    avers_labs_V=['aver_V_10','aver_V_25','avers_V_50']
    
    avers_df_V=pd.DataFrame(avers_arr_V,columns=avers_labs_V)
    
    avers_df_V.to_csv('MS2_pri11_avers_V.csv',sep=',',index=False)
    
    avers_arr_R=np.transpose([avers_R_10,avers_R_25,avers_R_50])
    
    avers_labs_R=['aver_R_10','aver_R_25','aver_R_50']
    
    avers_df_R=pd.DataFrame(avers_arr_R,columns=avers_labs_R)
    
    avers_df_R.to_csv('MS2_pri11_avers_R.csv',sep=',',index=False)
    
    
    avers_arr_Re=np.transpose([avers_Re_10,avers_Re_25,avers_Re_50])
    
    avers_labs_Re=['aver_Re_10','aver_Re_25','aver_Re_50']
    
    avers_df_Re=pd.DataFrame(avers_arr_Re,columns=avers_labs_Re)
    
    avers_df_Re.to_csv('MS2_pri11_avers_Re.csv',sep=',',index=False)
   
    ####write replicate correlations to file
    cor_labs=['corAB_10','corBC_10','corAC_10','corAB_25','corBC_25','corAC_25','corAB_50','corBC_50','corAC_50']
    
    cor_labs_V=['corAB_10','corAB_25','corAB_50']
    cor_arr_R=np.transpose([corAB_R_10,corBC_R_10,corAC_R_10,corAB_R_25,corBC_R_25,corAC_R_25,corAB_R_50,corBC_R_50,corAC_R_50])
    
    cor_arr_Re=np.transpose([corAB_Re_10,corAB_Re_25,corAB_Re_50])
    
    cor_arr_V=np.transpose([corAB_V_10,corAB_V_25,corAB_V_50])
    
    cor_arr_G=np.transpose([corAB_G_10,corBC_G_10,corAC_G_10,corAB_G_25,corBC_G_25,corAC_G_25,corAB_G_50,corBC_G_50,corAC_G_50])
    
    cor_df_R=pd.DataFrame(cor_arr_R,columns=cor_labs)
    
    cor_df_Re=pd.DataFrame(cor_arr_Re,columns=cor_labs_V)
    
    cor_df_V=pd.DataFrame(cor_arr_V,columns=cor_labs_V)
    
    cor_df_G=pd.DataFrame(cor_arr_G,columns=cor_labs)
    
    
    cor_df_R.to_csv('MS2_pri11_transcript_cor.csv',sep=',',index=False)
    
    cor_df_Re.to_csv('MS2_pri11_Reassembly_cor.csv',sep=',',index=False)
    
    cor_df_V.to_csv('MS2_pri11_virion_cor.csv',sep=',',index=False)
    
    cor_df_G.to_csv('MS2_pri11_genome_cor.csv',sep=',',index=False)
    
    ####write windowing correlations
    np.savetxt('cor_mat_D_25.csv',ca_D_25)
    np.savetxt('cor_mat_D_50.csv',ca_D_50)
    np.savetxt('cor_mat_D_10.csv',ca_D_10)
    
    mean_cor_arr=np.transpose([mean_ca_D_10,mean_ca_D_25,mean_ca_D_50])
    
    mean_cor_labs=['cor_D_10','cor_D_25','cov_D_50']
   
    mean_cor_df=pd.DataFrame(mean_cor_arr,columns=mean_cor_labs)
    
    mean_cor_df.to_csv('MS2_pri11_mean_cors.csv',sep=',',index=False)
    
    #calculate means of window profiles kept
    nad_R_25=np.mean(nad_R_25_arr[keep_ca_D_25,:],axis=0)
    nad_R_50=np.mean(nad_R_50_arr[keep_ca_D_50,:],axis=0)
    nad_R_10=np.mean(nad_R_10_arr[keep_ca_D_10,:],axis=0)
    
    nad_se_R_25=np.mean(nad_se_R_25_arr[keep_ca_D_25,:],axis=0)
    nad_se_R_50=np.mean(nad_se_R_50_arr[keep_ca_D_50,:],axis=0)
    nad_se_R_10=np.mean(nad_se_R_10_arr[keep_ca_D_10,:],axis=0)
    
    nad_Re_25=np.mean(nad_Re_25_arr[keep_ca_D_25,:],axis=0)
    nad_Re_50=np.mean(nad_Re_50_arr[keep_ca_D_50,:],axis=0)
    nad_Re_10=np.mean(nad_Re_10_arr[keep_ca_D_10,:],axis=0)
    
    nad_se_Re_25=np.mean(nad_se_Re_25_arr[keep_ca_D_25,:],axis=0)
    nad_se_Re_50=np.mean(nad_se_Re_50_arr[keep_ca_D_50,:],axis=0)
    nad_se_Re_10=np.mean(nad_se_Re_10_arr[keep_ca_D_10,:],axis=0)
    
    nad_V_25=np.mean(nad_V_25_arr[keep_ca_D_25,:],axis=0)
    nad_V_50=np.mean(nad_V_50_arr[keep_ca_D_50,:],axis=0)
    nad_V_10=np.mean(nad_V_10_arr[keep_ca_D_10,:],axis=0)
    
    nad_se_V_25=np.mean(nad_se_V_25_arr[keep_ca_D_25,:],axis=0)
    nad_se_V_50=np.mean(nad_se_V_50_arr[keep_ca_D_50,:],axis=0)
    nad_se_V_10=np.mean(nad_se_V_10_arr[keep_ca_D_10,:],axis=0)
    
    nad_G_25=np.mean(nad_G_25_arr[keep_ca_D_25,:],axis=0)
    nad_G_50=np.mean(nad_G_50_arr[keep_ca_D_50,:],axis=0)
    nad_G_10=np.mean(nad_G_10_arr[keep_ca_D_10,:],axis=0)
    
    nad_se_G_25=np.mean(nad_se_G_25_arr[keep_ca_D_25,:],axis=0)
    nad_se_G_50=np.mean(nad_se_G_50_arr[keep_ca_D_50,:],axis=0)
    nad_se_G_10=np.mean(nad_se_G_10_arr[keep_ca_D_10,:],axis=0)
    
    ad_R_25=np.mean(ad_R_25_arr[keep_ca_D_25,:],axis=0)
    ad_R_50=np.mean(ad_R_50_arr[keep_ca_D_50,:],axis=0)
    ad_R_10=np.mean(ad_R_10_arr[keep_ca_D_10,:],axis=0)
    
    ad_se_R_25=np.mean(ad_se_R_25_arr[keep_ca_D_25,:],axis=0)
    ad_se_R_50=np.mean(ad_se_R_50_arr[keep_ca_D_50,:],axis=0)
    ad_se_R_10=np.mean(ad_se_R_10_arr[keep_ca_D_10,:],axis=0)
    
    ad_Re_25=np.mean(ad_Re_25_arr[keep_ca_D_25,:],axis=0)
    ad_Re_50=np.mean(ad_Re_50_arr[keep_ca_D_50,:],axis=0)
    ad_Re_10=np.mean(ad_Re_10_arr[keep_ca_D_10,:],axis=0)
    
    ad_se_Re_25=np.mean(ad_se_Re_25_arr[keep_ca_D_25,:],axis=0)
    ad_se_Re_50=np.mean(ad_se_Re_50_arr[keep_ca_D_50,:],axis=0)
    ad_se_Re_10=np.mean(ad_se_Re_10_arr[keep_ca_D_10,:],axis=0)
    
    ad_V_25=np.mean(ad_V_25_arr[keep_ca_D_25,:],axis=0)
    ad_V_50=np.mean(ad_V_50_arr[keep_ca_D_50,:],axis=0)
    ad_V_10=np.mean(ad_V_10_arr[keep_ca_D_10,:],axis=0)
    
    ad_se_V_25=np.mean(ad_se_V_25_arr[keep_ca_D_25,:],axis=0)
    ad_se_V_50=np.mean(ad_se_V_50_arr[keep_ca_D_50,:],axis=0)
    ad_se_V_10=np.mean(ad_se_V_10_arr[keep_ca_D_10,:],axis=0)
    
    ad_G_25=np.mean(ad_G_25_arr[keep_ca_D_25,:],axis=0)
    ad_G_50=np.mean(ad_G_50_arr[keep_ca_D_50,:],axis=0)
    ad_G_10=np.mean(ad_G_10_arr[keep_ca_D_10,:],axis=0)
    
    ad_se_G_25=np.mean(ad_se_G_25_arr[keep_ca_D_25,:],axis=0)
    ad_se_G_50=np.mean(ad_se_G_50_arr[keep_ca_D_50,:],axis=0)
    ad_se_G_10=np.mean(ad_se_G_10_arr[keep_ca_D_10,:],axis=0)
    
    nad_D_25=np.mean(nad_D_25_arr[keep_ca_D_25,:],axis=0)
    nad_D_50=np.mean(nad_D_50_arr[keep_ca_D_50,:],axis=0)
    nad_D_10=np.mean(nad_D_10_arr[keep_ca_D_10,:],axis=0)
    
    nad_se_D_25=np.mean(nad_se_D_25_arr[keep_ca_D_25,:],axis=0)
    nad_se_D_50=np.mean(nad_se_D_50_arr[keep_ca_D_50,:],axis=0)
    nad_se_D_10=np.mean(nad_se_D_10_arr[keep_ca_D_10,:],axis=0)
    
    
    
    #specify nucleotide position
    nuc_pos=np.arange(350)+1573
    
    ####write normalised data to files
    full_arr_R=np.transpose([nuc_pos,nad_R_10[::-1],nad_se_R_10[::-1],nad_R_25[::-1],nad_se_R_25[::-1],nad_R_50[::-1],nad_se_R_50[::-1]])
    
    col_names_R=['genome_position','norm_react_R_10','norm_react_se_R_10','norm_react_R_25','norm_react_se_R_25','norm_react_R_50','norm_react_se_R_50']
    
    full_df_R=pd.DataFrame(full_arr_R,columns=col_names_R)
    
   
    
    full_df_R.to_csv('MS2_pri11_transcript_data.csv',sep=',',index=False)
    
    
    full_arr_Re=np.transpose([nuc_pos,nad_Re_10[::-1],nad_se_Re_10[::-1],nad_Re_25[::-1],nad_se_Re_25[::-1],nad_Re_50[::-1],nad_se_Re_50[::-1]])
    
    col_names_Re=['genome_position','norm_react_Re_10','norm_react_se_Re_10','norm_react_Re_25','norm_react_se_Re_25','norm_react_Re_50','norm_react_se_Re_50']
    
    full_df_Re=pd.DataFrame(full_arr_Re,columns=col_names_Re)
    full_df_Re.to_csv('MS2_pri11_Reassembly_data.csv',sep=',',index=False)
    
    full_arr_G=np.transpose([nuc_pos,nad_G_10[::-1],nad_se_G_10[::-1],nad_G_25[::-1],nad_se_G_25[::-1],nad_G_50[::-1],nad_se_G_50[::-1]])
    
    col_names_G=['genome_position','norm_react_G_10','norm_react_se_G_10','norm_react_G_25','norm_react_se_G_25','norm_react_G_50','norm_react_se_G_50']
    
    full_df_G=pd.DataFrame(full_arr_G,columns=col_names_G)
    
    
    full_df_G.to_csv('MS2_pri11_genome_data.csv',sep=',',index=False)
   
    full_arr_V=np.transpose([nuc_pos,nad_V_10[::-1],nad_se_V_10[::-1],nad_V_25[::-1],nad_se_V_25[::-1],nad_V_50[::-1],nad_se_V_50[::-1]])
    
    col_names_V=['genome_position','norm_react_V_10','norm_react_se_V_10','norm_react_R_25','norm_react_se_V_25','norm_react_V_50','norm_react_se_V_50']
    
    full_df_V=pd.DataFrame(full_arr_V,columns=col_names_V)
    
    full_df_V.to_csv('MS2_pri11_virion_data.csv',sep=',',index=False)
    
    ####write non normalised data to files
    full_arr_R=np.transpose([nuc_pos,ad_R_10[::-1],ad_se_R_10[::-1],ad_R_25[::-1],ad_se_R_25[::-1],ad_R_50[::-1],ad_se_R_50[::-1]])
    
    col_names_R=['genome_position','ad_react_R_10','ad_react_se_R_10','ad_react_R_25','ad_react_se_R_25','ad_react_R_50','ad_react_se_R_50']
    
    full_df_R=pd.DataFrame(full_arr_R,columns=col_names_R)
    full_df_R.to_csv('MS2_pri11_transcript_data_ad.csv',sep=',',index=False)
    full_arr_Re=np.transpose([nuc_pos,ad_Re_10[::-1],ad_se_Re_10[::-1],ad_Re_25[::-1],ad_se_Re_25[::-1],ad_Re_50[::-1],ad_se_Re_50[::-1]])
    
    col_names_Re=['genome_position','ad_react_Re_10','ad_react_se_Re_10','ad_react_Re_25','ad_react_se_Re_25','ad_react_Re_50','ad_react_se_Re_50']
    
    full_df_Re=pd.DataFrame(full_arr_Re,columns=col_names_Re)
    
    full_df_Re.to_csv('MS2_pri11_Reassembly_data_ad.csv',sep=',',index=False)
    full_arr_G=np.transpose([nuc_pos,ad_G_10[::-1],ad_se_G_10[::-1],ad_G_25[::-1],ad_se_G_25[::-1],ad_G_50[::-1],ad_se_G_50[::-1]])
    
    col_names_G=['genome_position','ad_react_G_10','ad_react_se_G_10','ad_react_G_25','ad_react_se_G_25','ad_react_G_50','ad_react_se_G_50']
    
    full_df_G=pd.DataFrame(full_arr_G,columns=col_names_G)
    
    print full_df_G
    
    full_df_G.to_csv('MS2_pri11_genome_data_ad.csv',sep=',',index=False)
   
    full_arr_V=np.transpose([nuc_pos,ad_V_10[::-1],ad_se_V_10[::-1],ad_V_25[::-1],ad_se_V_25[::-1],ad_V_50[::-1],ad_se_V_50[::-1]])
    
    col_names_V=['genome_position','ad_react_V_10','ad_react_se_V_10','ad_react_R_25','ad_react_se_V_25','ad_react_V_50','ad_react_se_V_50']
    
    full_df_V=pd.DataFrame(full_arr_V,columns=col_names_V)
    
    full_df_V.to_csv('MS2_pri11_virion_data_ad.csv',sep=',',index=False)
    
    ####write out difference maps
    full_arr_D=np.transpose([nuc_pos,nad_D_10[::-1],nad_se_D_10[::-1],nad_D_25[::-1],nad_se_D_25[::-1],nad_D_50[::-1],nad_se_D_50[::-1]])
    
    col_names_D=['`genome_position','norm_react_V_T_10','norm_react_se_V_T_10','norm_react_V_T_25','norm_react_se_V_T_25','norm_react_V_T_50','norm_react_se_V_T_50']
    
    full_df_D=pd.DataFrame(full_arr_D,columns=col_names_D)
    
    full_df_D.to_csv('MS2_pri11_difference_data.csv',sep=',',index=False)

    ####produce normalised reactivity profile snapshots 
    sam_funcs.sequence_snapshots(nuc_pos,nad_R_25,nad_se_R_25,col='g',virus='MS2',primer='pri11',condition='R',treatment='25')
    
    sam_funcs.sequence_snapshots(nuc_pos,nad_R_50,nad_se_R_50,col='g',virus='MS2',primer='pri11',condition='R',treatment='50')
    sam_funcs.sequence_snapshots(nuc_pos,nad_R_10,nad_se_R_10,col='g',virus='MS2',primer='pri11',condition='R',treatment='10')
    
    
    
    sam_funcs.sequence_snapshots(nuc_pos,nad_G_25,nad_se_G_25,col='r',virus='MS2',primer='pri11',condition='G',treatment='25')
    
    sam_funcs.sequence_snapshots(nuc_pos,nad_G_50,nad_se_G_50,col='r',virus='MS2',primer='pri11',condition='G',treatment='50')
    sam_funcs.sequence_snapshots(nuc_pos,nad_G_10,nad_se_G_10,col='r',virus='MS2',primer='pri11',condition='G',treatment='10')
    
    sam_funcs.sequence_snapshots(nuc_pos,nad_V_25,nad_se_V_25,col='b',virus='MS2',primer='pri11',condition='V',treatment='25')
    
    sam_funcs.sequence_snapshots(nuc_pos,nad_V_50,nad_se_V_50,col='b',virus='MS2',primer='pri11',condition='V',treatment='50')
    sam_funcs.sequence_snapshots(nuc_pos,nad_V_10,nad_se_V_10,col='b',virus='MS2',primer='pri11',condition='V',treatment='10')

    
    sam_funcs.sequence_snapshots(nuc_pos,nad_D_25,nad_se_D_25,col='m',virus='MS2',primer='pri11',condition='V_T',treatment='25')
    
    sam_funcs.sequence_snapshots(nuc_pos,nad_D_50,nad_se_D_50,col='m',virus='MS2',primer='pri11',condition='V_T',treatment='50')
    sam_funcs.sequence_snapshots(nuc_pos,nad_D_10,nad_se_D_10,col='m',virus='MS2',primer='pri11',condition='V_T',treatment='10')
    

