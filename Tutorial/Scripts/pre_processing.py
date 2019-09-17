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
    