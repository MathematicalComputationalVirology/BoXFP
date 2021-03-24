import os
import sys
import numpy as np
import math
from copy import deepcopy
import matplotlib.pyplot as plt
import pandas as pd
import pickle
import funcFile
import funcPeakAlign
import funcSeqAll
import funcToolsAll
import funcByRef
import funcGeneral
import funcTimeWarp
import shutil
from random import sample
from Qushape import *
import matplotlib.ticker as plticker

from Bio import SeqIO
import Bio.pairwise2 as pwise
from scipy import stats

from scipy.stats.stats import pearsonr

from scipy.stats.stats import spearmanr
from scipy.optimize import curve_fit
import glob

def RX_preprocess(primer,start,end,wlab,Top=999,wins=10,inc=5):

    """
    preprocess chromatographs
    
    Args:
        primer (str): Name of primer under investigation
        start (int): First position in chromatographs to consider in preprocessing
        end (int): Last position in chromatographs to consider in preprocessing
        wlab (str): Name of files into which preprocessed data is deposited
        Top (int): Define the end the position from which the top boundary of the windowing is taken (999 indicates the reference point is the position of the final size marker, n indicates the reference point is the nth size marker and None indicates the end of the chromatograph is the reference point)
        wins (int): Number of different windows used for preprocessing
        inc (int): Increment of sizes between windows in elution points	
    Returns:
        (None): Preprocessed data for each window is exported as a pickle .obj file
    """


    file_list=glob.glob('*.fsa')

    file_list=[i for i in file_list if 'pri'+str(primer) in i]
    print file_list
 

    data_arr0 = data_reader(file_list,end,start)

    data_arr = preprocess(data_arr0)
    n=len(data_arr0)

    data_arr1=mobility_shift(data_arr)
    
    #find peaks in TM traces
    peaksTM=peak_finder(data_arr1,4,.12,cap=3000,TM=1)

    
    #list all those peask that have a disproportionate number of SM peaks
    for i in range(len(peaksTM)):
        
        lenTM=len(peaksTM[i])
        
        if lenTM!=21:
            print file_list[i]+' '+str(i)
    
    #plot the SM traces with the peak positions marked to make sure the peak finder function has found the right peaks.
    sm_plotter(data_arr1,peaksTM,file_list)
    
    #run the data reader version two that carries out the windowing and stores the windows in a pickle .obj file
    DR_windowing(file_list,peaksTM,wlab,top=Top,increment=inc,windows=wins) 

    
def preprocess(data_arr,smooth=4,TM_smooth=1):

    """
    preprocess chromatographs

    Args:
        data_arr (list): List containing the raw chromotograph data
        smooth (int): Width of triangular smoothing of main data traces
        TM_smooth (int): Width of triangular smoothing for size marker
    Returns:
        (list): List containing the processed chromatograph data
    """

    data_arr1 = []

    for i,data in enumerate(data_arr):
        
        dataset_list = [] 
        #add elution time to data array
        dataset_list.append(np.array(data['Position']))
        
        #smoothing
        for ind in range(1,len(data.columns)-1):
            dataset_list.append(funcToolsAll.fSmooth( data.values[:,ind], smooth, 'triangle' ))
            
        dataset_list.append(funcToolsAll.fSmooth( data.values[:,4], TM_smooth, 'triangle' ))
        
        #baseline adjustment
        for ind in range(1,len(dataset_list)):
            dataset_list[ind] = funcToolsAll.baselineAdjust(dataset_list[ind],60, 'triangle')
        
        #autodecay (not performed on tape measure trace)   
        for ind in range(1,len(dataset_list)-1):
            dataset_list[ind]=funcToolsAll.autoDecaySum(dataset_list[ind], 200)
        
        data_arr1.append(dataset_list)
        
    return data_arr1


def mobility_shift(data_arr):
    
    """
    Mobility shift between footprinting trace and the ddA sequence trace

    Args:
        data_arr (list): List containing the preprocessed chromatograph data

    Return:
        (list): List containg the mobility shift corrected chromatograph data
    """
    
    data_arr1 = []
    
    for data in data_arr:
        data_RX= funcToolsAll.fMobilityShift(data[2],data[1],'HEX','5-FAM','posSim')
    
        data_out = [data[0],data_RX,data[2],data[3],data[4]]
        data_arr1.append(data_out)
    
    return data_arr1
    

def signal_alignment(data,align_data,move1,move2,fifcorr = True):

    """
    Align traces in chromatographs to each other based on trace alignment sequence outlined in trace_align

    Args:
        data (array): data array containing chromatograph data
        align_data (array): data array containing chromotograph data used for alignment
        move1 (array): Sequence of warping path used on the non align dataset
        move2 (array): Sequence of warping path used on the align dataset
        fifcorr (bool): Use correlations in alignment process

    Returns:
       (array): array of data for first dataset aligned to the second
    """
    
    #align RX data
    data_RX = funcPeakAlign.splineSampleData(data[1],align_data[1], move1, move2,isCorr=ifcorr)

    #align S1 data
    data_S1 = funcPeakAlign.splineSampleData(data[2],align_data[2], move1, move2,isCorr=ifcorr)

    #align S2 data
    data_S2 = funcPeakAlign.splineSampleData(data[3],align_data[3], move1, move2,isCorr=ifcorr)

    #align TM data 
    data_TM = funcPeakAlign.splineSampleData(data[4],align_data[4], move1, move2,isCorr=ifcorr)

    new_data = [np.array(data[0]),data_RX,data_S1,data_S2,data_TM]
        
    return new_data
    

def trace_align(data_arr,ind,ref,ifmob=True,ifdtw=True,ifpeak=True,gap1=-0.2):
  

    """
    Aligns traces to each other according to a reference data_set (align_data) and specific trace [ind]

    Args:
        data_arr (list): list containing all of the datasets in the ensemble
        ind (int): integer specifying the trace used for the alignment process
        ref (int): reference dataset used for alignment
        ifmob (bool): specify whether mobility shift should be performed
        ifdtw (bool): specify whether dynamic time warp should be performed
        ifpeak (bool): specify whether peak aliggment should be performed
        gap1 (float): gap penalty

    Returns:
        (list): List containing the ensemble with all datasets aligned to each other

    """ 
    #set up arrays to put data in at various stages
    data_out = []
    data_arr2 = []
    data_arr3 = []
     
        
    data_arr0 = data_arr
    
    if ifmob:
        data_arr1 = mobility_shift(data_arr0)
    else:
        data_arr1 = data_arr0
    
    #isolate the alignment data
    align_data = np.array(data_arr1[ref])
    
    #align dataset extraction
    if ifdtw:
        for i,data in enumerate(data_arr1):

            #find alignments based on ind trace
            move1,move2 = find_DTW_match(align_data[ind],data[ind])

            #align signals
            data2 = signal_alignment(data,align_data,move1,move2)

            data_arr2.append(data2)
    else:
        data_arr2 = data_arr1
    
    #isolate alignment data
    align_data = np.array(data_arr2[ref])

    if ifpeak:
    #peak alignment process
        for i,data in enumerate(data_arr2):

            #find alignments based on ind trace
            move1,move2 = find_peak_match_X(align_data[ind],data[ind],gap=gap1)

            #align signals
            data3 = signal_alignment(data,align_data,move1,move2)

            data_arr3.append(data3)
    else:
        data_arr3 = data_arr2
    
    return data_arr3


def normalise(data):
    
    """
    Normalise data

    Args:
        data (array): Array containing unnormalised  data profile


    Returns:
        (array): Array containing normalised data profile
    """
   
    data_deep  = deepcopy(data)
    
    mean = np.mean(data_deep)
    std = np.std(data_deep)
    
    data_out = (data_deep-mean)/std
    
    return data_out
    
    
def find_DTW_match(align_data,data, ifwrt = False):

    """
    Find optimal dynamic time warp (DTW) path

    Args:
        align_data (array): Dataset of chromatograph data to align the second dataset too.
        data (array): Second chromatograph dataset

    Returns:
        (tuple):
            linkX0 (array): Array containing the warping path for the alignment dataset
            linkX1 (array): Array containing the warping path for the secondary dataset


    """

    #if normalisation is required ask for it
    if ifwrt:
        data1,align_data1 = normalise_wrt(data,align_data)
    else:
        data1 = normalise(data)
        align_data1 = normalise(align_data)
        
    #set up sakoe-chiba window length
    r1=len(align_data1)*0.05
    #find optimal warping paths
    pathX,pathY= funcTimeWarp.myDTW(align_data1, data1, derivative=True, costMType='1',D=0.05, bandType = "noband", r=r1, gap=True)
    
    #perform post peak matching?
    linkX0,linkX1 = funcToolsAll.postpeakMatch0(pathX,pathY,step= 100)
    
    return linkX0,linkX1

    

def find_peak_match_X(align_data,data,ifwrt = True,gap=-0.2):

   """
    Alignment of two datasets based on peak matching

    Args:
        align_data (array): Chromatograph data used for alignment
        data (array): Chromatograph data to be realigned

    Returns:
        (tuple):
            linkX0 (array): Array containing the warping path for the alignment dataset
            linkX1 (array): Array containing the warping path for the secondary dataset


    """

    #normalise if requested
    if ifwrt:
        data1,align_data1 = normalise_wrt(data,align_data)
    else:
        data1 = normalise(data)
        align_data1 = normalise(align_data)
      
    #obtain parameters object
    dParameters= funcPeakAlign.DPeakAlignParams()
    
    #change similarity function
    #dParameters['simFunc'] = 'Derivative'
    
    dParameters['gap'] = gap
    
    #obtain peaklist for alignment data
    peakList0=funcPeakAlign.fPeakList(align_data1, isDel=False, isAdd=False,repType='Cubic')
    
    #obtain peaklist for data
    peakList1=funcPeakAlign.fPeakList(data1, isDel=False, isAdd=False,repType='Cubic')
    
    #calculate peak alignments 
    aligned0,aligned1=funcPeakAlign.myPeakAlignment(peakList0,peakList1,dParameters)
    
    #select peaks that only have a peak associated in other trace
    linkX0,linkX1=funcToolsAll.findAlignedPeakX(aligned0,aligned1)
    
    #perform post peak matching
    linkX0,linkX1=funcToolsAll.postpeakMatch0(linkX0,linkX1)
    
    return np.array(linkX0,int),np.array(linkX1,int)
    

def data_reader(file_list,top,bottom):
    
    """
    Generic data reader

    Args:
        file_list (list): List of the file names for datasets in the ensembles
        top (int): The highest elution time point in the region of interest
        bottom: the lowest elution time point in the region of interest
    returns:
        (list): List containing the chromatograph datasets from the ensemble

    """
    
    #intialise data array
    data_arr = []
    
    #iteratively read files and extract data
    for i,file in enumerate(file_list):
        data = pd.read_csv(file.strip('.fsa\n')+'_raw.csv')

        #tidy data
	
        data_arr.append(data_tidy(data,top,bottom))
    
    return data_arr

def DR_windowing(file_list,TM_peaks,date,top=999,bottom=0,increment=5,windows=10):
    
    """
    Performs preprocessing over several ROI windows

    Args:
        file_list (list):List containg the file names for the datasets in the ensemble
        TM_peaks (list): List containing all of the size marker peak positions in the datasets
        name (str): String indicating the common name for the .obj pickle files
        top (int): The highest peak in the peak list used
        bottom (int): The lowest peak in the peak lists used
        increment (int): Integer specifying the number of elution points by which the ROI windows change on either end
        windows (int): Integer specifying the number of windows to be used.

        Returns:
            None

    """
    
    for j in range(windows):
        
        
        #intialise data array
        data_arr = []
        increm=None
        #iteratively read files and extract data
        for i,file in enumerate(file_list):
            data = pd.read_csv(file.strip('.fsa\n')+'_raw.csv')
            if top!=999 and top!=None:
                top1=TM_peaks[i][top-1]+50+(j)*increment
	    elif top==None:
		top1=None
		increm=(j)*increment
            else:
                top1=TM_peaks[i][-1]+(j)*increment+30
            #smooth TM trace
            bottom1=TM_peaks[i][bottom]-j*increment-30
            #tidy data and add to final array
            data_arr.append(data_tidy(data,top1,bottom1,step=increm))
        data_arr1=preprocess(data_arr)
        data_arr2=mobility_shift(data_arr1)
        file_path=date+'_'+str(j)+'.obj'
        print str(bottom1)+'_'+str(top1)
        
        file1=open(file_path,'wb')
        
        pickle.dump(data_arr2,file1)
        
        file1.close()
    
def data_tidy(data,top,bottom,step=0):

    """
    Remove all data not in the region of interest (ROI)

    Args:
        data (array): Data array containing the electropherogram data
        top (int): Integer specifying the highest elution time point for the ROI
        bottom (int): Integer specifying the lowest elution time point for the ROI
        step (int): Integer specifying the top when full trace is used.
    Returns:
        (array): data array containing the
        electropherogram data for the ROI only

    """

    #remove data above the ROI
    if top==None:
    	data_out = data[data['Position']<(np.max(data['Position'])-10-step)]
    else:
    	data_out = data[data['Position']<top]
    #print data_out 
    #remove data below the ROI
    data_out = data_out[data_out['Position']>bottom]

    #set highly negative data to 0
    #data_out[data_out<=0] = 0
    
    return data_out


def remove_outliers(x,y, outlierConstant=1.5):
    
    """
    Remove outliers from trace based on amplitude

    Args:
        x (array): Array containing x coordinates of trace
        y (array): Array containing y coordinates of trace
        outlierConstant (float): IQR multiplier used for outlier determination
    Returns:
        (tuple):
            resultListx (array): Array of x coordinates in trace with outliers removed

            resultListy (array): Array of y coordinates in trace with outliers removed
    """

    ax = np.array(x)
    ay = np.array(y)
    upper_quartile = np.percentile(ay, 75)
    lower_quartile = np.percentile(ay, 25)
    IQR = (upper_quartile - lower_quartile) * outlierConstant
    quartileSet = upper_quartile + IQR
    resultListx = []
    resultListy = []
    for i in range(len(ay)):
        if ay[i] <= quartileSet:
            resultListx=np.append(resultListx,ax[i])
            resultListy=np.append(resultListy,ay[i])
    return resultListx,resultListy

def peak_finder(data_arr,ind,perc,TM=0,pn=None,cap=None,lower_limit=0,Marker_set='ROX'):
    
    """
    Peak finding function

    Args:
        data_arr (list): List containing the chromatograph data
        ind (int): Integer specifying the trace channel
        perc (float): Ratio of max peak intensity used as cutoff criterion
        TM (int): Set to find the specific number of peaks
        pn (int): Specific number of peaks to find
        cap (int): Cap specifying the highest peak amplitude considered in peak finding
        lower_limit (int): Lowest elution time point considered in peak finding
        Marker_set (str): Name of the size marker set used
    Returns:
        (list): List of recorded peak values for each dataset

    """
    
    if pn==None:
        if Marker_set=='ROX':
            pn=21
            close_limit=15
        elif Marker_set=='LIZ':
            pn=34
            close_limit=80
        
    #initialise peak array
    peak_arr = []
    peak_rec = []
    #iterate through datasets
    for i,data in enumerate(data_arr):
        
        #initialise peak position and value arrays
        peak_pos = []
        peak_val = []
        
        
        #deal with pandas objects
        if isinstance(data,pd.DataFrame):
            
            #extract values
            data = data.values
            
            #iterate through datapoints
            for k in range(1,len(data[:,0])-1):
                #workout differences between data point and those before and after
                diff1 = data[k+1,ind]-data[k,ind]
                diff2 = data[k,ind]-data[k-1,ind]
                
                #if derivative crosses bin value
                if (diff2>0) and (diff1<=0):
                    peak_pos = np.append(peak_pos,data[k,0])
                    peak_val = np.append(peak_val,data[k,4])
                    
        #contingency for objects that are not pandas
        else:
            for k in range(1,len(data[0])-1):
                diff1 = data[ind][k+1]-data[ind][k]
                diff2 = data[ind][k]-data[ind][k-1]
                if (diff2>0) and (diff1<=0):
                    peak_pos = np.append(peak_pos,data[0][k])
                    peak_val = np.append(peak_val,data[4][k])
 
        #remove all of those peaks below a certain threshold
        #peak_pos,peak_val=remove_outliers(peak_pos,peak_val)
        if cap!=None:
            peak_pos=peak_pos[peak_val<cap]
            peak_val=peak_val[peak_val<cap]
            
        peak_val=peak_val[peak_pos>lower_limit]    
        peak_pos=peak_pos[peak_pos>lower_limit]
        
    
        
        
        if TM==1:
            peak_inds = peak_val.argsort()[-pn:][::-1]
            peaks=peak_pos[peak_inds]
            peak_val=peak_val[peak_inds]
            
            peaks_inds = np.argsort(peaks)
            
            peaks1=peaks[peaks_inds]
         
            
        else:
            peaks1=peak_pos[peak_val>np.max(peak_val)*perc]   
        
        j=0
        while j<(len(peaks1)-1):
            diff2=peaks1[j+1]-peaks1[j]

            if diff2<close_limit:

                av_pos = (peaks1[j]+peaks1[j+1])/2

                peaks1[j] = av_pos

                peaks1 = np.delete(peaks1,j+1)
                j+=1
            j+=1
        
	if Marker_set=='LIZ':
		peaks1=peaks1[-pn:]      
       
        peak_arr.append(peaks1)
          
    return peak_arr

def peak_diffs(peak_arr, plot = 0):    
    
    """
    Calculate differences between adjacent peak positions

    Args:
        peak_arr (list): Peak lists for the datasets in ensemble
    Returns:
        (tuple):
            peak_diff_av (array): Average distances between peaks in ensemble
            peak_diff_arr (list): All the peak differences for the entire ensemble
    """
    
    #initialise arrays
    peak_diff_arr = []
    peak_rec = []
    peak_diff_av = []
    
    #read thorough peak lists
    for peaks in peak_arr:
        
        #initialise arrays
        peak_diff = []
        peak_av= []
        
        #read through peak lists and calculate differences and averages
        for i in range(1,len(peaks)):
            pdiff = peaks[i]-peaks[i-1]
            pave=(peaks[i]+peaks[i-1])/2
            
            peak_diff.append(pdiff)
            peak_av.append(pave)
            peak_rec = np.append(peak_rec,pdiff)
        
        peak_diff_av.append(peak_av)
        peak_diff_arr.append(peak_diff)
    
    #plot histogram of differences
    if plot == 1:
        plt.hist(peak_rec,bins=20)
        plt.xlabel('Elution time difference')
        plt.show()
    
    #return differences and averages
    return peak_diff_av,peak_diff_arr


def find_first(a,b):
    
    """
    Find the first occurence of sub array on array
    Args:
        a (array): Larger array
        b (array): sub array
    Returns:
       (int): first postion of array where sub array occurs
    """

    #calculate length of sub array
    len_b = len(b)
    
    #if lengths are equal return first index
    if len(a)==len(b):
        return 0
    
    #if arrays are input the wrong way round return warning
    elif len(a)<len(b):
        print 'wrong arrangement of arrays!'
    
    #scan through sections and compare
    else:
        for i in range(len(a)-len(b)+1):
            
            #if section and subarray are equal return start index of section
            if np.array_equal(a[i:i+len_b],b):
                return i

def RX_position(dfile,refSeq,seqInd=2,dinds=[],clip=None,searchStart=0,searchSpan=None,corr=0.7,voteInc=20):

    """
	Identify first position o fsize markers in based on the ddA trace and a reference sequence

    Args:
	dfile (str): Name of data file
	refSeq (str): Name of fasta file containing the reference sequence
	seqInd (int): Channel containg the ddA ladder trace
	dinds (list): Indices of datasets to use in the ensemble
	clip (int): Number of nucleotides in the size marker trace to be used, taken from the start of the size marker
	searchStart (int): Start position of search if a specific section of the reference sequence is required 
	searchSpan (int): Number of nucleotides from the start position of region of interest to search
	corr (float): Correlation cutoff for poor sequence traces in the ensemble
	voteInc (int): Number of intervals over which to assess voting score

    Returns:
	None: Print most likely nucleotide position of start of size marker set  
    """    

    sfile=open(dfile,'rb')
    data_arr=pickle.load(sfile)
    data_arr1=[]
    if dinds>0:
	for i in dinds:
		data_arr1.append(data_arr[i])
    else:
	data_arr1=data_arr

    partition_seq=seq_partitioning(data_arr1,seqInd)

    col=iter(plt.cm.rainbow(np.linspace(0,1,100)))

    for i in range(len(partition_seq)):
        c=next(col)
        plt.plot(np.arange(1,351),partition_seq[i][:,1],c=c)
    plt.show()    

    S_cov_matrix=correl_assessor(partition_seq,1)

    sum_array=np.nanmean(S_cov_matrix,axis=1)

    #print sum_arrayB
    #find those sequences with a high mean correlation
    keep_seqs=np.where(sum_array>corr)[0]


    plt.hist(sum_array)
    plt.show()

    partition_seq=[partition_seq[i] for i in keep_seqs]

    pos_test_mat = np.empty([20,20])
    corr_test_mat = np.empty([20,20])
    nuc_count_mat = np.empty([20,20])
    signif_mat = np.empty([20,20])

    x = np.linspace(0,1,voteInc)

    for i,x1 in enumerate(x):
        for j,x2 in enumerate(x):

            #generate consensus sequence
            corel_val,vote_arr = position_vote(partition_seq,x1,x2)

            #get nucleotide count
            nuc_count_mat[j] = np.count_nonzero(vote_arr)/float(len(vote_arr))*100

            #perform sequence search
	    if searchSpan==None:

            	signif_mat[i,j],pos_test_mat[i,j],corr_test_mat[i,j]= sequence_search(refSeq,vote_arr)

	    else:

		signif_mat[i,j],pos_test_mat[i,j],corr_test_mat[i,j]= sequence_search_area(refSeq,vote_arr,searchStart,searchSpan)


    print corr_test_mat
    print signif_mat
    print pos_test_mat.astype(int)

    #get max correlation values of cutoff1 and cutoff2
    max_corr = np.where(corr_test_mat == np.nanmax(corr_test_mat))

    #print the position of this highest correlation 
    print pos_test_mat[max_corr[0],max_corr[1]]
    print np.nanmax(corr_test_mat)


def seq_partitioning(data_arr,ind,tm_cutoff=21,marker_set='ROX'):
    
    """
    Partitioning function for the sequencing traces

    Args:
        data_arr (list): All of the datasets in the ensemble
        ind (int): Sequence channel in data for consideration
        tm_cutoff (int): The number of size marker peaks to consider in these calculations
        Marker_set (str): Name of the size marker set used

    Return:
        (list): Partitioned sequence trace peaks for all the datasets in the ensemble

    """

   
    #vectors dictationg the marker sizes and differences between markers in the tape measure
    if marker_set.upper()=='ROX':
        
        marker_sizes = np.array([50, 60, 90, 100, 120, 150, 160, 180, 190, 200, 220, 240, 260, 280, 290, 300, 320, 340, 360, 380,400])
    
    elif marker_set.upper()=='LIZ':
        
        marker_sizes = np.array([20, 40, 60, 80, 100, 114, 120, 140, 160, 180, 200, 214, 220, 240, 250, 260, 280, 300, 314, 320, 340, 360, 380, 400, 414, 420, 440, 460, 480, 500, 514, 520, 540, 560, 580, 600])
    else:
        print 'unknown marker set!'
        quit()
        
    marker_diffs=np.diff(marker_sizes)
    
    #find peaks in tape measure
    peaksTM = peak_finder(data_arr,4,.25)
    #sm_plotter(data_arr,peaksTM,fl) 
    #array of nucleotide position
    nuc_pos = np.arange(351)
    
    #initialise peak binning array
    reduced_peaks = []
    
    #iterate through data 
    for i in range(len(data_arr)):
        nuc_pos = range(350)
        
        #extract data index
        data = data_arr[i]
        
        #extract peaks
        peaks = peaksTM[i]
        
        if tm_cutoff<21:
            peaks = peaksTM[i][:tm_cutoff]
        else:
            peaks = peaksTM[i]
        
        
        
        
        
        #initials peak lists
        peak_list = []
        pos_list = []
        bin_list = []
        
        if isinstance(data,pd.DataFrame):
                
            #extract data values
            data1 = data.values
            
            
        else:

            #transpose data - fudge factor can be improved
            data1=np.transpose(data)
 
    

        peaks_trace1 = funcPeakAlign.fPeakList(data1[:,ind],repType = 'Cubic')

        p_ind = peaks_trace1['pos']
        
        peak_ind=np.unique(p_ind)
        
        p_pos = data1[peak_ind,0]
        
    
        p_amp =data1[peak_ind,ind]
        p_width = peaks_trace1['averW']
        p_std=peaks_trace1['stdW']
        
        peak_list=[p_pos,p_amp]
        peak_list=np.transpose(peak_list)
        
        shoulder_data = shoulder_finder(p_pos,data1,ind,i)
        #if shoulders found add them to peak list
        if len(shoulder_data)>0:
            peak_list1 = np.vstack((peak_list,shoulder_data))
            ind1=np.argsort(peak_list1[:,0])
            peak_list1=peak_list1[ind1]
        else:
            peak_list1 = peak_list
        
                
        #delete peaks if they are too close to each other
        x = 0
        while x<len(peak_list1)-1:
            
            #find the position and amplitudes of two neighbouring peaks
            p_pos1 = peak_list1[x,0]
            p_pos2 = peak_list1[x+1,0]
            
            p_amp1 = peak_list1[x,1]
            p_amp2 = peak_list1[x+1,1]
            
            
            #calculate difference in position between the two peaks
            diff1 = p_pos2-p_pos1
            
            #if the difference is less than 4
            if diff1<4:
                
                #if peak1 is in original peak list but not peak2 delete peak2     
                if (p_pos1 in peak_list[:,0]) and (p_pos2 not in peak_list[:,0]):
                    del_x = x+1
                    peak_list1 = np.delete(peak_list1,del_x,0)
                    x+=1
                    
                #if peak2 is in original peak list but not peak1 delete peak1    
                elif(p_pos2 in peak_list[:,0]) and (p_pos1 not in peak_list[:,0]):
                    del_x = x
                    peak_list1 = np.delete(peak_list1,del_x,0)
                else:
                    x+=1
                     
                
                
            else:
                x+=1
        
        x=0
        
        
        miss_pos = []
        miss_amp = []
        
        #fill in any positions that have missing peaks due to poor data quality
        
        #read through peak list
        for x in range(len(peak_list1)-1):
            
            #look at the difference in position between neighbouring peaks 
            diff1=peak_list1[x+1,0]-peak_list1[x,0]

            #determine spacing
            spacing_we=2*p_width-4*p_std
            
            #calculate number of missing peaks expected in a space
            miss_peaks=np.floor(diff1/(p_width-2*p_std))
            
            
            #if more that one missing peak space suspected
            if miss_peaks>1:
                

                #fill the missing space with peak data
                for z in range(int(miss_peaks)-1):
                
                
                    #determine position of new peak
                    av_pos=peak_list1[x,0]+(z+1)*int(p_width)
                    
                    #add new position to missing position array
                    miss_pos=np.append(miss_pos,av_pos)

                    #find amplitude of new peak
                    av_amp=data1[data1[:,0]==av_pos,ind]
                    
                    #add to missing amplitude data
                    miss_amp = np.append(miss_amp,av_amp)

        #create array for missing peak information
        miss_arr = np.transpose([miss_pos,miss_amp])
        
        
        #add the missing peaks to the new peak array
        if len(miss_arr)>0:
            peak_list2=np.vstack((peak_list1,miss_arr))
            ind2=np.argsort(peak_list2[:,0])
            peak_list2=peak_list2[ind2]
        else:
            peak_list2=peak_list1

       
        
        peak_list5=np.zeros((350,2))
    
        pos_ind=0
        #iterate through all peaks
        for j in range(len(peaks)-1):
            #extract TM peaks that define the trace region of interest
            start = peaks[j]
            end = peaks[j+1]
            
            #calculate difference between peak postions
            diff = end-start 
            
            
            #extract number of bins for specific region
            bins = marker_diffs[j]
        
            #set up the bin width
            nuc_sep = diff/bins
                
            #extract all of those trace values above TM start position
            peak_list3=peak_list2[peak_list2[:,0]>int(start),:]
            
            peak_list4=peak_list3[peak_list3[:,0]<int(end),:]
           
            #iterate through bins
            for t in range(bins):
                
                #extract all trace values within bin
                S_bin = peak_list4[(peak_list4[:,0]>(start+(t)*nuc_sep)),:]
                
            
                
                S_bin=S_bin[(S_bin[:,0]<(start+(t+1)*nuc_sep)),:]
                
                peak_list5[pos_ind,0]=start+(t+0.5)*nuc_sep
                #initialise binning array
                if len(S_bin)==0:
                    pos_ind+=1  
                    continue
                elif len(S_bin)==1:
                    peak_list5[pos_ind,1]=S_bin[:,1]
                
                elif len(S_bin)>1:
                    
                    peak_list5[pos_ind,1]=np.max(S_bin[:,1])
                    
                    
                pos_ind+=1    
                
        
             
        
        reduced_peaks.append(peak_list5)
                
        
    return reduced_peaks

        
def RX_partitioning_single(data_arr,ind,file_list,Issues=[],Skip=False,ll=0,perc=0.25,tm=0,extension=0,marker_set='ROX'):

    """
    partition footprinted trace for single replicates


    Partitioning function for the sequencing traces

    Args:
        data_arr (list): All of the datasets in the ensemble
        ind (int): Footprint channel in datasets
        file_list (list): File names for the datasets in the ensemble
        Issues (list): datasets with size marker trace irregularity
        Skip (bool): indicates whether datasets should be skipped
        ll (int): Lower limit of elution time points to be considered for size marker trace
        perc (float): ratio of maximum peak intensity for cutoff
        tm (int): Set to find specific number of peaks
        extension (int): How many peaks to extend the size marker trace by. Negative values indicate  a reduced number of size marker peaks is being considered
        Marker_set (str): Name of the size marker set used


    Return:

        data_out2 (list): The partitioned footprinted data traces peakList format

    """    
   
    #vectors dictationg the marker sizes and differences between markers in the tape measure
    
    if marker_set.upper()=='ROX':
        
        marker_sizes = np.array([50, 60, 90, 100, 120, 150, 160, 180, 190, 200, 220, 240, 260, 280, 290, 300, 320, 340, 360, 380,400])
    
    elif marker_set.upper()=='LIZ':
        
        marker_sizes = np.array([20, 40, 60, 80, 100, 114, 120, 140, 160, 180, 200, 214, 220, 240, 250, 260, 280, 300, 314, 320, 340, 360, 380, 400, 414, 420, 440, 460, 480, 500, 514, 520, 540, 560, 580, 600])
    else:
        print 'unknown marker set!'
        quit()
        
    marker_diffs=np.diff(marker_sizes)
    
    #find peaks in tape measure
    peaksTM = peak_finder(data_arr,4,perc,TM=tm,lower_limit=ll)
    
    sam_funcs.sm_plotter(data_arr,peaksTM,file_list)
    #array of nucleotide position
    if extension>0:
	    for k in range(len(peaksTM)):
		if k in Issues:
			if Skip:
				print 'horse'
			else:
				continue
		TM_diffs=np.diff(peaksTM[k])
		values=np.divide(TM_diffs,marker_diffs/10)
		
		positions=np.arange(1,21)
		if k==0:
			pos_out=positions
			vals_out=values
		else:
			pos_out=np.append(pos_out,positions)
			vals_out=np.append(vals_out,values)
	    popt, pcov = curve_fit(TM_func, pos_out,vals_out,maxfev=8000)
	    extend=np.arange(22,22+extension)
	    added_peak_diffs=TM_func(extend,*popt) 
	    #print added_peak_diffs
	   
	    #array of nucleotide position
	    for pp in range(len(extend)):
		marker_diffs=np.append(marker_diffs,10)
		marker_sizes=np.append(marker_sizes,marker_sizes[-1]+10)	   


    
    #initialise peak binning array
    reduced_peaks = []
    
    data_out2 = []
    #iterate through data 
    for i in range(len(data_arr)):
        
        
        
        #extract data index
        data = data_arr[i]
        
        #extract peaks
        if extension<0:
                peaks=peaksTM[i][:abs(extension)]

        else:
                peaks = peaksTM[i]
                if extension>0:
                        for pt in added_peak_diffs:
                                peaks=np.append(peaks,peaks[-1]+pt)
        
    

        #initials peak lists
        peak_list = []
        pos_list = []
        bin_list = []
        
        
        if isinstance(data,pd.DataFrame):
                
            #extract data values
            data1 = data.values
            
            
        else:

            #transpose data - fudge factor can be improved
            data1=np.transpose(data)
       
        #extract peak list for thetarget trace
        peaks_trace1 = funcPeakAlign.fPeakList(data1[:,ind],repType = 'Cubic')
        
        #extract average and standard deviations of widths
        p_width = peaks_trace1['averW']
        p_std=peaks_trace1['stdW']
        
        
        #extract indices of peak positions
        p_ind = peaks_trace1['pos']
        
        #remove duplicates of indices
        p_ind =np.unique(p_ind)
        
        #determine peak amp and position
        p_pos = data1[p_ind,0]
    
        p_amp = data1[p_ind,ind]
      

        #find shoulders 
        shoulder_data = shoulder_finder(p_pos,data1,ind,i)
        
        # transpose peak list
        peak_list=np.transpose([p_pos,p_amp])

        #if shoulders found add them to peak list
        if len(shoulder_data)>0:
            peak_list1 = np.vstack((peak_list,shoulder_data))
            ind1=np.argsort(peak_list1[:,0])
            peak_list1=peak_list1[ind1]
        else:
            peak_list1 = peak_list
        
                
        #delete peaks if they are too close to each other
        x = 0
        while x<len(peak_list1)-1:
            
            #find the position and amplitudes of two neighbouring peaks
            p_pos1 = peak_list1[x,0]
            p_pos2 = peak_list1[x+1,0]
            
            p_amp1 = peak_list1[x,1]
            p_amp2 = peak_list1[x+1,1]
            
            
            #calculate difference in position between the two peaks
            diff1 = p_pos2-p_pos1
            
            #if the difference is less than 4
            if diff1<4:
                
                #if peak1 is in original peak list but not peak2 delete peak2     
                if (p_pos1 in peak_list[:,0]) and (p_pos2 not in peak_list[:,0]):
                    del_x = x+1
                    peak_list1 = np.delete(peak_list1,del_x,0)
                    x+=1
                    
                #if peak2 is in original peak list but not peak1 delete peak1    
                elif(p_pos2 in peak_list[:,0]) and (p_pos1 not in peak_list[:,0]):
                    del_x = x
                    peak_list1 = np.delete(peak_list1,del_x,0)
                else:
                    x+=1
                     
                
                
            else:
                x+=1
        
        x=0
        
        
        miss_pos = []
        miss_amp = []
        
        #fill in any positions that have missing peaks due to poor data quality
        
        #read through peak list
        for x in range(len(peak_list1)-1):
            
            #look at the difference in position between neighbouring peaks 
            diff1=peak_list1[x+1,0]-peak_list1[x,0]

            #determine spacing
            spacing_we=2*p_width-4*p_std
            
            #calculate number of missing peaks expected in a space
            miss_peaks=np.floor(diff1/(p_width-2*p_std))
            
            
            #if more that one missing peak space suspected
            if miss_peaks>1:
                

                #fill the missing space with peak data
                for z in range(int(miss_peaks)-1):
                
                
                    #determine position of new peak
                    av_pos=peak_list1[x,0]+(z+1)*int(p_width)
                    
                    #add new position to missing position array
                    miss_pos=np.append(miss_pos,av_pos)

                    #find amplitude of new peak
                    av_amp=data1[data1[:,0]==av_pos,ind]
                    
                    #add to missing amplitude data
                    miss_amp = np.append(miss_amp,av_amp)

        #create array for missing peak information
        miss_arr = np.transpose([miss_pos,miss_amp])
        
        
        
        #add the missing peaks to the new peak array
        if len(miss_arr)>0:
            peak_list2=np.vstack((peak_list1,miss_arr))
            ind2=np.argsort(peak_list2[:,0])
            peak_list2=peak_list2[ind2]
        else:
            peak_list2=peak_list1
      
        
       p_pos2=[]
        p_amp2=[]
        
        did_bin=[]
        bin_counter=[]
        bin_pos=[]
        
           
        for j in range(len(peaks)-1):
            #extract TM peaks that define the trace region of interest
            tm_bins=0
            
            #extract adjacent peaks in TM
            start = peaks[j]
            end = peaks[j+1]
            
            #calculate distance between TM peaks
            diff=end-start
            
            # find the nearest target peak downwind of start TM peak
            start_ind = find_nearest_ind(peak_list2[:,0],start)
            # find the nearest target peak downwind of start TM peak
            end_ind = find_nearest_ind(peak_list2[:,0],end)
            
            #extract target peak list for space between TM peaks
            peak_list3=peak_list2[start_ind:end_ind,:]
            nuc_sep=diff/marker_diffs[j]
           
            #add clipped target peaks to list
            did_bin.append(peak_list3)
                
              
        #work through clipped target peak list  
        for j in range(len(did_bin)):
            did=did_bin[j]
            did2=np.reshape(did,(-1,2))
          
                
            #set TM peaks at the beginning and end of the clipped target peak list    
            start = peaks[j]
            end=peaks[j+1]

            #determine the difference between end and start
            diff5=end-start

            #calculate seperation of nucleotides
            nuc_sep=diff5/marker_diffs[j]
            
            
            bin_alloc=[]
            
            did_bin2=[]
            
            
            #determine peaks in different marker bins
            for kk in range(marker_diffs[j]):
                
                did3=did2[did2[:,0]>=start+(kk-1)*nuc_sep,:]
                did3=did3[did3[:,0]<start+(kk)*nuc_sep,:]

                #convert bins with 4 peaks in to 2 peaks and remove surpluses
                if len(did3)>2:
                    
                    args2=np.argsort(did3[:,1][::-1])
                    did3=did3[args2[:1],:]
                    
                bin_alloc=np.append(bin_alloc,len(did3))
                
                did_bin2.append(did3)
                
                
            skip=0
            
            #reassignment of data to evenly spread between bins 
            
            #if all bins have one nucleotide in assign each to respective bin 
            if  np.all(bin_alloc==1):
                for kk in range(marker_diffs[j]):
                    if j == 0 and kk==0 :
                        did_arr1=did_bin2[kk]
                    else:
                        did_arr1=np.vstack((did_arr1,did_bin2[kk]))
                        
                    
            #if the total number of bins in the system is equal to the predicted amount 
            elif np.sum(bin_alloc)==marker_diffs[j]:
                
                #if a bin has no peaks in it skip
                for kk in range(marker_diffs[j]):
                    if bin_alloc[kk]==0:
                        if j == 0 and kk==0 :
                            skip=1
                        continue
                    #else add the data to the new list    
                    else:
                        if j == 0 and kk==0 :
                            did_arr1=did_bin2[kk]
                        elif skip==1:
                            did_arr1=did_bin2[kk]
                        else:
                            did_arr1=np.vstack((did_arr1,did_bin2[kk]))
                    
                    
            #if number of peaks is greater than expected
            elif np.sum(bin_alloc)>marker_diffs[j]:
                diff4=np.sum(bin_alloc)-marker_diffs[j]
                
                #if any bins exist with no peaks in them
                if np.any(bin_alloc==0):
                    
                    #extract bin positions with two peaks in them
                    args2=np.where(bin_alloc==2)
                    
                    args2 = np.array(args2)[0]
                    
                    #extract bin positions with no peaks in them 
                    args0=np.where(bin_alloc==0)
                    
                    args0 = np.array(args0)[0]

                    #calculate the sum of distances between args2 and args0
                    dist_arr=distance_determiner(args2,args0)

                    #sort the arguments based on descending sum off distance
                    dist_args=np.argsort(dist_arr[::-1])

                    #order args2 based on Sum of Distance
                    del_args=args2[dist_args]
                    
                    #calculate difference of numbers of bins containing 2 and 0 peaks
                    num2=np.count_nonzero(bin_alloc==2)
                    
                    num0=np.count_nonzero(bin_alloc==0)
                    
              
                    diff5=num2-num0
                    
                    #set the indices of the bins with 2 peaks that should be averaged
                    del_args=del_args[:(diff5)]
                    
                    
                    for kk in range(marker_diffs[j]):
                        #if no peaks in bin continue
                        if bin_alloc[kk]==0:
                            if j == 0 and kk==0 :
                                skip=1
                            continue

                        #if 1 peak in bin place into array    
                        elif bin_alloc[kk]==1:

                            if j == 0 and kk==0 :
                                did_arr1=did_bin2[kk]
                            elif skip==1:
                                did_arr1=did_bin2[kk]
                                skip=0
                            else:
                                did_arr1=np.vstack((did_arr1,did_bin2[kk]))
                                
                        #if 2 peaks in bin
                        elif bin_alloc[kk]==2:
                            
                            #if kk in del args choose maximum amplitude peak in bin.  
                            if (kk in del_args):
                                if j == 0 and kk==0 :
                                    select_arg=np.argmax(did_bin2[kk][:,1])

                                    did_arr1=did_bin2[kk][select_arg,:]
                                elif skip==1:
                                    select_arg=np.argmax(did_bin2[kk][:,1])

                                    did_arr1=did_bin2[kk][select_arg,:]
                                    skip=0
                                else:
                                    select_arg=np.argmax(did_bin2[kk][:,1])
                                    did_arr1=np.vstack((did_arr1,did_bin2[kk][select_arg,:]))
                            #if bins not in del args add both peaks to final array
                            else:
                                if j == 0 and kk==0 :
                                    did_arr1=did_bin2[kk]
                                    
                                elif skip==1:
                                    did_arr1=did_bin2[kk]
                                    skip=0 
                                else:
                                    did_arr1=np.vstack((did_arr1,did_bin2[kk]))

                 
                
                # no bins are empty
                if np.all(bin_alloc>0):
                    for kk in range(marker_diffs[j]):
                    
                        if bin_alloc[kk]==1:
                            if j == 0 and kk==0 :
                                did_arr1=did_bin2[kk]
                            else:
                                did_arr1=np.vstack((did_arr1,did_bin2[kk]))
                                
                        #if peaks in bin take one with maximum peak amplitude                        
                        elif bin_alloc[kk]==2:
                            if j == 0 and kk==0 :
                                select_arg=np.argmax(did_bin2[kk][:,1])
                        
                                did_arr1=did_bin2[kk][select_arg,:]
                            else:
                                select_arg=np.argmax(did_bin2[kk][:,1])
                                did_arr1=np.vstack((did_arr1,did_bin2[kk][select_arg,:]))
                    
    
            #if less peaks in space than expected
            elif np.sum(bin_alloc)<marker_diffs[j]:
            
                #determine which empty bins to keep
                if np.any(bin_alloc==2):
                
                    #calculate bins with two peaks in
                    args2=np.where(bin_alloc==2)
                    
                    args2 = args2[0]
                    
                    
                    #calculate bins with no peaks in
                    args0=np.where(bin_alloc==0)
                    
                    args0 = args0[0]

                    #calculate the sum of distances for zero bins
                    dist_arr=distance_determiner(args0,args2)
                    
                    
                    dist_args=np.argsort(dist_arr[::-1])
                    
                    #organise args of zer bins based on descending distances
                    del_args=args0[dist_args]
                    
                    
                    #calculate difference between number of bins containing 2 and zero peaks
                    num2=np.count_nonzero(bin_alloc==2)
                    
                    num0=np.count_nonzero(bin_alloc==0)
                    
                    
                    
                    diff5=num0-num2
                    del_args=del_args[:diff5]

                    
                    
                    for kk in range(marker_diffs[j]):
                        
                        #if no peaks in bin and argument not in del args continue
                        if bin_alloc[kk]==0:
                            if (kk not in del_args):
                                if j == 0 and kk==0 :
                                    skip=1
                                continue
                        #else if no bins in peak add zero peak to list
                            else:
                                if j == 0 and kk==0 :
                                    did_arr1=np.array(((start+nuc_sep*(kk-0.5)),0))
                                elif skip==1:
                                    did_arr1=did_bin2[kk]
                                    skip=0
                                else:
                                    did_arr1=np.vstack((did_arr1,np.array(((start+nuc_sep*(kk-0.5)),0))))
                        else:
                            if j == 0 and kk==0 :
                                did_arr1=did_bin2[kk]
                            elif skip==1:
                                    did_arr1=did_bin2[kk]
                                    skip=0
                            else:
                                did_arr1=np.vstack((did_arr1,did_bin2[kk]))
                #if no bins have 2 peaks in them
                else:
                    for kk in range(marker_diffs[j]):
                        
                        if bin_alloc[kk]==0:
                                
                            if j == 0 and kk==0 :
                              
                            
                               did_arr1=np.array(((start+nuc_sep*(kk-0.5)),0))
                                
                            
                            else:
                                did_arr1=np.vstack((did_arr1,np.array(((start+nuc_sep*(kk-0.5)),0))))
                        
                        else:
                            if j == 0 and kk==0 :
                                did_arr1=did_bin2[kk]
                      
                                
                            else:
                                did_arr1=np.vstack((did_arr1,did_bin2[kk]))
                            
            
            #add peak list to larger peak list
            if j==0:
                did_arr2=did_arr1
                
            else:
                did_arr2=np.vstack((did_arr2,did_arr1))            
                

            
            
            if j==0:
                peak_list4=did2
            else:
                peak_list4=np.vstack((peak_list4,did2))
                    
        
           
        #add final target peak data to original peak list
        pos_ind=did_arr1[:,0]-data[0][0]
        peaks_trace1['pos']=pos_ind   
        peaks_trace1['amp']=did_arr1[:,1]
        
        
        
        peaks_trace1['NPeak']=len(did_arr1)
                    
       
        
        data_out2.append(peaks_trace1)
    
    return data_out2

def TM_func(x, a,b):

    """
    Function for extrapolation of size marker
    Args:
	x (int): Size marker position
	a (float): Offset variable
	b (float): Scale variable

    Returns:
	(float): Extrapolated size marker seperation
    """

    return  a-b*np.log(x)

def DpartList():

    """
    PartList object generator

    Args:
	None

    Returns:
	(dict): Partlist object

    """

    dPartList={}
    dPartList['binAlloc']=np.array([],dtype='i4')
    dPartList['partRX']=np.array([],dtype='f4')
    dPartList['peakInfo']=np.array([],dtype='f4')
    dPartList['SM']=np.array([],dtype='f4')
    dPartList['MD']=np.array([],dtype='f4')

    return dPartList



def RX_partitioning_replicates_extended(data_arr,ind,perc,Issues=[],Skip=False,extension=0,Cap=None,ll=0):
    
    """
    partition sequence trace using tape measure 
    processes involved:
    
    find tape measure peaks (peak_finder)
    extract sequencing traces between tape measure peaks
    calculate bin widths between peaks
    find peaks within bins - if no peak found assign max value in bin
    return peak amplitude, position and relative nucleotide number dependent on tape measure
    """
    
    #vectors dictationg the marker sizes and differences between markers in the tape measure
    marker_diffs1 = np.array([10, 30, 10, 20, 30, 10, 20, 10, 10, 20, 20, 20, 20, 10, 10, 20, 20, 20, 20, 20])
    marker_diffs = [10,10, 30, 10, 20, 30, 10, 20, 10, 10, 20, 20, 20, 20, 10, 10, 20, 20, 20, 20, 20,10]
    marker_sizes = [50, 60, 90, 100, 120, 150, 160, 180, 190, 200, 220, 240, 260, 280, 290, 300, 320, 340, 360, 380,400]
    marker_diffs=np.diff(marker_sizes)
    #find peaks in tape measure
    peaksTM = peak_finder(data_arr,4,perc,TM=1)
    sm_plotter(data_arr,peaksTM)
    if extension>0:
	    for k in range(len(peaksTM)):
		if k in Issues:
			continue
		TM_diffs=np.diff(peaksTM[k])
		values=np.divide(TM_diffs,marker_diffs/10)
		
		positions=np.arange(1,21)
		if k==0:
			pos_out=positions
			vals_out=values
		else:
			pos_out=np.append(pos_out,positions)
			vals_out=np.append(vals_out,values)
	    popt, pcov = curve_fit(TM_func, pos_out,vals_out,maxfev=8000)

	    extend=np.arange(22,22+extension)
            added_peak_diffs=TM_func(extend,*popt)
	    if len(Issues)>0 and not Skip:
	    	extend2=np.arange(2,22+extension)
	    	added_peak_diffs2=TM_func(extend2,*popt) 
	    #print added_peak_diffs
	   
	    #array of nucleotide position
	    for pp in range(len(extend)):
		marker_diffs=np.append(marker_diffs,10)
		marker_sizes=np.append(marker_sizes,marker_sizes[-1]+10)	   

    #initialise peak binning array
    reduced_peaks = []
    data_out=[] 
    data_out2=[]
    data_out3=[]
    peaks_TM2=[]
    #iterate through data 
    for i in range(len(data_arr)):
        
        
        
        #extract data index
        data = data_arr[i]
        
        #extract peaks
        if extension<0:
		peaks=peaksTM[i][:abs(extension)]
	
	else:
		if i in Issues and not Skip:
        		peaks = peaksTM[i][:1]
			if extension>0:
				for pt in added_peak_diffs2:
					peaks=np.append(peaks,peaks[-1]+pt)
    		else: 
			peaks = peaksTM[i]
                        if extension>0:
                                for pt in added_peak_diffs:
                                        peaks=np.append(peaks,peaks[-1]+pt)
        #initials peak lists
        peak_list = []
        pos_list = []
        bin_list = []
        
        
        if isinstance(data,pd.DataFrame):
                
            #extract data values
            data1 = data.values
            
            
        else:

            #transpose data - fudge factor can be improved
            data1=np.transpose(data)
       
        #extract peak list for thetarget trace
        peaks_trace1 = funcPeakAlign.fPeakList(data1[:,ind],repType = 'Cubic')
        
        #extract average and standard deviations of widths
        p_width = peaks_trace1['averW']
        p_std=peaks_trace1['stdW']
        
        
        #extract indices of peak positions
        p_ind = peaks_trace1['pos']
        
        #remove duplicates of indices
        p_ind =np.unique(p_ind)
        
        #determine peak amp and position
        p_pos = data1[p_ind,0]
    
        p_amp = data1[p_ind,ind]
        
        #find shoulders 
        shoulder_data = shoulder_finder(p_pos,data1,ind,i)
        
        # transpose peak list
        peak_list=np.transpose([p_pos,p_amp])

        #if shoulders found add them to peak list
        if len(shoulder_data)>0:
            peak_list1 = np.vstack((peak_list,shoulder_data))
            ind1=np.argsort(peak_list1[:,0])
            peak_list1=peak_list1[ind1]
        else:
            peak_list1 = peak_list
        
                
        #delete peaks if they are too close to each other
        x = 0
        while x<len(peak_list1)-1:
            
            #find the position and amplitudes of two neighbouring peaks
            p_pos1 = peak_list1[x,0]
            p_pos2 = peak_list1[x+1,0]
            
            p_amp1 = peak_list1[x,1]
            p_amp2 = peak_list1[x+1,1]
            
            
            #calculate difference in position between the two peaks
            diff1 = p_pos2-p_pos1
            
            #if the difference is less than 4
            if diff1<4:
                
                #if peak1 is in original peak list but not peak2 delete peak2     
                if (p_pos1 in peak_list[:,0]) and (p_pos2 not in peak_list[:,0]):
                    del_x = x+1
                    peak_list1 = np.delete(peak_list1,del_x,0)
                    x+=1
                    
                #if peak2 is in original peak list but not peak1 delete peak1    
                elif(p_pos2 in peak_list[:,0]) and (p_pos1 not in peak_list[:,0]):
                    del_x = x
                    peak_list1 = np.delete(peak_list1,del_x,0)
                else:
                    x+=1
                     
                
                
            else:
                x+=1
        
        x=0
        
        
        miss_pos = []
        miss_amp = []
        
        #fill in any positions that have missing peaks due to poor data quality
        
        #read through peak list
        for x in range(len(peak_list1)-1):
            
            #look at the difference in position between neighbouring peaks 
            diff1=peak_list1[x+1,0]-peak_list1[x,0]

            #determine spacing
            spacing_we=2*p_width-4*p_std
            
            #calculate number of missing peaks expected in a space
            miss_peaks=np.floor(diff1/(p_width-2*p_std))
            
            
            #if more that one missing peak space suspected
            if miss_peaks>1:
                

                #fill the missing space with peak data
                for z in range(int(miss_peaks)-1):
                
                
                    #determine position of new peak
                    av_pos=peak_list1[x,0]+(z+1)*int(p_width)
                    
                    #add new position to missing position array
                    miss_pos=np.append(miss_pos,av_pos)

                    #find amplitude of new peak
                    av_amp=data1[data1[:,0]==av_pos,ind]
                    
                    #add to missing amplitude data
                    miss_amp = np.append(miss_amp,av_amp)

        #create array for missing peak information
        miss_arr = np.transpose([miss_pos,miss_amp])
        
        
        #add the missing peaks to the new peak array
        if len(miss_arr)>0:
            peak_list2=np.vstack((peak_list1,miss_arr))
            ind2=np.argsort(peak_list2[:,0])
            peak_list2=peak_list2[ind2]
        else:
            peak_list2=peak_list1

        
       
        
        p_pos2=[]
        p_amp2=[]
        
        did_bin=[]
        bin_counter=[]
        bin_pos=[]

        did_bin3=[]
        bin_alloc2=[]
        for j in range(len(peaks)-1):
            #extract TM peaks that define the trace region of interest
            tm_bins=0
            
            #extract adjacent peaks in TM
            start = peaks[j]
            end = peaks[j+1]
            
            #calculate distance between TM peaks
            diff=end-start
            #print 'cow'
	    #print start
	    #print end
            # find the nearest target peak downwind of start TM peak
            start_ind = find_nearest_ind(peak_list2[:,0],start)
            # find the nearest target peak downwind of start TM peak
            end_ind = find_nearest_ind(peak_list2[:,0],end)
            
            #extract target peak list for space between TM peaks
            peak_list3=peak_list2[start_ind:end_ind,:]
            
            
            #add clipped target peaks to list
            did_bin.append(peak_list3)
                
              
        #work through clipped target peak list  
        for j in range(len(did_bin)):
            did=did_bin[j]
            did2=np.reshape(did,(-1,2))
          
                
            #set TM peaks at the beginning and end of the clipped target peak list    
            start = peaks[j]
            end=peaks[j+1]

            #determine the difference between end and start
            diff5=end-start

            #calculate seperation of nucleotides
            nuc_sep=diff5/marker_diffs[j]
            
            
            bin_alloc=[]
            
            did_bin2=[]
            
            
            #determine peaks in different marker bins
            for kk in range(marker_diffs[j]):
                
                did3=did2[did2[:,0]>=start+(kk-1)*nuc_sep,:]
                did3=did3[did3[:,0]<start+(kk)*nuc_sep,:]

                #convert bins with 4 peaks in to 2 peaks and remove surpluses
                if len(did3)>2:
                    
                    args2=np.argsort(did3[:,1][::-1])
                    did3=did3[args2[:1],:]
                    
                bin_alloc=np.append(bin_alloc,len(did3))
                
                did_bin2.append(did3)
            bin_alloc2.append(bin_alloc)    
                
            skip=0
            
            #reassignment of data to evenly spread between bins 
            
            #if all bins have one nucleotide in assign each to respective bin 
            if  np.all(bin_alloc==1):
                for kk in range(marker_diffs[j]):
                    if kk==0 :
                        did_arr1=did_bin2[kk]
                    else:
                        did_arr1=np.vstack((did_arr1,did_bin2[kk]))
                        
                               
            #if the total number of bins in the system is equal to the predicted amount 
            elif np.sum(bin_alloc)==marker_diffs[j]:
                
                #if a bin has no peaks in it skip
                for kk in range(marker_diffs[j]):
                    if bin_alloc[kk]==0:
                        if  kk==0 :
                            skip=1
                        continue
                    #else add the data to the new list    
                    else:
                        if  kk==0 :
                            did_arr1=did_bin2[kk]
                        elif skip==1:
                            did_arr1=did_bin2[kk]
                            skip=0
                        else:
                            did_arr1=np.vstack((did_arr1,did_bin2[kk]))
                    
                    
            #if number of peaks is greater than expected
            elif np.sum(bin_alloc)>marker_diffs[j]:
                diff4=np.sum(bin_alloc)-marker_diffs[j]
                
                #if any bins exist with no peaks in them
                if np.any(bin_alloc==0):
                    
                    #extract bin positions with two peaks in them
                    args2=np.where(bin_alloc==2)
                    
                    args2 = np.array(args2)[0]
                    
                    #extract bin positions with no peaks in them 
                    args0=np.where(bin_alloc==0)
                    
                    args0 = np.array(args0)[0]

                    #calculate the sum of distances between args2 and args0
                    dist_arr=distance_determiner(args2,args0)

                    #sort the arguments based on descending sum off distance
                    dist_args=np.argsort(dist_arr[::-1])

                    #order args2 based on Sum of Distance
                    del_args=args2[dist_args]
                    
                    #calculate difference of numbers of bins containing 2 and 0 peaks
                    num2=np.count_nonzero(bin_alloc==2)
                    
                    num0=np.count_nonzero(bin_alloc==0)
                    
              
                    diff5=num2-num0
                    
                    #set the indices of the bins with 2 peaks that should be averaged
                    del_args=del_args[:(diff5)]
                    
                    
                    for kk in range(marker_diffs[j]):
                        #if no peaks in bin continue
                        if bin_alloc[kk]==0:
                            if kk==0 :
                                skip=1
                            continue

                        #if 1 peak in bin place into array    
                        elif bin_alloc[kk]==1:

                            if kk==0 :
                                did_arr1=did_bin2[kk]
                            elif skip==1:
                                did_arr1=did_bin2[kk]
                                skip=0
                            else:
                                did_arr1=np.vstack((did_arr1,did_bin2[kk]))
                                
                        #if 2 peaks in bin
                        elif bin_alloc[kk]==2:
                            
                            #if kk in del args choose maximum amplitude peak in bin.  
                            if (kk in del_args):
                                if kk==0 :
                                    select_arg=np.argmax(did_bin2[kk][:,1])

                                    did_arr1=did_bin2[kk][select_arg,:]
                                elif skip==1:
                                    select_arg=np.argmax(did_bin2[kk][:,1])

                                    did_arr1=did_bin2[kk][select_arg,:]
                                    skip=0
                                else:
                                    select_arg=np.argmax(did_bin2[kk][:,1])
                                    did_arr1=np.vstack((did_arr1,did_bin2[kk][select_arg,:]))
                            #if bins not in del args add both peaks to final array
                            else:
                                if kk==0 :
                                    did_arr1=did_bin2[kk]
                                    
                                elif skip==1:
                                    did_arr1=did_bin2[kk]
                                    skip=0 
                                else:
                                    did_arr1=np.vstack((did_arr1,did_bin2[kk]))

                 
                
                # no bins are empty
                if np.all(bin_alloc>0):
                    for kk in range(marker_diffs[j]):
                    
                        if bin_alloc[kk]==1:
                            if kk==0 :
                                did_arr1=did_bin2[kk]
                            else:
                                did_arr1=np.vstack((did_arr1,did_bin2[kk]))
                                
                        #if peaks in bin take one with maximum peak amplitude                        
                        elif bin_alloc[kk]==2:
                            if kk==0 :
                                select_arg=np.argmax(did_bin2[kk][:,1])
                        
                                did_arr1=did_bin2[kk][select_arg,:]
                            else:
                                select_arg=np.argmax(did_bin2[kk][:,1])
                                did_arr1=np.vstack((did_arr1,did_bin2[kk][select_arg,:]))
            
    
            #if less peaks in space than expected
            elif np.sum(bin_alloc)<marker_diffs[j]:
            
                #determine which empty bins to keep
                if np.any(bin_alloc==2):
                
                    #calculate bins with two peaks in
                    args2=np.where(bin_alloc==2)
                    
                    args2 = args2[0]
                    
                    
                    #calculate bins with no peaks in
                    args0=np.where(bin_alloc==0)
                    
                    args0 = args0[0]

                    #calculate the sum of distances for zero bins
                    dist_arr=distance_determiner(args0,args2)
                    
                    
                    dist_args=np.argsort(dist_arr[::-1])
                    
                    #organise args of zer bins based on descending distances
                    del_args=args0[dist_args]
                    
                    
                    #calculate difference between number of bins containing 2 and zero peaks
                    num2=np.count_nonzero(bin_alloc==2)
                    
                    num0=np.count_nonzero(bin_alloc==0)
                    
                    
                    
                    diff5=num0-num2
                    del_args=del_args[:diff5]

                    
                    
                    for kk in range(marker_diffs[j]):
                        
                        #if no peaks in bin and argument not in del args continue
                        if bin_alloc[kk]==0:
                            if (kk not in del_args):
                                if kk==0 :
                                    skip=1
                                continue
                        #else if no bins in peak add zero peak to list
                            else:
                                if kk==0 :
                                    did_arr1=np.array(((start+nuc_sep*(kk-0.5)),np.nan))
                                elif skip==1:
                                    did_arr1=did_bin2[kk]
                                    skip=0
                                else:
                                    did_arr1=np.vstack((did_arr1,np.array(((start+nuc_sep*(kk-0.5)),np.nan))))
                        else:
                            if kk==0 :
                                did_arr1=did_bin2[kk]
                            elif skip==1:
                                    did_arr1=did_bin2[kk]
                                    skip=0
                            else:
                                did_arr1=np.vstack((did_arr1,did_bin2[kk]))
                #if no bins have 2 peaks in them
                else:
                    for kk in range(marker_diffs[j]):
                        
                        if bin_alloc[kk]==0:
                                
                            if kk==0 :
                              
                            
                               did_arr1=np.array(((start+nuc_sep*(kk-0.5)),np.nan))
                                
                            
                            else:
                                did_arr1=np.vstack((did_arr1,np.array(((start+nuc_sep*(kk-0.5)),np.nan))))
                        
                        else:
                            if kk==0 :
                                did_arr1=did_bin2[kk]
                      
                                
                            else:
                                did_arr1=np.vstack((did_arr1,did_bin2[kk]))
                            
            did_bin3.append(did_arr1)
            #add peak list to larger peak list
            if j==0:
                did_arr2=did_arr1
                
            else:
                did_arr2=np.vstack((did_arr2,did_arr1))            
                

            
            
            if j==0:
                peak_list4=did2
            else:
                peak_list4=np.vstack((peak_list4,did2))
                    
        
           
        #add final target peak data to original peak list
        
        
        peaks_trace1['pos']=did_arr1[:,0].astype(int)
        peaks_trace1['amp']=did_arr1[:,1]
        
        peaks_trace1['NPeak']=len(did_arr1)
        
        data_out3.append(peaks_trace1)
                    
        data_out.append(bin_alloc2)
        
        data_out2.append(did_bin3)
        peaks_TM2.append(peaks)

    partlist=DpartList()
    
    partlist['binAlloc']=data_out

    partlist['partRX']=data_out2
    partlist['peakInfo']=data_out3
    partlist['SM']=peaks_TM2
    partlist['MD']=marker_diffs
       
    return partlist



def RX_partition_realignment_extended(partition,inds,data_arr1,fl=None,corr_b=0.7,inspect=1000,Cap=None,perc=0.25,tm=0,ll=0,marker_set='ROX'):
    

    """
    Realign partitioned footprinting data based on up to 3 replicates with size marker region extended by 10 nucleotides on either side

    Args:
        partition (list): Partitioned footprinting data produced by RX_partitioning_replicates
        data_arr (list): All of the datasets in the ensemble
        peak_info1 (list): Peaklists of footprinting data for all datasets in the ensemble
        inds (int): Indices in the ensemble indicating the location of the replicates
        data_arr1 (list): All of the datasets in the ensemble
        Peaks_TM1 (list): peak positions for size marker traces of the datasets in the ensemble
        corr_b (float): Correlation value indicating the minimum pairwise pearson correlation considered in realignment process
        Cap: Maximum peak intensity considered for peak finding function
        perc (float): Ratio of maximum peak intensity for cutoff
        inspect (int): Specify a particular dataset that you want to inspect



    Return:
        (list): The aligned partitioned footprinted stored as a peakList object

    """
   
    #initialise marker diffs
    if marker_set.upper()=='ROX':
        
        marker_sizes = np.array([50, 60, 90, 100, 120, 150, 160, 180, 190, 200, 220, 240, 260, 280, 290, 300, 320, 340, 360, 380,400])
    
    elif marker_set.upper()=='LIZ':
         
        marker_sizes = np.array([20, 40, 60, 80, 100, 114, 120, 140, 160, 180, 200, 214, 220, 240, 250, 260, 280, 300, 314, 320, 340, 360, 380, 400, 414, 420, 440, 460, 480, 500, 514, 520, 540, 560, 580, 600])
        marker_sizes=marker_sizes[2:]
    else:
        print 'unknown marker set!'
        quit()
        
    marker_diffs=deepcopy(partition['MD'])
    
    parts=[]
    bins=[]
    TM=[]    
    new_arr=[]
    #get letters for replicates
    labs=map(chr, range(65, (65+len(inds))))
    
    partition2=deepcopy(partition['partRX'])
    peak_info=deepcopy(partition['peakInfo'])
    data_arr=deepcopy(data_arr1)
    bin_alloc=deepcopy(partition['binAlloc'])
    peaksTM=deepcopy(partition['SM'])
    
    
    if fl!=None:
        sam_funcs.sm_plotter(data_arr,peaksTM,fl)    

    #extract partitions and bins from lists
    for i in range(len(inds)):
        parts.append(deepcopy(partition2[inds[i]]))
        bins.append(deepcopy(bin_alloc[inds[i]]))
    	TM.append(peaksTM[inds[i]])
    
    for i in range(len(parts[0])):
       
        bins=marker_diffs[i]
       
        spaces=[]
        spaces_c=[]
        #print i
	#print parts
        
        #remove intensity nan data
        for j in range(len(inds)):
            if fl!=None:
		print fl[inds[j]]    
           
            
            spaces.append(parts[j][i])
            spaces_c.append(parts[j][i][~np.isnan(parts[j][i][:,1]),:])
        
        
        bars=[]
        lens=[]
        
        #barcode generation
        for j in range(len(inds)):
            
            bar = barcode_generator(spaces_c[j][:,1])
            lens=np.append(lens,len(bar))
            bars.append(bar)
            
       
            
        #if less peaks are observed than expected
        if np.any(lens<bins):     
            
            #if two replicates
            if len(inds)==2:
                
                #perform alignment
                new_arrays,lens1=nw_align(bars,bins,labs)
		#print lens1
		bar_lens=[]

                #if no realigned sequences have correct length
                if np.all(np.array(lens1)==0):
                    
                    #set nans in original sublists to zero 
                    spaceA_new=np.nan_to_num(spaces[0])
                    spaces_c[0]=spaceA_new
                    
                    #re barcode using new data
                    bars[0]=barcode_generator(spaceA_new[:,1])

                    #nw_align on new barcodes
                    new_arrays,lens2=nw_align(bars,bins,labs)
                    
                    #determine correct insertionsi
                    new_traces=trace_align2(spaces_c,new_arrays,bins)

                # if aligned barcodes with correct length observed.    
                else:
		    for fg in range(len(new_arrays)):

                        bar_lens.append(len(new_arrays[fg][0][0]))
                    bar_lens=np.array(bar_lens)
		    if np.any(bar_lens<bins):

                    	#determine correct alignments
                    	new_traces=trace_align2(spaces_c,new_arrays,bins,bar_lens,all_short=True)
		    else: 
		    	new_traces=trace_align2(spaces_c,new_arrays,bins)
            #if 3 replicates used
            else:
		#print bars
                #perform alignment of barcodes
                new_arrays,lens1=nw_align(bars,bins,labs)
		#print 'new_array'
                #print new_arrays 
                
		bar_lens=[]
		for fg in range(len(new_arrays)):

                        bar_lens.append(len(new_arrays[fg][0][0]))
		#if no alignments have correct length
                if np.all(np.array(bar_lens)!=bins):
                    
                    #convert first replicates original sublists nan to 0
                    spaceA_new=np.nan_to_num(spaces[0])
                    spaces_c[0]=spaceA_new
                    
                    #reproduce barcodes
                    bars[0]=barcode_generator(spaceA_new[:,1])
                    
                    #realign barcodes
                    new_arrays,lens2=nw_align(bars,bins,labs)
                    bar_lens1=[]
                    for fg in range(len(new_arrays)):

                    	bar_lens1.append(len(new_arrays[fg][0][0]))
		    bar_lens1=np.array(bar_lens1)
		    #print 'caraboa'
		    #print bar_lens1
		    #print lens2
                    #determine corret insertions
       		    zero_len = np.where(lens2==0)[0]

		    if np.any(np.array(bar_lens1)<bins):
                        new_traces=trace_align3(spaces_c,new_arrays,bins,bar_lens1,some_short=True)

                    elif zero_len>0:

			zero_len = np.where(lens2==0)[0]
                        spaceA_new=np.nan_to_num(spaces[zero_len[0]])
                        bars_new=barcode_generator(spaceA_new[:,1])
                        new_arrs=np.transpose([bars_new,labs[zero_len[0]],5.0])
                        
                        spaces_c[zero_len[0]]=spaceA_new
                        new_arrays[zero_len[0]]=np.array([new_arrs])
                    
                    else:
                    	new_traces=trace_align3(spaces_c,new_arrays,bins)
                    
                #if aligned sequences have correct length
                elif np.all(np.array(bar_lens)==bins):

                    #determine correct insertion 
                    new_traces=trace_align3(spaces_c,new_arrays,bins)
		elif np.any(np.array(bar_lens)==bins): 
                    new_traces=trace_align3(spaces_c,new_arrays,bins,bar_lens,some_short=True)
                else:    
                    for q,arr in enumerate(new_arrays):
                        #if one of the barcodes doesn't have correct length after alignment
                        if lens1[q]==0 and np.all(np.delete(lens1,q)>0):
                            
                            #remove data for replicate with short barcode
                            new_arrays1=deepcopy(new_arrays)
                            del new_arrays1[q]
                            spaces1=deepcopy(spaces_c)
                            del spaces1[q]
                            
                            #determine correct insertions for reduced 
                            new_traces1=trace_align2(spaces1,new_arrays1,bins)
                            #add short data back in
                            new_traces1.insert(q,spaces_c[q])
                            new_bars=[]
                            
                            #recalculate barcodes
                            for traces in new_traces1:
                                new_bars = np.append(new_bars,barcode_generator(traces[:,1]))
                            
                            new_arrays2,lens2=nw_align(new_bars,bins,labs)
                            
                            if np.all(np.array(lens2)==0):
                                new_traces=deepcopy(spaces)
                            else:
                                new_traces=trace_align3(new_traces1,new_arrays2,bins)
                    
        else:
            new_traces=deepcopy(spaces)
        
        if len(inds)==3:
            starts=[TM[0][i],TM[1][i],TM[2][i]]
            ends=[TM[0][i+1],TM[1][i+1],TM[2][i+1]]
            #print 'bull'
	    #print new_traces
	    #print i
	    #print bins
	    #print marker_diffs
            if i==inspect:
                new_traces1=final_realignment(new_traces,data_arr,starts,ends,bins,inds,labs,corr=0.95,corr_b=0.7,inspect=True)
            
            else:
                new_traces1=final_realignment(new_traces,data_arr,starts,ends,bins,inds,labs,corr=0.95,corr_b=corr_b)
            
        else:
            new_traces1=deepcopy(new_traces)        
        if i==0:
            for j in range(len(inds)):
                new_arr.append(new_traces1[j])
                
        else:
            for j in range(len(inds)):

                new_arr[j]=np.vstack((new_arr[j],new_traces1[j]))

    peak_infos=[]
    
    #add data to peeak list. 
    for j in range(len(inds)):
        dataA=data_arr[inds[j]]
        
        
        posA=nan_remover(new_arr[j][:,0],TM[j][0],TM[j][-1])
        pos_indA=np.in1d(dataA[0],posA.astype(int)).nonzero()[0]
        peak_infoA=peak_info[inds[j]]
        peak_infoA['amp']=np.nan_to_num(new_arr[j][:,1])
        peak_infoA['pos']=posA-dataA[0][0]
        peak_infoA['pos']=peak_infoA['pos'].astype(int)
        peak_infoA['NPeak']=len(new_arr[j])
        peak_infos.append(peak_infoA)
    
    
    
    return peak_infos



def final_realignment(new_traces2,data_arr1,starts,ends,bins1,inds1,labs,corr=0.95,corr_b=0.7,inspect=False):
    
    """
    Final realignment wrapper function (only called by RX_realignment functions)

    Args:
        new_traces2 (list): List containing partitioned footprinting data package between size marker positions
        data_arr1 (list): All of the datasets in the ensemble
        starts (list): Size marker peaks at the start of the packaged region for each replicate
        ends (list): Size marker peaks at the end of the packaged region for each replicate
        inds1 (int): Indices in the ensemble indicating the location of the replicates
        data_arr1 (list): All of the datasets in the ensemble
        labs (list): Replicate labels
        corr(float): Starting correlation value for alignment process
        corr_b (float): Correlation value indicating the minimum pairwise pearson correlation considered in realignment process
        inspect (bool): Specify whether inspection of data is required


    Return:
        (list): The aligned partitioned footprinted data for each replicate in the specific part of size marker region.

    """

    #deepcopy data    
    new_traces=deepcopy(new_traces2)
    data_arr=deepcopy(data_arr1)
    bins=deepcopy(bins1)
    inds=deepcopy(inds1)
    #print 'fox'
    #print bins


    #split bins large than 10 nulceotides into 10 nucleotide bins
    if bins>19:
        split=float(bins)/10
        
        

        for q in range(int(np.ceil(split))):
            nt_split=[]
	    #print bins
            #print starts
	    #print ends
	    #print new_traces
            new_starts=[]
            new_ends=[]
	    #print q
            if q==np.ceil(split)-1 and split<np.ceil(split):
   		#print 'mouse'                
                for j in range(len(inds)):
                    #create new starts and ends
                    diff=ends[j]-starts[j]
                    new_starts.append(starts[j]+q*10*diff/float(bins))
                    new_ends.append(ends[j])
                    nt_split.append(new_traces[j][(q*10):,:])
                    new_bins=bins-q*10
                    
            else:
                #print 'rat'
                for j in range(len(inds)):
                    #create new starts and ends
                    diff=ends[j]-starts[j]
                    new_starts.append(starts[j]+q*10*diff/float(bins))
                    new_ends.append(starts[j]+(q+1)*10*diff/float(bins))
		    #print 'cat'
		    #print len(new_traces[j])
                    nt_split.append(new_traces[j][(q*10):(q+1)*10,:])
                    new_bins=10
                    #print nt_split  
	     
            
            #carry out subprocess
            if  inspect and q==2:
                new_traces2=fr_subprocess(nt_split,data_arr1,new_starts,new_ends,new_bins,inds,labs,corr=0.95,corr_b=0.7,inspect=True)
            else:
		#print 'nt_split'
		#print nt_split
                new_traces2=fr_subprocess(nt_split,data_arr1,new_starts,new_ends,new_bins,inds,labs,corr=0.95,corr_b=0.7,inspect=False)
            if q==0:
                new_traces1=deepcopy(new_traces2)
            else:
                for j in range(len(inds)):
                    new_traces1[j]=np.vstack((new_traces1[j],deepcopy(new_traces2[j])))
            
    else: 
        new_traces1=fr_subprocess(new_traces,data_arr,starts,ends,bins,inds,labs,corr=0.95,corr_b=0.7,inspect=False)
    return new_traces1

def fr_subprocess(new_traces2,data_arr1,starts,ends,bins1,inds1,labs,corr=0.95,corr_b=0.7,inspect=False):

    """
    Final realignment sub process (only called by final_realignment wrapper function)

    Args:
        new_traces2 (list): List containing partitioned footprinting data package between size marker positions
        data_arr1 (list): All of the datasets in the ensemble
        starts (list): Size marker peaks at the start of the packaged region for each replicate
        ends (list): Size marker peaks at the end of the packaged region for each replicate
        bins1 (array): number of bins (i.e. number of nucleotides between size markers) expected in each part of the size marker set.
        inds1 (int): Indices in the ensemble  indicating the location of the replicates
        data_arr1 (list): All of the datasets in the ensemble
        labs (list): Replicate labels
        corr(float): Starting correlation value for alignment process
        corr_b (float): Correlation value indicating the minimum pairwise pearson correlation considered in realignment process
        inspect (bool): Specify whether inspection of data is required


    Return:
        (list): The aligned partitioned footprinted data for each replicate in the specific part of size marker region.

    """

    
    hispace=[]
    p_inds=[]
    bars1=[]
    new_traces=deepcopy(new_traces2)
    data_arr=deepcopy(data_arr1)
    bins=deepcopy(bins1)
    inds=deepcopy(inds1)
    for j in range(len(new_traces)):
	#print 'goose'
	#print new_traces[j]
        bars1.append(coarse_grainer(new_traces[j][:,1]))
    cg_corr=np.zeros([len(bars1),len(bars1)])
    for t in range(len(bars1)):
        for j in range(len(bars1)):
            cg_corr[t,j]=pearsonr(bars1[t],bars1[j])[0]


    mat=correl_assessor(new_traces,1)
    if inspect:
        plt.plot(np.arange(len(new_traces[0])),new_traces[0][:,1])
        plt.plot(np.arange(len(new_traces[0])),new_traces[1][:,1])
        plt.plot(np.arange(len(new_traces[0])),new_traces[2][:,1])
        plt.show()
        
    space_cors=np.transpose(np.where(cg_corr>corr))


    space_cors2=space_cors[space_cors[:,0]>space_cors[:,1],:]



    if len(space_cors2)==1:

        for j in range(len(inds)):
            if j not in space_cors2:
                hispace.append(new_traces[j][new_traces[j][:,1]>np.max(new_traces[j][:,1])*0.33,:]) 
                p_inds=np.transpose(np.where(new_traces[j][:,1]>np.max(new_traces[j][:,1])*0.33))
            else:
                hispace.append(new_traces[j])
        
        bars=[]
        for j in range(len(hispace)):
            bars.append(barcode_generator(hispace[j][:,1]))
        
        new_arrays1,lens1=nw_align(bars,bins,labs)
        
        if np.all(np.array(lens1)==0):
            return new_traces
        else:
            
            new_traces1=trace_align3_v2(hispace,new_arrays1,data_arr,inds,bins,starts,ends,inspect=inspect,peak_inds=p_inds)
            
            if new_traces1=='error':
                return new_traces
    elif corr>=corr_b:
        new_traces1=fr_subprocess(new_traces,data_arr,starts,ends,bins,inds,labs,corr=corr-0.05,corr_b=0.7,inspect=inspect)
    else:
        new_traces1=deepcopy(new_traces)
    
 
    return new_traces1
    

def trace_align3(spaces,new_array,bins,bar_lens=None,all_short=False,some_short=False):
    

    """
    Partitioned data alignment process for 3 replicates (secondary function)

    Args:
        spaces (list): unaligned partitioned footprinting data for each replicate
        new_array (list): List containing new NW alignments for each replicate
        bar_lens (list): lens of barcodes for the different replicates
        all_short (bool): if all aligned are still short
        some_short (bool): if some of the aligned barcodes are short

    Return:
        (list): The aligned partitioned footprinted data for each replicate in the specific part of size marker region.
    """


    cov_arr=[]
    insert_indsA=[]
    insert_indsB=[]
    insert_indsC=[]
    
    new_array1=deepcopy(new_array)
    spaces1=deepcopy(spaces)
     
    #barcodes extract
    new_arrayA=new_array1[0]
    new_arrayB=new_array1[1]
    new_arrayC=new_array1[2]
    einsA=False
    einsB=False
    einsC=False
    diffA=0
    diffB=0
    diffC=0
    #print 'rooster'
    #print bar_lens
    if some_short:
    	if bar_lens[0]<bins:
		einsA=True
		diffA=bar_lens[0]
	if bar_lens[1]<bins:
		einsB=True
		diffB=bar_lens[1]
	if bar_lens[2]<bins:
		einsC=True
		diffC=bar_lens[2]

    #peak data extract
    spaceA_c=spaces1[0]
    spaceB_c=spaces1[1]
    spaceC_c=spaces1[2]
    
    #iterate thorugh combinations
    for j in range(len(new_arrayA)):
        #find insertion points
        insertA=indexes(new_arrayA[j,0],bins,diffA,eins=einsA)
        
        #find correction length
        correctA=np.arange(len(insertA))

        #correct insertions
        insertA_f=np.subtract(insertA,correctA)
        
        #create new intensity profile
        spaceA_n=np.insert(spaceA_c,insertA_f.astype(int),np.array((np.nan,0.1)),axis=0)
	#print 'spaceA'
	#print len(spaceA_n)
        for k in range(len(new_arrayB)):
	    #print new_arrayB[k,0]
            insertB=indexes(new_arrayB[k,0],bins,diffB,eins=einsB)
	    #print insertB
            correctB=np.arange(len(insertB))

            insertB_f=np.subtract(insertB,correctB)
            spaceB_n=np.insert(spaceB_c,insertB_f.astype(int),np.array((np.nan,0.1)),axis=0)
	    #print 'spaceB'
	    #print len(spaceB_n)
            for l in range(len(new_arrayC)):
                insertC=indexes(new_arrayC[l,0],bins,diffC,eins=einsC)
                correctC=np.arange(len(insertC))

                insertC_f=np.subtract(insertC,correctC)
                
                spaceC_n=np.insert(spaceC_c,insertC_f.astype(int),np.array((np.nan,0.1)),axis=0)
                #print 'spaceC'
		#print len(spaceC_n)
                #calculate PCCs between possible pairwise combinations 
                
                cgA=coarse_grainer(spaceA_n[:,1])
                cgB=coarse_grainer(spaceB_n[:,1])
                cgC=coarse_grainer(spaceC_n[:,1])

                
                cor_arrA=np.nan_to_num(spaceA_n[:,1])
                cor_arrB=np.nan_to_num(spaceB_n[:,1])
                cor_arrC=np.nan_to_num(spaceC_n[:,1])
		#print 'ABC'
                #print cor_arrA
		#print cor_arrB
		#print cor_arrC 
                normA=spaceA_n[:,1]/np.max(spaceA_n[:,1])
                normB=spaceB_n[:,1]/np.max(spaceB_n[:,1])
                normC=spaceC_n[:,1]/np.max(spaceC_n[:,1])
                #print len(spaceA_n[:,1])
		#print len(spaceB_n[:,1])
                cov_AB=get_cov(spaceA_n[:,1],spaceB_n[:,1])
                
                cov_BC=get_cov(spaceB_n[:,1],spaceC_n[:,1])
                cov_AC=get_cov(spaceA_n[:,1],spaceC_n[:,1])
                         
                #summate PCCs
                cov_arr=np.append(cov_arr,cov_AB+cov_BC+cov_AC)
                #print cov_arr
                #bin new profiles
                insert_indsA.append(spaceA_n)
                insert_indsB.append(spaceB_n)
                insert_indsC.append(spaceC_n)
    
    #find max summate PCC
    #print cov_arr
    max_cor_ind=np.argmax(cov_arr)
    

    #return intensity profiles for max sum PCC
     
    newA = insert_indsA[max_cor_ind]
    newB = insert_indsB[max_cor_ind]
    newC = insert_indsC[max_cor_ind]
 
    return [newA,newB,newC]



def trace_align3_v2(spaces,new_array,data_arr,inds,bins, starts, ends,peak_inds,inspect=False):
    
    """
    Partitioned data alignment process for 3 replicates (secondary function)

    Args:
        spaces (list): Unaligned partitioned footprinting data for each replicate
        new_array (list): List containing new NW alignments for each replicate
        data_arr (list): All the preprocessed data in the ensemble
        inds (list): Indices indicating the positions of the replicates in the ensemble list
        bins (array): number of bins (i.e. number of nucleotides between size markers) expected in each part of the size marker set
        starts (list): Size marker peaks at the start of the packaged region for each replicate
        ends (list): Size marker peaks at the end of the packaged region for each replicate
        peak_inds (int): Position of bins with assigned peaks before realignment

    Returns:
        (list): The aligned partitioned footprinted data for each replicate in the specific part of size marker region.
    """
    
    cov_arr=[]
    insert_indsA=[]
    insert_indsB=[]
    insert_indsC=[]
    cg_arrA=[]
    cg_arrB=[]
    cg_arrC=[]
    peak_travel=[]
    
    new_array1=deepcopy(new_array)
    spaces1=deepcopy(spaces)
    
    #barcodes extract
    new_arrayA=new_array1[0]
    new_arrayB=new_array1[1]
    new_arrayC=new_array1[2]

    #peak data extract
    spaceA_c=spaces1[0]
    spaceB_c=spaces1[1]
    spaceC_c=spaces1[2]
    
 
    startA=starts[0]
    startB=starts[1]
    startC=starts[2]
    
    endA=ends[0]
    endB=ends[1]
    endC=ends[2]
    
    dataA=data_arr[inds[0]]
    dataB=data_arr[inds[1]]
    dataC=data_arr[inds[2]]
    
    #iterate thorugh combinations
    for j in range(len(new_arrayA)):
        #find insertion points
        insertA=notindexes(new_arrayA[j,0])
      
        if len(insertA)<1:
            spaceA_n=spaceA_c
            
        else:
            if len(insertA)<bins:
                peak_distances=travel_determiner(insertA,peak_inds)
                
            spaceA_n=np.zeros([bins,2])
            for t in range(len(insertA)):
                spaceA_n[insertA[t],:]=spaceA_c[t,:]
            if insertA[0]!=0:
                diff=insertA[0]
                
                space=spaceA_c[0,0]-startA
             
                
                new_ns=space/diff
                
                for q in range(int(diff)):
                        new_pos=startA+(q)*new_ns
                        datanew=dataA[1][(new_pos+0.5*new_ns)>dataA[0]]
                        
                        datanp=dataA[0][(new_pos+0.5*new_ns)>dataA[0]]
                        
                        datanew=datanew[datanp>(new_pos-0.5*new_ns)]
                        new_value=np.nanmean(datanew)
                        
                    
                        spaceA_n[q,0]=new_pos
                        spaceA_n[q,1]=new_value
                           
            if insertA[-1]!=bins-1:
                diff=bins-insertA[-1]
                
                space=endA-spaceA_c[-1,0]
                
                new_ns=space/diff
                
                for q in range(1,int(diff)):
                        new_pos=spaceA_c[-1,0]+(q)*new_ns
                        datanew=dataA[1][(new_pos+0.5*new_ns)>dataA[0]]
                        
                        datanp=dataA[0][(new_pos+0.5*new_ns)>dataA[0]]
                        
                        datanew=datanew[datanp>(new_pos-0.5*new_ns)]
                        new_value=np.nanmean(datanew)
                        
                        spaceA_n[insertA[-1]+q,0]=new_pos
                        spaceA_n[insertA[-1]+q,1]=new_value
                        
            for t in range(len(insertA)-1):
                diff=insertA[t+1]-insertA[t]
                space=spaceA_c[t+1,0]-spaceA_c[t,0]
                new_ns=space/diff
                if diff>1:
                    for q in range(1,int(diff)-1):
                        new_pos=spaceA_c[t,0]+(q)*new_ns
                        datanew=dataA[1][(new_pos+0.5*new_ns)>dataA[0]]
                        
                        datanp=dataA[0][(new_pos+0.5*new_ns)>dataA[0]]
                        
                        datanew=datanew[datanp>(new_pos-0.5*new_ns)]
                        
                        new_value=np.nanmean(datanew)
                        
                        spaceA_n[int(insertA[t])+q,0]=new_pos
                        spaceA_n[int(insertA[t])+q,1]=new_value
               

        for k in range(len(new_arrayB)):
            
            
            insertB=notindexes(new_arrayB[k,0])

        
            if len(insertB)<1:
                spaceB_n=spaceB_c
            else:
                if len(insertB)<bins:
                    peak_distances=travel_determiner(insertB,peak_inds)
                spaceB_n=np.zeros([bins,2])
                for t in range(len(insertB)):
                    spaceB_n[insertB[t],:]=spaceB_c[t,:]
                    
                if insertB[0]!=0:
                    diff=insertB[0]
                
                    space=spaceB_c[0,0]-startB
                
                    new_ns=space/diff
                
                    for q in range(int(diff)):
                        new_pos=startB+(q+0.5)*new_ns
                        datanew=dataB[1][(new_pos+0.5*new_ns)>dataB[0]]
                        
                        
                        
                        datanp=dataB[0][(new_pos+0.5*new_ns)>dataB[0]]
                        
                        datanew=datanew[datanp>(new_pos-0.5*new_ns)]
                        
                        new_value=np.nanmean(datanew)
                        
                        
                        spaceB_n[q,0]=new_pos
                        spaceB_n[q,1]=new_value
                        
                if insertB[-1]!=bins-1:
                    diff=bins-insertB[-1]
                
                    space=endB-spaceB_c[-1,0]
                
                    new_ns=space/diff
                
                    for q in range(1,int(diff)):
                        new_pos=spaceB_c[-1,0]+(q)*new_ns
                        datanew=dataB[1][(new_pos+0.5*new_ns)>dataB[0]]
                        
                        datanp=dataB[0][(new_pos+0.5*new_ns)>dataB[0]]
                        
                        datanew=datanew[datanp>(new_pos-0.5*new_ns)]
                        new_value=np.nanmean(datanew)
                        
                        spaceB_n[insertB[-1]+q,0]=new_pos
                        spaceB_n[insertB[-1]+q,1]=new_value


                for t in range(len(insertB)-1):
                    diff=insertB[t+1]-insertB[t]
                    
                    space=spaceB_c[t+1,0]-spaceB_c[t,0]
                    new_ns=space/diff
                    if diff>1:
                        for q in range(1,int(diff)):
                            new_pos=spaceB_c[t,0]+(q)*new_ns
                            datanew=dataB[1][(new_pos+0.5*new_ns)>dataB[0]]
                            
                            datanp=dataB[0][(new_pos+0.5*new_ns)>dataB[0]]
                        
                            datanew=datanew[datanp>(new_pos-0.5*new_ns)]
                            
                            
                            
                            new_value=np.nanmean(datanew)

                            spaceB_n[int(insertB[t])+q,0]=new_pos
                            spaceB_n[int(insertB[t])+q,1]=new_value
               
               
                
               
                            
                
                
            for l in range(len(new_arrayC)):
                insertC=notindexes(new_arrayC[l,0])
                
                if len(insertC)<1:
                    spaceC_n=spaceC_c
                else:
                    if len(insertC)<bins:
                        peak_distances=travel_determiner(insertC,peak_inds)
                    spaceC_n=np.zeros([bins,2])
                    for t in range(len(insertC)):
                        spaceC_n[insertC[t],:]=spaceC_c[t,:]
                    if insertC[0]!=0:
                        diff=insertC[0]
                
                        space=spaceC_c[0,0]-startC
                
                        new_ns=space/diff
                
                        for q in range(int(diff)):
                            new_pos=startC+(q+1)*new_ns
                            datanew=dataC[1][(new_pos+0.5*new_ns)>dataC[0]]

                            datanp=dataC[0][(new_pos+0.5*new_ns)>dataC[0]]

                            datanew=datanew[datanp>(new_pos-0.5*new_ns)]
                            new_value=np.nanmean(datanew)

                            spaceC_n[q,0]=new_pos
                            spaceC_n[q,1]=new_value
                        
                    if insertC[-1]!=bins-1:
                        diff=bins-insertC[-1]
                
                        space=endC-spaceC_c[-1,0]
                
                        new_ns=space/diff
                
                        for q in range(1,int(diff)):
                            new_pos=spaceC_c[-1,0]+(q)*new_ns
                            datanew=dataC[1][(new_pos+0.5*new_ns)>dataC[0]]

                            datanp=dataC[0][(new_pos+0.5*new_ns)>dataC[0]]

                            datanew=datanew[datanp>(new_pos-0.5*new_ns)]
                            new_value=np.nanmean(datanew)

                            spaceC_n[insertC[-1]+q,0]=new_pos
                            spaceC_n[insertC[-1]+q,1]=new_value


                    for t in range(len(insertC)-1):
                        diff=insertC[t+1]-insertC[t]
                        space=spaceC_c[t+1,0]-spaceC_c[t,0]
                        new_ns=space/diff
                        if diff>1:
                            for q in range(1,int(diff)-1):
                                new_pos=spaceC_c[t,0]+(q)*new_ns
                                datanew=dataC[1][(new_pos+0.5*new_ns)>dataC[0]]
                                
                                datanp=dataC[0][(new_pos+0.5*new_ns)>dataC[0]]
                        
                                datanew=datanew[datanp>(new_pos-0.5*new_ns)]
                                new_value=np.nanmean(datanew)

                                spaceC_n[int(insertC[t])+q,0]=new_pos
                                spaceC_n[int(insertC[t])+q,1]=new_value
                #calculate PCCs between possible pairwise combinations 
                
                
                
                cor_arrA=np.nan_to_num(spaceA_n[:,1])
                cor_arrB=np.nan_to_num(spaceB_n[:,1])
                cor_arrC=np.nan_to_num(spaceC_n[:,1])
                
                cgA=coarse_grainer(cor_arrA)
                cgB=coarse_grainer(cor_arrB)
                cgC=coarse_grainer(cor_arrC)
                cg_arrA.append(cgA)
                cg_arrB.append(cgB)
                cg_arrC.append(cgC)
                normA=spaceA_n[:,1]/np.max(spaceA_n[:,1])
                normB=spaceB_n[:,1]/np.max(spaceB_n[:,1])
                normC=spaceC_n[:,1]/np.max(spaceC_n[:,1])

                normA=np.nan_to_num(normA)
                normB=np.nan_to_num(normB)
                normC=np.nan_to_num(normC)
               
                
                
                
                cov_AB=pearsonr(normA,normB)[0]
               
                cov_BC=pearsonr(normB,normC)[0]
                cov_AC=pearsonr(normA,normC)[0]
                
                
                peak_travel=np.append(peak_travel,np.max(peak_distances))
                
                 
                #summate PCCs
                cov_arr=np.append(cov_arr,cov_AB+cov_BC+cov_AC)
                
                #bin new profiles
                insert_indsA.append(spaceA_n)
                insert_indsB.append(spaceB_n)
                insert_indsC.append(spaceC_n)
    
    #find max summate PCC
    
    if np.min(peak_travel)>4:
        return 'error'
    if len(peak_travel)<2:
        cov_arr2=cov_arr
    else:
        cov_arr2=cov_arr[np.where(peak_travel<5)]
    max_cor_ind=np.argmax(np.nan_to_num(cov_arr2))
    
    if np.all(np.isnan(cov_arr2)):
        return 'error'
    newA= insert_indsA[max_cor_ind]
    newB = insert_indsB[max_cor_ind]
    newC = insert_indsC[max_cor_ind]
   
    return [newA,newB,newC]

def nan_remover(arr,start,end):

    """
    Remove nans and zero values from peak position lists.

    Args:
        arr (arr): Array of values
        start (int): Size marker peak at the start of the packaged region
        end (int): Size marker peak at the end of the packaged region


    Returns:
        (arr): Peak positions with nan and zero values removed
    """    
   
    new_arr=arr
    
    num_inds=np.transpose(np.where(~np.isnan(new_arr)))
    num_inds=num_inds.astype(int)
    
    if num_inds[0]!=0:
        
        diff=num_inds[0]
        
        sep=new_arr[num_inds[0]]-start
        
        new_ns=sep/diff
        
        for i in range(int(diff)):
            new_arr[i]=start+(i+1)*new_ns
    
        
        
    for i in range(len(num_inds)-1):
        if (num_inds[i+1]-num_inds[i])>1:
            diff=num_inds[i+1]-num_inds[i]
            sep=new_arr[num_inds[i+1]]-new_arr[num_inds[i]]
            new_ns=sep/diff
            
            for t in range(int(diff)-1):
                new_arr[num_inds[i]+t+1]=new_arr[num_inds[i]]+(t+1)*new_ns
            
    if num_inds[-1]!=(len(new_arr)-1):
        
        diff=len(new_arr)-num_inds[-1]
        
        sep=end-new_arr[num_inds[-1]]
        
        new_ns=sep/diff
        
        for i in range(int(diff)-1):
            
            
            new_arr[num_inds[-1]+i+1]=new_arr[num_inds[-1]]+(i+1)*new_ns 
    
    
    num_inds=np.transpose(np.nonzero(new_arr))
    
    num_inds=num_inds.astype(int)
    if len(num_inds)>0:
        if num_inds[0]!=0:

            diff=num_inds[0]

            sep=new_arr[num_inds[0]]-start

            new_ns=sep/diff

            for i in range(int(diff)):
                new_arr[i]=start+(i+1)*new_ns



        for i in range(len(num_inds)-1):
            if (num_inds[i+1]-num_inds[i])>1:
                diff=num_inds[i+1]-num_inds[i]
                sep=new_arr[num_inds[i+1]]-new_arr[num_inds[i]]
                new_ns=sep/diff

                for t in range(int(diff)-1):
                    new_arr[num_inds[i]+t+1]=new_arr[num_inds[i]]+(t+1)*new_ns

        if num_inds[-1]!=(len(new_arr)-1):

            diff=len(new_arr)-num_inds[-1]-1

            sep=end-new_arr[num_inds[-1]]

            new_ns=sep/diff

            for i in range(int(diff)):
               
               
                new_arr[num_inds[-1]+i]=new_arr[num_inds[-1]]+(i+1)*new_ns 

 
    
            
                                                      
    return new_arr

def trace_align2(spaces,new_array,bins,bar_lens=None,all_short=False):
    
    """
    Partitioned data alignment process for 2 replicates (secondary function)

    Args:
        spaces (list): Unaligned partitioned footprinting data for each replicate
        new_array (list): List containing new NW alignments for each replicate
        bins (int): Number of position bins in size marker space currently being investigated
        bar_lens (list): List of aligned lengths for different replicates
        all_short (bool): indicate whether there are any alignments which fall short
    Return:
        (list): The aligned partitioned footprinted data for each replicate in the specific part of size marker region.
    """
    
    
    cov_arr=[]
    insert_indsA=[]
    insert_indsB=[]
    
    new_array1=deepcopy(new_array)
    spaces1=deepcopy(spaces)
    
    #extract barcodes
    new_arrayA=new_array[0]
    new_arrayB=new_array[1]

    #extract data
    spaceA_c=spaces[0]
    spaceB_c=spaces[1]

    einsA=False
    einsB=False
    diffA=0
    diffB=0
    if all_short==True:
	    if bar_lens[0]<bins:
		einsA=True
		diffA=bar_lens[0]
	    if bar_lens[1]<bins:
		einsB=True
		diffB=bar_lens[1]

    #iterate through combinations
    for j in range(len(new_arrayA)):
        
        #find insertion indices
        insertA=indexes(new_arrayA[j,0],bins,diffA,eins=einsA)
        correctA=np.arange(len(insertA))

        #correct for insertion indices
        insertA_f=np.subtract(insertA,correctA)

        #create new intensity profile
        spaceA_n=np.insert(spaceA_c,insertA_f.astype(int),np.array((np.nan,0.1)),axis=0)


        for k in range(len(new_arrayB)):
            insertB=indexes(new_arrayB[k,0],bins,diffB,eins=einsB)
            correctB=np.arange(len(insertB))

            insertB_f=np.subtract(insertB,correctB)
            spaceB_n=np.insert(spaceB_c,insertB_f.astype(int),np.array((np.nan,0.1)),axis=0)

            #calculate PCC
            cov_AB=get_cov(spaceA_n[:,1],spaceB_n[:,1])
            cov_arr=np.append(cov_arr,cov_AB)
            insert_indsA.append(spaceA_n)
            insert_indsB.append(spaceB_n)


    #find max PCC
    
    
    max_cor_ind=np.argmax(cov_arr)

    #return profiles associated with max PCC
    newA = insert_indsA[max_cor_ind]
    newB = insert_indsB[max_cor_ind]
    return [newA,newB]
            
            
def nw_align(bars,bins,labels):
    
    """
    Needleman-Wunsch alignment of barcodes

    Args:
        bars (list): Barcodes for partitioned footprint data before alignment
        bins (array): Number of bins (i.e. number of nucleotides between size markers) expected in each part of the size marker set
        labels (list): Replicate labels

    Returns:
        (tuple):
            output (list): Aligned barcodes
            lens_out (list): Number of aligned barcodes generated for each replicate

    """

    align_bin=[]
    score_bin=[]
    lab_bin=[]

    output=[]
    lens_out=[]
    lens_in=[]
    len_bin=[]
    for i in range(len(bars)):
	lens_in.append(len(bars[i]))
 
    #print 'lens_in'
    #print lens_in
    #print bins    
    #iterate through barcodes 
    for j in range(len(bars)):
        for i in range(len(bars)): 
            #if aligning barcode with itself skip
            if j==i:
                continue
            #if one barcode is the right length severely punish insertions on that sequences
            elif (len(bars[j])==bins) and (len(bars[i])!=bins):
                nw=pwise.align.globalxd(bars[j],bars[i],-5,-1,-1,-1)
            
            #if the other barcode is the right length severely punish insertions on that sequences
            elif (len(bars[j])!=bins) and (len(bars[i])==bins):
                nw=pwise.align.globalxd(bars[j],bars[i],-1,-1,-5,-1)
            #if both barcodes have less than the expected length equal gap penalties.     
            elif (len(bars[j])!=bins) and (len(bars[i])!=bins):
                nw=pwise.align.globalxx(bars[j],bars[i])
            else:
                continue
            
            #iterate through aligned data
            for align in nw:
		if len(align[0])<=bins and len(align[1])<=bins:
                    #bin aligned sequences
                    align_bin=np.append(align_bin,align[0])
                    align_bin=np.append(align_bin,align[1])
                    len_bin=np.append(len_bin,len(align[0]))
                    len_bin=np.append(len_bin,len(align[1]))

                    #bin alignment scores
                    score_bin=np.append(score_bin,[align[2]]*2)
                    
                    #bin replicate labels
                    lab_bin=np.append(lab_bin,labels[j])
                    lab_bin=np.append(lab_bin,labels[i])

   
    #create array with all info
    new_array=np.transpose([align_bin,lab_bin,score_bin,len_bin])
    
#if nothing in new array return empty array
    if len(new_array)==0:
        output=[[],[],[]]
        lens_out=[0,0,0]
        
    #split array based on replicate    
    else:
	max_len=np.max(len_bin)
        for t in range(len(labels)):
            new_array1=new_array[new_array[:,1]==labels[t],:]
	    
            if len(new_array1)>0:
		max_len=np.max(new_array1[:,-1].astype(float)) 
                
                #filter out sequences with highest score
                new_array2=new_array1[new_array1[:,-1].astype(float)==max_len]


                #extract unique alignments
                u,u_ind1=np.unique(new_array2[:,0],return_index=True)

                new_array2=new_array2[u_ind1,:]

                #if possible alignments are too numerous, prune
                if len(new_array2)>100:
                    new_array2=new_array2[:100,:]
            else:
                new_array2=np.array([[bars[t],labels[t],'8',len(bars[t])]])
                
            lens_out=np.append(lens_out,len(new_array2))
            
            #add aligned sequences to output array
            output.append(new_array2)
    
            #print output
    #print output
    return output,lens_out
            
            
def indexes(string,bins,bl,character='-',eins=False):

    """
    Determine the insertion points from aligned sequences

    Args:
        string (str): The aligned barcodes
        bins (int): Number of bins in size marker space under consideration
        bl (list): Lengths for aligned traces in the different replicates
        character (chr): The character to find
        eins (bool): Indicate whether extra indices are required

    Returns:
        (array): indices where the chosen characters occurs
    """

    
    output=[]
    for i,c in enumerate(string):
        if c==character:
            output = np.append(output,i)
    
    if eins:
	output=np.append(output,np.arange(bl,bins))
    return output

def notindexes(string,character='-'):

    """
    Determine the non insertion points from aligned sequences

    Args:
        string (str): The aligned barcodes
        character (chr): the character to find
    Returns:
        (array): indices where the chosen characters does not occur
    """
    
    output=[]
    for i,c in enumerate(string):
        if c!=character:
            output = np.append(output,i)
    return output
                
        
def barcode_generator(array):
    
    """
    Barcode generator

    Args:
        array (array): Profile of the partition footprinting data
    Returns:
        (str): Barcode based on profile

    """
    
    max_val=np.max(array)
    
    string=''
    
    for i in range(len(array)):
        
        if array[i]>max_val*0.66:
            string=string+'H'
            
        elif array[i]>max_val*0.33:
            string=string+'M'
            
        else:
            string=string+'L'
    
    return string

def coarse_grainer(array):
    
    """
    Coarse grained profile generator

    Args:
        array (array): Profile of the partition footprinting data
    Returns:
        (array): Coarse graining of profile
    """

    max_val=np.max(array)
    
    cg_arr=[]
    
    for i in range(len(array)):
        
        if array[i]>max_val*0.66:
            cg_arr=np.append(cg_arr,2)
             
        elif array[i]>max_val*0.33:
            cg_arr=np.append(cg_arr,1)
        else:
            cg_arr=np.append(cg_arr,0)
    
    return cg_arr
            
def RX_calculator_single(partition_RX,data_arr2,RX):
    
    """
    Calculate reactivities of single replicate

    Args:
        Partition_RX (list): Partitioned footprinting data profile in peakList object
        data_arr (list): All of the preprocessed datasets in the ensemble
        RX (int): Index of dataset under investigation

    Returns:
        (list): PeakList object containing calculated peak areas and widths

    """
  
    new_peak_list1=fit_shape_gauss(partition_RX,data_arr2[RX])
    
    
    return new_peak_list1




def error_propagation(arr1,arr2):
    
    """
    Trigonometric error propogation

    Args:
        arr1 (array): First array of errors
        arr2 (array): Second array of errors

    Returns:
        (array): Combined errors
    """

    sq_sum_err=np.add(np.square(arr1),np.square(arr2))
    
    new_errs=np.sqrt(sq_sum_err)

    return new_errs

def distance_determiner(arr1,arr2):
    
    """
    Distance calculator between different nucleotide positions

    Args:
        arr1 (array):  first array of indices
        arr2 (array): Second array of indices

    Returns:
        (array): distances between the selected indices in the two arrays

    """
   
    arr_out=[]
    
    
    #iterate through arr1
    for i in range(len(arr1)):
        
        #extract each value in arr1
        value=arr1[i]

        #calculate absolute differences 
        arr3=np.absolute(arr2-value)

        #sum and bin absolute differences
        arr_out=np.append(arr_out,np.sum(arr3))

    return arr_out


def travel_determiner(arr1,arr2):
    
    """
    Calculate distance that a peak has shifted following realignment

    Args:
        arr1 (array): Initial positions of peaks
        arr2 (array): Positions of peaks post alignment

    Returns:
        (array): Distances travelled in nucleotides by peaks

    """

    arr_out=[]
    
    #iterate through arr1
    for i in range(len(arr1)):
        

        #sum and bin absolute differences
        arr_out=np.append(arr_out,np.absolute(arr2[i]-arr1[i]))

    return arr_out


def fit_shape_gauss(dPeakList,data,isOptPos=True,controlA=None):    

    """
    Calculate the reactivities for each peak in the footprinting data

    Args:

        dPeakList (list): peakList object containing footprinting data
        data (array): Corresponding preprocesed
        dataset
        isOptPos (bool): Specify whether position should be optimised
        controlA (float): control values for optimizePosition function
    Returns:
        (list): peak information with reactivities (peak areas) calculated based on the amplitude and the width
    """    

    #deepcopy peak list
    peak_list=deepcopy(dPeakList)
    #deepcopy data
    data_dc=deepcopy(data)
    
    
    #calculate first sigma 
    sigma = funcSeqAll.optimizeOneSigma(data_dc[1],peak_list['pos'],peak_list['amp'])
    
    #if requested optimise position
    if isOptPos:                
        peak_list1=funcSeqAll.optimizePosition(data_dc[1],peak_list,sigma,controlA)
    else:
        peak_list1=peak_list
    #optimise all sigmas  
    peak_list1['wid']=funcSeqAll.optimizeAllSigma(data_dc[1],peak_list1,sigma)
    
    #optimise amplitudes
    peak_list1['amp']=funcSeqAll.optimizeAmp(data_dc[1],peak_list1)
    
    #calculate the areas
    peak_list1['area']=np.abs(peak_list1['amp']*peak_list1['wid'])
    
    
    return peak_list1
        
   
def find_nearest_ind(array, value):
   
    """
    Find the nearest peak to a particular position

    Args: 
        array (array): Array of peak positions
        value (float): Particular position

    Returns:
        (int): The position of the nearest peak
    """

    array1 = np.asarray(array)
    #calculate difference of elements to the chosen value
    array2= array1 - value
  
    #find the absolute values of the negative differences
    array3 = [i for i in range(len(array2)) if array2[i]>0]
    
    
    #find argument of data with minimum difference
    ind =array3[0]
    
    #return index
    return ind            
            
    
def shoulder_finder(peak_arr,data,ind,i):
    
    """
    Find shoulders in a trace 

    Args:
        peak_arr (list): peakList object containing peaks in investigated trace
        data (array): Dataset under investigation
        ind (int): index in the ensemble for the dataset under investigation
    Returns:
        (array): shoulder positions and amplitudes
    """
    
    shoulder_pos_arr = []
    shoulder_amp_arr = []
    
    #read through peak array
    for k in range(len(peak_arr)-1):

        #extract a peak pair
        peak1 = peak_arr[k]
        peak2 = peak_arr[k+1]
        
        #extract the data between the peaks
        trace_between_peaks=deepcopy(data[data[:,0]>peak1])

        trace_between_peaks1 = deepcopy(trace_between_peaks[trace_between_peaks[:,0]<peak2])

        #calculate distance between peaks
        bp_len = len(trace_between_peaks1)

        #if there is no data between peaks continue
        if bp_len == 0:
            continue
        #find trough in data subset
        trough = np.argmin(trace_between_peaks1[:,ind])

        #split data according to the trough position
        first_half=deepcopy(trace_between_peaks1[0:(trough),:])

        second_half = deepcopy(trace_between_peaks1[(trough):(bp_len),:])

        #if the len of the first half is greater than two data points continue
        if len(first_half)>2:
            
            #calculate derivative
            fh_trace_dx = funcGeneral.deriv1(first_half[:,ind])

            #find peaks in derivative
            fh_dx_peaks = funcPeakAlign.peakDetection_v2(fh_trace_dx,isY=False)
            
        
            #if peaks have been found
            if len(fh_dx_peaks)>0:

             
                #bin the shoulder position and amplitude
                shoulder_pos = first_half[fh_dx_peaks,0]
                shoulder_amp = first_half[fh_dx_peaks,ind]
                
                #remove those data points with zero amplitude
                shoulder_pos = shoulder_pos[shoulder_amp>0]

                shoulder_amp = shoulder_amp[shoulder_amp>0]

                #add data to shoulder arrays
                shoulder_pos_arr = np.append(shoulder_pos_arr,shoulder_pos)

                shoulder_amp_arr = np.append(shoulder_amp_arr,shoulder_amp)

               

        if len(second_half)>2:
            
            #calculate derivative
            sh_trace_dx = funcGeneral.deriv1(second_half[:,ind])

            #invert the derivative trace
            mod_sh_trace_dx = -np.absolute(sh_trace_dx)


            #find peaks in the derivative trace
            sh_dx_peaks = funcPeakAlign.peakDetection_v2(mod_sh_trace_dx,isY=False)
       
            
            #check that peaks have been found
            if len(sh_dx_peaks)>0:

                #bin the data 
                shoulder_pos = second_half[sh_dx_peaks,0]

                shoulder_amp = second_half[sh_dx_peaks,ind]

                #remove all those points that have zero amplitude
                shoulder_pos = shoulder_pos[shoulder_amp>0]

                shoulder_amp = shoulder_amp[shoulder_amp>0]


                shoulder_pos_arr = np.append(shoulder_pos_arr,shoulder_pos)

                shoulder_amp_arr = np.append(shoulder_amp_arr,shoulder_amp)

                
        #merge shoulder position and amplitudes
        shoulder_data=[shoulder_pos_arr,shoulder_amp_arr]
        
       
    #transpose data and return         
    return np.transpose(shoulder_data)
                  
    
def position_vote(partition_data,cut1,cut2,plot=0,clip=350):
    
    """
    position voting over the ensemble 


    Args: 
        partition_data (list): Partitioned sequence traces across ensemble
        cut1 (float): Value of first cutoff used to convert each partitioned trace into a binary sequence
        cut2 (float): Threshold ratio for voting
        clip (int): Point at which to clip the binned sequences 
    Returns:
        (tuple):
            correl_return (float): mean correlation between binary arrays in ensemble
            ballot_box (array): Consensus binary array
    """
   
    vote_arr = []
    
    #convert  partition data array to binary
    for i,data in enumerate(partition_data):
        if isinstance(data,np.ndarray):
            data_=deepcopy(data[:,1])
            
        else:
            data_ = deepcopy(data[2])
        
        #convert all points below cutoff to 0 and all above to 1
        data_1=data_[:clip]
        data_1[data_1<np.max(data_1)*cut1]=0
        data_1[data_1>0]=1
        
              
        #add binary arrays to voting matrix
        if i == 0:
            vote_arr = data_1
        else:
            vote_arr = np.vstack((vote_arr,data_1))
    #sum matrix over all entries
    ballot_box = np.sum(vote_arr,axis=0)
    
    #convert ballot box into 
    ballot_box[ballot_box<cut2*len(partition_data)]=0
    ballot_box[ballot_box>0]=1
    
    
    #create correlation matrix
    correl_mat=np.empty([len(partition_data),len(partition_data)])
    for i in range(len(partition_data)):
        for j in range(len(partition_data)):
            correl_mat[i,j] = get_cov(vote_arr[i,:],vote_arr[j,:])
            
    if plot == 1: 
        correl_list = correl_mat.flatten()
        plt.hist(correl_list)
        plt.show()
        correl_return = correl_mat
        
    else:
        correl_return = correl_mat.mean()

            
    return correl_return,ballot_box
    

    
    
def sequence_content(seq_file, nuc = 'T'):

    """
    determine content of a specific nucleotide

    Args:
        seq_file (str): Name of reference sequence fasta file
        nuc (chr): The nucleotide under investigation
    Returns:
        (float): Percentage of the target nucleotide in sequence
    """ 
    
    #read the fasta file 
    for record in SeqIO.parse(seq_file,'fasta'):
        seq_arr = record.seq
        count  = 0
        
        #count the number of occurences of the selected nucleotide
        for i in range(len(seq_arr)):
            if seq_arr[i] == nuc:
                count = count+1
    #calculate percentage            
    return float(count)/float(len(seq_arr))*100


def seq_to_bin(seq_arr,nuc = 'T'):
    
    """
    sequence searching function 

    Args:
        seq_file (str): Name of reference sequence fasta file
        ballot_box (array): Consensus binary representation of sequence from data
        top (int): Highest nucleotide position considered
        bottom (int): Lowest nucleotide position considered
        Nuc (chr): The nucleotide under investigation
    Returns:
        (tuple):
            signif (float): Z score for maximum correlation value
            corr_argmax (int): Index for max correlation
            corr_max (float): Max correlation

    """

    #convert sequence to binary
    for i in range(len(seq_arr)):
            
        if i == 0:
            if seq_arr[i].upper() == nuc:
                bin_seq = 1
            else:
                bin_seq = 0
        else:
            if seq_arr[i].upper() == nuc:
                bin_seq = np.append(bin_seq,1)
            else:
                bin_seq = np.append(bin_seq,0)
    return bin_seq
        

def sequence_search(seq_file, ballot_box,top=10000,bottom=0,Nuc='T'):

    """
    sequence searching function 

    Args:
        seq_file (str): Name of reference sequence fasta file
        ballot_box (array): Consensus binary representation of sequence from data
        top (int): Highest nucleotide position considered
        bottom (int): Lowest nucleotide position considered
        Nuc (chr): The nucleotide under investigation
    Returns:
        (tuple):
            signif (float): Z score for maximum correlation value
            corr_argmax (int): Index for max correlation
            corr_max (float): Max correlation

    """

    #initialise correlation array
    correl_rec = []
    accuracy_rec = []
    
    #extract reference sequence data 
    for record in SeqIO.parse(seq_file,'fasta'):
        seq_arr = record.seq
        
        
        #convert sequence to binary
        bin_seq = seq_to_bin(seq_arr,nuc=Nuc)
        
        bin_seq=bin_seq[bottom:top]
        

        
        
        #search through sequence binary array 
        for i in range(len(seq_arr)-len(ballot_box)):
            
            #extract subarray 
            seq_sec = bin_seq[(i):i+len(ballot_box)]
            
            #calculate correlation between sub array and voting matrix
            
            
            
            if len(ballot_box)==len(seq_sec):
                correl = get_cov(seq_sec,ballot_box[::-1])

                accuracy = accuracy_measure(seq_sec,ballot_box[::-1])
           
            else:
                correl=0
                accuracy=0
            
            correl_rec.append(correl)
            
            accuracy_rec.append(accuracy)
            
        correl_rec = np.array(correl_rec)
        
        accuracy_rec = np.array(accuracy_rec)
    
        signif = signif_assessor(correl_rec)
        
        #return maximum correlation and coordinates
        return signif, np.argmax(correl_rec),np.max(correl_rec)
    
    
def sequence_search_area(seq_file, ballot_box,start,window,Nuc='T'):

    """
    sequence searching function over specific area 

    Args:
        seq_file (str): Name of reference sequence fasta file
        ballot_box (array): Consensus binary representation of sequence from data
        start (int): start of rejoin to scan over
        window (int): size of area to search over
        Nuc (chr): The nucleotide under investigation
    Returns:
        (tuple):
            signif (float): Z score for maximum correlation value
            corr_argmax (int): Index for max correlation
            corr_max (float): Max correlation

    """

    #initialise correlation array
    correl_rec = []
    accuracy_rec = []
    
    #extract reference sequence data 
    for record in SeqIO.parse(seq_file,'fasta'):
        seq_arr = record.seq
        
        
        #convert sequence to binary
        bin_seq = seq_to_bin(seq_arr,nuc=Nuc)        
        
        
        
        #search through sequence binary array 
        for i in range(start-1,start+window-1):
            
            #extract subarray 
            seq_sec = bin_seq[i:i+len(ballot_box)]
            
            #calculate correlation between sub array and voting matrix
            
            
            
            if len(ballot_box)==len(seq_sec):
                correl = get_cov(seq_sec,ballot_box[::-1])

                accuracy = accuracy_measure(seq_sec,ballot_box[::-1])
            
            else:
                correl=0
                accuracy=0
            
            correl_rec.append(correl)
            
            accuracy_rec.append(accuracy)
            
        correl_rec = np.array(correl_rec)
        
        accuracy_rec = np.array(accuracy_rec)
        
        signif = signif_assessor(correl_rec)
        
        #plt.hist(accuracy_rec)
        #plt.show()
            
        #return maximum correlation and coordinates
        return signif, start+np.argmax(correl_rec),np.max(correl_rec)
        
def renormalisation(area_arr):
    
    """
    Renormalises data for two or more reactivity datasets

    Args:
        area_arr (list): area arrays to be combined

    Returns:
        (tuple):
            new_area_arr (list): Normalised reactivities for datasets supplied
            new_aver (array): Normalisation factors produced
    """  
    
    comb_data=[]
    #read over all data
    for i in range(len(area_arr)):
        
        #combine the data into one array
        comb_data=np.append(comb_data,area_arr[i])
    
    #find outliers and average for cobined set
    Pout,Pav=funcSeqAll.findPOutlierBox(comb_data)
    
    new_area_arr=[]
    new_aver=[]
    
    #normalise each data set based on new outlier and average parameters
    for t in range(len(area_arr)):
        
        new_norm,aver=funcSeqAll.normSimple(area_arr[t],Pout,Pav)
        
        #set negative values to zero
        new_norm[new_norm<0]=0
        
        new_area_arr.append(new_norm)
        
        new_aver=np.append(new_aver,aver)
        
    return new_area_arr,aver
              
    
def accuracy_measure(seq_arr,vote_arr):
    
    """
    calculates the accuracy of the sequence alignment 
    processes involved:
    
    if elements are equal add to count
    divide the count by the length of the vote array
    """
    
    count=0
    for i in range(len(vote_arr)):
        if seq_arr[i]==1 and vote_arr[i] == 1:
            count = count + 1
        elif seq_arr[i]==0 and vote_arr[i] == 0:
            count = count + 1
    return float(count)/float(len(vote_arr))
            

def correl_assessor(data_arr,ind):
    
    """
    calculates the accuracy of the sequence alignment 

    Args: 
        seq_arr (array): Binary representation of
        sequence
        vote_arr (array): Binary consensus sequence from data

    Returns: 
        (float): accuracy of consensus sequence
    """

    
    correl_arr = np.empty([len(data_arr),len(data_arr)])
    
    correl_list = []
    for i,data1 in enumerate(data_arr):
        if isinstance(data1,pd.DataFrame):
            data1 = data1.values
        else:
            data1=data1
        for j,data2 in enumerate(data_arr):
            
            if isinstance(data2,pd.DataFrame):
                data2 = data2.values
            
                correl_arr[i,j]=pearsonr(data1[:,ind]/np.max(data1[:,ind]),data2[:,ind]/np.max(data2[:,ind]))[0]
                
            elif isinstance(data2,np.ndarray):
                data2 = data2
                
                
                correl_arr[i,j]=pearsonr(data1[:,ind]/np.max(data1[:,ind]),data2[:,ind]/np.max(data2[:,ind]))[0]
                
            else:
                
                correl_arr[i,j] = pearsonr(data1[ind]/np.max(data1[ind]),data2[ind]/np.max(data2[ind]))[0]
                

            
    return correl_arr


def count_correl_above(correl_mat,limit):
    
    """
    count numbers of correlation matrix elements above a certain threshold

    Args:
        correl_mat (array): Matrix correlation values
        limit: Threshold for counting

    Returns: 
        (float): Percentage of correlations above the limit

    """ 
    
    correl_list = correl_mat.flatten()
    
    full_len = len(correl_list)
    
    
    above_len = len([p for p in correl_list if p>limit])
    
    
    return float(above_len)/float(full_len)*100


def get_cov(array1,array2):
    
    """
    Calculate correlation between two arrays 

    Args:
        array1 (array): First array
        array2 (array): Second array

    Returns:
        (float): correlation value
    """



    return np.cov(array1,array2)[0][1]/(np.std(array1)*np.std(array2))


def signif_assessor(data_arr):
    
    """
    Calculate Z-score of max value

    Args:
        data_arr (array): array of correlation values
    
    Returns:
        (float): Z-score of max correlation value
    """
    
    max_val = np.max(data_arr)
    
    data_arr1 = np.delete(data_arr,np.argmax(data_arr))
    
    mean = np.mean(data_arr1)
    std = np.std(data_arr1)
    
    sigma_2 = (max_val-mean)/std
    
    return sigma_2


def scanned_correl(partition_arr,ind_arr,window=10):
    
    """
    scanning correlation assessor:
    
    extract the datasets you want to work with (max 3)
    calculate moving correlations between different datasets.
    bin correlations    
    """
    
    #extract datasets
    partA=partition_arr[ind_arr[0]]
    
    partB=partition_arr[ind_arr[1]]
    
    partC=partition_arr[ind_arr[2]]
    
    #initialise correl bins
    correlsAB=[]
    correlsAC=[]
    correlsBC=[]
    
    #iterate through the data
    for i in range(len(partC['amp'])-window):
        
        #calculate correlations across windows
        correlAB=get_cov(partA['amp'][i:((i+1)*window)],partB['amp'][i:((i+1)*window)])
        correlAC=get_cov(partA['amp'][i:((i+1)*window)],partC['amp'][i:((i+1)*window)])
        correlBC=get_cov(partC['amp'][i:((i+1)*window)],partB['amp'][i:((i+1)*window)])
        
        #bin correlations
        correlsAB=np.append(correlsAB,correlAB)
        correlsAC=np.append(correlsAC,correlAC)
        correlsBC=np.append(correlsBC,correlBC)
    
    #return binned correlations
    return correlsAB,correlsAC,correlsBC
        
def sequence_snapshots(xcoord,ycoord,yerr,col,window=100,virus='',primer='',condition='',treatment='',diff=False):
    
    """
    Produces snapshots reactivity profile in primer region 

    Args:
        xcoord (array): xcoordinates of reactivity data
        ycoord (array): ycoordinates of reactivity data
        yerr (array): Errors on reactivities
        col (chr): Colour used for
        window (int): Window size for snapshots
        virus (str): Virus under investigation
        primer (str): Primer used
        condition (str): Treatment condition
        treatment (str): treatment exposure time
        diff (bool): specify whether difference map is being plotted
    Returns:
        None
    """
    
    intervals=5
    loc=plticker.MultipleLocator(base=intervals)
    #calculate bottom and top lim
    b_lim=np.floor(xcoord[0]/window)*window
    t_lim=np.ceil(xcoord[-1]/window)*window
    i=0
    #iterate through windows
    while (b_lim+(i-1)*window)<t_lim:
        #calculate edges of image
        start=(b_lim+(i)*window)
        end=(b_lim+(i+1)*window)
        
        
        #generate filename
        file_name=virus+'_'+primer+'_'+condition+'_'+treatment+'_'+str(int(start))+'_'+str(int(end))+'.png'
        fig,ax1=plt.subplots()
        
        
        #plot and save data
        ax1.bar(xcoord,ycoord[::-1],color=col)
        
        if yerr!=None:
            ax1.errorbar(xcoord,ycoord[::-1],yerr=yerr[::-1],color='k',fmt=None,elinewidth=1,capsize=None)
        ax1.set_xlim(start,end)
        if diff:
            ax1.set_ylim(-5,5)
        else:    
            ax1.set_ylim(0,5)
        ax1.xaxis.set_minor_locator(loc)
        ax1.set_xlabel('genome position')
        ax1.set_ylabel('normalised reactivity')
        ax1.grid(which='both', axis='x', linestyle='--', color='#b9babd')
        ax1.xaxis.set_tick_params(labeltop='on')
        plt.title(virus+'_'+primer+'_'+condition+'_'+treatment+'_'+str(int(start))+'_'+str(int(end)), y=1.08)         
        plt.savefig(file_name)
        plt.close()
        i+=1


def RX_calculator_replicates(partition,data_arr,inds,single_0=False,get_correls=False):

    """
    Calculate average and standard error of areas and average of signal amplitude

    Args: 
        partition (list): Partitioned footprinting data
        data_arr (list): Preprocessed datasets in the ensemble
        inds (list): Indices indicating the locations of the replicates in the ensemble data list
        get_correls (bool): Specify whether correlations should be calculated

    Returns:
        (tuple):
            amp_av (array): Average amplitudes of partitioned data
            area_av (array): Average peak areas of partitioned data
            area_sd (array): Standard deviations on areas
    """
    
    areas=[]
    
    #calculate areas 
    for k in range(len(partition)):
        
        peak_list=RX_calculator_single(partition[k],data_arr,inds[k])
        
        areas.append(peak_list)
    
    #if requested calculate correls
    if get_correls:
        print correl_assessor(areas,'area')
    
    #gather areas and amplitudes together
    for t in range(len(partition)):
        
        if t==0:
            
            area_arr=areas[t]['area']
            amp_arr=areas[t]['amp']
            
        else:
            area_arr=np.vstack((area_arr,areas[t]['area']))
            amp_arr=np.vstack((amp_arr,areas[t]['amp']))
    
    #calculate averages and standard errors
    area_av=np.nanmean(area_arr,axis=0)
    amp_av=np.nanmean(amp_arr,axis=0)
    area_sd=np.std(area_arr,axis=0)
    
    
    return amp_av,area_av,area_sd

def RX_correction(area_RX,area_BG,scaling):
    
    """
    Background correction of reactivities and normalisation

    Args:
        area_RX (array): Array of areas for treatment 
        area_BG (array): Array of areas for background
        scaling (array): scaling factors between RX and BG

    Returns:
        (tuple):
            area_diff (array): Unnormalised reactivites
            norm_area_diff (array): Normalised reactivities
            aver (float): Normalisation factor

    """    
    #calculate corrected reactivities
    area_diff=area_RX-scaling*area_BG
    
    #calculate outliers and average percentage
    Pout,Paver=funcSeqAll.findPOutlierBox(area_diff)
    
    #normalise reactivities
    norm_area_diff,aver=funcSeqAll.normSimple(area_diff,Pout,Paver)
    
    #remove negative values
    norm_area_diff[norm_area_diff<0]=0
    
    return area_diff,norm_area_diff,aver

def raw_trace_plotter(file_list):
    
    """
    plots the raw electropherograms to pngs

    Args:
        file_list (list): list of all the data files in the ensemble

    Returns:
        None
    """

    for i,file in enumerate(file_list):
        data = pd.read_csv(file.strip('\n')+'_raw.csv')
        
        
        plt.plot(data['Position'],data['ReactionChannel#1'],'b',label='RX')
        plt.plot(data['Position'],data['SequenceChannel#1'],'r',label='ddA')
        plt.plot(data['Position'],data['SequenceChannel#2'],'g',label='ddC')
        plt.plot(data['Position'],data['SizeMarker'],'k',label='SM')
        plt.legend()
        plt.savefig(file.strip('\n')+'_raw.png')
        plt.close()
        
def sm_plotter(data_arr,TM_peaks,file_list=None):

    """
    Plot size marker traces and positions of size markers determined by peak_finder

    Args: 
        data_arr (list): Datasets in the ensemble
        TM_peaks (list): Size marker peak positions
        file_list (list): file names of datasets in the ensemble

    Returns:
        None
    """
    
    exists=os.path.isdir('./sm_plots')
    
    if exists:
        shutil.rmtree('./sm_plots')
        
    os.mkdir('./sm_plots')
    
    for i in range(len(data_arr)):
        peaks=TM_peaks[i]
        data=data_arr[i]
	if file_list!=None:
        	file=file_list[i]
	else:
		file='this1_'+str(i)
        fig,ax=plt.subplots()
        ax.plot(data[0],data[4],'g',lw=2)
        for t in range(len(peaks)):
            ax.plot([peaks[t]]*2,[-100,10000],'k',linestyle='--',alpha=0.7,lw=2)
      
        plt.savefig('./sm_plots/'+file.strip('\n')+'_TM.png')
        plt.close()
            
       

def RX_analyse(wfiles,indsBG,indsRX,virus,primer,start_pos,condition,exposure,skip=[],wrange=[0,10],sm_extend=0,wcut=0.7):

    """
    Wrapper algorithm for reactivity analysis 

    Args: 
        wfiles (str): String giving the prefix of the windowing files 
        indsBG (list): Indices of the background (0 ms) datasets
        indsRX (list): Indices of the reactivity datasets
        virus (str): Virus name
        primer (str): Primer name
        start_pos (int): Position of first nucleotide in the reactivity profile
        condition (str): Treatment name
        exposure (int): X-ray exposure time
        skip (list): Indicate which datasets should be ignored
        wrange (list): List of preprocessed windowed datasets for examination
        sm_extension (int): Number of size marker points (each representing +10 nts) by which to extend the size marker trace. a negaitve value indicates that less than the normal number of markers are being used (-n indicates that only the first n size marker peaks are to be used)
	wcut (float): Cutoff correlation for windowing 

    Returns:
        None: data deposited in various .csv files
    """

    cor_AB=[]
    cor_BC=[]
    cor_AC=[]

    avers=[] 

    
    if len(wrange)==0:
	wlist=glob.glob(wfiles+'*')
    else:
	wlist=[]
	for w in wrange:
		wlist.append(wfiles+'_'str(w)+'.obj'
		
    for i in range(len(wlist)):
	wpath=wlist[i]
	#open data file
        file_1= open(wpath,'rb')

        #load data
        data_arr2=pickle.load(file_1)

	"""
	for kk in indsBG:
	    plt.plot(data_arr2[kk][0],data_arr2[kk][1]) 
	plt.show()

        for kk in indsRX:
            plt.plot(data_arr2[kk][0],data_arr2[kk][1])    
        plt.show()
        """
    	if len(indsBG)==1 and len(indsRX)==1:
 
		partition_RX_single=RX_partitioning_single(data_arr2,1,Issues=issues,Skip=skip,extension=sm_extend)
	
		peak_list_BG=RX_calculator_single(partition_RX_single[indsBG[0]],data_arr2,indsBG[0])

        	peak_list_RX=RX_calculator_single(partition_RX_single[indsRX[0]],data_arr2,indsRX[0])
	    	sf=funcSeqAll.scaleShapeDataWindow(peak_list_BG['amp'],peak_listRX['amp'],deg=40,rate=0.25,step=10,fit='linear')

   	     	ad_RX,nad_RX,aver_RX=RX_correction(peak_list_RX['area'],peak_list_BG['area'],sf)

    	elif len(indsBG)==1 and len(indsRX)>1:


             	partition_RX_single=RX_partitioning_single(data_arr2,1,Issues=issues,Skip=skip,extension=sm_extend)
	     	partition_RX_replicates=RX_partitioning_replicates_extended(data_arr2,1,.12,Issues=issues,Skip=skip,extension=sm_extend)

             	peak_list_BG=RX_calculator_single(partition_RX_single[indsBG[0]],data_arr2,indsBG[0])

             	part_RX=RX_partition_realignment_extended(partition_RX_replicates,indsRX,data_arr2)


		cor_mat=correl_assessor(part_RX,'amp')

                if len(indsRX)==2:
                        cor_AB=np.append(cor_AB,cor_mat[0,1])
                else:
                        cor_AB=np.append(cor_AB,cor_mat[0,1])
                        cor_AC=np.append(cor_AC,cor_mat[0,2])
                        cor_BC=np.append(cor_BC,cor_mat[2,1])

                print 'Background replicate correlations'
                print correl_assessor(part_BG,'amp')

                print 'Treatment replicate correlations'
                print correl_assessor(part_RX,'amp')

	     	amp_av_RX,area_av_RX,area_sd_RX=RX_calculator_replicates(part_RX,data_arr2,indsRX)

		nad_se_RX=area_sd_RX/(aver_RX*len(indsRX))
		ad_se_RX=area_sd_RX/(len(indsRX))

             	sf=funcSeqAll.scaleShapeDataWindow(amp_av_RX,peak_list_BG['amp'],deg=40,rate=0.25,step=10,fit='linear')

		ad_RX,nad_RX,aver_RX=BoXFP.RX_correction(area_av_RX,peak_list_BG['area'],sf)
	else:
	 	print 'cat'

	        partition_RX_replicates=RX_partitioning_replicates_extended(data_arr2,1,.12,Issues=issues,Skip=skip,extension=sm_extend)

	     	part_BG=RX_partition_realignment_extended(partition_RX_replicates,indsBG,data_arr2)

	     	part_RX=RX_partition_realignment_extended(partition_RX_replicates,indsRX,data_arr2)


	     	cor_mat=correl_assessor(part_RX,'amp')

	     	if len(indsRX)==2:
			cor_AB=np.append(cor_AB,cor_mat[0,1])
	     	else:
			cor_AB=np.append(cor_AB,cor_mat[0,1])
        		cor_AC=np.append(cor_AC,cor_mat[0,2])
        		cor_BC=np.append(cor_BC,cor_mat[2,1])

	     	print 'Background replicate correlations'
             	print correl_assessor(part_BG,'amp')

             	print 'Treatment replicate correlations'
             	print correl_assessor(part_RX,'amp')
            
             	amp_av_BG,area_av_BG,area_sd_BG=RX_calculator_replicates(part_BG,data_arr2,indsBG)

             	amp_av_RX,area_av_RX,area_sd_RX=RX_calculator_replicates(part_RX,data_arr2,indsRX)

             	sf=funcSeqAll.scaleShapeDataWindow(amp_av_RX,amp_av_BG,deg=40,rate=0.12,step=10,fit='linear')


             	ad_RX,nad_RX,aver_RX=RX_correction(area_av_RX,area_av_BG,sf)
	     	ad_se_RX=error_propagation((area_sd_RX/len(indsRX)),sf*(area_sd_BG/len(indsBG)))

		nad_se_RX=area_sd_RX/(aver_RX*len(indsRX))
                ad_se_RX=area_sd_RX/(len(indsRX))


	avers=np.append(avers,aver_RX)

  
	print 'Normalisation factors: '+str(aver_RX)


        if i==wrange[0]:
			
             	nad_RX_arr=deepcopy(nad_RX)
		ad_RX_arr=deepcopy(ad_RX)

		if len(indsRX)>1:
			 nad_se_RX_arr=deepcopy(nad_se_RX)
                         ad_se_RX_arr=deepcopy(ad_se_RX)

	else:
		nad_RX_arr=np.vstack((nad_RX_arr,nad_RX))
                ad_RX_arr=np.vstack((ad_RX_arr,ad_RX))
		
		if len(indsRX)>1: 	
                 	nad_se_RX_arr=np.vstack((nad_se_RX_arr,nad_se_RX))
                        ad_se_RX_arr=np.vstack((ad_se_RX_arr,ad_se_RX))


    ca_RX=np.zeros([len(nad_RX_arr),len(nad_RX_arr)])
 
    for i in range(len(nad_RX_arr)):
	for j in range(len(nad_RX_arr)):
		ca_RX[i,j]=pearsonr(nad_RX_arr[i,:],nad_RX_arr[j,:])[0]
    
    mean_ca_RX=np.mean(ca_RX,axis=1)
    print mean_ca_RX
    keep_ca_RX=np.where(mean_ca_RX>wcut)[0]
    print keep_ca_RX

    nad_RX=np.mean(nad_RX_arr[keep_ca_RX,:],axis=0)
    ad_RX=np.mean(ad_RX_arr[keep_ca_RX,:],axis=0)
    nuc_pos=np.arange(len(nad_RX))+start_pos
    if len(indsRX)>1 or len(indsBG)>1:
        nad_se_RX=np.mean(nad_se_RX_arr[keep_ca_RX,:],axis=0)
        ad_se_RX=np.mean(ad_se_RX_arr[keep_ca_RX,:],axis=0)	

	full_arr=np.transpose([nuc_pos,nad_RX[::-1],nad_se_RX[::-1]])

    	col_names=['genome_position','norm_react_RX','norm_react_se_RX']

        full_arr_ad=np.transpose([nuc_pos,ad_RX[::-1],ad_se_RX[::-1]])

        col_names_ad=['genome_position','ad_react_RX','ad_react_se_RX']
	if len(indsRX)==3:

        	cor_labs=['corAB','corBC','corAC']

    		cor_arr=np.transpose([cor_AB,cor_BC,cor_AC])

	elif len(indsRX)==2:

                cor_labs=['corAB']

                cor_arr=np.transpose([cor_AB])
	
   	cor_df=pd.DataFrame(cor_arr,columns=cor_labs)
	repCorrDir=virus+'_pri'+primer+'_repCorrs'
	if not os.path.exists(repCorrDir):
    		os.makedirs(repCorrDir)
    	cor_df.to_csv(repCorrDir+'/'+virus+'_pri'+primer+'_'+condition+'_'+str(exposure)+'ms_repCorrs.csv',sep=',',index=False)

    else:
	

    	full_arr=np.transpose([nuc_pos,nad_RX[::-1]])

        col_names=['genome_position','norm_react_RX']

        full_arr_ad=np.transpose([nuc_pos,ad_RX[::-1]])

        col_names_ad=['genome_position','ad_react_RX']

    full_df=pd.DataFrame(full_arr,columns=col_names)

    reacDir=virus+'_pri'+primer+'_reacs'
    if not os.path.exists(reacDir):
    	os.makedirs(reacDir)

    full_df.to_csv(reacDir+'/'+virus+'_pri'+primer+'_'+condition+'_'+str(exposure)+'ms_reacsNorm.csv',sep=',',index=False)


    full_df_ad=pd.DataFrame(full_arr_ad,columns=col_names_ad)

    full_df_ad.to_csv(reacDir+'/'+virus+'_pri'+primer+'_'+condition+'_'+str(exposure)+'ms_reacs.csv',sep=',',index=False)


    avers_arr=np.transpose([avers])

    avers_labs=['avers']

    avers_df=pd.DataFrame(avers_arr,columns=avers_labs)

    averDir=virus+'_pri'+primer+'_normFactors'
    if not os.path.exists(averDir):
    	os.makedirs(averDir)

    avers_df.to_csv(averDir+'/'+virus+'_pri'+primer+'_'+condition+'_'+str(exposure)+'ms_normFactors.csv',sep=',',index=False)

    windDir=virus+'_pri'+primer+'_windowingCorrs'
    if not os.path.exists(windDir):
        os.makedirs(windDir)

    np.savetxt(windDir+'/'+virus+'_pri'+primer+'_'+condition+'_'+str(exposure)+'ms_wcorrMat.csv',ca_RX)

    #### write mean window correlations to file
    mean_cor_arr=np.transpose([mean_ca_RX])

    mean_cor_labs=['cor_RX']

    mean_cor_df=pd.DataFrame(mean_cor_arr,columns=mean_cor_labs)

    mean_cor_df.to_csv(windDir+'/'+virus+'_pri'+primer+'_'+condition+'_'+str(exposure)+'ms_windCorrs.csv',sep=',',index=False)

def DM_generator(data1,data2):
    
    """
    Difference map calculations
    
    Args:
        data1 (array): First data array
        data2 (array): Second data array 
    Returns:
        (tuple):
            data_out_arr (array): Data array of difference maps
            avers_arr (array): Array of normalisation factors 
    """
    
    data_v1=data1.values
    
    data_v2=data2.values
    
    data_out_arr=data_v1[:,0].reshape((350,1))
    
    
    for i in range(1,data1.shape[1],2):
        
        
        comb_data=np.append(data_v1[:,i],data_v2[:,i])
        
        Pout,Pav=funcSeqAll.findPOutlierBox(comb_data)
        
        data_n1,aver1=funcSeqAll.normSimple(data_v1[:,i],Pout,Pav)
        
        data_n1[data_n1<0]=0
        
        data_n2,aver2=funcSeqAll.normSimple(data_v2[:,i],Pout,Pav)
        
        data_n2[data_n2<0]=0
        
        avers=np.array([aver1,aver2])
        
        data_o=data_n2-data_n1
        
        err_o=error_propagation(data_v1[:,i+1]/aver1,data_v2[:,i+1]/aver2)
        
        
        data_out=np.transpose([data_o,err_o])
        #print data_out
        if i==1:
            avers_arr=avers
        else:
            avers_arr=np.append(avers_arr,avers)
        
        
        data_out_arr=np.append(data_out_arr,data_out,axis=1)
    return data_out_arr,avers_arr
			     
			     
def write_out_raw_csv(data_bloc, data_list):
    """
    Writes out the CSV data from the raw FSA file
    
    Args:
        data_bloc (array): Array of unpacked data from ABIF files
        data_list (list): List of data file names
    Returns:
        None
    """
    for i,data_file in enumerate(data_list):
            data = data_bloc[i]
	    f = open('%s' % data_file.replace('.fsa', '_raw.csv'), 'w')
	    f.write('Position,ReactionChannel#1,SequenceChannel#1,SequenceChannel#2,SizeMarker\n')
	    for position in range(len(data)):
		f.write('%d,%d,%d,%d,%d\n' % (position+1,
		                              data[position][0],    # reaction channel #1
		                              data[position][1],    # sequencing channel #1
		                              data[position][2],    # reaction channel #2
		                              data[position][3]))   # sequencing channel #2
	    f.close()

def readABI(dir_list):
    """
    Read ABIF files for output into other forms
    
    Args:
        dir_list (list): list of datafiles
    Returns:
        (array): array of unpacked data
        
    """
    data_bloc = []
    for fsa_file in dir_list:
            reader = ABIFReader.ABIFReader(fsa_file)
	    col0 = reader.getData('DATA',1)
	    col1 = reader.getData('DATA',2)
	    col2 = reader.getData('DATA',3)
	    col3 = reader.getData('DATA',4)
	    
	    data=np.zeros([len(col0),4],dtype='f4')
	    data[:,0]=np.array(col0)
	    data[:,1]=np.array(col1)
	    data[:,2]=np.array(col2)
	    data[:,3]=np.array(col3)
	    data_bloc.append(data)
    return data_bloc

def signal_assessor(dir_list):
    
    """
    Signal assessment function
    
    Args:
        dir_list (list): list of datafiles
    
    Returns: 
        None
    """
    dir_name = os.path.basename(os.getcwd())

    f = open(dir_name+'_signal_assess.csv','w+')
    f.write('sample,av_RC,sd_RC,av_SC1,sd_SC1,av_SC2,sd_SC2\n')

    for fsa_file in dir_list:

        sample_id = fsa_file.strip('.fsa')

        print sample_id

        csv_file = sample_id + '_raw.csv'

        data_temp = open(csv_file,'r').readlines()

        data_1 = []

        labels = ['Position','ReactionChannel#1','SequenceChannel#1','SequenceChannel#2','SizeMarker']

        for line in data_temp[1:]:
            data_1.append([int(i) for i in line.strip('\n').split(',')])



        data = [[i[j] for i in data_1] for j in range(5)] 

        if np.argmax(data[1]) > 5000:

            quit()

        data = [i[np.argmax(data[1]):] for i in data] 

        sub_samples = [[data[k][j:j+1000] for j in range(0,len(data[k]),1000)] for k in range(5)] 


        avg_arr = [[],[],[],[]]

        for k in range(1,4+1):
            for sample in sub_samples[k]:
                tmp_arr = sample 
                tmp_arr.sort()
                avg_min = round(np.mean(tmp_arr[:250]),3)
                avg_max = round(np.mean(tmp_arr[-250:]),3)

                avg_arr[k-1].append([avg_min,avg_max,sub_samples[0][sub_samples[k].index(sample)][0],sub_samples[0][sub_samples[k].index(sample)][-1]])


        for k in range(0,4):
            avg_arr[k].pop(np.argmax([i[1] for i in avg_arr[k]])) # - highest high
            avg_arr[k].pop(np.argmin([i[1] for i in avg_arr[k]])) # - lowest high ?


        sig_strength = np.mean([avg_arr[0][i][1]-avg_arr[0][i][0] for i in range(len(avg_arr[0]))])
        seq_strength = np.mean([avg_arr[1][i][1]-avg_arr[1][i][0] for i in range(len(avg_arr[0]))])
        seq2_strength = np.mean([avg_arr[2][i][1]-avg_arr[2][i][0] for i in range(len(avg_arr[0]))])

        sig_strength2 = np.std([avg_arr[0][i][1]-avg_arr[0][i][0] for i in range(len(avg_arr[0]))])
        seq_strength2 = np.std([avg_arr[1][i][1]-avg_arr[1][i][0] for i in range(len(avg_arr[0]))])
        seq2_strength2 = np.std([avg_arr[2][i][1]-avg_arr[2][i][0] for i in range(len(avg_arr[0]))])



        f.write('%s,%f,%f,%f,%f,%f,%f\n'% (sample_id,round(sig_strength,3),round(sig_strength2,3),round(seq_strength,3),round(seq_strength2,3),round(seq2_strength,3),round(seq2_strength2,3)))
		
		
    f.close()

