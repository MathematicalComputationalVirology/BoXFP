import os
import sys

sys.path.append(os.path.abspath('/Users/mqbppsc3/Desktop/externalLibraries'))

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
from random import sample
import struct

import matplotlib.ticker as plticker

from mpl_toolkits.mplot3d import Axes3D

from mpl_toolkits.mplot3d import proj3d

from difflib import SequenceMatcher

from Bio import SeqIO
import Bio.pairwise2 as pwise
from scipy import stats

from scipy.stats.stats import pearsonr

from scipy.stats.stats import spearmanr


#convert sequence into gaussian peaks
def sequence_converter_gaussians(seq_file,letter1='T',letter2='C',height1=100.0,height2=50.0,width=2.0,low_height=10.0,sep=10,plot = True):
    
    """
    Convert nucleotide sequence into profile of gaussians
    
    Args: 
    
        seq_file(str):The path of the reference genome fasta file
        letter1(chr):The nucleotide associated with the higher gaussian
        letter2(chr):The nucleotide associated with the medium gaussian
        height1 (float): The amplitude of the highest gaussian
        height2 (float): The amplitude of the medium gaussian
        width (float): The width of the gaussians 
        low_height (float): The amplitude of gaussians associated with the non specific nucleotides
        sep (float): The distances between the centres of gaussians in the profile
    
    Returns:
        (tuple):
            seq_arr (np.ndarray): guassian profile
            nt_arr (np.ndarray): locations of nucleotides across profile 
        
    """
    
    #extract sequence
    for record in SeqIO.parse(seq_file,'fasta'):
        seq = record.seq
    
    #sequence array for gaussian peaks to bin into
    seq_arr = [0 for i in range(sep*(len(seq)+2))]

    #peak width as nucleotide number
    p_width = int(round(width))

    #array of nucleotide positions 
    nt_arr = [i/10 for i in range(sep*(len(seq)+2),0,-1)]

    #gaussians
    #gaussian for T nucleotides
    g_h1_arr = [height1*math.exp(-(i**2)/(2.0*width**2)) for i in range(-p_width*4,p_width*4)]
    
    #gaussian for C nucleotides
    g_h2_arr = [height2*math.exp(-(i**2)/(2.0*width**2)) for i in range(-p_width*4,p_width*4)]
    
    #gaussian for A and G nucleotides
    g_l_arr = [low_height*math.exp(-(i**2)/(2.0*width**2)) for i in range(-p_width*4,p_width*4)]

    #go through nucleotide sequence
    for i in range(len(seq)):
        #add gaussian for T values
        if (seq[i]==letter1):
            for k in range(len(g_h1_arr)):
                seq_arr[i*sep-int(len(g_h1_arr)/2.0)+k] += g_h1_arr[k]
                
        #add gaussian for C values
        elif (seq[i]==letter2):
            for k in range(len(g_h2_arr)):
                seq_arr[i*sep-int(len(g_h2_arr)/2.0)+k] += g_h2_arr[k]
                
        # add gaussian for anything else
        else:
            for k in range(len(g_l_arr)):
                seq_arr[i*sep-int(len(g_l_arr)/2.0)+k] += g_l_arr[k]

    #plot the gaussian trace
    if plot:
        plt.plot(nt_arr,seq_arr,'r')
        plt.savefig('combined_sequence_trace.png')
    
    #return nucleotide position list and trace
    return seq_arr,nt_arr

    
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
    

def trace_align(data_arr,ind,ref,samp =False,ifmob=True,ifdtw=True,ifpeak=True,gap1=-0.2):
    """
    Aligns traces to each other according to a reference data_set (align_data) and specific trace [ind]

    Args:
        data_arr (list): list containing all of the datasets in the ensemble
        ind (int): integer specifying the trace used for the alignment process
        ref: reference dataset used for alignment
        samp (bool): specify whether a sample of the ensemble should be used instead of the whole ensemble
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
    
    if samp:
        
        #restriction of size of array for testing 
        data_arr0 = data_arr[10:20]
        
        ##print data_arr0
    else:
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
    
    
def find_DTW_match(align_data,data):
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

    data1 = normalise(data)
    align_data1 = normalise(align_data)
        
    #set up sakoe-chiba window length
    r1=len(align_data1)*0.05
    #find optimal warping paths
    pathX,pathY= funcTimeWarp.myDTW(align_data1, data1, derivative=True, costMType='1',D=0.05, bandType = "noband", r=r1, gap=True)
    
    #perform post peak matching?
    linkX0,linkX1 = funcToolsAll.postpeakMatch0(pathX,pathY,step= 100)
    
    return linkX0,linkX1
    

def find_peak_match_X(align_data,data,gap=-0.2):
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
        data = pd.read_csv(file.strip('\n')+'_raw.csv')

        #tidy data
        data_arr.append(data_tidy(data,top,bottom))
    
    return data_arr

def DR_windowing(file_list,TM_peaks,name,top=20,bottom=0,increment=5,windows=10):
    
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

        #iteratively read files and extract data
        for i,file in enumerate(file_list):
            data = pd.read_csv(file.strip('\n')+'_raw.csv')
            if top!=20 and len(TM_peaks[i])==21:
                top1=TM_peaks[i][top]-10-(j)*increment
            else:
                top1=TM_peaks[i][-1]+50+(j)*increment
            #smooth TM trace
            bottom1=TM_peaks[i][bottom]-50-j*increment
            #tidy data and add to final array
            data_arr.append(data_tidy(data,top1,bottom1))
        data_arr1=preprocess(data_arr)
        data_arr2=mobility_shift(data_arr1)
        file_path=name+'_'+str(j)+'.obj'
        print str(bottom1)+'_'+str(top1)
        
        file1=open(file_path,'wb')
        
        pickle.dump(data_arr2,file1)
        
        file1.close()
    
def data_tidy(data,top,bottom):

    """
    Remove all data not in the region of interest (ROI)
    
    Args:
        data (array): Data array containing the electropherogram data
        top (int): Integer specifying the highest elution time point for the ROI
        bottom (int): Integer specifying the lowest elution time point for the ROI
    Returns:
        (array): data array containing the 
        electropherogram data for the ROI only
        
    """
    #remove data above the ROI
    data_out = data[data['Position']<top]
    
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

def peak_finder(data_arr,ind,perc,TM=0,pn=21,cap=None,lower_limit=0):
    
    """
    Peak finding function 
    
    Args: 
        data_arr (list): List containing the chromatograph data
        ind (int): Integer specifying the trace channel 
        perc (float): Ratio of max peak intensity used as cutoff criterion
        TM (int): Set to find the specific number of peaks
        pn (int): Specific number of peaks to find
        cap (int): Cap specifying the highest peak amplitude considered in peak finding 
        lower limit (int): Lowest elution time point considered in peak finding
    Returns:
        (list): List of recorded peak values for each dataset
        
    """
    
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
                if (diff2>0) and (diff1<0):
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
            peak_val2=peak_val[peaks_inds]
            
        else:
            peaks1=peak_pos[peak_val>np.max(peak_val)*perc]
            peak_val2=peak_val[peak_val>np.max(peak_val)*perc]
            
            
            
            j=0
            while j<(len(peaks1)-1):
                diff2=peaks1[j+1]-peaks1[j]

                if diff2<15:

                    av_pos = (peaks1[j]+peaks1[j+1])/2

                    peaks1[j] = av_pos

                    peaks1 = np.delete(peaks1,j+1)
                    j+=1
                j+=1
               
        peak_rec.append(peak_val2)
        peak_arr.append(peaks1)
        ##print peaks1
        ##print i
        
    return peak_arr

def peak_diffs(peak_arr):    
    
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




def S1_partitioning(data_arr,ind):
    
    """
    Partitioning function for the sequencing traces
    
    Args: 
        data_arr (list): All of the datasets in the ensemble
        ind (int): Sequence channel in data for consideration
    Return:
        (list): Partitioned sequence trace peaks for all the datasets in the ensemble
         
    """
    
    #vectors dictationg the marker sizes and differences between markers in the tape measure
    marker_diffs = [10, 30, 10, 20, 30, 10, 20, 10, 10, 20, 20, 20, 20, 10, 10, 20, 20, 20, 20, 20]
    marker_sizes = [50, 60, 90, 100, 120, 150, 160, 180, 190, 200, 220, 240, 260, 280, 290, 300, 320, 340, 360, 380]
    
    #find peaks in tape measure
    peaksTM = peak_finder(data_arr,4,.25,TM=1)
    
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
        
        
        
        ##print number of peaks found
        ##print peaks
        
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
        p_std=findStdW(peaks_trace1['pos'])
        
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

        
def RX_partitioning_single(data_arr,ind,file_list,ll=0,perc=0.25,tm=False,tm_cutoff=21):
    
    """
    partition footprinted trace for single replicates
    
    
    Partitioning function for the sequencing traces
    
    Args: 
        data_arr (list): All of the datasets in the ensemble
        ind (int): Footprint channel in datasets
        file_list (list): File names for the datasets in the ensemble
        ll (int): Lower limit of elution time points to be considered for size marker trace
        perc (float): ratio of maximum peak intensity for cutoff
        tm (bool): Set to find specific number of peaks
        tm_cutoff (int): Highest peak marker to consider 
        
    Return:
        (tuple):
            data_out (list): The partitioned footprinted data traces in numpy.ndarray format
            data_out2 (list): The partitioned footprinted data traces peakList format
            
    """
    
    #vectors dictationg the marker sizes and differences between markers in the tape measure
    marker_diffs = [10, 30, 10, 20, 30, 10, 20, 10, 10, 20, 20, 20, 20, 10, 10, 20, 20, 20, 20, 20]
    marker_sizes = [50, 60, 90, 100, 120, 150, 160, 180, 190, 200, 220, 240, 260, 280, 290, 300, 320, 340, 360, 380]
    
    #find peaks in tape measure
    peaksTM = peak_finder(data_arr,4,perc,TM=tm,lower_limit=ll)
    sam_funcs.sm_plotter(data_arr,peaksTM,file_list)
    #array of nucleotide position
    nuc_pos = np.arange(351)
    
    #initialise peak binning array
    reduced_peaks = []
    data_out = []
    data_out2 = []
    #iterate through data 
    for i in range(len(data_arr)):
        
        
        
        #extract data index
        data = data_arr[i]
        
        #extract peaks
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
       
        #extract peak list for thetarget trace
        peaks_trace1 = funcPeakAlign.fPeakList(data1[:,ind],repType = 'Cubic')
        
        #extract average and standard deviations of widths
        p_width = peaks_trace1['averW']
        p_std=findStdW(peaks_trace1['pos'])
        
        
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

        
        ##print peak_list2
        
        if i ==100:
            fig,ax = plt.subplots(1)
            ax.plot(data[0],data[1],'r',lw=2)
            ax.scatter(peak_list[:,0],peak_list[:,1],color='k',s=50)
            
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xlim(1550,1850)
            ax.set_ylim(-50,2050)


            plt.show()
            
            fig,ax = plt.subplots(1)
            ax.plot(data[0],data[1],'r',lw=2) 
      
            ax.scatter(shoulder_data[:,0],shoulder_data[:,1],color='k',marker='s',s=50)
            ax.set_xlim(1550,1850)
            ax.set_ylim(-50,2050)
            ax.set_xticks([])
            ax.set_yticks([])

            plt.show()
            
            fig,ax = plt.subplots(1)
            ax.plot(data[0],data[1],'r',lw=2) 
            ax.scatter(miss_arr[:,0],miss_arr[:,1],color='k',marker='^',s=50)
            
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xlim(1550,1850)
            ax.set_ylim(-50,2050)



            plt.show()
        if i ==100:
            fig,ax = plt.subplots(1)

            #ax.plot(data[0],data[4],'g')
            ax.plot(data[0],data[1],'r')
            for t in range(len(peaks)):
                ax.plot([peaks[t]]*2,[-100,10000],'k',linestyle='--',alpha=0.7,lw=2)
            ax.set_ylim(-10,2500)
            ax.scatter(peak_list2[:,0],peak_list2[:,1])
            plt.show()
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
            ##print start
            # find the nearest target peak downwind of start TM peak
            start_ind = find_nearest_ind(peak_list2[:,0],start)
            # find the nearest target peak downwind of start TM peak
            end_ind = find_nearest_ind(peak_list2[:,0],end)
            
            #extract target peak list for space between TM peaks
            peak_list3=peak_list2[start_ind:end_ind,:]
            nuc_sep=diff/marker_diffs[j]
            if i ==100:
                fig,ax = plt.subplots(1)
                ax.plot(data[0],data[1],'r',label='Foot#printed Sample',lw=2)
                ax.plot(data[0]-nuc_sep,data[4],'g',label='Size Marker',lw=2)
                ax.scatter(peak_list3[:,0],peak_list3[:,1],color='k',s=50)
                for q in range(marker_diffs[j]):
                    plt.plot([start+(q-1)*nuc_sep]*2,[-100,10000],'k',linestyle='--',alpha=0.7,lw=2)
                
                ax.plot([start+(marker_diffs[j]-1)*nuc_sep]*2,[-100,10000],'k',linestyle='--',alpha=0.7,label='Bin Edge',lw=2)
                ax.legend( prop={'size': 20})
             
                ax.set_xticks([])
                ax.set_yticks([])
    
                plt.show()

            
            
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
        
        #print did_arr1[:,1]
        
        peaks_trace1['NPeak']=len(did_arr1)
                    
        data_out.append(peak_list2)
        
        data_out2.append(peaks_trace1)
    
    return data_out,data_out2





def RX_partitioning_single_500(data_arr0,ind,file_list,ll=0,perc=0.25,tm=0,tm_cutoff=21,Pn=21,Cap=None):
    
    """
    partition footprinted trace for single replicates using ROX 500 size marker set
    
    
    Partitioning function for the sequencing traces
    
    Args: 
        data_arr (list): All of the datasets in the ensemble
        ind (int): Footprint channel in datasets
        file_list (list): File names for the datasets in the ensemble
        ll (int): Lower limit of elution time points to be considered for size marker trace
        perc (float): ratio of maximum peak intensity for cutoff
        tm (bool): Set to find specific number of peaks
        tm_cutoff (int): Highest peak marker to consider 
        
    Return:
        (tuple):
            data_out (list): The partitioned footprinted data traces in numpy.ndarray format
            data_out2 (list): The partitioned footprinted data traces peakList format
            
    """
    
    data_arr=deepcopy(data_arr0)
    #vectors dictationg the marker sizes and differences between markers in the tape measure
    marker_diffs = [25,25,40,10,10,40,50,50,40,10,50,50,40,10]
    marker_sizes = [50, 75, 100, 139, 150, 160, 200, 250,300, 340, 350, 400,450,490,500]
    
    #find peaks in tape measure
    peaksTM = peak_finder(data_arr,4,perc,TM=tm,lower_limit=ll,cap=Cap,pn=Pn)
    sam_funcs.sm_plotter(data_arr,peaksTM,file_list)
    #array of nucleotide position
    nuc_pos = np.arange(351)
    
    #initialise peak binning array
    reduced_peaks = []
    data_out = []
    data_out2 = []
    #iterate through data 
    for i in range(len(data_arr)):
        
        
        
        #extract data index
        data = data_arr[i]
        
        #extract peaks
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
       
        #extract peak list for thetarget trace
        peaks_trace1 = funcPeakAlign.fPeakList(data1[:,ind],repType = 'Cubic')
        
        #extract average and standard deviations of widths
        p_width = peaks_trace1['averW']
        p_std=findStdW(peaks_trace1['pos'])
        
        
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

        
        ##print peak_list2
        
        if i ==100:
            fig,ax = plt.subplots(1)
            ax.plot(data[0],data[1],'r',lw=2)
            ax.scatter(peak_list[:,0],peak_list[:,1],color='k',s=50)
            
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xlim(1550,1850)
            ax.set_ylim(-50,2050)


            plt.show()
            
            fig,ax = plt.subplots(1)
            ax.plot(data[0],data[1],'r',lw=2) 
      
            ax.scatter(shoulder_data[:,0],shoulder_data[:,1],color='k',marker='s',s=50)
            ax.set_xlim(1550,1850)
            ax.set_ylim(-50,2050)
            ax.set_xticks([])
            ax.set_yticks([])

            plt.show()
            
            fig,ax = plt.subplots(1)
            ax.plot(data[0],data[1],'r',lw=2) 
            ax.scatter(miss_arr[:,0],miss_arr[:,1],color='k',marker='^',s=50)
            
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xlim(1550,1850)
            ax.set_ylim(-50,2050)



            plt.show()
        if i ==100:
            fig,ax = plt.subplots(1)

            #ax.plot(data[0],data[4],'g')
            ax.plot(data[0],data[1],'r')
            for t in range(len(peaks)):
                ax.plot([peaks[t]]*2,[-100,10000],'k',linestyle='--',alpha=0.7,lw=2)
            ax.set_ylim(-10,2500)
            ax.scatter(peak_list2[:,0],peak_list2[:,1])
            plt.show()
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
            ##print start
            # find the nearest target peak downwind of start TM peak
            start_ind = find_nearest_ind(peak_list2[:,0],start)
            # find the nearest target peak downwind of start TM peak
            end_ind = find_nearest_ind(peak_list2[:,0],end)
            
            #extract target peak list for space between TM peaks
            peak_list3=peak_list2[start_ind:end_ind,:]
            nuc_sep=diff/marker_diffs[j]
            if i ==100:
                fig,ax = plt.subplots(1)
                ax.plot(data[0],data[1],'r',label='Foot#printed Sample',lw=2)
                ax.plot(data[0]-nuc_sep,data[4],'g',label='Size Marker',lw=2)
                ax.scatter(peak_list3[:,0],peak_list3[:,1],color='k',s=50)
                for q in range(marker_diffs[j]):
                    plt.plot([start+(q-1)*nuc_sep]*2,[-100,10000],'k',linestyle='--',alpha=0.7,lw=2)
                
                ax.plot([start+(marker_diffs[j]-1)*nuc_sep]*2,[-100,10000],'k',linestyle='--',alpha=0.7,label='Bin Edge',lw=2)
                ax.legend( prop={'size': 20})
             
                ax.set_xticks([])
                ax.set_yticks([])
    
                plt.show()

            
            
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
        
        #print did_arr1[:,1]
        
        peaks_trace1['NPeak']=len(did_arr1)
                    
        data_out.append(peak_list2)
        
        data_out2.append(peaks_trace1)
    
    return data_out,data_out2

def RX_partitioning_replicates(data_arr,ind,perc,Cap=None,tm=0,ll=0,tm_cutoff=21,fl=None,Pn=21):
    
    """
    partition footprinted trace for up to 3 replicates
    
    
    Partitioning function for the sequencing traces
    
    Args: 
        data_arr (list): All of the datasets in the ensemble
        ind (int): Footprint channel in datasets
        fl (list): File names for the datasets in the ensemble
        ll (int): Lower limit of elution time points to be considered for size marker trace
        perc (float): ratio of maximum peak intensity for cutoff
        tm (bool): Set to find specific number of peaks
        Pn (int): specific number of peaks to find in peak_finder function
        tm_cutoff (int): Highest peak marker to consider 
        
    Return:
        (tuple):
            data_out (list): The initial bin allocations for footprinting partition
            data_out2 (list): The partitioned footprinting data in discrete packets dictated by the size marker regions
            data_out3 (list): Peaklists of footprinting data for al datasets in the ensemble. 
            
    """
    
    #vectors dictationg the marker sizes and differences between markers in the tape measure
    marker_diffs = [10, 30, 10, 20, 30, 10, 20, 10, 10, 20, 20, 20, 20, 10, 10, 20, 20, 20, 20, 20]
    marker_sizes = [50, 60, 90, 100, 120, 150, 160, 180, 190, 200, 220, 240, 260, 280, 290, 300, 320, 340, 360, 380]
    
    #find peaks in tape measure
    peaksTM = peak_finder(data_arr,4,perc,TM=tm,cap=Cap,lower_limit=ll,pn=Pn)
    
    if fl!=None:
        sam_funcs.sm_plotter(data_arr,peaksTM,fl)
    
    
    ##print peaksTM
    #array of nucleotide position
    nuc_pos = np.arange(351)
    
    #initialise peak binning array
    reduced_peaks = []
    data_out = []
    data_out2 = []
    data_out3=[]
    #iterate through data 
    for i in range(len(data_arr)):
        
        
        
        #extract data index
        data = data_arr[i]
        #print data[0][0]
        
        #extract peaks
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
       
        #extract peak list for thetarget trace
        peaks_trace1 = funcPeakAlign.fPeakList(data1[:,ind],repType = 'Cubic')
        
        #extract average and standard deviations of widths
        p_width = peaks_trace1['averW']
        p_std=findStdW(peaks_trace1['pos'])
        
        
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

        
        ##print peak_list2
        
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
            ##print start
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
    
    return data_out,data_out2,data_out3


def RX_partitioning_replicates_500(data_arr,ind,perc,Cap=None,tm=0,ll=0,tm_cutoff=21,fl=None,Pn=21):
    
    """
    Partition footprinted trace for up to 3 replicates given ROX 500 size marker set
    
    
    Partitioning function for the sequencing traces
    
    Args: 
        data_arr (list): All of the datasets in the ensemble
        ind (int): Footprint channel in datasets
        fl (list): File names for the datasets in the ensemble
        ll (int): Lower limit of elution time points to be considered for size marker trace
        perc (float): ratio of maximum peak intensity for cutoff
        tm (bool): Set to find specific number of peaks
        Pn (int): specific number of peaks to find in peak_finder function
        tm_cutoff (int): Highest peak marker to consider 
        
    Return:
        (tuple):
            data_out (list): The initial bin allocations for footprinting partition
            data_out2 (list): The partitioned footprinting data in discrete packets dictated by the size marker regions
            data_out3 (list): Peaklists of footprinting data for al datasets in the ensemble. 
            
    """
    
    #vectors dictationg the marker sizes and differences between markers in the tape measure
    marker_diffs = [25,25,40,10,10,40,50,50,40,10,50,50,40,10]
    marker_sizes = [50, 75, 100, 139, 150, 160, 200, 250,300, 340, 350, 400,450,490,500]
    
    #find peaks in tape measure
    peaksTM = peak_finder(data_arr,4,perc,TM=tm,cap=Cap,lower_limit=ll,pn=Pn)
    
    if fl!=None:
        sam_funcs.sm_plotter(data_arr,peaksTM,fl)
    
    
    ##print peaksTM
    #array of nucleotide position
    nuc_pos = np.arange(351)
    
    #initialise peak binning array
    reduced_peaks = []
    data_out = []
    data_out2 = []
    data_out3=[]
    #iterate through data 
    for i in range(len(data_arr)):
        
        
        
        #extract data index
        data = data_arr[i]
        #print data[0][0]
        
        #extract peaks
        if tm_cutoff<21:
            peaks = peaksTM[i][:tm_cutoff]
        else:
            peaks = peaksTM[i]
        print i
        print len(peaks)
        if fl!=None:
            print fl[i]
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
        p_std=findStdW(peaks_trace1['pos'])
        
        
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

        
        ##print peak_list2
        
        p_pos2=[]
        p_amp2=[]
        
        did_bin=[]
        bin_counter=[]
        bin_pos=[]

        did_bin3=[]
        bin_alloc2=[]
        for j in range(len(marker_diffs)):
            #extract TM peaks that define the trace region of interest
            tm_bins=0
            
            #extract adjacent peaks in TM
            start = peaks[j]
            end = peaks[j+1]
            
            #calculate distance between TM peaks
            diff=end-start
            ##print start
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
    
    return data_out,data_out2,data_out3



def RX_partitioning_replicates_extended(data_arr,ind,perc,Cap=None,ll=0,tm=0):
    
    """
    partition footprinted trace for up to 3 replicates with extension of 10 nucleotides either side of the trace.
    
    
    Partitioning function for the sequencing traces
    
    Args: 
        data_arr (list): All of the datasets in the ensemble
        ind (int): Footprint channel in datasets
        fl (list): File names for the datasets in the ensemble
        ll (int): Lower limit of elution time points to be considered for size marker trace
        perc (float): ratio of maximum peak intensity for cutoff
        tm (bool): Set to find specific number of peaks
     
        
    Return:
        (tuple):
            data_out (list): The initial bin allocations for footprinting partition
            data_out2 (list): The partitioned footprinting data in discrete packets dictated by the size marker regions
            data_out3 (list): Peaklists of footprinting data for al datasets in the ensemble. 
            
    """
    
    #vectors dictationg the marker sizes and differences between markers in the tape measure
    marker_diffs1 = np.array([10, 30, 10, 20, 30, 10, 20, 10, 10, 20, 20, 20, 20, 10, 10, 20, 20, 20, 20, 20])
    marker_diffs = [10,10, 30, 10, 20, 30, 10, 20, 10, 10, 20, 20, 20, 20, 10, 10, 20, 20, 20, 20, 20,10]
    marker_sizes = [50, 60, 90, 100, 120, 150, 160, 180, 190, 200, 220, 240, 260, 280, 290, 300, 320, 340, 360, 380]
    
    #find peaks in tape measure
    peaksTM = peak_finder(data_arr,4,perc,TM=1,cap=Cap,lower_limit=ll)
    
    ##print peaksTM
    #array of nucleotide position
    nuc_pos = np.arange(351)
    
    #initialise peak binning array
    reduced_peaks = []
    data_out = []
    data_out2 = []
    data_out3=[]
    peaks_TM2=[]
    #iterate through data 
    for i in range(len(data_arr)):
        
        
        
        #extract data index
        data = data_arr[i]
        
        #extract peaks
        peaks = peaksTM[i]
        peak_diffs=np.diff(peaks)
        
        nuc_seps=np.divide(peak_diffs,marker_diffs1.astype(float))
        
        mean_ns=np.mean(nuc_seps)
        #print peaks
        peaks=np.insert(peaks,0,peaks[0]-mean_ns*10)
        peaks=np.append(peaks,peaks[-1]+mean_ns*10)
        
        #print peaks


        
    

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
        p_std=findStdW(peaks_trace1['pos'])
        
        
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

        
        ##print peak_list2
        
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
            ##print start
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
       
    
    return data_out,data_out2,data_out3



def RX_partition_realignment(partition, bin_alloc1,peak_info1,inds,data_arr1,fl=None,corr_b=0.7,inspect=1000,Cap=None,perc=0.25,tm=0,tm_cutoff=21):
    
    """
    Realign partitioned footprinting data based on up to 3 replicates
        
    Args:
        partition (list): Partitioned footprinting data produced by RX_partitioning_replicates
        data_arr (list): All of the datasets in the ensemble
        peak_info1 (list): Peaklists of footprinting data for all datasets in the ensemble
        inds (int): Indices in the ensemble indicating the location of the replicates
        data_arr1 (list): All of the datasets in the ensemble
        fl (list): File names for the datasets in the ensemble
        corr_b (float): Correlation value indicating the minimum pairwise pearson correlation considered in realignment process
        ll (int): Lower limit of elution time points to be considered for size marker trace
        Cap: Maximum peak intensity considered for peak finding function
        perc (float): Ratio of maximum peak intensity for cutoff
        inspect (int): Specify a particular dataset that you want to inspect
        tm (bool): Set to find specific number of peaks
        tm_cutoff (int): Highest peak marker to consider 
        
    Return:
        (list): The aligned partitioned footprinted stored as a peakList object
            
    """
    
    #initialise marker diffs
    marker_diffs = [10, 30, 10, 20, 30, 10, 20, 10, 10, 20, 20, 20, 20, 10, 10, 20, 20, 20, 20, 20]
    marker_sizes = [50, 60, 90, 100, 120, 150, 160, 180, 190, 200, 220, 240, 260, 280, 290, 300, 320, 340, 360, 380]
    parts=[]
    bins=[]
    
    new_arr=[]
    #get letters for replicates
    labs=map(chr, range(65, (65+len(inds))))
    
    partition2=deepcopy(partition)
    peak_info=deepcopy(peak_info1)
    data_arr=deepcopy(data_arr1)
    bin_alloc=deepcopy(bin_alloc1)
    
    peaksTM = peak_finder(data_arr,4,perc,TM=tm,cap=Cap)
    
    
    if fl!=None:
        sam_funcs.sm_plotter(data_arr,peaksTM,fl)
    
    

    
    TM=[]
    for i in range(len(inds)):
        if tm_cutoff<21:
            TM.append(peaksTM[inds[i]][:tm_cutoff])
        else:
            TM.append(peaksTM[inds[i]])
    #extract partitions and bins from lists
    for i in range(len(inds)):
        parts.append(deepcopy(partition2[inds[i]]))
        bins.append(deepcopy(bin_alloc[inds[i]]))
    

    for i in range(len(parts[0])):
        
        bins=marker_diffs[i]
       
        spaces=[]
        spaces_c=[]
        
        #remove intensity nan data
        for j in range(len(inds)):
            spaces.append(parts[j][i])
            spaces_c.append(parts[j][i][~np.isnan(parts[j][i][:,1]),:])
        
        
        bars=[]
        lens=[]
        
        #barcode generation
        for j in range(len(inds)):
            bar = barcode_generator(spaces_c[j][:,1])
            lens=np.append(lens,len(bar))
            bars.append(bar)
            
        if i==inspect:
            print lens
        #if less peaks are observed than expected
        if np.any(lens<bins):     
            
            #if two replicates
            if len(inds)==2:
                
                #perform alignment
                new_arrays,lens1=nw_align(bars,bins,labs)
                
                #if no realigned sequences have correct length
                if np.all(np.array(lens1)==0):
                    
                    #set nans in original sublists to zero 
                    spaceA_new=np.nan_to_num(spaces[0])
                    spaces_c[0]=spaceA_new
                    
                    #re barcode using new data
                    bars[0]=barcode_generator(spaceA_new[:,1])

                    #nw_align on new barcodes
                    new_arrays,lens2=nw_align(bars,bins,labs)
                    
                    #determine correct insertions
                    new_traces=trace_align2(spaces_c,new_arrays)

                # if aligned barcodes with correct length observed.    
                elif np.all(lens1>0):
                    #determine correct alignments
                    new_traces=trace_align2(spaces_c,new_arrays)

            #if 3 replicates used
            else:
                #perform alignment of barcodes
                new_arrays,lens1=nw_align(bars,bins,labs)
                
                #if no alignments have correct length
                if np.all(np.array(lens1)==0):
                    
                    #convert first replicates original sublists nan to 0
                    spaceA_new=np.nan_to_num(spaces[0])
                    spaces_c[0]=spaceA_new
                    
                    #reproduce barcodes
                    bars[0]=barcode_generator(spaceA_new[:,1])
                    
                    #realign barcodes
                    new_arrays,lens2=nw_align(bars,bins,labs)
                    
                    #determine correct insertions
                    zero_len = np.where(lens2==0)[0]
                    
                    if zero_len>0:

                        
                        spaceA_new=np.nan_to_num(spaces[zero_len[0]])
                        bars_new=barcode_generator(spaceA_new[:,1])
                        new_arrs=np.transpose([bars_new,labs[zero_len[0]],5.0])
                        
                        spaces_c[zero_len[0]]=spaceA_new
                        new_arrays[zero_len[0]]=np.array([new_arrs])
                    
                   
                    new_traces=trace_align3(spaces_c,new_arrays)
                    
                #if aligned sequences have correct length
                elif np.all(lens1>0):

                    #determine correct insertion 
                    new_traces=trace_align3(spaces_c,new_arrays)

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
                            new_traces1=trace_align2(spaces1,new_arrays1)
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
                                new_traces=trace_align3(new_traces1,new_arrays2)

            
            
        
                    
        else:
            new_traces=deepcopy(spaces)
        
        if len(inds)==3:
            starts=[TM[0][i],TM[1][i],TM[2][i]]
            ends=[TM[0][i+1],TM[1][i+1],TM[2][i+1]]
            
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


def RX_partition_realignment_500(partition, bin_alloc1,peak_info1,inds,data_arr1,fl=None,corr_b=0.7,inspect=1000,Cap=None,perc=0.25,tm=0,tm_cutoff=21,Pn=21):
    
    """
    Realign partitioned footprinting data based on up to 3 replicates given ROX 500 size marker set
        
    Args:
        partition (list): Partitioned footprinting data produced by RX_partitioning_replicates
        data_arr (list): All of the datasets in the ensemble
        peak_info1 (list): Peaklists of footprinting data for all datasets in the ensemble
        inds (int): Indices in the ensemble indicating the location of the replicates
        data_arr1 (list): All of the datasets in the ensemble
        fl (list): File names for the datasets in the ensemble
        corr_b (float): Correlation value indicating the minimum pairwise pearson correlation considered in realignment process
        ll (int): Lower limit of elution time points to be considered for size marker trace
        Cap: Maximum peak intensity considered for peak finding function
        perc (float): Ratio of maximum peak intensity for cutoff
        inspect (int): Specify a particular dataset that you want to inspect
        tm (bool): Set to find specific number of peaks
        tm_cutoff (int): Highest peak marker to consider
        Pn (int): specific number of peaks to find
        
    Return:
        (list): The aligned partitioned footprinted stored as a peakList object
            
    """
    
    #initialise marker diffs
    marker_diffs = [25,25,40,10,10,40,50,50,40,10,50,50,40,10]
    marker_sizes = [50, 75, 100, 139, 150, 160, 200, 250,300, 340, 350, 400,450,490,500]
    parts=[]
    bins=[]
    
    new_arr=[]
    #get letters for replicates
    labs=map(chr, range(65, (65+len(inds))))
    
    partition2=deepcopy(partition)
    peak_info=deepcopy(peak_info1)
    data_arr=deepcopy(data_arr1)
    bin_alloc=deepcopy(bin_alloc1)
    
    peaksTM = peak_finder(data_arr,4,perc,TM=tm,cap=Cap,pn=Pn)
    
    
    if fl!=None:
        sam_funcs.sm_plotter(data_arr,peaksTM,fl)
    
    

    
    TM=[]
    for i in range(len(inds)):
        if tm_cutoff<21:
            TM.append(peaksTM[inds[i]][:tm_cutoff])
        else:
            TM.append(peaksTM[inds[i]])
    #extract partitions and bins from lists
    for i in range(len(inds)):
        parts.append(deepcopy(partition2[inds[i]]))
        bins.append(deepcopy(bin_alloc[inds[i]]))
    

    for i in range(len(parts[0])):
        
        bins=marker_diffs[i]
       
        spaces=[]
        spaces_c=[]
        
        #remove intensity nan data
        for j in range(len(inds)):
            spaces.append(parts[j][i])
            spaces_c.append(parts[j][i][~np.isnan(parts[j][i][:,1]),:])
        
        
        bars=[]
        lens=[]
        
        #barcode generation
        for j in range(len(inds)):
            bar = barcode_generator(spaces_c[j][:,1])
            lens=np.append(lens,len(bar))
            bars.append(bar)
            
        if i==inspect:
            print lens
        #if less peaks are observed than expected
        if np.any(lens<bins):     
            
            #if two replicates
            if len(inds)==2:
                
                #perform alignment
                new_arrays,lens1=nw_align(bars,bins,labs)
                
                #if no realigned sequences have correct length
                if np.all(np.array(lens1)==0):
                    
                    #set nans in original sublists to zero 
                    spaceA_new=np.nan_to_num(spaces[0])
                    spaces_c[0]=spaceA_new
                    
                    #re barcode using new data
                    bars[0]=barcode_generator(spaceA_new[:,1])

                    #nw_align on new barcodes
                    new_arrays,lens2=nw_align(bars,bins,labs)
                    
                    #determine correct insertions
                    new_traces=trace_align2(spaces_c,new_arrays)

                # if aligned barcodes with correct length observed.    
                elif np.all(lens1>0):
                    #determine correct alignments
                    new_traces=trace_align2(spaces_c,new_arrays)

            #if 3 replicates used
            else:
                #perform alignment of barcodes
                new_arrays,lens1=nw_align(bars,bins,labs)
                
                #if no alignments have correct length
                if np.all(np.array(lens1)==0):
                    
                    #convert first replicates original sublists nan to 0
                    spaceA_new=np.nan_to_num(spaces[0])
                    spaces_c[0]=spaceA_new
                    
                    #reproduce barcodes
                    bars[0]=barcode_generator(spaceA_new[:,1])
                    
                    #realign barcodes
                    new_arrays,lens2=nw_align(bars,bins,labs)
                    
                    #determine correct insertions
                    zero_len = np.where(lens2==0)[0]
                    
                    if zero_len>0:

                        
                        spaceA_new=np.nan_to_num(spaces[zero_len[0]])
                        bars_new=barcode_generator(spaceA_new[:,1])
                        new_arrs=np.transpose([bars_new,labs[zero_len[0]],5.0])
                        
                        spaces_c[zero_len[0]]=spaceA_new
                        new_arrays[zero_len[0]]=np.array([new_arrs])
                    
                   
                    new_traces=trace_align3(spaces_c,new_arrays)
                    
                #if aligned sequences have correct length
                elif np.all(lens1>0):

                    #determine correct insertion 
                    new_traces=trace_align3(spaces_c,new_arrays)

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
                            new_traces1=trace_align2(spaces1,new_arrays1)
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
                                new_traces=trace_align3(new_traces1,new_arrays2)

            
            
        
                    
        else:
            new_traces=deepcopy(spaces)
        
        if len(inds)==3:
            starts=[TM[0][i],TM[1][i],TM[2][i]]
            ends=[TM[0][i+1],TM[1][i+1],TM[2][i+1]]
            
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

def RX_partition_realignment_extended(partition, bin_alloc1,peak_info1,inds,data_arr1,peaks_TM1,corr_b=0.7,inspect=1000,Cap=None,perc=0.25):
    
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
    marker_diffs = [10,10, 30, 10, 20, 30, 10, 20, 10, 10, 20, 20, 20, 20, 10, 10, 20, 20, 20, 20, 20,10]
    marker_sizes = [50, 60, 90, 100, 120, 150, 160, 180, 190, 200, 220, 240, 260, 280, 290, 300, 320, 340, 360, 380]
    parts=[]
    bins=[]
    
    new_arr=[]
    #get letters for replicates
    labs=map(chr, range(65, (65+len(inds))))
    
    partition2=deepcopy(partition)
    peak_info=deepcopy(peak_info1)
    data_arr=deepcopy(data_arr1)
    bin_alloc=deepcopy(bin_alloc1)
    peaksTM=deepcopy(peaks_TM1)
    
    
    
    TM=[]
    for i in range(len(inds)):
        TM.append(peaksTM[inds[i]])

    #extract partitions and bins from lists
    for i in range(len(inds)):
        parts.append(deepcopy(partition2[inds[i]]))
        bins.append(deepcopy(bin_alloc[inds[i]]))
    
    
    for i in range(len(parts[0])):
        
        bins=marker_diffs[i]
       
        spaces=[]
        spaces_c=[]
        
        #remove intensity nan data
        for j in range(len(inds)):
            spaces.append(parts[j][i])
            spaces_c.append(parts[j][i][~np.isnan(parts[j][i][:,1]),:])
        
        
        bars=[]
        lens=[]
        
        #barcode generation
        for j in range(len(inds)):
     
            bar = barcode_generator(spaces_c[j][:,1])
            lens=np.append(lens,len(bar))
            bars.append(bar)
            
        if i==inspect:
            print lens
        #if less peaks are observed than expected
        if np.any(lens<bins):     
            
            #if two replicates
            if len(inds)==2:
                
                #perform alignment
                new_arrays,lens1=nw_align(bars,bins,labs)
                
                #if no realigned sequences have correct length
                if np.all(np.array(lens1)==0):
                    
                    #set nans in original sublists to zero 
                    spaceA_new=np.nan_to_num(spaces[0])
                    spaces_c[0]=spaceA_new
                    
                    #re barcode using new data
                    bars[0]=barcode_generator(spaceA_new[:,1])

                    #nw_align on new barcodes
                    new_arrays,lens2=nw_align(bars,bins,labs)
                    
                    #determine correct insertions
                    new_traces=trace_align2(spaces_c,new_arrays)

                # if aligned barcodes with correct length observed.    
                elif np.all(lens1>0):
                    #determine correct alignments
                    new_traces=trace_align2(spaces_c,new_arrays)

            #if 3 replicates used
            else:
                #perform alignment of barcodes
                new_arrays,lens1=nw_align(bars,bins,labs)
                
                #if no alignments have correct length
                if np.all(np.array(lens1)==0):
                    
                    #convert first replicates original sublists nan to 0
                    spaceA_new=np.nan_to_num(spaces[0])
                    spaces_c[0]=spaceA_new
                    
                    #reproduce barcodes
                    bars[0]=barcode_generator(spaceA_new[:,1])
                    
                    #realign barcodes
                    new_arrays,lens2=nw_align(bars,bins,labs)
                    #print lens2
                    
                    #determine correct insertions
                    #print spaces_c
                    #print i
                    new_traces=trace_align3(spaces_c,new_arrays)
                    
                #if aligned sequences have correct length
                elif np.all(lens1>0):

                    #determine correct insertion 
                    new_traces=trace_align3(spaces_c,new_arrays)

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
                            new_traces1=trace_align2(spaces1,new_arrays1)
                            #add short data back in
                            new_traces1.insert(q,spaces_c[q])
                            new_bars=[]
                            
                            #recalculate barcodes
                            for traces in new_traces1:
                                new_bars = np.append(new_bars,barcode_generator(traces[:,1]))
                            
                            new_arrays2,lens2=nw_align(new_bars,bins,labs)

                            new_traces=trace_align3(new_traces1,new_arrays2)

            
            
        
                    
        else:
            new_traces=deepcopy(spaces)
        
        if len(inds)==3:
            starts=[TM[0][i],TM[1][i],TM[2][i]]
            ends=[TM[0][i+1],TM[1][i+1],TM[2][i+1]]
            
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

    #split bins large than 10 nulceotides into 10 nucleotide bins
    if bins>19:
        split=float(bins)/10
        
        

        for q in range(int(np.ceil(split))):
            nt_split=[]
            new_starts=[]
            new_ends=[]
            if q==np.ceil(split)-1 and split<np.ceil(split):
                #print 'horse'
                for j in range(len(inds)):
                    #create new starts and ends
                    diff=ends[j]-starts[j]
                    new_starts.append(starts[j]+q*10*diff/float(bins))
                    new_ends.append(ends[j])
                    nt_split.append(new_traces[j][(q*10):,:])
                    new_bins=bins-q*10
                    #print new_bins
            else:
                for j in range(len(inds)):
                    #create new starts and ends
                    diff=ends[j]-starts[j]
                    new_starts.append(starts[j]+q*10*diff/float(bins))
                    new_ends.append(starts[j]+(q+1)*10*diff/float(bins))
                    nt_split.append(new_traces[j][(q*10):(q+1)*10,:])
                    new_bins=10
                

            
            #carry out subprocess
            if  inspect and q==2:
                new_traces2=fr_subprocess(nt_split,data_arr1,new_starts,new_ends,new_bins,inds,labs,corr=0.95,corr_b=0.7,inspect=True)
            else:
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
        bars1.append(coarse_grainer(new_traces[j][:,1]))
    cg_corr=np.zeros([len(bars1),len(bars1)])
    for t in range(len(bars1)):
        for j in range(len(bars1)):
            cg_corr[t,j]=pearsonr(bars1[t],bars1[j])[0]


    gggiv,mat=correl_assessor(new_traces,1)
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
        #print hispace
        bars=[]
        for j in range(len(hispace)):
            bars.append(barcode_generator(hispace[j][:,1]))
        #print bars
        new_arrays1,lens1=nw_align(bars,bins,labs)
        #print new_arrays1
        if np.all(np.array(lens1)==0):
            return new_traces
        else:
            #print p_inds
            new_traces1=trace_align3_v2(hispace,new_arrays1,data_arr,inds,bins,starts,ends,inspect=inspect,peak_inds=p_inds)
            
            if new_traces1=='error':
                return new_traces
    elif corr>=corr_b:
        new_traces1=fr_subprocess(new_traces,data_arr,starts,ends,bins,inds,labs,corr=corr-0.05,corr_b=0.7,inspect=inspect)
    else:
        new_traces1=deepcopy(new_traces)
    
 
    return new_traces1
    

def trace_align3(spaces,new_array):
    
    """
    Partitioned data alignment process for 3 replicates (secondary function) 
        
    Args:
        spaces (list): unaligned partitioned footprinting data for each replicate
        new_array (list): List containing new NW alignments for each replicate
        
    Return:
        (list): The aligned partitioned footprinted data for each replicate in the specific part of size marker region.   
    """
    
    cov_arr=[]
    insert_indsA=[]
    insert_indsB=[]
    insert_indsC=[]
    
    new_array1=deepcopy(new_array)
    spaces1=deepcopy(spaces)
    
    #print spaces1
    #print new_array1
    
    #barcodes extract
    new_arrayA=new_array1[0]
    new_arrayB=new_array1[1]
    new_arrayC=new_array1[2]
    
    #print new_array1

    #peak data extract
    spaceA_c=spaces1[0]
    spaceB_c=spaces1[1]
    spaceC_c=spaces1[2]
    
    #iterate thorugh combinations
    for j in range(len(new_arrayA)):
        #find insertion points
        insertA=indexes(new_arrayA[j,0])
        
        #find correction length
        correctA=np.arange(len(insertA))

        #correct insertions
        insertA_f=np.subtract(insertA,correctA)
        
        #create new intensity profile
        spaceA_n=np.insert(spaceA_c,insertA_f.astype(int),np.array((np.nan,0.1)),axis=0)

        for k in range(len(new_arrayB)):
            insertB=indexes(new_arrayB[k,0])
            correctB=np.arange(len(insertB))

            insertB_f=np.subtract(insertB,correctB)
            spaceB_n=np.insert(spaceB_c,insertB_f.astype(int),np.array((np.nan,0.1)),axis=0)

            for l in range(len(new_arrayC)):
                insertC=indexes(new_arrayC[l,0])
                correctC=np.arange(len(insertC))

                insertC_f=np.subtract(insertC,correctC)
                
                spaceC_n=np.insert(spaceC_c,insertC_f.astype(int),np.array((np.nan,0.1)),axis=0)
                
                #calculate PCCs between possible pairwise combinations 
                
                cgA=coarse_grainer(spaceA_n[:,1])
                cgB=coarse_grainer(spaceB_n[:,1])
                cgC=coarse_grainer(spaceC_n[:,1])
                #cov_cg_AB=pearsonr(cgA,cgB)[0]
                #cov_cg_BC=pearsonr(cgC,cgB)[0]
                #cov_cg_AC=pearsonr(cgA,cgC)[0]
                
                cor_arrA=np.nan_to_num(spaceA_n[:,1])
                cor_arrB=np.nan_to_num(spaceB_n[:,1])
                cor_arrC=np.nan_to_num(spaceC_n[:,1])
                
                normA=spaceA_n[:,1]/np.max(spaceA_n[:,1])
                normB=spaceB_n[:,1]/np.max(spaceB_n[:,1])
                normC=spaceC_n[:,1]/np.max(spaceC_n[:,1])
                
                cov_AB=get_cov(spaceA_n[:,1],spaceB_n[:,1])
                
                cov_BC=get_cov(spaceB_n[:,1],spaceC_n[:,1])
                cov_AC=get_cov(spaceA_n[:,1],spaceC_n[:,1])
                
                
                """
                cov_AB=get_cov(normA,normB)
                
                cov_BC=get_cov(normB,normC)
                cov_AC=get_cov(normA,normC)
                """
                
                #summate PCCs
                cov_arr=np.append(cov_arr,cov_AB+cov_BC+cov_AC)
                
                #bin new profiles
                insert_indsA.append(spaceA_n)
                insert_indsB.append(spaceB_n)
                insert_indsC.append(spaceC_n)
    
    #find max summate PCC
    #print cov_arr
    max_cor_ind=np.argmax(cov_arr)
    

    #return intensity profiles for max sum PCC
    newA= insert_indsA[max_cor_ind]
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
        bins (array): number of bins (i.e. number of nucleotides between size markers) expected in each part of the size marker set.
        starts (list): Size marker peaks at the start of the packaged region for each replicate
        ends (list): Size marker peaks at the end of the packaged region for each replicate
        peak_inds (int): Position of bins with assigned peaks before realignment. 
        inspect (bool): Specify whether inspection of data is required
        
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
                if inspect:
                    print peak_distances
               
               

                #summate PCCs
                cov_arr=np.append(cov_arr,cov_AB+cov_BC+cov_AC)
                
                
                
                
                
                #bin new profiles
                insert_indsA.append(spaceA_n)
                insert_indsB.append(spaceB_n)
                insert_indsC.append(spaceC_n)

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
    if inspect:
        
        plt.plot(np.arange(len(newA)),newA[:,1])
        plt.plot(np.arange(len(newB)),newB[:,1])
        plt.plot(np.arange(len(newC)),newC[:,1])
        plt.show()


    
    
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

def trace_align2(spaces,new_array):
    
    """
    Partitioned data alignment process for 2 replicates (secondary function) 
        
    Args:
        spaces (list): unaligned partitioned footprinting data for each replicate
        new_array (list): List containing new NW alignments for each replicate
        
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

    #iterate through combinations
    for j in range(len(new_arrayA)):
        
        #find insertion indices
        insertA=indexes(new_arrayA[j,0])
        correctA=np.arange(len(insertA))

        #correct for insertion indices
        insertA_f=np.subtract(insertA,correctA)

        #create new intensity profile
        spaceA_n=np.insert(spaceA_c,insertA_f.astype(int),np.array((np.nan,0.1)),axis=0)


        for k in range(len(new_arrayB)):
            insertB=indexes(new_arrayB[k,0])
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
                if len(align[0])==bins and len(align[1])==bins:
                    #bin aligned sequences
                    align_bin=np.append(align_bin,align[0])
                    align_bin=np.append(align_bin,align[1])
                    
                    #bin alignment scores
                    score_bin=np.append(score_bin,[align[2]]*2)
                    
                    #bin replicate labels
                    lab_bin=np.append(lab_bin,labels[j])
                    lab_bin=np.append(lab_bin,labels[i])

   
    #create array with all info
    new_array=np.transpose([align_bin,lab_bin,score_bin])
    
    #if nothing in new array return empty array
    if len(new_array)==0:
        output=[[],[],[]]
        lens_out=[0,0,0]
        
    #split array based on replicate    
    else:
        for t in range(len(labels)):
            new_array1=new_array[new_array[:,1]==labels[t],:]
    
            
            if len(new_array1)>0:
                max_score = np.max(new_array1[:,2].astype(float)) 
                
                #filter out sequences with highest score
                new_array2=new_array1[new_array1[:,2].astype(float)==max_score,:]


                #extract unique alignments
                u,u_ind1=np.unique(new_array2[:,0],return_index=True)

                new_array2=new_array2[u_ind1,:]

                #if possible alignments are too numerous, prune
                if len(new_array2)>100:
                    new_array2=new_array2[:100,:]
            else:
                new_array2=new_array1
                
            lens_out=np.append(lens_out,len(new_array2))
            
            #add aligned sequences to output array
            output.append(new_array2)
    
            
    return output,lens_out
            
            
def indexes(string,character='-'):
    """
    Determine the insertion points from aligned sequences
    
    Args:
        string (str): The aligned barcodes
        character (chr): the character to find
    Returns:
        (array): indices where the chosen characters occurs
    """
    
    output=[]
    for i,c in enumerate(string):
        if c==character:
            output = np.append(output,i)
    return np.array(output).astype(int)

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
    return np.array(output).astype(int)
                
        
def barcode_generator(array):
    """
    Barcode generator
    
    Args: 
        array (array): Profile of the partition footprinting data
    Returns:
        (str): Barcode based on profile
        
    """
    #print array
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
    coarse grained profile generator
    
    Args: 
        array (array): Profile of the partition footprinting data
    Returns:
        (array): coarse graining of profile
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
            
def RX_calculator_single(partition_RX,data_arr,RX):
    
    """
    Calculate reactivities of single replicate
    
    Args: 
        Partition_RX (list): Partitioned footprinting data profile in peakList object
        data_arr (list): All of the preprocessed datasets in the ensemble
        RX (int): Index of dataset under investigation
        
    Returns:
        (list): PeakList object containing calculated peak areas and widths
    
    """
   
    new_peak_list=fit_shape_gauss(partition_RX,data_arr[RX])
    
    
    return new_peak_list




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
    #print arr1
    #print arr2
    
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
        (list): peak information with reactivities (peak areas) calculated
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
    array3 = np.absolute(array2[array2<0])
    
    #find argument of data with minimum difference
    ind =array3.argmin()
    
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
            fh_dx_peaks = funcPeakAlign.peakDetection(fh_trace_dx,isY=False)
            
        
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
            sh_dx_peaks = funcPeakAlign.peakDetection(mod_sh_trace_dx,isY=False)
       
            
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
    
    
        
def gaussian_sequence_trace(bin_data_arr):
    
    """
    convert binary array into gaussian trace 
  
    Args:
        bin_data_arr (array): binary array
    Returns:
        (tuple):
            x (array): x positions of trace
            trace_arr (array): Gaussian profile generated
    """
    
    #set position points
    x = np.arange(1500,7501,1)
    
    #set bins 
    bins = np.arange(6,6006,12)
    
    #set gaussian parameters
    height = 100
    p_width = 2
    
    trace_arr = []
    
    #set up gaussian function
    g_h_arr = [height*math.exp(-(i**2)/(2.0*p_width**2)) for i in range(-p_width*3,p_width*3)]
    
    #iterate through binned data
    for bin_data in bin_data_arr:
        
        #set empty array to add gaussian traces to
        trace = [0 for i in range(len(x))]
        
        #read through binary array
        for i in range(len(bin_data[0])):
            ##print int(bin_data[0][i])
            
            #if value of array element is 1 add gaussian to trace at position
            if bin_data[0][i]==1:
                for k in range(len(g_h_arr)):
                    trace[int(bins[i])-int(len(g_h_arr)/2.0)+k] += g_h_arr[k]
        trace_arr.append(trace)        
                
    return x, trace_arr    
                    
    
def position_vote(partition_data,cut1,cut2,plot=0,clip=350):
    
    """
    position voting over the ensemble 
    
    
    Args: 
        partition_data (list): Partitioned sequence traces across ensemble
        cut1 (float): Value of first cutoff used to convert each partitioned trace into a binary sequence
        cut2 (float): Threshold ratio for voting 
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
    convert sequence array to binary based on nucleotide 
    
    Args: 
        seq_arr (str): Nucleotide sequence
        nuc (chr): The nucleotide under investigation
    Returns:
        (array): binary representation of sequence
    """ 
    #convert sequence to binary
    for i in range(len(seq_arr)):
            
        if i == 0:
            if seq_arr[i] == nuc:
                bin_seq = 1
            else:
                bin_seq = 0
        else:
            if seq_arr[i] == nuc:
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
            
            ##print len(ballot_box),len(seq_sec)
            
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
        
    return new_area_arr,new_aver
              
    
def accuracy_measure(seq_arr,vote_arr):
    
    """
    calculates the accuracy of the sequence alignment 
    
    Args: 
        seq_arr (array): Binary representation of 
        sequence
        vote_arr (array): Binary consensus sequence from data
    
    Returns: 
        (float): accuracy of consensus sequence
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
    correlation matrix determination 
   
    Args: 
        data_arr (list): Datasets in the ensemble
        ind (int): Index specifying the trace under investigation
    
    Returns:
        (tuple):
            correl_list (array): list of correlation values
            correl_mat (array): matrix correlation values
            
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
                
            correl_list.append(correl_arr[i,j])
            
    return correl_list,correl_arr


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
    import numpy as np

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


def RX_calculator_replicates(partition,data_arr,inds,get_correls=False):
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
        data = pd.read_csv(file.strip('.fsa')+'_raw.csv')
        plt.plot(data['Position'],data['ReactionChannel#1'],'b',label='RX')
        plt.plot(data['Position'],data['SequenceChannel#1'],'r',label='ddA')
        plt.plot(data['Position'],data['SequenceChannel#2'],'g',label='ddC')
        plt.plot(data['Position'],data['SizeMarker'],'k',label='SM')
        plt.legend()
        plt.savefig(file.strip('\n')+'_raw.png')
        plt.close()
        
def sm_plotter(data_arr,TM_peaks,file_list):
    
    """Plot size marker traces and positions of size markers determined by peak_finder
    
    Args: 
        data_arr (list): Datasets in the ensemble
        TM_peaks (list): Size marker peak positions
        file_list (list): file names of datasets in the ensemble
    
    Returns:
        None
    """
    
    for i in range(len(data_arr)):
        peaks=TM_peaks[i]
        data=data_arr[i]
        file=file_list[i]
        fig,ax=plt.subplots()
        ax.plot(data[0],data[4],'g',lw=2)
        for t in range(len(peaks)):
            ax.plot([peaks[t]]*2,[-100,10000],'k',linestyle='--',alpha=0.7,lw=2)
        
        plt.savefig(file.strip('\n')+'_TM.png')
        plt.close()
            

def findStdW(peakX,rate=0.33):
    """
    Find standard deviation of peak widths
    
    Args:
        peakX (array): array of peak positions 
        rate (float): central percentage of ordered peak widths to calculate standard deviations over
    Returns:
        (float): standard deviations of widths
    
    """
    NPeak=len(peakX)
    diffW=peakX[1:]-peakX[:-1]
    diffW=np.sort(diffW)
    s=int(NPeak*rate)
    e=int(NPeak*(1-rate))
    stdW = np.std(diffW[s:e])
    
    return stdW


ABIF_TYPES = {1: 'byte', 2: 'char', 3: 'word', 4: 'short', 5: 'long', 7: 'float', 8: 'double',\
        10: 'date', 11: 'time', 12: 'thumb', 13: 'bool', 18: 'pString', 19: 'cString'}

class ABIFReader:
    def __init__(self, fn):
        self.filename = fn
        self.file = open(fn, 'rb')
        self.type = self.readNextString(4)
        if self.type != 'ABIF':
            self.close()
            raise SystemExit("error: No ABIF file '%s'" % fn)
        self.version = self.readNextShort()
        dir = DirEntry(self)
        self.seek(dir.dataoffset)
        self.entries = [DirEntry(self) for i in range(dir.numelements)]

    def getData(self, name, num = 1):
        entry = self.getEntry(name, num)
        if not entry:
            raise SystemExit("error: Entry '%s (%i)' not found in '%s'" % (name, num, self.filename))
        self.seek(entry.mydataoffset())
        data = self.readData(entry.elementtype, entry.numelements)
        if data != NotImplemented and len(data) == 1:
            return data[0]
        else:
            return data

    def showEntries(self):
        for e in self.entries:
            print e

    def getEntry(self, name, num):
        for e in self.entries:
            if e.name == name and e.number == num:
                return e
        return None

    def readData(self, type, num):
        if type == 1:
            return [self.readNextByte() for i in range(num)]
        elif type == 2:
            return self.readNextString(num)
        elif type == 3:
            return [self.readNextUnsignedInt() for i in range(num)]
        elif type == 4:
            return [self.readNextShort() for i in range(num)]
        elif type == 5:
            return [self.readNextLong() for i in range(num)]
        elif type == 7:
            return [self.readNextFloat() for i in range(num)]
        elif type == 8:
            return [self.readNextDouble() for i in range(num)]
        elif type == 10:
            return [self.readNextDate() for i in range(num)]
        elif type == 11:
            return [self.readNextTime() for i in range(num)]
        elif type == 12:
            return [self.readNextThumb() for i in range(num)]
        elif type == 13:
            return [self.readNextBool() for i in range(num)]
        elif type == 18:
            return self.readNextpString()
        elif type == 19:
            return self.readNextcString()
        elif type >= 1024:
            return self.readNextUserData(type, num)
        else:
            return NotImplemented

    def readNextBool(self):
        return self.readNextByte() == 1

    def readNextByte(self):
        return self.primUnpack('B', 1)

    def readNextChar(self):
        return self.primUnpack('c', 1)

    def readNextcString(self):
        chars = []
        while True:
            c = self.readNextChar()
            if ord(c) == 0:
                return ''.join(chars)
            else:
                chars.append(c)

    def readNextDate(self):
        return datetime.date(self.readNextShort(), self.readNextByte(), self.readNextByte())

    def readNextDouble(self):
        return self.primUnpack('>d', 8)

    def readNextInt(self):
        return self.primUnpack('>i', 4)

    def readNextFloat(self):
        return self.primUnpack('>f', 4)

    def readNextLong(self):
        return self.primUnpack('>l', 4)

    def readNextpString(self):
        nb = self.readNextByte()
        chars = [self.readNextChar() for i in range(nb)]
        return ''.join(chars)

    def readNextShort(self):
        return self.primUnpack('>h', 2)

    def readNextString(self, size):
        chars = [self.readNextChar() for i in range(size)]
        return ''.join(chars)
    
    def readNextThumb(self):
        return (self.readNextLong(), self.readNextLong(), self.readNextByte(), self.readNextByte())

    def readNextTime(self):
        return datetime.time(self.readNextByte(), self.readNextByte(), self.readNextByte(), self.readNextByte())

    def readNextUnsignedInt(self):
        return self.primUnpack('>I', 4)
    
    def readNextUserData(self, type, num):

        return NotImplemented

    def primUnpack(self, format, nb):
        val=self.file.read(nb)
        x = struct.unpack(format, val )
        return x[0]
    
    def close(self):
        self.file.close()

    def seek(self, pos):
        self.file.seek(pos)

    def tell(self):
        return self.file.tell()

class DirEntry:
    def __init__(self, reader):
        self.name = reader.readNextString(4)
        self.number = reader.readNextInt()
        self.elementtype = reader.readNextShort()
        self.elementsize = reader.readNextShort()
        self.numelements = reader.readNextInt()
        self.datasize = reader.readNextInt()
        self.dataoffsetpos = reader.tell()
        self.dataoffset = reader.readNextInt()
        self.datahandle = reader.readNextInt()

    def __str__(self):
        return "%s (%i) / %s (%i)" % (self.name, self.number, self.mytype(), self.numelements)

    def mydataoffset(self):
        if self.datasize <= 4:
            return self.dataoffsetpos
        else:
            return self.dataoffset

    def mytype(self):
        if self.elementtype < 1024:
            return ABIF_TYPES.get(self.elementtype, 'unknown')
        else:
            return 'user'
        

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
            reader = ABIFReader(fsa_file)
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