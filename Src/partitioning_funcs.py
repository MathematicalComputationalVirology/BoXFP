from __future__ import print_function
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
import sam_funcs
from random import sample

from mpl_toolkits.mplot3d import Axes3D


from Bio import SeqIO
import Bio.pairwise2 as pwise
from scipy import stats


def TM_partitioning_v1(data_arr,ind):
    
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
    marker_diffs = [10, 30, 10, 20, 30, 10, 20, 10, 10, 20, 20, 20, 20, 10, 10, 20, 20, 20, 20, 20]
    marker_sizes = [50, 60, 90, 100, 120, 150, 160, 180, 190, 200, 220, 240, 260, 280, 290, 300, 320, 340, 360, 380]
    
    #find peaks in tape measure
    peka,peaksTM = peak_finder(data_arr,4,.25,TM=1)
    
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
        
        
        
        #print number of peaks found
        #print peaks
        
        #initials peak lists
        peak_list = []
        pos_list = []
        bin_list = []
               
    
        #iterate through all peaks
        for j in range(len(peaks)-1):
            #extract TM peaks that define the trace region of interest
            start = peaks[j]
            end = peaks[j+1]
            
            #calculate difference between peak postions
            diff = end-start 
            #print diff
            
            #extract number of bins for specific region
            bins = marker_diffs[j]
        
            #contingency for pandas object
            if isinstance(data,pd.DataFrame):
                
                #extract data values
                data = data.values
                
                #extract all of those trace values above TM start position
                S_TM=[p for p in data if p[0]>int(start)]
            
            #contingency for all other objects
            else:
                if j == 0:
                    #transpose data - fudge factor can be improved
                    data=np.transpose(data)
                
                #extract all of those trace values above TM start position
                S_TM=data[data[:,0]>int(start),:]
                
                #print S_TM
            S_TM=np.array(S_TM)
            
            #remove extract those trace values below TM end position
            S_TM=S_TM[S_TM[:,0]<int(end),:]
            
            #print S_TM[:,0]
            S_TM=np.array(S_TM)
            
            #set up the bin width
            nuc_sep = diff/bins
            
            #iterate through bins
            for t in range(bins):
                
                #extract all trace values within bin
                S_bin = S_TM[(S_TM[:,0]>(start+(t)*nuc_sep)),:]
                
                #print S_bin
                #print nuc_sep
                #print start
                #print t
                
                S_bin=S_bin[(S_bin[:,0]<(start+(t+1)*nuc_sep)),:]
                
                #print S_bin
                
                #initialise binning array
                bin_peaks = []
                
                #calculate derivative for trace in bin
                for k in range(1,(len(S_bin)-1)):
                    diff1 = S_bin[k,ind]-S_bin[k-1,ind]
                    diff2 = S_bin[k+1,ind]-S_bin[k,ind]
                    
                    #if derivative crosses 0 with -ve slope bin peak amp and position 
                    if (diff1>0) and (diff2<0):
                        bin_peaks.append(S_bin[k,ind])
                
                #if derivative does not cross 0 with -ve slope take max value in bin and average postion in bin 
                if len(bin_peaks)==0:
                    peak_val=np.max(S_bin[:,ind])
                    peak_pos = start+(t-0.5)*nuc_sep
                
                #take max peak in bin and its position
                else:
                    peak_val = np.max(bin_peaks)
                    peak_pos = S_bin[:,0][S_bin[:,ind]==peak_val]
                    
                
                #bin peak value
                peak_list = np.append(peak_list,peak_val)
                
                #bein peak position
                pos_list = np.append(pos_list,peak_pos)
                #print pos_list
        
        
        #locate postion of binned data in relation to bin list
        locator = find_first(marker_diffs,bin_list)
        
        #set start nucleotide of binned data 
        marker_start = marker_sizes[locator]-marker_sizes[0]
        nuc_pos = np.array(nuc_pos)+marker_start
        
        #remove nucleotide values larger than 350 
        nuc_pos = nuc_pos[nuc_pos<350]
    
        #assemble partitioned peak array
        
        if len(peak_list)>349:
            peak_arr = [nuc_pos,pos_list,peak_list]
        
        reduced_peaks.append(peak_arr)
                
        
    return reduced_peaks

def TM_partitioning_v2_1(data_arr,ind):
    
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
    marker_diffs = [10, 30, 10, 20, 30, 10, 20, 10, 10, 20, 20, 20, 20, 10, 10, 20, 20, 20, 20, 20]
    marker_sizes = [50, 60, 90, 100, 120, 150, 160, 180, 190, 200, 220, 240, 260, 280, 290, 300, 320, 340, 360, 380]
    
    #find peaks in tape measure
    peka,peaksTM = peak_finder(data_arr,4,.2,TM=1)
    
    #array of nucleotide position
    nuc_pos = np.arange(351)
    
    #initialise peak binning array
    reduced_peaks = []
    shoulder_arr = []
    #iterate through data 
    for i in range(len(data_arr)):
        
        #extract data index
        data = data_arr[i]
        
        #extract peaks
        peaks = peaksTM[i]
        

        #
        
        #print number of peaks found
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
 
    

        peaks_trace1 = funcPeakAlign.fPeakList(data1[:,ind],repType = 'Cubic')

        p_ind = peaks_trace1['pos']
        
        p_pos = np.unique(data1[p_ind,0])
        
        
    
        p_amp = peaks_trace1['amp']
        
        p_width = peaks_trace1['averW']
        
        shoulder_pos_arr = []
        shoulder_amp_arr = []
    
        
        for j in range(len(peaks)-1):
            #extract TM peaks that define the trace region of interest
            start = peaks[j]
            end = peaks[j+1]
            
            peaks_in_diff=p_pos[p_pos>start-p_width]
            
            pamps_in_diff=p_amp[p_pos>start-p_width]
            
            peaks_in_diff=peaks_in_diff[peaks_in_diff<end-p_width]
            
            pamps_in_diff=pamps_in_diff[peaks_in_diff<end-p_width]
            
            peak_arr = peaks_in_diff
    
            for k in range(len(peak_arr)-1):
                    #print peak_arr

                    peak1 = peak_arr[k]
                    peak2 = peak_arr[k+1]

                    trace_between_peaks=deepcopy(data1[data1[:,0]>peak1])

                    trace_between_peaks1 = deepcopy(trace_between_peaks[trace_between_peaks[:,0]<peak2])

                    bp_len = len(trace_between_peaks1)

                    if bp_len ==0:
                        break

                    trough = np.argmin(trace_between_peaks1[:,ind])
                    print(trough)

                    #print trace_between_peaks1

                    first_half=deepcopy(trace_between_peaks1[0:trough,:])

                    second_half = deepcopy(trace_between_peaks1[int(trough):int(bp_len),:])



                    if len(first_half)>2:
                        fh_trace_dx = funcGeneral.deriv1(first_half[:,ind])

                        fh_dx_peaks = funcPeakAlign.peakDetection_v2(fh_trace_dx,isY=False)

                        if len(fh_dx_peaks)>0:

                            #if fh_dx_peaks ==len(first_half):
                                #break

                            #print fh_dx_peaks

                            shoulder_pos = first_half[fh_dx_peaks,0]

                            shoulder_amp = first_half[fh_dx_peaks,ind]

                            shoulder_pos = shoulder_pos[shoulder_amp>0]

                            shoulder_amp = shoulder_amp[shoulder_amp>0]

                            shoulder_pos_arr = np.append(shoulder_pos_arr,shoulder_pos)

                            shoulder_amp_arr = np.append(shoulder_amp_arr,shoulder_amp)

                            print('shoulder')
                            if i == 38:
                                
                                print(peaks[0])
                                print(fh_trace_dx)
                                print(first_half)

                                fig,(plt1,plt2)=plt.subplots(2)

                                plt1.plot(first_half[:,0],fh_trace_dx)
                                plt2.plot(first_half[:,0],first_half[:,ind])

                                plt.show()


                    if len(second_half)>2:

                        sh_trace_dx = funcGeneral.deriv1(second_half[:,ind])

                        mod_sh_trace_dx = -np.absolute(sh_trace_dx)


                        sh_dx_peaks = funcPeakAlign.peakDetection_v2(mod_sh_trace_dx,isY=False)

                        if len(sh_dx_peaks)>0:

                            #if sh_dx_peaks ==len(second_half):
                                #break

                            shoulder_pos = second_half[sh_dx_peaks,0]

                            shoulder_amp = second_half[sh_dx_peaks,ind]

                            shoulder_pos = shoulder_pos[shoulder_amp>0]

                            shoulder_amp = shoulder_amp[shoulder_amp>0]


                            shoulder_pos_arr = np.append(shoulder_pos_arr,shoulder_pos)

                            shoulder_amp_arr = np.append(shoulder_amp_arr,shoulder_amp)

                            """
                            print sh_trace_dx.astype(int)
                            print second_half
                            if i==35:
                                fig,(plt1,plt2)=plt.subplots(2)

                                plt1.plot(second_half[:,0],sh_trace_dx)
                                plt2.plot(second_half[:,0],second_half[:,ind])

                                plt.show()

                            """

                            print('shoulder')

            shoulder_data=[shoulder_pos_arr,shoulder_amp_arr]
        shoulder_arr.append(shoulder_data)

    return shoulder_arr

def RX_partitioning_single(data_arr,ind):
    
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
    marker_diffs = [10, 30, 10, 20, 30, 10, 20, 10, 10, 20, 20, 20, 20, 10, 10, 20, 20, 20, 20, 20]
    marker_sizes = [50, 60, 90, 100, 120, 150, 160, 180, 190, 200, 220, 240, 260, 280, 290, 300, 320, 340, 360, 380]
    
    #find peaks in tape measure
    peka,peaksTM = peak_finder(data_arr,4,.25,TM=1)
    
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

        
        #print peak_list2
        
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
            #print start
            # find the nearest target peak downwind of start TM peak
            start_ind = find_nearest_ind(peak_list2[:,0],start)
            # find the nearest target peak downwind of start TM peak
            end_ind = find_nearest_ind(peak_list2[:,0],end)
            
            #extract target peak list for space between TM peaks
            peak_list3=peak_list2[start_ind:end_ind,:]
            nuc_sep=diff/marker_diffs[j]
            if i ==100:
                plt.plot(data[0],data[1])
                plt.scatter(peak_list3[:,0],peak_list3[:,1],color='k')
                for q in range(marker_diffs[j]):
                    plt.plot([start+(q-1)*nuc_sep]*2,[-100,10000],'k',linestyle='--',alpha=0.7)
                
                plt.plot([start+(marker_diffs[j]-1)*nuc_sep]*2,[-100,10000],'k',linestyle='--',alpha=0.7,label='Bin Edge')
                plt.legend()
                plt.xlabel('Elution Time')
                plt.ylabel('Fluorescence Intensity')
                plt.xlim(1550,1850)
                plt.ylim(-100,4000)
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
        peaks_trace1['pos']=did_arr1[:,0].astype(int)
        peaks_trace1['amp']=did_arr1[:,1]
        
        peaks_trace1['NPeak']=len(did_arr1)
                    
        data_out.append(peak_list2)
        
        data_out2.append(peaks_trace1)
    
    return data_out,data_out2

def RX_partition_realignment(partition, bin_alloc,peak_info,inds,data_arr):
    marker_diffs = [10, 30, 10, 20, 30, 10, 20, 10, 10, 20, 20, 20, 20, 10, 10, 20, 20, 20, 20, 20]
    marker_sizes = [50, 60, 90, 100, 120, 150, 160, 180, 190, 200, 220, 240, 260, 280, 290, 300, 320, 340, 360, 380]
    
    partA=partition[inds[0]]
    
    binsA=bin_alloc[inds[0]]
    
    partB=partition[inds[1]]
    
    binsB=bin_alloc[inds[1]]
    
    partC=partition[inds[2]]
    
    binsC=bin_alloc[inds[2]]
    
    new_arrA=[]
    new_arrB=[]
    new_arrC=[]
    

    for i in range(len(partA)):
        
        spaceA=partA[i]
        spaceB=partB[i]
        spaceC=partC[i]
        
        spaceA_c=spaceA[~np.isnan(spaceA[:,1]),:]
        spaceB_c=spaceB[~np.isnan(spaceB[:,1]),:]
        spaceC_c=spaceC[~np.isnan(spaceC[:,1]),:]
        
        bins=marker_diffs[i]
        """
        plt.plot(np.arange(bins),spaceA[:,1])
        plt.plot(np.arange(bins),spaceB[:,1])
        plt.plot(np.arange(bins),spaceC[:,1])
        plt.show()
        """
    
        
        spaces=[spaceA_c,spaceB_c,spaceC_c]
        #print spaceB
        #print spaceB_c
        #print i
        barA=barcode_generator(spaceA_c[:,1])
        barB=barcode_generator(spaceB_c[:,1])
        barC=barcode_generator(spaceC_c[:,1])
        
        bars=[barA,barB,barC]
        #print bars
        lens=np.array([len(barA),len(barB),len(barC)])
        
        #print 'bins '+str(bins)
        #print 'lens '+str(lens)
        
        labs=['A','B','C']
        
        if np.any(lens<bins):
            new_arrays,lens1=nw_align(bars,bins,labs)
            #print lens1
            if np.all(np.array(lens1)==0):
                #print lens1
                spaceA_c=np.nan_to_num(spaceA)
                spaces[0]=spaceA_c
                bars[0]=barcode_generator(spaceA_c[:,1])
                new_arrays,lens2=nw_align(bars,bins,labs)
                #print lens2
                new_traces=trace_align3(spaces,new_arrays)
                
            
            elif np.all(lens1>0):
                new_traces=trace_align3(spaces,new_arrays)
            else:    
                for q,arr in enumerate(new_arrays):
                    
                    if lens1[q]==0 and np.all(np.delete(lens1,q)>0):
                        new_arrays1=deepcopy(new_arrays)
                        del new_arrays1[q]
                        spaces1=deepcopy(spaces)
                        del spaces1[q]
                        new_traces1=trace_align2(spaces1,new_arrays1)
                        #print new_traces1
                        new_traces1.insert(q,spaces[q])
                        #print spaces[q]
                        #print new_traces1
                        new_bars=[]
                        for traces in new_traces1:
                            #print traces
                            new_bars = np.append(new_bars,barcode_generator(traces[:,1]))
                        
                        new_arrays2,lens2=nw_align(new_bars,bins,labs)
                    
                        new_traces=trace_align3(new_traces1,new_arrays2)

            #print new_traces[0]
            #print 'i '+str(i)
            #print new_arrA
            
            if i==0:
                new_arrA=new_traces[0]
                new_arrB=new_traces[1]
                new_arrC=new_traces[2]
            else:
                
                new_arrA=np.vstack((new_arrA,new_traces[0]))
                new_arrB=np.vstack((new_arrB,new_traces[1]))
                new_arrC=np.vstack((new_arrC,new_traces[2]))
               
        else:
            if i==0:
                new_arrA=spaceA
                new_arrB=spaceB
                new_arrC=spaceC
            else:
                new_arrA=np.vstack((new_arrA,spaceA))
                new_arrB=np.vstack((new_arrB,spaceB))
                new_arrC=np.vstack((new_arrC,spaceC))
        #print 'new_arrA'
        #print new_arrA
    
    dataA=data_arr[inds[0]]
    dataB=data_arr[inds[1]]
    dataC=data_arr[inds[2]]
    
    print(new_arrA[:,0])
    
    posA=nan_remover(new_arrA[:,0])
    posB=nan_remover(new_arrB[:,0])
    posC=nan_remover(new_arrC[:,0])
    
    print(posA)
    print(len(posA))
    
    pos_indA=np.in1d(dataA[0],posA.astype(int)).nonzero()[0]
    pos_indB=np.in1d(dataB[0],posB.astype(int)).nonzero()[0]
    pos_indC=np.in1d(dataC[0],posC.astype(int)).nonzero()[0]
    
    print(pos_indA)
    print(len(pos_indA))
    
    
    peak_infoA=peak_info[inds[0]]
    peak_infoA['amp']=new_arrA[:,1]
    peak_infoA['pos']=pos_indA
    peak_infoA['NPeak']=len(new_arrA)
    
    peak_infoB=peak_info[inds[1]]
    peak_infoB['amp']=new_arrB[:,1]
    peak_infoB['pos']=pos_indB
    peak_infoB['NPeak']=len(new_arrB)
    
    peak_infoC=peak_info[inds[2]]
    peak_infoC['amp']=new_arrC[:,1]
    peak_infoC['pos']=pos_indC
    peak_infoC['NPeak']=len(new_arrC)
    
    return [peak_infoA,peak_infoB,peak_infoC]




