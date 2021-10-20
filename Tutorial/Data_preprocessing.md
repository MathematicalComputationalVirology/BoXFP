# Data Pre-processing 

## Script

The script given for this tutorial is `pre_processing_tcv_priB1.py`

The data files for this tutorial can be found in the folder Pre-processing/

## Input data

Input data should be provided as `_raw.csv` files for each dataset in the ensemble. 

## Pre-processing

Pre-processing involves, trace smoothing, baseline correction, decay correction, mobility shifting and trace alignment between foot printing profile and the sequencing ladder.
All these processes are outlined in. 

The most important factor in these pre-processing steps to consider is the region of the chromatograph over which the pre-processing is performed. Pre-processing should not be performed over the whole chromatograph as the high signal intensity of the entry peak (and exit peak if it exists) will disproportionately affect these processes, particularly the decay correction step. pre-processing should therefore be performed over what we define as the region of interest (ROI) after the entry peak (and before the exit peak if it exists). Despite these general rules of thumb some edge effects and signal anomalies may result in similar issues with pre-processing that requires very careful consideration of the ROI. in BoXFP the size marker peaks are used as a point of reference for the region of interest in each data set. To do so an initial searching region must be defined which encompass the size marker peaks. 


To reduce the likelihood of these pre-processing issues occurring pre-processing is performed over several regions of interest. Reactivity profiles are then generated for each of these pre-processing windows and correlations between profiles of windows are calculated pairwise. Those windows that have a poor correlation with the rest of windows on average are removed and an average over the remaining windows are taken. 

For each window the pre-processed data is deposited in to `.obj` pickle data files. 

All the different processes performed during pre-processing can be performed by implementing a single wrapper function `RX_preprocess`.

Pre-processing is performed for all datasets in which the same primer has been used. in the example given the primer B1 TCV samples are being pre-processed:

`xfp.RX_preprocess('B1',1350,None,'test',Top=None)`

+ The first argument specifies the primer to be under investigation. 

+ The second and third arguments specify the start and end of the initial search region. if the third argument is specified as None then the end of the initial search region is the end of the chromatograph. This is mainly to be used if extrapolation beyond the size marker is required. 

+ The second and third arguments specify the start and end of the initial search region. if the third argument is specified as None then the end of the initial search region is the end of the chromatograph. This is mainly to be used if extrapolation beyond the size marker is required. 

+ The fourth argument specifies the prefix to be used to identify the `.obj` data files. 

+ The `Top` argument specifies the size marker peak to use as reference point for the end of the ROI. it the exit peak is located within the set of size marker peaks then a peak below the exit peak should be chosen. if `Top` is not specified the last size marker peak is taken as the reference. if Top is set to None the end of the chromatograph is taken as the reference point. 

Other Arguments that can be specified in the function `RX_preprocess` are:

+ `wins`: The number of windows used in pre-processing. Default is 10.
+ `inc`: The increment in elution points of window size between the windows. Default is 5. 

Running on standard hardware and software the script should take 15-30 minutes to run.









