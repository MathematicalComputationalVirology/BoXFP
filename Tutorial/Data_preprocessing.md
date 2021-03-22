# Data Preprocessing 

## Script

The script given for this tutorial is `pre_processing_tcv_priB1.py`

## Input data

Input data should be provided as `_raw.csv` files for each dataset in the ensemble. 

## Preprocessing

Preprocessing involves, trace smoothing, baseline correction, decay correction, mobiility shifting and trace alignment between footprinting profile and the sequencing ladder.
All of these processes are outlined in. 

THe most important factor in these preprocessing steps to consider is the region of the chromatograph over which the preprocessing is performed. Preprocessing should not be performed over the whole chromatograph as the high signal intensity of the entry peak (and exit peak if it exists) will disproportionately effect these processes, particularly the decay correction step. preprocessing should therefore be performed over what we define as the region of interest (ROI) after the entry peak (and before the exit peak if it exists). Despite these general rules of thumb some edge effects and signal anomalies may result in similar issues with preprocessing that requires very careful consideration of the ROI. in BoXFP the size marker peaks are used as a point of reference for the region of interest in each data set. To do so an initial searching region has to be deifined which encompass the size marker peaks. 


To reduce the likehood of these preprocessing issues occuring preprocessing is performed over several regions of interest. Reactivity profiles are then generated for each of these preprocessing windows and correlations between profiles of windows are calculated pairwise. Those windows that have a poor correlation with the rest of windows on average are removed and an average over the remaining windows are taken. 

For each window the preprocessed data is deposited in to `.obj` pickle data files. 

All of the different processes performed during preprocessing can be performed by implementing a single wrapper function `RX_preprocess`.

Preprocessing is performed for all datasets in which the same primer has been used. in the example given the primer B1 TCV samples are being preprocessed:

`xfp.RX_preprocess('B1',1350,None,'test',Top=None)`

+ The first argument specifies the primer to be under investigation. 

+ The second and third arguments specify the start and end of the initial search region. if the third argument is specified as None then the end of the intial search region is the end of the chromatograph. This is mainly to be used if extrapolation beyond the size marker is required. 

+ The second and third arguments specify the start and end of the initial search region. if the third argument is specified as None then the end of the intial search region is the end of the chromatograph. This is mainly to be used if extrapolation beyond the size marker is required. 

+ The fourth argument specificies the prefix to be used to identify the `.obj` data files. 

+ The Top argument specifies the size marerk peak to use as reference point for the end of the ROI. it the exit peak is located within the set of size marker peaks then a peak below the exit peak should be chosen. if Top is not specified the last size marker peak is taken as the reference. if Top is set to None the enbd of the chromatograph is taken as the reference point. 

Other Arguments that can be specified in the function `RX_preprocess` are:

+ wins: The number of windows used in preprocessing
+ inc: the increment in elution points of window size between the windows





