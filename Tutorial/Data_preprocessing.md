# Data Preprocessing 

## Script

The script given for this tutorial is `pre_processing.py`

## Input data

Input data should be provided as `_raw.csv` files for each dataset in the ensemble. 

## Initial data read-in 

In the first instance a text file listing all of the datasets in the ensemble. this can be created by typing the following command in the  command line:

`for i in *raw.csv; do echo $i | sed -e "s/.fsa//g" ; done > data_files.fl`

this data file is then read in as an array which is used to open all of the dataset files intially using the function `data_reader`. This function requires the specification of the start and end of the initial region of interest (ROI), which should cover the peaks of the size marker trace for all datasets in the ensemble: 

`data_arr0 = sam_funcs.data_reader(file_list,6000,1400)`

correct placement of the ROI in the script provided is determined by plotting all of the size marker traces onto a sinlge plot, to ensure that size marker regions for each of the datasets is covered by the initial ROI:

```
n=len(data_arr0)

#plot the size marker traces to ensure that all peaks have been captured
color=iter(plt.cm.rainbow(np.linspace(0,1,n)))
for i in range(n):
    col=next(color)
    plt.plot(data_arr0[i]['Position'],data_arr0[i]['SizeMarker'],c=col)
plt.show()
```
preprocessing is then carried out over this initial ROI for each dataset using the `preprocess` function. 

## ROI determination for individual datasets 

Peak positions in the size marker trace are determined using the function `peak_finder`. The plotting function `sm_plotter` is then used to plot the size marker traces with the peak position annotated on it for reference. These plots should be consulted to make sure that the `peak_finder` function has correctly located all the peaks in the size marker trace. 

## Preprocessing over windows

With the positions of the size marker peaks determined for each dataset in the ensemble, the function `data_reader_v2` to generate preprocessed versions of the datasets in the ensemble over several windows in the data. This is performed to ensure that error due to preprocessing has been reduced. by defualt 10 windows are used that increase in size from 5 elution points either side of the size marker region to 50 elution points either side of the size marker region. The preprocessed data is then stored in several pickle `.obj` files, each with a specified name for the entire batch and an index t indicate the window. 



