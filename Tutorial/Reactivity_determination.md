# Reactivity determination

## Scripts and data

The scripts given for this tutorial are `final_calcs_tcv_priB1.py`
The data files for this tutorial can be found in the folder Reactivity_and_postioning/

## Input data

Input data should be provided as the `.obj` pickle data files created during preprocessing.

## Output data

Output data is deposited in several sub folders in the following manner:

- reacs/: Folder containing the calculated reactivity profiles. 
  - reacs.csv: Files containing unnormalised reactivities per nucleotide.
  - reacsNorm.csv: Files containing normalised reactivities per nucleotide.

- repCorrs/: Folder contain replicate correlations.
  - repCorrs.csv: File containing pairwise correlations between replicates.
 
- windowingCorrs/: Folder containing correlation information for windowing.
  - windCorrs.csv: File containing mean pairwise correlations of different preprocessing windows. 
  - wcorrMat.csv: File containing all the pairwise correlations of data from different windows. 

- normFactors/: Folder containing normalisation factors.
  - normFactors.csv: File containing normalisation factors. 

##  Reactivity trace partitioning

Each `.obj` data file is unpacked and the partitioned using the `RX_partitioning` functions. there are two functions created for partitioning of the reactivity traces; `RX_partitioning_replicates` and  `RX_partitioning_single`. 

- `RX_partitioning_single` is designed for data with only a single replicate. 
- `RX_partitioning_replicates` is designed for data with upto and including 3 replicates. 

### Reactivity profile alignment

`RX_partitioning_replicates` outputs the partitioned reativity in several 'packets'. In this form, the data is input into the function `RX_partition_realignment`. This function aligns the partitioned reactivity profiles from the different replicates with each other.

After this process the pearson correlation coefficients between the replicates are calculated and stored in array which will later be output to a csv file. 

## Reactivity calculations 

Reactivities for each nucleus are calculated as the area of the peak associated with that nucleotide. This is acheived using two different functions depending on the number of replicates used. `RX_calculator_single` is used to calculate areas for data with only a single replicate. `RX_calculator_replicates` is used to calculate average peak areas over upto 3 replicates. it also calculates the standard deviations of the peak areas and the average amplitudes of peaks associated with each nucleotide position. 

Following this calculated reactivities must me corrected for background. The scaling factor of the background treatment (0 ms) to the other treatments is calculated using the function `scaleShapeDataWindow` from the QuShape package `funcSeqAll`. The corrected reactivities are then calculated using the `RX_correction` functions, which produces a normalised reactivity profile and the normalisation factor used to produce it. The normalisation factors for are stored and later output into files. 

Errors from the background and the footprinted datasets are propogated using the function `error_propagation`. These errors are then normalised using the normalisation factors used for the reactivity profile. 

## Reactivity profiles averaged over windows

The above processes are performed for each of the preprocessed datasets stored in the `.obj` files. For each preprocessed set the output reactivity profiles are stored in an array.

Pearson correlations coefficients between the reactivity profiles are then calculated and the average correlation values are calculated for each window. Those preprocessed data with an average pearson correlation below 0.75 are discarded. The average over the remaining datasets for each profile is then taken as well as the standard error.

The unnormalised and normalised reactivities are then output as `.csv` files

## RX_analyse:
 
A wrapper function has been created that can perform all of the above processes sequentially using `RX_analyse`. In the example reactivity profiles are being generated for primer B1 extensions of TCV RNA extracted from virion in TCV buffer.

The first step required is to define which of the datasets in the ensemble are the background (0 ms) samples and the exposed reactivity samples (in this case the 25 ms exposure samples):

`
#list indices in data file that correspond to each dataset
A_0=[0,1,2]
A_25=[3,4,5]
`

As can be seen this involves listing the indicies for each, and as such it is often helpful to generate a list of the datasets under investigation. In this tutorial that list is provided in the file `tcv_priB1_data_files.fl`. 

The nucleotide position for the first position in the profile is then specified:

`
nuc_start=2097  
`

note that this value is based on the value generated using the python script file `sequence_alignment_tcv_priB1.py`. Given that the value produced by this script relates to the end of the size marker peaks and we are extrpolating 150 nts beyond the end of the size marker peaks, the `nuc_start` should be the number produced by `sequence_alignment_tcv_priB1.py` minus 150 nts. 


`RX_analyse` can then be called:

`
#Run realignment on partitioned data
xfp.RX_analyse('210316_tcv_B1',A_0,A_25,'TCV','B1',nuc_start,'Extracted_tb',25,sm_extend=15)
`

+ The first argument details the prefix of the `.obj` preprocessed data files. 
+ The second and third arguments are the lists detailing the background and reactivity datasets in the ensemble
+ The fourth and fifth arguments detail the virus and primer under investigation. 
+ The sixth argument gives the nucleotide position of the first position in the reactivity profile. 
+ The seventh and eighth arguments indicate the sample treatment condition and exposure time of the sample. 
+ The ninth and final argument called is `sm_extend` which specifies in 10s of nts how much to extrapolate beyond the size marker set. Note that the default for `sm_extend` is zero and if a negative value is specified, fewer than the given set of size marker peaks are used (if `sm_extend` is set to -n, the first n peaks in the size marker set are used. This should be used if an exit peak is located within or close to the end of the size marker set.

Other arguments that can be used in `RX_analyse` are:

+ `skip`: List of indices of datasets to be disregarded. Default is an empty list. 
+ `wrange`: Specifies the preprocessed data files to be used. Defaults to all datafiles used. 
+ `wcut`: cutoff correlation for windowing. Default is 0.7. 
