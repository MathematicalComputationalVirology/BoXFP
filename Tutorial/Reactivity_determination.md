# Reactivity determination

## Scripts

The scripts given for this tutorial are `final_calcs_v5_1_ms2_pr11.py`, `final_calcs_v5_1_tcv_pr1.py` and `final_calcs_v5_2_ms2_pr1.py`

## Input data

Input data should be provided as the `.obj` pickle data files created during preprocessing.

## Output data

- `data_ad.csv` files containing unnormalised reactivities per nucleotide.
- `data.csv` files containing normalised reactivities per nucleotide.
- `cor.csv` files containing pairwise correlations between replicates.
- `mean_cors.csv` files containing mean pairwise correlations of different preprocessing windows. 
- `cor_mat` files containing all the pairwise correlations of data from different windows. 
- `avers.csv` files containing normalisation factors. 
- `.png` files containing snapshots of the reactivities. 


##  Reactivity trace partitioning

Each `.obj` data file is unpacked and the partitioned using the `RX_partitioning` functions. there are two functions created for partitioning of the reactivity traces; `RX_partitioning_replicates` and  `RX_partitioning_single`. 

- `RX_partitioning_single` is designed for data with only a single replicate. 
- `RX_partitioning_replicates` is designed for data with upto and including 3 replicates. 

### Reactivity profile alignment

`RX_partitioning_replicates` outputs the partitioned reativity in several 'packets'. In this form, the data is input into the function `RX_partition_realignment`. This function aligns the partitioned reactivity profiles from the different replicates with each other.

After this process the pearson correlation coefficients between the replicates are calculated and stored in array which will later be output to a csv file. 

## Reactivity calculations 

Reactivities for each nucleus are calculated as the area of the peak associated with that nucleotide. This is acheived using two different functions depending on the number of replicates used. `RX_calculator_single` is used to calculate areas for data with only a single replicate. `RX_calculator_replicates` is used to calculate average peak areas over upto 3 replicates. it also calculates the standard deviations of the peak areas and the average amplitudes of peaks associated with each nucleotide position. 

Following this calculated reactivities must me corrected for background. The scaling factor of the background treatment (0 ms) to the other treatments is calculated using the function `scaleShapeDataWindow` from the QuShape package `funcSeqAll`. The corrected reactivities are then calculated using the function `reactivity_correction`, which produces a normalised reactivity profile and the normalisation factor used to produce it. The normalisation factors for each 

Errors from the background and the footprinted datasets are propogated using the function `error_propagation`. These errors are then normalised using the normalisation factors used for the reactivity profile. 

## Difference mapping 

to correctly calculate difference maps between the different conditions (in particular in this example between virion and transcript). The background corrected samples must be renormalised based on both profiles combined. this is acheived using the `findPoutlierBox` and `normSimple` functions in QuShapes `funcSeqAll` packages:

```
#Merge Virion and Transcript reactivities
diffs_25=np.append(ad_V_25,ad_R_25)


#calculate percentage of averages and outliers
PO_25,PA_25=funcSeqAll.findPOutlierBox(diffs_25)


#Renormalise Virion and transcript reactivities
nad_V_25_2,aver_V_25_2=funcSeqAll.normSimple(ad_V_25,PO_25,PA_25)
nad_R_25_2,aver_R_25_2=funcSeqAll.normSimple(ad_R_25,PO_25,PA_25)

#calculate difference map
nad_D_25=nad_V_25_2-nad_R_25_2

#calculate new errors
nad_se_V_25_2=ad_se_V_25/aver_V_25_2
nad_se_R_25_2=ad_se_R_25/aver_R_25_2

#calculating errors on difference mapping
nad_se_D_25=BoXFP.error_propagation(nad_se_V_25_2,nad_se_R_25_2)
```

the corrected reactivities, both normalised and unnormalised are stored in arrays along with the difference maps for later use. 

## Reactivity profiles averaged over windows

The above processes are performed for each of the preprocessed datasets stored in the `.obj` files. For each preprocessed set the output reactivity profiles are stored in an array.

Pearson correlations coefficients between the reactivity profiles are then calculated and the average correlation values are calculated for each window. Those preprocessed data with an average pearson correlation below 0.75 are discarded. The average over the remaining datasets for each profile is then taken as well as the standard error.

The unnormalised and normalised reactivities are then output as `.csv` files, and snapshots of the reactivity profiles as `.png` files. 
