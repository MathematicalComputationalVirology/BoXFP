# Genome Position Determination 

## Script

The script given for this tutorial is `sequence_alignment_tcv_prB1.py`

The data files for this tutorial can be found in the folder Reactivity_and_positioning/

## Input data

Input data should be from one of the pickle `.obj` files created during preprocessing of the data. 
The reference genome sequence should be provided in a text file in FASTA format. 

## Data partitioning 

The ddA ladder trace in the data is partitioned using the `S1_partitioning` function. Correlations between partitioned ladder traces are assessed using the function `correl_assessor`. The mean correlation is then calculated for each dataset, and these means are plotted as a histogram,

```
sum_array=S_cov_matrix.mean(axis=1)

plt.hist(sum_array)
plt.show()
```

Datasets which have a poor correlation with the other datasets in the ensemble can be removed if requested,

```
keep_seqs=np.where(sum_array>cutoff)[0]

S1_partition1 = [S1_partition[i] for i in keep_seqs]

```

## Consensus sequence generation. 

A consensus binary is generated using the function `position_vote`. Two cut-off values, `x1` and `x2`, must be specified in this function. `x1` indicates the cut-off percentage which specifies whether a partitioned peak should be considered a labelled nucleotide in the trace, and then be designated as a 1 in the binary sequence. `x2` indicates the cut-off percentage of votes cast in each position which indicate whether a position is considered a labelled nucleotide or not. For best results these two values which should range between 0 and 1 should be varied iteratively, as in the tutorial script. 

## Genome scanning 

This consensus sequence is then scanned across a binary representation of the reference genome sequence, and the position with the greatest correlation is then isolated and printed out. There are two functions to perform this action. 

- `sequence_search` scans across the entire genome.
- `sequence_search_area` scans across a small section of the genome, with a start point specified and section size. 

## RX_position

All the above processes can be run in sequence using the wrapper function `RX_position`. In the example given we are using `RX_position` to find the first nucleotide position in the size marker trace. Using the information on the primer, we already know the start position of the primer (2647), and we know the length of the size marker trace based on the set being used (ROX-400HD), so we can therefore assume that the start of the size marker trace should be around the 2247 nt position. However due to mobility shifting the exact position may be up to 5 nts from that position. To get the exact position `RX_position` is used in the following manner:

`xfp.RX_position('210315_tcv_B1_0.obj','TCV_ref_genome.fasta',searchStart=2242,searchSpan=10)`

+ The first argument specifies the pre-processed data file to be used. Note that for position determination only one of the pre-processed data files is required. 
+ The second argument specifies the name of the reference genome file that the ddA ladders will be compared to.
+ The third argument specifies where in the reference genome to start the search. The default value for `searchStart` is 0.
+ The fourth argument specifies how many nts positions after the start position are to be analysed. Note that if `SearchSpan` is not specified the entire genome will be searched. 

Other Arguments that can be used in `RX_position` are:

+ `seqInd`: specify which channel in the chromatograph contains the sequence ladder to use. Default 2nd channel.
+ `dinds`: specify which datasets are to be used. Default all datasets are used. 
+ `clip`: specify how much of the size marker region is to be used for the position determination. Defaults to the entire size marker region. 
+ `corr`: Correlation cut-off for ignoring poorly correlating sequence traces. Default is 0.7.
+ `voteInc`: Specify how many intervals between 0 and 100% the `x1` and `x2` cut-offs should be iterated through. Default value is 20. 

Note that if you are extrapolating beyond the given size marker the start position produced needs to be offset to reflect this. For instance, in the case of the primer B1 example, the final reactivity profile will extend 150 nts beyond the size marker region, so the position determined must be shifted -150 nts. 

On standard hardware and software, it should take 15 minutes to run the script.


