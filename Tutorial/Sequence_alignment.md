# Genome Position Determination 

## Script

The script given for this tutorial is `sequence_alignment_v1_TCV_pr1.py`

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

A consensus binary is generated using the function `position_vote`. Two cutoffs values, `x1` and `x2`, must be specified in this function. `x1` indicates the cutoff percentage which specifies whether a partitioned peak should be considered a labelled nucleotide in the trace, and then be designated as a 1 in the binary sequence. `x2` indicates the cutoff percentage of votes cast in each position which indicate whether a position is considered a labelled nucleotide or not. For best results these two values which should range between 0 and 1 should be varied iteratively, as in the tutorial script. 

## Genome scanning 

This consensus sequence is then scanned across a binary representation of the reference genome sequence, and the position with the greatest correlation is then isolated and printed out. There are two functions to perform this action. 

- `sequence_search` scans across the entire genome.
- `sequence_search_area` scans across a small section of the genome, with a start point specified and section size. 


