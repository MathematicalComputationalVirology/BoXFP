# File conversion and Signal assessment

## Scripts

The scripts given for this tutorial are `file_converter_signal_assessor_v1.py`

## Input data

Input data should be provided as the `.fsa` chromatography files.

## Output data

`_raw.csv` files containing containing the unpacked chromatography data, `_raw.png` images of the chromatograms, and a `.csv` file containing the signal assessment metrics. 

## Script running

In the command line input the following command:

`python file_converter_signal_assessor_v1.py`

this is the only command required.

## Format of unpacked chromatography data

The `_raw.csv` files produced by this script contains 5 columns; *Position*, *ReactionChannel*, *SequenceChannel#1*, *SequenceChannel#2* and *SizeMarker*

## Signal assessment

The signal assessment file produced provides information on the signal quality of each trace in each dataset except the size marker traces. Two measures are used to assess signal quality;

- `av_`: defines the average signal strength of each trace along the chromatograph. 
- `sd_`: defines the standard deviation of the signal strength for each trace along the chromatograph. 