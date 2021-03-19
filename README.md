# BoXFP
 **Tool box for analysis of capillary electrophoresis footprinting data**

 ## What is BoXFP?

 BoxFP is a platform independent software package designed for quick, efficient and reliable analysis of capillary electrophoresis data of footprinted nucleic acid samples. 

## What can BoXFP do? 

BoXFP has been developed as an add-on package to the QuShape software packages developed by Dr Fethullah Karabiber of the Weeks Laboratory (https://weeks.chem.unc.edu/qushape/). It utilises a number of algorithms to effectively provide nucleotide specific reactivity information from both single experiments and  over several replicates. Capillary electrophoresis data from SHAPE, Hydroxyl-radical and enzyme based reactivity probing experiments can all be analysed using the package. 

# System Requirements

## Hardware requirements
`BoXFP` package requires only a standard computer with enough RAM to support the in-memory operations.

## Software requirements
### OS requirements
This package is supported for *macOS* and *Linux*. The package has been tested on the following systems:
+ macOS: Mojave (10.14.1)
+ Linux: Ubuntu 16.04

### Python Dependencies
`BoXFP` runs on python 2.7 and requires the following packages:

- `NumPy`: *1.16.3* 
- `Pandas`: *0.24.2*
- `Matplotlib`: *2.2.3*
- `Pickle`: *72223*
- `SciPy`: *1.2.1*
- `Bio`: *1.73*
- `pyqtgraph`: *0.10.0*

# Installation guide

BoXFP can be installed from source from the `Github Repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/PsamClark/BoXFP

Or download the `tarball`_:

.. code-block:: console

    $ curl  -OL https://github.com/PSamClark/BoXFP

Given the requirements outlined above are met installation time will take less than 30 seconds

## Usage

### Python Module

For greater user input options, and for processing multiple datasets at once, the `BoXFP` module can be imported.
FOr more information about scripting using the imported BoXFP module see the [Tutorial](https://github.com/PsamClark/BoXFP/tree/python3_update/Tutorial)

## Support and Contacts

 - The first point of contact for any issues: 
 https://github.com/PsamClark/BoXFP/issues
 
 - Further information can be obtained from the following email address:
 sam.clark@york.ac.uk
 
