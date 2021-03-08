# BoXFP
 **Tool box for analysis of capillary electrophoresis footprinting data**

 ## What is BoXFP?

 BoxFP is a platform independent software package designed for quick, efficient and reliable analysis of capillary electrophoresis data of footprinted nucleic acid samples. 

## What can BoXFP do? 

BoXFP has been developed as an add-on package to the QuShape software packages developed by Dr Fethullah Karabiber of the Weeks Laboratory (https://weeks.chem.unc.edu/qushape/). It utilises a number of algorithms to effectively provide nucleotide specific reactivity information from both single experiments and  over several replicates. Capillary electrophoresis data from SHAPE, Hydroxyl-radical and enzyme based reactivity probing experiments can all be analysed using the package. 

## Installation

### Requirements

  * Python 3.3+ or Python 2.7
  * macOS or Linux (Windows not officially supported, but might work)

BoXFP can be installed from source from the `Github Repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/PsamClark/BoXFP

Or download the `tarball`_:

.. code-block:: console

    $ curl  -OL https://github.com/PSamClark/BoXFP
    
## Usage

### Graphical User Interface (GUI)

GUIs for the Preprocessing, Reactivity determination, and Sequence analysis processes have been provided.
All that is required is to run the `Preprocessing.py`, `XFPAnalysis.py` and `SeqAnalysis.py` scripts respectively. 

### Python Module

For greater user input options, and for processing multiple datasets at once, the `BoXFP` module can be imported.
FOr more information about scripting using the imported BoXFP module see the [Tutorial](https://github.com/PSamClark/BoXFP/python3_update/Tutorial)

## Support and Contacts

 - The first point of contact for any issues: 
 https://github.com/PsamClark/BoXFP/issues
 
 - Further information can be obtained from the following email address:
 sam.clark@york.ac.uk
 
