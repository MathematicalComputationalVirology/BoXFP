import os
import sys
import numpy as np
import math
from copy import deepcopy
import matplotlib.pyplot as plt
import pandas as pd
import pickle


from scipy import stats
from random import sample
from Bio import SeqIO



if __name__ == '__main__':
    
    #sys.path.append(os.path.abspath('/Users/mqbppsc3/Desktop/externalLibraries'))
    sys.path.append(os.path.abspath('/home/samclark/sclark/Cap_elec/software/externalLibraries'))
    import funcFile
    import funcPeakAlign
    import funcSeqAll
    import funcToolsAll
    import funcByRef
    import funcGeneral
    import funcTimeWarp
    import BoXFP as xfp
    
    np.set_printoptions(threshold=sys.maxsize)
    
    xfp.RX_position('210316_tcv_B1_0.obj','TCV_ref_genome.fasta',searchStart=2245,searchSpan=10)
