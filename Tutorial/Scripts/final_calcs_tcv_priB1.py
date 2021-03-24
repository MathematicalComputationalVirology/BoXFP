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
from scipy.stats.stats import pearsonr




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
    import sam_funcs_v3 as sam_funcs
    
    np.set_printoptions(threshold=sys.maxsize)
      

        #list indices in data file that correspond to each dataset
    A_0=[0,1,2]
    A_25=[3,4,5]



    skip=[]

    nuc_start=2097  
    #Run realignment on partitioned data

    sam_funcs.RX_analyse('210316_tcv_B1',A_0,A_25,'TCV','B1',nuc_start,'Extracted_tb',25,skip,sm_extend=15)


