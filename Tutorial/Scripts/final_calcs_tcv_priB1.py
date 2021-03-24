import os
import sys
import numpy as np


if __name__ == '__main__':
    
    sys.path.append(os.path.abspath('BoXFP directory'))

    import BoXFP as xfp
    
    np.set_printoptions(threshold=sys.maxsize)
      

    #list indices in data file that correspond to each dataset
    A_0=[0,1,2]
    A_25=[3,4,5]

    skip=[]

    nuc_start=2097  
    #Run realignment on partitioned data

    xfp.RX_analyse('210315_tcv_B1',A_0,A_25,'TCV','B1',nuc_start,'Virion_water',25,skip,sm_extend=15)


