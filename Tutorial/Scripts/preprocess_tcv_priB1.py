import os
import sys
import numpy as np

if __name__ == '__main__':
    
    sys.path.append(os.path.abspath('BoXFP directory'))

    import BoXFP as xfp
    
    np.set_printoptions(threshold=sys.maxsize)

    xfp.RX_preprocess('B1',1350,None,'test',Top=None)
      

