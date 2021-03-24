import os
import sys
import numpy as np
import datetime

if __name__ == '__main__':
    

    sys.path.append(os.path.abspath('BoXFP directory'))

    import BoXFP as xfp
    
    start = time.time()
    
    np.set_printoptions(threshold=sys.maxsize)
    
    xfp.RX_position('210315_tcv_B1_0.obj','TCV_ref_genome.fasta',searchStart=2245,searchSpan=10)

    end = time.time()

    print 'Runtime is: '+str(end-start)+' seconds')
    
    
