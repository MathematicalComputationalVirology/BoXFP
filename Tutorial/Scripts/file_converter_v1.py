import os
import sys
import datetime
import numpy as np


if __name__ == '__main__':
    
    sys.path.append(os.path.abspath('BoXFP directory')) 

    import BoXFP as xfp


    start = time.time()

    
    reference_data_bloc = xfp.readABI(reference_data_files)

    xfp.write_out_raw_csv(reference_data_bloc, reference_data_files)

    xfp.raw_trace_plotter(reference_data_files)

    end = time.time()

    print 'Runtime is: '+str(end-start)+' seconds')

