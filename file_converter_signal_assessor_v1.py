import os
import sys
import struct
import datetime
import numpy as np
import glob
import matplotlib.pyplot as plt
import pandas as pd
import time
import struct




if __name__ == '__main__':
    
    sys.path.append(os.path.abspath('/Users/mqbppsc3/Desktop/Functions')) #sys.path.append(os.path.abspath('/home/samclark/sclark/Cap_elec/software/externalLibraries'))
    import funcFile
    import funcPeakAlign
    import funcSeqAll
    import funcToolsAll
    import funcByRef
    import funcGeneral
    import funcTimeWarp
    import BoXFP


    start = time.time()

    reference_data_files = glob.glob('*.fsa')
    #print reference_data_files
    
    reference_data_bloc = BoXFP.readABI(reference_data_files)

    BoXFP.write_out_raw_csv(reference_data_bloc, reference_data_files)

    BoXFP.signal_assessor(reference_data_files)

    BoXFP.raw_trace_plotter(reference_data_files)

    end = time.time()

    print (end-start)

