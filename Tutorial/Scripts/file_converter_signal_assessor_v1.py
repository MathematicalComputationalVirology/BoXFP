import os
import sys
import struct
import datetime
import numpy as np
import glob
import matplotlib.pyplot as plt
import pandas as pd
import time
import pyqtgraph as pq
import pyqtgraph.exporters


ABIF_TYPES = {1: 'byte', 2: 'char', 3: 'word', 4: 'short', 5: 'long', 7: 'float', 8: 'double',\
        10: 'date', 11: 'time', 12: 'thumb', 13: 'bool', 18: 'pString', 19: 'cString'}

class ABIFReader:
    def __init__(self, fn):
        self.filename = fn
        self.file = open(fn, 'rb')
        self.type = self.readNextString(4)
        if self.type != 'ABIF':
            self.close()
            raise SystemExit("error: No ABIF file '%s'" % fn)
        self.version = self.readNextShort()
        dir = DirEntry(self)
        self.seek(dir.dataoffset)
        self.entries = [DirEntry(self) for i in range(dir.numelements)]

    def getData(self, name, num = 1):
        entry = self.getEntry(name, num)
        if not entry:
            raise SystemExit("error: Entry '%s (%i)' not found in '%s'" % (name, num, self.filename))
        self.seek(entry.mydataoffset())
        data = self.readData(entry.elementtype, entry.numelements)
        if data != NotImplemented and len(data) == 1:
            return data[0]
        else:
            return data

    def showEntries(self):
        for e in self.entries:
            print e

    def getEntry(self, name, num):
        for e in self.entries:
            if e.name == name and e.number == num:
                return e
        return None

    def readData(self, type, num):
        if type == 1:
            return [self.readNextByte() for i in range(num)]
        elif type == 2:
            return self.readNextString(num)
        elif type == 3:
            return [self.readNextUnsignedInt() for i in range(num)]
        elif type == 4:
            return [self.readNextShort() for i in range(num)]
        elif type == 5:
            return [self.readNextLong() for i in range(num)]
        elif type == 7:
            return [self.readNextFloat() for i in range(num)]
        elif type == 8:
            return [self.readNextDouble() for i in range(num)]
        elif type == 10:
            return [self.readNextDate() for i in range(num)]
        elif type == 11:
            return [self.readNextTime() for i in range(num)]
        elif type == 12:
            return [self.readNextThumb() for i in range(num)]
        elif type == 13:
            return [self.readNextBool() for i in range(num)]
        elif type == 18:
            return self.readNextpString()
        elif type == 19:
            return self.readNextcString()
        elif type >= 1024:
            return self.readNextUserData(type, num)
        else:
            return NotImplemented

    def readNextBool(self):
        return self.readNextByte() == 1

    def readNextByte(self):
        return self.primUnpack('B', 1)

    def readNextChar(self):
        return self.primUnpack('c', 1)

    def readNextcString(self):
        chars = []
        while True:
            c = self.readNextChar()
            if ord(c) == 0:
                return ''.join(chars)
            else:
                chars.append(c)

    def readNextDate(self):
        return datetime.date(self.readNextShort(), self.readNextByte(), self.readNextByte())

    def readNextDouble(self):
        return self.primUnpack('>d', 8)

    def readNextInt(self):
        return self.primUnpack('>i', 4)

    def readNextFloat(self):
        return self.primUnpack('>f', 4)

    def readNextLong(self):
        return self.primUnpack('>l', 4)

    def readNextpString(self):
        nb = self.readNextByte()
        chars = [self.readNextChar() for i in range(nb)]
        return ''.join(chars)

    def readNextShort(self):
        return self.primUnpack('>h', 2)

    def readNextString(self, size):
        chars = [self.readNextChar() for i in range(size)]
        return ''.join(chars)
    
    def readNextThumb(self):
        return (self.readNextLong(), self.readNextLong(), self.readNextByte(), self.readNextByte())

    def readNextTime(self):
        return datetime.time(self.readNextByte(), self.readNextByte(), self.readNextByte(), self.readNextByte())

    def readNextUnsignedInt(self):
        return self.primUnpack('>I', 4)
    
    def readNextUserData(self, type, num):

        return NotImplemented

    def primUnpack(self, format, nb):
        val=self.file.read(nb)
        x = struct.unpack(format, val )
        return x[0]
    
    def close(self):
        self.file.close()

    def seek(self, pos):
        self.file.seek(pos)

    def tell(self):
        return self.file.tell()

class DirEntry:
    def __init__(self, reader):
        self.name = reader.readNextString(4)
        self.number = reader.readNextInt()
        self.elementtype = reader.readNextShort()
        self.elementsize = reader.readNextShort()
        self.numelements = reader.readNextInt()
        self.datasize = reader.readNextInt()
        self.dataoffsetpos = reader.tell()
        self.dataoffset = reader.readNextInt()
        self.datahandle = reader.readNextInt()

    def __str__(self):
        return "%s (%i) / %s (%i)" % (self.name, self.number, self.mytype(), self.numelements)

    def mydataoffset(self):
        if self.datasize <= 4:
            return self.dataoffsetpos
        else:
            return self.dataoffset

    def mytype(self):
        if self.elementtype < 1024:
            return ABIF_TYPES.get(self.elementtype, 'unknown')
        else:
            return 'user'
        

def write_out_raw_csv(data_bloc, data_list):
    """
    Writes out the CSV data from the raw FSA file.
    """
    for i,data_file in enumerate(data_list):
            data = data_bloc[i]
	    f = open('%s' % data_file.replace('.fsa', '_raw.csv'), 'w')
	    f.write('Position,ReactionChannel#1,SequenceChannel#1,SequenceChannel#2,SizeMarker\n')
	    for position in range(len(data)):
		f.write('%d,%d,%d,%d,%d\n' % (position+1,
		                              data[position][0],    # reaction channel #1
		                              data[position][1],    # sequencing channel #1
		                              data[position][2],    # reaction channel #2
		                              data[position][3]))   # sequencing channel #2
	    f.close()

def readABI(dir_list):
    data_bloc = []
    for fsa_file in dir_list:
            reader = ABIFReader(fsa_file)
	    col0 = reader.getData('DATA',1)
	    col1 = reader.getData('DATA',2)
	    col2 = reader.getData('DATA',3)
	    col3 = reader.getData('DATA',4)
	    
	    data=np.zeros([len(col0),4],dtype='f4')
	    data[:,0]=np.array(col0)
	    data[:,1]=np.array(col1)
	    data[:,2]=np.array(col2)
	    data[:,3]=np.array(col3)
	    data_bloc.append(data)
    return data_bloc

def signal_assessor(dir_list,virus_name):
	dir_name = os.path.basename(os.getcwd())

	f = open(dir_name+'_signal_assess.csv','w+')
	f.write('sample,virus,av_RC,sd_RC,av_SC1,sd_SC1,av_SC2,sd_SC2\n')

	for fsa_file in dir_list:

		sample_id = fsa_file.strip('.fsa')

		csv_file = sample_id + '_raw.csv'



		data_temp = open(csv_file,'r').readlines()

		data_1 = []

		labels = ['Position','ReactionChannel#1','SequenceChannel#1','SequenceChannel#2','SizeMarker']

		for line in data_temp[1:]:
		    data_1.append([int(i) for i in line.strip('\n').split(',')])



		data = [[i[j] for i in data_1] for j in range(5)] 

		if np.argmax(data[1]) > 5000:

		    quit()

		data = [i[np.argmax(data[1]):] for i in data] 

		sub_samples = [[data[k][j:j+1000] for j in range(0,len(data[k]),1000)] for k in range(5)] 


		avg_arr = [[],[],[],[]]

		for k in range(1,4+1):
		    for sample in sub_samples[k]:
			tmp_arr = sample 
			tmp_arr.sort()
			avg_min = round(np.mean(tmp_arr[:250]),3)
			avg_max = round(np.mean(tmp_arr[-250:]),3)

			avg_arr[k-1].append([avg_min,avg_max,sub_samples[0][sub_samples[k].index(sample)][0],sub_samples[0][sub_samples[k].index(sample)][-1]])

		




		for k in range(0,4):
		    avg_arr[k].pop(np.argmax([i[1] for i in avg_arr[k]])) # - highest high
		    avg_arr[k].pop(np.argmin([i[1] for i in avg_arr[k]])) # - lowest high ?




		sig_strength = np.mean([avg_arr[0][i][1]-avg_arr[0][i][0] for i in range(len(avg_arr[0]))])
		seq_strength = np.mean([avg_arr[1][i][1]-avg_arr[1][i][0] for i in range(len(avg_arr[0]))])
		seq2_strength = np.mean([avg_arr[2][i][1]-avg_arr[2][i][0] for i in range(len(avg_arr[0]))])

		sig_strength2 = np.std([avg_arr[0][i][1]-avg_arr[0][i][0] for i in range(len(avg_arr[0]))])
		seq_strength2 = np.std([avg_arr[1][i][1]-avg_arr[1][i][0] for i in range(len(avg_arr[0]))])
		seq2_strength2 = np.std([avg_arr[2][i][1]-avg_arr[2][i][0] for i in range(len(avg_arr[0]))])

			

		f.write('%s,%s,%f,%f,%f,%f,%f,%f\n'% (sample_id.strip('.fsa'),virus_name,round(sig_strength,3),round(sig_strength2,3),round(seq_strength,3),round(seq_strength2,3),round(seq2_strength,3),round(seq2_strength2,3)))
		
		
	f.close()

def plotter(dir_list):
    for fsa_file in dir_list:

	data = pd.read_csv('%s' % fsa_file.replace('.fsa', '_raw.csv'))

	pos = data['Position']
	RC1 = data['ReactionChannel#1']
	SC1 = data['SequenceChannel#1']
	RC2 = data['SequenceChannel#2']
	SC2 = data['SizeMarker']

	pq.setConfigOption('background', 'w')

	trace_plot = pq.plot()
	trace_plot.addLegend()
	trace_plot.plot(pos,RC1,pen = 'b', name = 'RC1')
	trace_plot.plot(pos,SC1,pen = 'r', name = 'SC1')
	trace_plot.plot(pos,RC2,pen = 'g', name = 'SC2')
	trace_plot.plot(pos,SC2,pen = 'k', name = 'SM')
	


	exporter = pq.exporters.ImageExporter(trace_plot.plotItem)
	exporter.parameters()['width'] = 1800
	exporter.export('%s' % fsa_file.replace('.fsa', '_raw.png'))
     



if __name__ == '__main__':


    start = time.time()

    reference_data_files = glob.glob('*.fsa')
    
    reference_data_bloc = readABI(reference_data_files)

    write_out_raw_csv(reference_data_bloc, reference_data_files)

    signal_assessor(reference_data_files,sys.argv[1])

    plotter(reference_data_files)

    end = time.time()

    print (end-start)

