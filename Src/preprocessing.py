import os
import sys
from tkinter import *
import matplotlib
import matplotlib.pyplot as plt
import glob
import numpy as np
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,NavigationToolbar2Tk)
from matplotlib.figure import Figure
import sam_funcs_v3 as sam_funcs


	
def set_cwd():
    file_path=fentry.get()
    foutput.delete(0.0,END)
    if os.path.exists(file_path):
        os.chdir(file_path)
        out_text='Path found!'
    else:
        out_text='Path NOT found!'
    foutput.insert(END,out_text)

def preprocess_wrap():
    file_list=glob.glob('*.fsa')
    print(file_list)
    start=sentry.get()
    end=eentry.get()
    wlab=wentry.get()

    data_arr0 = sam_funcs.data_reader(file_list,end,start)

    data_arr = sam_funcs.preprocess(data_arr0)
    n=len(data_arr0)

    data_arr1=sam_funcs.mobility_shift(data_arr)
    
    #find peaks in TM traces
    peaksTM=sam_funcs.peak_finder(data_arr1,4,.12,cap=3000,TM=1)

    
    #list all those peask that have a disproportionate number of SM peaks
    for i in range(len(peaksTM)):
        
        lenTM=len(peaksTM[i])
        
        if lenTM!=21:
            print (file_list[i]+' '+str(i))
    
    #plot the SM traces with the peak positions marked to make sure the peak finder function has found the right peaks.
    sam_funcs.sm_plotter(data_arr1,peaksTM,file_list)
    
    #run the data reader version two that carries out the windowing and stores the windows in a pickle .obj file
    sam_funcs.DR_windowing(file_list,peaksTM,wlab,top=None) 

def plot_check():
    file_list=glob.glob('*.fsa')
    start=sentry.get()
    end=eentry.get()
    data_arr0 = sam_funcs.data_reader(file_list,end,start)

    n=len(data_arr0)
    #plot the size marker traces to ensure that all peaks have been captured
    fig=Figure()
    a=fig.add_subplot(111)
    color=iter(plt.cm.rainbow(np.linspace(0,1,n)))
    for i in range(n):
        col=next(color)
        a.plot(data_arr0[i]['Position'],data_arr0[i]['SequenceChannel#2'],c=col)
    canvas=FigureCanvasTkAgg(fig,master=window)
    canvas.draw()

    canvas.get_tk_widget().grid(row=4,column=0,columnspan=4)

    toolbarFrame=Frame(master=window)
    toolbarFrame.grid(row=5,column=1,columnspan=2)
    
    toolbar=NavigationToolbar2Tk(canvas,toolbarFrame)
window=Tk()
window.configure(bg='black')
window.title('Preprocessing')
window.geometry('800x800')
Label(window,text='File path:',bg='white',fg='black',width=20).grid(row=0,column=0)
fentry=Entry(fg='black',bg='white',width=20)
fentry.grid(row=0,column=1)
fbutton = Button(text='find path',bg='red',fg='white',command=set_cwd).grid(row=0,column=2)
foutput=Text(window,width=20,height=1,wrap=WORD,background='white')
foutput.grid(row=1,column=0)
Label(window,text='Start:',bg='white',fg='black',width=20).grid(row=2,column=0)
sentry=Entry(window,fg='black',bg='white',width=20)
sentry.grid(row=2,column=1)
eentry=Entry(window,fg='black',bg='white',width=20)
eentry.grid(row=2,column=3)
Label(window,text='End:',bg='white',fg='black',width=20).grid(row=2,column=2)
Label(window,text='wfile label:',bg='white',fg='black',width=20).grid(row=3,column=0)
wentry=Entry(window,fg='black',bg='white',width=20)
wentry.grid(row=3,column=1)
wbutton = Button(text='Check',bg='red',fg='white',command=plot_check).grid(row=3,column=2)

pbutton=Button(text='Preprocess',bg='green',fg='black',command=preprocess_wrap).grid(row=3,column=3)
window.mainloop()
