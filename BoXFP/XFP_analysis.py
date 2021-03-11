import os
import sys
from tkinter import *
import tkinter.font as TKfont
import matplotlib
import matplotlib.pyplot as plt
import glob
import numpy as np
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,NavigationToolbar2Tk)
from matplotlib.figure import Figure
import sam_funcs_v3 as sam_funcs

class StdOutConsole(object):
    def __init__(self,text_area):
        self.text_out=text_area

    def write(self,text):
        self.text_out.insert(END,text)

    def flush(self):
        pass

def set_cwd():
    file_path=fentry.get().strip('\n')
    print(os.stat(file_path))
    foutput.delete(0.0,END)
    if os.path.exists(file_path):
        os.chdir(file_path)
        out_text='Path found!'
        pri_pick()
        pp_list()
    else:
        out_text='Path NOT found!'
    foutput.insert(END,out_text)


def pp_list():
    lb0.delete(0,END)
    file_list=glob.glob('*.obj')

    clabs=[]

    for i in file_list:
        cond=i.split('_')[:-1]

        clabs.append('_'.join(cond))

    uclabs=list(set(clabs))
    for i in uclabs:
        lb0.insert(END,i)

def pri_pick():
    lb1.delete(0,END)
    file_list=glob.glob('*.fsa')

    clabs=[]
    plabs=[]
    elabs=[]
    rlabs=[]
    for i in file_list:
        cond=i.split('_')[0]

        if cond!='Blank':
            pri=i.split('_')[1]
            clabs.append(cond)
            plabs.append(pri)
            elabs.append(i.split('_')[2])
            rlabs.append(i.split('_')[4])

    uclabs=list(set(plabs))
    for i in uclabs:
        lb1.insert(END,i)



def cond_pick(event):
    lb2.delete(0,END)
    file_list=glob.glob('*.fsa')

    clabs=[]
    plabs=[]
    elabs=[]
    rlabs=[]
    pri=lb1.get(lb1.curselection())
    for i in file_list:
        cond=i.split('_')[0]

        if cond!='Blank' and i.split('_')[1]==pri:
            clabs.append(cond)

    uclabs=list(set(clabs))
    for i in uclabs:
        lb2.insert(END,i)


def exp_pick(event):
    lb3.delete(0,END)
    file_list=glob.glob('*.fsa')
    clabs=[]
    plabs=[]
    elabs=[]
    rlabs=[]
    pri=lb1.get(lb1.curselection())
    condition=lb2.get(lb2.curselection())
    for i in file_list:
        cond=i.split('_')[0]

        if (cond==condition and i.split('_')[1]==pri) and i.split('_')[2]!='0':
            clabs.append(i.split('_')[2])

    uclabs=list(set(clabs))
    for i in uclabs:
        lb3.insert(END,i)

def exclusion_info():


    if var1.get()==1:
        window.geometry('830x560')
        bgentry.grid(row=3,column=3,sticky='NESW')
        rxentry.grid(row=3,column=5,sticky='NESW')
        labbg.grid(row=3,column=2,sticky='NESW')
        labrx.grid(row=3,column=4,sticky='NESW')

    elif var1.get()==0:
        window.geometry('720x560')
        bgentry.grid_forget()
        rxentry.grid_forget()
        labbg.grid_forget()
        labrx.grid_forget()


def analyse_wrap():
    file_list=glob.glob('*.fsa')

    if var1.get()==1:
        bgskip=[char for char in bgentry.get()]
        rxskip=[char for char in rxentry.get()]
    elif var1.get()==0:
        bgskip=[]
        rxskip=[]

    bgind=[]
    rxind=[]
    skip=[]
    pp=lb0.get(lb0.curselection())
    pri=lb1.get(lb1.curselection())
    condition=lb2.get(lb2.curselection())
    exposure=lb3.get(lb3.curselection())
    for t,i in enumerate(file_list):
        cond=i.split('_')[0]

        if (cond==condition and i.split('_')[1]==pri):
            if i.split('_')[2]=='0' and i.split('_')[4] not in bgskip:
                bgind.append(t)
            elif i.split('_')[2]=='0' and i.split('_')[4] in bgskip:
                skip.append(t)

            elif i.split('_')[2]==exposure and i.split('_')[4] not in rxskip:
                rxind.append(t)

            elif i.split('_')[2]==exposure and i.split('_')[4] in rxskip:
                skip.append(t)


    virus=ventry.get()
    condName=centry.get()
    extension=int(eentry.get())
    nuc_start=int(nentry.get())
       

    sam_funcs.RX_analyse(pp,bgind,rxind,virus,pri,nuc_start,condName,int(exposure),issues=skip,sm_extend=extension)


print(os.getcwd())
window=Tk()
window.configure(bg='#b5c9c5')
window.title('XFP analysis')
window.geometry('720x560')

myfont=TKfont.Font(family='Helvetica CY',size=11)

lb0=Listbox(window,exportselection=False)
lb0.grid(row=1,column=0,rowspan=5)

lb1=Listbox(window,exportselection=False)
lb1.grid(row=7,column=0,rowspan=5)
lb1.bind('<<ListboxSelect>>',cond_pick)


lb2=Listbox(window,exportselection=False)
lb2.grid(row=14,column=0,rowspan=5)
lb2.bind('<<ListboxSelect>>',exp_pick)


lb3=Listbox(window,exportselection=False)
lb3.grid(row=21,column=0,rowspan=5)


Label(window,text='File list:',bg='#999999',fg='black',width=10,font='helvetica 10').grid(row=0,column=1,sticky='NESW')
fentry=Entry(fg='black',bg='white',width=20)
fentry.grid(row=0,column=2,sticky='NESW')
fbutton = Button(text='find path',bg='red',fg='white',command=set_cwd,width=10,font='helvetica 10').grid(row=0,column=3,sticky='NESW')
foutput=Text(window,width=20,height=1,wrap=WORD,background='white',font='helvetica 10')
foutput.grid(row=0,column=4,sticky='NESW')

Label(window,text='Preprocess files',bg='#999999',fg='black',width=20,font='helvetica 10').grid(row=0,column=0,sticky='NESW')

Label(window,text='Primer',bg='#999999',fg='black',width=20,font='helvetica 10').grid(row=6,column=0,sticky='NESW')

Label(window,text='Condition',bg='#999999',fg='black',width=20,font='helvetica 10').grid(row=13,column=0,sticky='NESW')

Label(window,text='Exposure (ms):',bg='#999999',fg='black',width=20,font='helvetica 10').grid(row=20,column=0,sticky='NESW')


Label(window,text='Virus:',bg='#999999',fg='black',width=10,font='helvetica 10').grid(row=1,column=1,sticky='NESW')
ventry=Entry(window,fg='black',bg='white',width=20)
ventry.grid(row=1,column=2,sticky='NESW')

Label(window,text='Condition name:',bg='#999999',fg='black',width=10,font='helvetica 10').grid(row=1,column=3,sticky='NESW')
centry=Entry(window,fg='black',bg='white',width=20)
centry.grid(row=1,column=4,sticky='NESW')

Label(window,text='Extension:',bg='#999999',fg='black',width=10,font='helvetica 10').grid(row=2,column=1,sticky='NESW')
eentry=Entry(window,fg='black',bg='white',width=20)
eentry.grid(row=2,column=2,sticky='NESW')

Label(window,text='First position:',bg='#999999',fg='black',width=20,font='helvetica 10').grid(row=2,column=3,sticky='NESW')
nentry=Entry(window,fg='black',bg='white',width=20)
nentry.grid(row=2,column=4,sticky='NESW')

labbg=Label(window,text='BG exclusions:',bg='#999999',fg='black',width=20,font='helvetica 10')

labrx=Label(window,text='RX exclusions:',bg='#999999',fg='black',width=20,font='helvetica 10')

bgentry=Entry(fg='black',bg='white',width=20)

rxentry=Entry(fg='black',bg='white',width=20)


var1=IntVar()
var1.set(0)
exclusion_box=Checkbutton(window,text='Exclusions',variable=var1,onvalue=1,offvalue=0,command=exclusion_info,font='helvetica 10')

exclusion_box.grid(row=3,column=1,sticky='NESW')
wbutton = Button(text='RUN',bg='green',fg='Black',command=analyse_wrap,font='helvetica 10').grid(row=5,column=1,sticky='NESW')
out_text=Text(window,background='white',wrap=WORD)
out_text.grid(row=6,column=1,rowspan=20,columnspan=5,sticky='NESW')
scrollb=Scrollbar(window,command=out_text.yview)
scrollb.grid(row=6,column=7,rowspan=20,sticky='NESW')
out_text['yscrollcommand']=scrollb.set
out1=StdOutConsole(out_text)

sys.stdout=out1


window.mainloop()
