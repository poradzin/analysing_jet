#!/bin/env python3
#fitz 2015
# Script to look at KK3 data to find sawteeth
# usage: ./toother.py <shotNumber>

import ppf
import numpy as np
import matplotlib.widgets as wid
import matplotlib.pyplot as plt
import sys

pulse=int(sys.argv[1])
seq=0
dda="KK3"

ierr=ppf.ppfgo(pulse=pulse,seq=seq)
ppf.ppfsetdevice("JET")
ppf.ppfuid("JETPPF","r")


get=lambda dtyp:ppf.ppfget(pulse,dda,dtyp)

initial=59
initial_time = 40.

traceNums=np.arange(96)+1
#get all the traces
traces=[]
times=[]
Rs=[]
Rt=[]
maxT=0
maxTime=0
data = []
for i in traceNums:
    [_,_,trace_data,_,trace_time,ierr]=get("TE"+str(i).zfill(2))
    [_,_,R_data,_,R_time,ierr]=get("RC"+str(i).zfill(2))
    if(len(R_data>0)):
        maxT=max(maxT,trace_data.max())
        maxTime=max(maxTime,trace_time.max())
    traces.append(trace_data)
    times.append(trace_time)
    Rs.append(R_data)
    Rt.append(R_time)
    if len(R_data)==0:
        data.append([i,0])
    else:
        data.append([i,np.amax(R_data)])
data = np.array(data)
data = data[data[:,1].argsort()]
ind=data[:,0].astype(int)-1


fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)
ax.plot(times[initial],traces[initial])
ax.set_ylim([0,maxT])
ax.set_xlim([initial_time,maxTime])
ax.set_ylabel('Te (eV)')
ax.set_xlabel('t (s)')
ax.set_title(f'Channel: {initial}')
ax2=ax.twinx()
ax2.plot(Rt[initial],Rs[initial],'r.')
ax2.set_xlim([initial_time,maxTime])
ax2.set_ylim([2.0,4.0])
ax2.set_ylabel('$R_{mid}$ (m)',color='r')
sliderAxis=plt.axes([0.25, 0.1, 0.65, 0.03])
slider=wid.Slider(sliderAxis,'Index',1,len(traces),valinit=initial,valfmt='%0.0f')

def update(val):
    channel=ind[int(slider.val)]
    ylim=ax.get_ylim()
    xlim=ax.get_xlim()
    ax.clear()
    ax.plot(times[channel],traces[channel])
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)
    ax.set_title(f'Channel {channel+1}')
    ax2.clear()
    ax2.plot(Rt[channel],Rs[channel],'r.')
    ax2.set_xlim(xlim)
    ax.set_ylabel('Te (eV)')
    ax2.set_ylabel('$R_{mid}$ (m)',color='r')



slider.on_changed(update)

plt.show()

