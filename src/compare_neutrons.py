import plotWindow as pw
import matplotlib.pyplot as plt
import profiles as ps
import sys
import numpy as np

def getind(val,data):
    return np.argmin(np.abs(data-val))

def compare_outputs():
    '''
    Relevant for regression tests.
    Compares all p1 output with p2 and prints the ones that do not agree within 
    some relative error.
    '''
    rel = 0.05
    out={}
    return None


def plot_profile(signal,t):
    '''
    '''
    
    fig = plt.figure() 
    fig.suptitle(f'{p1.pulse} {p1.runid} vs. {p2.runid} at t= {t:.2f}', fontsize=13) 
    ax = fig.add_subplot(111)

    color = iter(plt.cm.rainbow(np.linspace(0,1,len(p))))
    for run in p:
        c=next(color)
        ax.plot(
            p[run].x[gi(p[run].t,t),:],
            p[run].transp(signal)[gi(p[run].t,t),:],
            color=c,
            linewidth=2,
            label=f'{p[run].runid}: {signal}'
            )
        
    ax.set_xlabel('X')
    ax.set_ylabel(r'$units$')
    cornernote(ax)
    ax.set_xlim(0.0,1.0)
    leg=ax.legend()
    leg.set_draggable(True)
    win.addPlot(f'{signal} profiles',fig)

def plot_time_trace(signal,x_pos):
    '''
    '''

    color = iter(plt.cm.rainbow(np.linspace(0,1,len(p))))
    fig = plt.figure() 
    fig.suptitle(f'{p1.pulse} {p1.runid} vs. {p2.runid} at x= {x_pos:.2f}', fontsize=13) 
    ax = fig.add_subplot(111)
    for run in p:
        c=next(color)
        ax.plot(
            p[run].t,
            p[run].transp(signal)[:,gi(p[run].x[0],x_pos)],
            color=c,
            linewidth=2,
            label=f'{p[run].runid}: {signal}',
            )
            
    ax.set_xlabel('X')
    ax.set_ylabel(r'$units$')
    cornernote(ax)
    #ax.set_xlim(0.0,1.0)
    leg=ax.legend()
    leg.set_draggable(True)
    win.addPlot(f'{signal} time trace',fig)

def cornernote(axis,text):                                                                           
    fontsize=9
    axis.annotate(
            text, 
            xy=(0.98, 0.02), 
            xycoords='figure fraction', 
            fontsize=fontsize, 
            ha='right', 
            va='bottom', 
            label='cornernote'
            )
    return None

############################################################

gi = getind
pulse=int(sys.argv[1])
p={}
runs = sys.argv[2:]

for run in runs:
    p[run]=ps.Neutrons(pulse,str(run))
    p[run].get_transp_neutrons()

p[runs[0]].get_exp()


rnt_time, rnt = p[runs[0]].rnt
rnt_unc = 0.1

win=pw.plotWindow()

linestyles = ('solid','dashed','dotted','dashdot','loosely dotted')
lnstyle = iter(linestyles)
fig = plt.figure() 
fig.suptitle(f'{pulse}', fontsize=13) 
ax = fig.add_subplot(111)

ax.plot(rnt_time,rnt,color='r',linewidth=2,label='TIN/RNT')
ax.fill_between(rnt_time, 
        (1 - rnt_unc) * rnt, 
        (1 + rnt_unc) * rnt, 
        alpha=0.4, 
        facecolor='r', 
        edgecolor='none'
        )

for run in runs:
    ln = next(lnstyle)
    ax.plot(p[run].transp_time+40,
            p[run].transp('NEUTT'),
            color='k',
            linewidth=2,
            linestyle=ln,
            label=f'{run} NEUTT'
            )




ax.set_xlabel('time [s]')
ax.set_ylabel(r'$s^{-1}$')
cornernote(ax,pulse)
ax.set_xlim(p[runs[0]].transp_time[0]+40,p[runs[0]].transp_time[-1]+40)
leg=ax.legend()
leg.set_draggable(True)

win.addPlot('RNT',fig)

win.show()

#ax.plot(p1.transp_time+40,
#        p1.transp('NEUTT'),
#        color='k',
#        linewidth=2,
#        label=f'{runid} total'
#        )
#ax.plot(p2.transp_time+40,
#        p2.transp('NEUTT'),
#        color='k',
#        linestyle = 'dashed',
#        linewidth=2,
#        label=f'{runid2} total'
#        )
#
#ax.plot(p1.transp_time+40, 
#        p1.transp('BTNTS_DT'),
#        color='magenta',
#        #linestyle = (0, (5, 1)),
#        linewidth=2, 
#        label=f'{runid}: BTNTS_DT'
#        )
#
#ax.plot(p2.transp_time+40, 
#        p2.transp('BTNTS_DT'),
#        color='magenta',
#        #linestyle = (0, (5, 1)),
#        linewidth=2,
#        linestyle = 'dashed',
#        label=f'{runid2}: BTNTS_DT'
#
#        )
#
#if p1.signal('BTNTS_TD'):
#    ax.plot(p1.transp_time+40, 
#        p1.transp('BTNTS_TD'),
#        color='darkred',
#        #linestyle= (0, (5, 1)),
#        linewidth=2,
#        label=f'{runid}: BTNTS_TD'
#        )
#
##V04
#if p1.signal('BTNTS_TD'):
#    ax.plot(p2.transp_time+40, 
#        p2.transp('BTNTS_TD'),
#        color='purple',
#        #linestyle= (0, (5, 1)),
#        linewidth=2,
#        linestyle = 'dashed',
#        label = f'{runid2}: BTNTS_TD',
#        )
#
#ax.set_xlabel('time [s]')
#ax.set_ylabel(r'$s^{-1}$')
#cornernote(ax)
#ax.set_xlim(p1.transp_time[0]+40,p1.transp_time[-1]+40)
#leg=ax.legend()
#leg.set_draggable(True)
#win.addPlot(f'{runid} vs. {runid2}',fig)



#signal= 'NEUT'
#p={}
#for run in sys.argv[2:]:
#    p[run]=ps.Neutrons(pulse,str(run))
#    p[run].add_data(signal)


#t1=1.9
#plot_profile(signal,t1)
#t2 = 4.0
#plot_profile(signal,t2)


#xpos=0.45
#plot_time_trace(signal,xpos)
#xpos = 0.2
#plot_time_trace(signal,xpos)





