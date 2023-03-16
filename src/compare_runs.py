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

    p1.add_data(signal)
    p2.add_data(signal)
    
    fig = plt.figure() 
    fig.suptitle(f'{p1.pulse} {p1.runid} vs. {p2.runid} at t= {t:.2f}', fontsize=13) 
    ax = fig.add_subplot(111)

    it1 = gi(p1.t,t)
    it2 = gi(p2.t,t)

    x1 =  p1.x[it1,:]
    x2 =  p2.x[it2,:]

    ax.plot(
            x1,
            p1.transp(signal)[it1,:],
            color='r',
            linewidth=2,
            label=f'{p1.runid}: {signal}'
            )
    ax.plot(
            x2,
            p2.transp(signal)[it2,:],
            color='b',
            linestyle='dashed',
            linewidth=2,
            label=f'{p2.runid}: {signal}'
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
    p1.add_data(signal)
    p2.add_data(signal)

    x1 = p1.t    
    x2 = p2.t
    ix1 = gi(p1.x[0],x_pos)
    ix2 = gi(p2.x[0],x_pos)

    fig = plt.figure() 
    fig.suptitle(f'{p1.pulse} {p1.runid} vs. {p2.runid} at x= {x_pos:.2f}', fontsize=13) 
    ax = fig.add_subplot(111)
    ax.plot(
            x1,
            p1.transp(signal)[:,ix1],
            color='r',
            linewidth=2,
            label=f'{p1.runid}: {signal}'
            )
    ax.plot(
            x2,
            p2.transp(signal)[:,ix2],
            color='b',
            linestyle='dashed',
            linewidth=2,
            label=f'{p2.runid}: {signal}'
            )
    
    ax.set_xlabel('X')
    ax.set_ylabel(r'$units$')
    cornernote(ax)
    #ax.set_xlim(0.0,1.0)
    leg=ax.legend()
    leg.set_draggable(True)
    win.addPlot(f'{signal} time trace',fig)

def cornernote(axis):                                                                           
    fontsize=9
    text = p1.pulse
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
runid=str(sys.argv[2])
runid2=str(sys.argv[3])

p1 = ps.Neutrons(pulse,runid)
#p1.get_transp_neutrons()
#p1.get_exp()

print(f'Times from {p1.t[0]:.3f}s to {p1.t[-1]:.3f}s.')

rnt_time, rnt = p1.rnt
rnt_unc = 0.1
win=pw.plotWindow()

#p2.Densities(pulse,rnid2)
p2= ps.Neutrons(pulse,runid2)
#p2.get_transp_neutrons()

fig = plt.figure() 
fig.suptitle(f'{pulse} {runid} vs. {runid2}', fontsize=13) 
#ax = fig.add_subplot(111)
#
#ax.plot(rnt_time,rnt,color='r',linewidth=2,label='TIN/RNT')
#ax.fill_between(rnt_time, 
#        (1 - rnt_unc) * rnt, 
#        (1 + rnt_unc) * rnt, 
#        alpha=0.4, 
#        facecolor='r', 
#        edgecolor='none'
#        )
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



t1=2.5
xpos = 0.2
signal = 'Q'

plot_profile(signal,t1)
plot_time_trace(signal,xpos)


win.show()



