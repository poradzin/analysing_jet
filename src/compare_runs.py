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
    fig.suptitle(f'{signal} at t= {t:.2f}', fontsize=13) 
    ax = fig.add_subplot(111)

    color = iter(plt.cm.rainbow(np.linspace(0,1,len(p))))
    for run in p:
        c=next(color)
        ax.plot(
            p[run].x[gi(p[run].t,t),:],
            p[run].transp(signal)[gi(p[run].t,t),:],
            color=c,
            linewidth=2,
            label=f'{p[run].runid}'
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
    fig.suptitle(f'{signal} at x= {x_pos:.2f}', fontsize=13) 
    ax = fig.add_subplot(111)
    for run in p:
        c=next(color)
        ax.plot(
            p[run].t,
            p[run].transp(signal)[:,gi(p[run].x[0],x_pos)],
            color=c,
            linewidth=2,
            label=f'{p[run].runid}',
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
    text = pulse
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
if __name__=="__main__":
    gi = getind
    pulse=int(sys.argv[1])
    signal='Q'

    p={}
    for run in sys.argv[2:]:
        p[run]=ps.Neutrons(pulse,str(run))
        p[run].add_data(signal)
    win=pw.plotWindow()
    t1=1.9
    plot_profile(signal,t1)
    t2 = 4.0
    plot_profile(signal,t2)
    
    xpos=0.45
    plot_time_trace(signal,xpos)
    xpos = 0.2
    plot_time_trace(signal,xpos)
    
    
    win.show()



