import plotWindow as pw
import matplotlib.pyplot as plt
import profiles as ps
import sys
import numpy as np

def gi(val,data):
    '''
    find index of data closesest to the required value
    '''
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

    #color = iter(plt.cm.rainbow(np.linspace(0,1,len(p))))
    name = 'tab10'
    cmap = plt.cm.get_cmap(name)
    color = iter(cmap.colors)
    units=[]
    for run in p:
        units.append(p[run].units(signal))
        c=next(color)
        ax.plot(
            p[run].x[gi(p[run].t,t),:],
            p[run].transp(signal)[gi(p[run].t,t),:],
            color=c,
            linewidth=1,
            label=f'{p[run].runid}'
            )
    if all(unit==units[0] for unit in units):
        ax.set_ylabel(units[0])
    else:
        ax.set_ylabel(r'$units$')
    ax.set_xlabel('X')
    cornernote(ax)
    ax.set_xlim(0.0,1.0)
    leg=ax.legend()
    leg.set_draggable(True)
    win.addPlot(f'{signal} profiles',fig)

def plot_time_trace(signal,x_pos):
    '''
    '''
    #color = iter(plt.cm.rainbow(np.linspace(0.1,0.9,len(p))))
    name = 'tab10'
    cmap = plt.cm.get_cmap(name)
    color = iter(cmap.colors)

    fig = plt.figure() 
    fig.suptitle(f'{signal} at x= {x_pos:.2f}', fontsize=13) 
    ax = fig.add_subplot(111)
    for run in p:
        c=next(color)
        ax.plot(
            p[run].t,
            p[run].transp(signal)[:,gi(p[run].x[0],x_pos)],
            color=c,
            linewidth=1,
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
    pulse=int(sys.argv[1])

    signal='BDENS'
    signal2='BDENS_D'
    times=[6.0,7.0]
    xpos=[0.1]

    p={}
    for run in sys.argv[2:]:
        p[run]=ps.Transp(pulse,str(run))
        p[run].add_data(signal)
        p[run].add_data(signal2)

    win=pw.plotWindow()
    for t in times: 
        plot_profile(signal,t)
        plot_profile(signal2,t)
    for x in xpos:
        plot_time_trace(signal,x)
    win.show()



