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

def convert_units(units):
    conversion = {'CM**2/SEC':('$m^2/s$', 1e-4)}
    if units in conversion:
        return conversion[units]
    else:
        return (units,1.0)

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
    if all(unit==units[0] for unit in units):
        unit,factor = convert_units(units[0])
    else:
        unit = r'$units$'
        factor = 1.

    for run in p:
        c=next(color)
        ax.plot(
            p[run].x[gi(p[run].t,t),:],
            p[run].transp(signal)[gi(p[run].t,t),:]*factor,
            color=c,
            linewidth=1,
            label=f'{p[run].runid} at {p[run].t[gi(p[run].t,t)]:.3f}s'
           )
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-3,3))
    ax.set_xlabel('X')
    ax.set_ylabel(unit)
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
    
    units = []
    for run in p:
        units.append(p[run].units(signal))
    if all(unit==units[0] for unit in units):
        unit,factor = convert_units(units[0])
    else:
        unit = r'$units$'
        factor = 1.

    for run in p:
        c=next(color)
        ax.plot(
            p[run].t,
            p[run].transp(signal)[:,gi(p[run].x[0],x_pos)]*factor,
            color=c,
            linewidth=1,
            label=f'{p[run].runid} at x={p[run].x[0][gi(p[run].x[0],x_pos)]:.4f}',
            )
    ax.set_xlabel('X')
    ax.set_ylabel(unit)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-3,3))
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
    #signals1 = ['ETA_SP', 'ETA_USE', 'ETA_NC', 'ETA_TSC']
    signals1 = ['NE']
    signals2 = ['TE','TI']
    signals = [signals1,signals2]
    times=[9.00]
    xpos=[0.5]

    p={}
    for run in sys.argv[2:]:
        p[run]=ps.Transp(pulse,str(run))
        for signal in signals:
            p[run].add_data(*signal)

    win=pw.plotWindow()
    for t in times: 
        for sig in signals:
            plot_profile(sig,t)
        #plot_profile(signal2,t)
    for x in xpos:
        for sig in signals:
            plot_time_trace(sig,x)
    win.show()





