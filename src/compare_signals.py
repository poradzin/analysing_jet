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

def plot_profile(signals,times):
    '''
    '''
    fig = plt.figure() 
    fig.suptitle(f'{signals[0]} profiles', fontsize=13) 
    ax = fig.add_subplot(111)

    #color = iter(plt.cm.rainbow(np.linspace(0,1,len(p))))
    name = 'tab10'
    cmap = plt.cm.get_cmap(name)
    color = iter(cmap.colors)
    units=[]
    for run in p:
        for signal in signals:
            units.append(p[run].units(signal))
    if all(unit==units[0] for unit in units):
        unit,factor = convert_units(units[0])
    else:
        unit = r'$units$'
        factor = 1.

    for run in p:
        for signal in signals:
            for t in times:
                c=next(color)
                ax.plot(
                    p[run].x[gi(p[run].t,t),:],
                    p[run].transp(signal)[gi(p[run].t,t),:]*factor,
                    color=c,
                    linewidth=2,
                    label=f'{p[run].runid}/{signal} at {p[run].t[gi(p[run].t,t)]:.3f}s'
                    )
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-3,3))
    ax.tick_params(axis='both', labelsize=14)
    ax.set_xlabel('X')
    ax.set_ylabel(unit)
    cornernote(ax)
    ax.set_xlim(0.0,1.0)
    leg=ax.legend()
    leg.set_draggable(True)
    win.addPlot(f'{signals[0]} profiles',fig)

def plot_time_trace(signals,x_pos):
    '''
    '''
    print(f'Signals: {signals}')
    #color = iter(plt.cm.rainbow(np.linspace(0.1,0.9,len(p))))
    name = 'tab10'
    cmap = plt.cm.get_cmap(name)
    color = iter(cmap.colors)

    fig = plt.figure() 
    fig.suptitle(f'{signals[0]} at x= {x_pos:.2f}', fontsize=13) 
    ax = fig.add_subplot(111)
    
    units = []
    for run in p:
        for signal in signals:
            print(f'Signal: {signal}')
            units.append(p[run].units(signal))
    if all(unit==units[0] for unit in units):
        unit,factor = convert_units(units[0])
    else:
        unit = r'$units$'
        factor = 1.

    for run in p:
        for signal in signals:
            c=next(color)
            ax.plot(
                p[run].t,
                p[run].transp(signal)[:,gi(p[run].x[0],x_pos)]*factor,
                color=c,
                linewidth=2,
                label=f'{p[run].runid}/{signal} at x={p[run].x[0][gi(p[run].x[0],x_pos)]:.4f}',
                )
    ax.set_xlabel('X')
    ax.set_ylabel(unit)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-3,3))
    ax.tick_params(axis='both', labelsize=14)
    cornernote(ax)
    #ax.set_xlim(0.0,1.0)
    leg=ax.legend()
    leg.set_draggable(True)
    win.addPlot(f'{signals[0]} time trace',fig)

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
    #signals1 = ['ETA_SP', 'ETA_SPS', 'ETA_WNC', 'ETA_TSC','ETA_USE']
    #signals2 = ['NUSTE','NUSTI']
    signals1 = ['BDENS']
    signals2 = ['UFASTPA']
    #signals = [signals1, signals2]
    signals = [signals1,signals2]
    times1=[7.1,7.15,7.2,7.25,7.285,7.295]
    times = [times1]
    xpos=[0.99]

    p={}
    for run in sys.argv[2:]:
        p[run]=ps.Transp(pulse,str(run))
        for signal in list(np.concatenate(signals)):
            p[run].add_data(signal)

    win=pw.plotWindow()
    for t in times: 
        for sig in signals:
            plot_profile(sig,t)
        #plot_profile(signal2,t)
    for x in xpos:
        for sig in signals:
            plot_time_trace(sig,x)
    win.show()





