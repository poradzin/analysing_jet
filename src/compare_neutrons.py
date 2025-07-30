import plotWindow as pw
import matplotlib.pyplot as plt
import profiles as ps
import sys
import numpy as np
import argparse

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

# Function to parse command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description='Compares neutron rates of several TRANSP runs')
    parser.add_argument('pulse', type=int, help='JET discharge number')
    parser.add_argument('--label', nargs='+', help='List of labels for runs')
    parser.add_argument('--color',nargs='+', help='List of colors for runs')
    parser.add_argument('--linestyle',nargs='+', help='List of linestyles for runs')
    parser.add_argument('--total',action='store_true', 
            help='Plots only total neutron rate.')
    parser.add_argument('--thermal',action='store_true', 
            help='Adds thermal neutron rate component.')
    parser.add_argument('--beam',action='store_true', 
            help='Adds beam target neutron rate component.')
    parser.add_argument('--time_trace',action='store_true', 
            help='Creates time-trace plot comparing different components.')
    parser.add_argument('runs', nargs='+', help='List of runs')
    return parser.parse_args()

args = parse_arguments()

pulse = args.pulse
runs = args.runs

# Check if the number of labels and colors provided matches the number of runs
if args.label and len(args.label) != len(runs):
    print("Error: Number of labels provided does not match the number of TRASP runs.")
    sys.exit(1)
if args.color and len(args.color) != len(runs):
    print("Error: Number of colors provided does not match the number of TRANSP runs.")
    sys.exit(1)
if args.linestyle and len(args.linestyle) != len(runs):
    print("Error: Number of linestyles provided does not match the number of TRANSP runs.")
    sys.exit(1)

# Convert linestyle strings or tuples to appropriate format
linestyles = []
if args.linestyle:
    for style in args.linestyle:
        # Check if the style is a tuple
        if style.startswith('(') and style.endswith(')'):
            try:
                linestyle_tuple = eval(style)  # Safely evaluate the tuple
                linestyles.append(linestyle_tuple)
            except SyntaxError:
                print(f"Error: Invalid linestyle tuple '{style}'.")
                sys.exit(1)
        else:
            linestyles.append(style)  # Otherwise, assume it's a string


p = {}
for run in runs:
    p[run] = ps.Neutrons(pulse, str(run))
    p[run].get_transp_neutrons()

p[runs[0]].get_exp()

rnt_time, rnt = p[runs[0]].rnt
rnt_unc = 0.1

win = pw.plotWindow()

if args.total:
    linestyles = ('solid', (0,(5,1)), 'dotted', 'dashdot')
    #color = iter(plt.cm.rainbow(np.linspace(0,1,len(p))))
else:
    linestyles = ('solid', 'dotted', (0,(5,1)), 'dotted', 'dashdot', 'loosely dotted')
lnstyle = iter(linestyles)
    
fig = plt.figure()
fig.suptitle(f'{pulse} neutron rate', fontsize=13)
ax = fig.add_subplot(111)

ax.plot(rnt_time, rnt, color='r', linewidth=2, label='TIN/RNT')
ax.fill_between(rnt_time,
                (1 - rnt_unc) * rnt,
                (1 + rnt_unc) * rnt,
                alpha=0.4,
                facecolor='r',
                edgecolor='none')
if args.total:
    ln = next(lnstyle)
    for ind, run in enumerate(runs):
        lbt = 'total\n' if ind==0 else ''
        label = lbt + args.label[runs.index(run)] if args.label else f'{run} NEUTT'
        ax.plot(p[run].transp_time + 40,
                p[run].transp('NEUTT'),
                color=args.color[runs.index(run)] if args.color else 'k',
                linewidth=2,
                linestyle=args.linestyle[runs.index(run)] if args.linestyle else ln,
                label=label
               )


if args.beam: 
    ln = next(lnstyle)
    for ind, run in enumerate(runs):
        bt = 'beam-target' if ind==0 else ''
        label = bt+ args.label[runs.index(run)] if args.label and not args.total else bt 
        ax.plot(p[run].transp_time + 40,
                p[run].transp('BTNTS'),
                color=args.color[runs.index(run)] if args.color else 'b',
                linewidth=2,
                linestyle=args.linestyle[runs.index(run)] if args.linestyle else ln,
                #linestyle = 'dashdot',
                label=label
               )
if args.thermal:
    ln = next(lnstyle) 
    for ind, run in enumerate(runs):
        lth = 'thermal\n'if ind==0 and not (args.total or args.beam) else 'thermal' if ind==0 else ''
        label = lth+ args.label[runs.index(run)] if args.label and not (args.total or args.beam) else lth
        ax.plot(p[run].transp_time + 40,
                p[run].transp('NEUTX'),
                color=args.color[runs.index(run)] if args.color else 'orange',
                linewidth=2,
                linestyle=args.linestyle[runs.index(run)] if args.linestyle else ln,
                #linestyle = 'dashed',
                label=label
               )
 
if not (args.total or args.thermal or args.beam):
    lnstyle = iter(linestyles)
    for ind, run in enumerate(runs):
        ln = next(lnstyle)
        lbt = 'beam-target'if ind==0 else ''
        ax.plot(p[run].transp_time + 40,
                p[run].transp('BTNTS'),
                color=args.color[runs.index(run)] if args.color else 'b',
                linewidth=2,
                linestyle=args.linestyle[runs.index(run)] if args.linestyle else ln,
                label=lbt
               )
    
        lth = 'thermal'if ind==0 else ''
        ax.plot(p[run].transp_time + 40,
                p[run].transp('NEUTX'),
                color=args.color[runs.index(run)] if args.color else 'orange',
                linewidth=2,
                linestyle=args.linestyle[runs.index(run)] if args.linestyle else ln,
                label=lth
               )
    
        lth = 'beam-beam'if ind==0 else ''
        ax.plot(p[run].transp_time + 40,
                p[run].transp('BBNTS'),
                color=args.color[runs.index(run)] if args.color else 'green',
                linewidth=2,
                linestyle=args.linestyle[runs.index(run)] if args.linestyle else ln,
                label=lth
               )



ax.set_xlabel('time [s]')
ax.set_ylabel(r'$s^{-1}$')
# cornernote(ax,pulse) # This function seems to be custom, you can uncomment and use if needed
ax.set_xlim(p[runs[0]].transp_time[0] + 40, p[runs[0]].transp_time[-1] + 40)
leg = ax.legend(loc='best')
leg.set_draggable(True, use_blit=False, update='loc')

win.addPlot('RNT', fig)


if args.total:
    linestyles = ('solid', (0,(5,1)), 'solid', 'dotted', 'dashdot')
    #color = iter(plt.cm.rainbow(np.linspace(0,1,len(p))))
else:
    linestyles = ('solid', 'dotted', (0,(5,1)), 'dotted', 'dashdot', 'loosely dotted')
lnstyle = iter(linestyles)
    
fig = plt.figure()
fig.suptitle(f'{pulse} neutron ratios', fontsize=13)
ax = fig.add_subplot(111)

lnstyle = iter(linestyles)
ln = next(lnstyle)
for ind, run in enumerate(runs):
    bt = 'beam-target\n' if ind==0 else ''
    label = bt+ args.label[runs.index(run)] if args.label else bt 

    lbt = 'beam-target'if ind==0 else ''
    ax.plot(p[run].transp_time + 40,
            p[run].transp('BTNTS')/p[run].transp('NEUTT'),
            color=args.color[runs.index(run)] if args.color else 'b',
            linewidth=2,
            linestyle=args.linestyle[runs.index(run)] if args.linestyle else ln,
            label=label
           )

ln = next(lnstyle)
for ind, run in enumerate(runs):
    th = 'thermal/tot' if ind==0 else ''
    label = th
    ax.plot(p[run].transp_time + 40,
            p[run].transp('NEUTX')/p[run].transp('NEUTT'),
            color=args.color[runs.index(run)] if args.color else 'orange',
            linewidth=2,
            linestyle=args.linestyle[runs.index(run)] if args.linestyle else ln,
            label=label
           )

ax.set_xlabel('time [s]')
#ax.set_ylabel(r'$s^{-1}$')
# cornernote(ax,pulse) # This function seems to be custom, you can uncomment and use if needed
ax.set_xlim(p[runs[0]].transp_time[0] + 40, p[runs[0]].transp_time[-1] + 40)
leg = ax.legend(loc='best')
leg.set_draggable(True, use_blit=False, update='loc')

win.addPlot('Ratios', fig)

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





