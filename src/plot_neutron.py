#!/usr/bin/env python
import plotWindow as pw
import matplotlib.pyplot as plt
import numpy as np
import profiles as ps
import sys

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('pulse', type=int)
parser.add_argument('runid', type=str)
parser.add_argument('--add_strap', nargs=2, type=float, metavar=('XMIN', 'XMAX'),
                    help='Add a shaded vertical strap between XMIN and XMAX (in seconds)')
parser.add_argument('--strap_color', default='lightgrey', help='Color of the strap (default: lightgrey)')
parser.add_argument('--strap_alpha', type=float, default=0.5, help='Transparency of the strap (default: 0.5)')

args = parser.parse_args()

pulse = args.pulse
runid = args.runid

#pulse=int(sys.argv[1])
#runid=str(sys.argv[2])


neutrons= ps.Neutrons(pulse,runid)
neutrons.get_transp_neutrons()
neutrons.get_exp()

print(f'Times from {neutrons.t[0]:.3f}s to {neutrons.t[-1]:.3f}s.')
print(f"Total TRANSP neutrons: {neutrons.total('NEUTT')}")
print(f"Total TIN/RNT neutrons: {neutrons.tot_exp()}")  
def cornernote(axis):                                                                           
    fontsize=9
    text = neutrons.transpcid
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

import numpy as np

def safe_divide(x, y):
    """Returns (x, x/y) where y > 0, and 0 elsewhere"""
    result = np.zeros_like(x)
    mask = y > 0
    result[mask] = x[mask] / y[mask]
    return result

rnt_time, rnt = neutrons.rnt
rnt_unc = 0.1


win=pw.plotWindow()
###################################################################
fig = plt.figure() 
fig.suptitle(f'Thermal vs. beam', fontsize=13) 
ax = fig.add_subplot(111)

# Optionally add strap
if args.add_strap:
    xmin, xmax = args.add_strap
    ax.axvspan(xmin, xmax, color=args.strap_color, alpha=args.strap_alpha)


ax.plot(rnt_time,rnt,color='r',linewidth=2,label='Measured')
ax.fill_between(rnt_time, 
        (1 - rnt_unc) * rnt, 
        (1 + rnt_unc) * rnt, 
        alpha=0.4, 
        facecolor='r', 
        edgecolor='none'
        )
ax.plot(neutrons.transp_time+40,
        neutrons.transp('NEUTT'),
        color='k',
        linewidth=2,
        label='total'
        )
ax.plot(neutrons.transp_time+40,
        neutrons.transp('NEUTX'),
        color = 'orange',
        linewidth=2,
        label='thermal')
ax.plot(neutrons.transp_time+40,
        neutrons.transp('BTNTS'),
        color = 'blue',
        linewidth=2,
        label='beam-target'
        )
ax.plot(neutrons.transp_time+40,
        neutrons.transp('BBNTS'),
        color = 'green', 
        linewidth=2,
        label='beam-beam')

#ax.plot(neutrons.transp_time+40, neutrons_th,color = 'gold',linewidth=2, label='RNT - beam')
#ax.plot(time_av,data_av,color='purple', linewidth=3)

ax.set_xlabel('Time [s]')
ax.set_ylabel('Neutron rate [#/s]')
ax.set_xlim(neutrons.transp_time[0]+40,neutrons.transp_time[-1]+40)
cornernote(ax)
leg = ax.legend(fontsize='x-small', title='Neutron Rate')
leg.set(draggable=True)

#win.cornernote(neutrons.transpcid)
win.addPlot('th vs. beam',fig)

##################################################################################################
#################################################################################################
fig = plt.figure() 
#fig.suptitle(f'neutron_rate - DT split up', fontsize=13) 
ax = fig.add_subplot(111)

# Optionally add strap
if args.add_strap:
    xmin, xmax = args.add_strap
    ax.axvspan(xmin, xmax, color=args.strap_color, alpha=args.strap_alpha)


ax.set_title(f'{neutrons.transpcid} Thermal vs. beam-target neutrons')

ax.plot(rnt_time,rnt,color='r',linewidth=2,label='TIN/RNT')
#ax.plot(R14_time,R14,color='orange',linewidth=2,label='TIN/R14')

ax.plot(neutrons.transp_time+40,
        neutrons.transp('NEUTT'),
        color='k',
        linewidth=2,
        label='Total'
        )

if neutrons.transp('NEUTX_DD') is not None: 
    ax.plot(neutrons.transp_time+40, 
            neutrons.transp('NEUTX_DD'),
            color='green',
            linewidth=2, 
            label='Thermal DD'
            )
if neutrons.transp('BTNTS_DD') is not None: 
    ax.plot(neutrons.transp_time+40, 
            neutrons.transp('BTNTS_DD'),
            color='blue',
            linewidth=2, 
            label='Beam-target DD'
            )
#print(f"neutrons.transp_dt: {neutrons.transp_dt}")
#print(f"neutrons.transp('EMUTX_DT'): {neutrons.transp('NEUTX_DT')}")
#print(f"isinstance(neutrons.transp_dt, (list,np.ndarray)) : {isinstance(neutrons.transp_dt, (np.ndarray,list))}")
#print(f"type(neutrons.transp_dt)): {type(neutrons.transp_dt)}") 
if isinstance(neutrons.transp_dt, (list,np.ndarray)):
    ax.plot(neutrons.transp_time+40,
        neutrons.transp_dt,
	color='k',
	linestyle= (0, (5, 1)),
	linewidth=2,
	label='Total DT'
	)
if isinstance(neutrons.transp('NEUTX_DT'), (list,np.ndarray)):
    ax.plot(neutrons.transp_time+40, 
        neutrons.transp('NEUTX_DT'),
        color='green',
        linestyle= (0, (5, 1)),
        linewidth=2, 
        label='Thermal DT'
        )
#print(f"isinstance(neutrons.transp('BTNTS_DT'), list): {isinstance(neutrons.transp('BTNTS_DT'),list)}")
#print(f"type(neutrons.transp('BTNTS_DT')): {type(neutrons.transp('BTNTS_DT'))}")
if isinstance(neutrons.transp('BTNTS_DT'), (list,np.ndarray)):
    print('BTNTS_DT went through')
    ax.plot(neutrons.transp_time+40, 
        neutrons.transp('BTNTS_DT'),
        color='blue',
        linestyle = (0, (5, 1)),
        linewidth=2, 
        label='Beam-target DT'
        )
print(f"neutrons.signal('BTNTS_TD'):{neutrons.signal('BTNTS_TD')}")
if isinstance(neutrons.signal('BTNTS_TD'),(list,np.ndarray)):
    ax.plot(neutrons.transp_time+40, 
        neutrons.transp('BTNTS_TD'),
        color='rebeccapurple',
        linestyle= (0, (5, 1)),
        linewidth=2,
        label='Beam-target TD'
        )
#################################################
ax.set_xlabel('time [s]')
ax.set_ylabel(r'$s^{-1}$')
cornernote(ax)
ax.set_xlim(neutrons.transp_time[0]+40,neutrons.transp_time[-1]+40)
leg=ax.legend()
leg.set_draggable(True)
win.addPlot('thermal vs. DT',fig)

################################################################################################
#################################################################################################
fig = plt.figure() 
#fig.suptitle(f'neutron_rate - DT split up', fontsize=13) 
ax = fig.add_subplot(111)

# Optionally add strap
if args.add_strap:
    xmin, xmax = args.add_strap
    ax.axvspan(xmin, xmax, color=args.strap_color, alpha=args.strap_alpha)


ax.set_title(f'{neutrons.transpcid} Thermal, beam-target ratio to total. ')

ax.plot(neutrons.transp_time+40,safe_divide(neutrons.transp('NEUTX'),neutrons.transp('NEUTT')),color='orange',linewidth=2, label="TH/TOT")
ax.plot(neutrons.transp_time+40,safe_divide(neutrons.transp('BTNTS'),neutrons.transp('NEUTT')),color='blue',linewidth=2, label="BT/TOT")

ax.set_xlabel('time [s]')
#ax.set_ylabel(r'$s^{-1}$')
cornernote(ax)
ax.set_xlim(neutrons.transp_time[0]+40,neutrons.transp_time[-1]+40)
leg=ax.legend()
leg.set_draggable(True)
win.addPlot('thermal/DT ratio',fig)

#

win.show()



