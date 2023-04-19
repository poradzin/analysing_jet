#!/usr/bin/env python
import plotWindow as pw
import matplotlib.pyplot as plt
import profiles as ps
import sys

pulse=int(sys.argv[1])
runid=str(sys.argv[2])


neutrons= ps.Neutrons(pulse,runid)
neutrons.get_transp_neutrons()
neutrons.get_exp()

print(f'Times from {neutrons.t[0]:.3f}s to {neutrons.t[-1]:.3f}s.')

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
rnt_time, rnt = neutrons.rnt
rnt_unc = 0.1
win=pw.plotWindow()

fig = plt.figure() 
fig.suptitle(f'Thermal vs. beam', fontsize=13) 
ax = fig.add_subplot(111)

ax.plot(rnt_time,rnt,color='r',linewidth=2,label='TIN/RNT')
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

fig = plt.figure() 
#fig.suptitle(f'neutron_rate - DT split up', fontsize=13) 
ax = fig.add_subplot(111)

ax.set_title(f'{neutrons.transpcid} Thermal vs. beam-target neutrons')

ax.plot(rnt_time,rnt,color='r',linewidth=2,label='TIN/RNT')
#ax.plot(R14_time,R14,color='orange',linewidth=2,label='TIN/R14')
ax.plot(neutrons.transp_time+40,
        neutrons.transp('NEUTT'),
        color='k',
        linewidth=2,
        label='Total'
        )
ax.plot(neutrons.transp_time+40, 
        neutrons.transp('NEUTX_DD'),
        color='green',
        linewidth=2, 
        label='Thermal DD'
        )
ax.plot(neutrons.transp_time+40, 
        neutrons.transp('BTNTS_DD'),
        color='blue',
        linewidth=2, 
        label='Beam-target DD'
        )
ax.plot(neutrons.transp_time+40,
        neutrons.transp_dt,
        color='k',
        linestyle= (0, (5, 1)),
        linewidth=2,
        label='Total DT'
        )
ax.plot(neutrons.transp_time+40, 
        neutrons.transp('NEUTX_DT'),
        color='green',
        linestyle= (0, (5, 1)),
        linewidth=2, 
        label='Thermal DT'
        )
ax.plot(neutrons.transp_time+40, 
        neutrons.transp('BTNTS_DT'),
        color='blue',
        linestyle = (0, (5, 1)),
        linewidth=2, 
        label='Beam-target DT'
        )
if neutrons.signal('BTNTS_TD'):
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



win.show()



