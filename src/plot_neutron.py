import plotWindow as pw
import matplotlib.pyplot as plt
import profiles as ps
import sys

pulse=int(sys.argv[1])
runid=str(sys.argv[2])


neutrons= ps.Neutrons(pulse,runid)
neutrons.get_transp_neutrons()
neutrons.get_exp()

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
ax.legend(fontsize='x-small', title='Neutron Rate')

win.addPlot('th vs. beam',fig)


win.show()



