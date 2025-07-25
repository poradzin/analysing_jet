#!/usr/bin/env python
import plotWindow as pw
import numpy as np
import matplotlib.pyplot as plt
import profiles as ps
import sys
    
pulse=int(sys.argv[1])                                                                        
runid=str(sys.argv[2])                                                                        
profs = ps.Neutrons(pulse, runid)                                                                
profs.get_transp_neutrons()

#Ion_fraction_dda = 'KS3B'
Ion_fraction_dda = 'KT5P'
kt5p = ps.EXP(pulse,runid,dda=Ion_fraction_dda)
#pulse averaged Wiesen coefficient
if all(profs.Rdtdd>0.0): 
    ratio = np.average(profs.ntnd/profs.Rdtdd,weights=profs.dt)
else: 
    ratio = 0.0
ratio_W = np.average(profs.wiesen,weights=profs.dt)

Meff = np.average(profs.meff,weights=profs.dt)
print(f'np.any(profs.ntnd_av>0): {np.any(profs.ntnd_av>0)}')
NDoverNTND = f'{1/(1+profs.ntnd_av)*100:.6f}%\n' if np.any(profs.ntnd_av>0) else '0\n'
print(NDoverNTND)
text = (
        r'$\left<\frac{n_T}{n_H+n_D+n_T}\right>=$'+f'{profs.tttd_av*100:.4f}%\n'
        r'$\left<\frac{n_D}{n_D+n_T}\right>=$'+NDoverNTND+
        #r'$ \left<n_T/n_D\right>=$'+f'{profs.ntnd_av*100:.4f}%\n'
        r'$ \left<n_T/n_e\right>=$'+f'{profs.ntne_av*100:.4f}%'
       )

print(f'Meff = {Meff}')
print(f'nT/nD = {profs.ntnd_av:6f}')
props = dict(boxstyle='round', facecolor='white', alpha=0.5)

win=pw.plotWindow()

#fn = FigureNotebook(0, 'Tritium profiles')
fig = plt.figure()
fig.suptitle(f'T concentration', fontsize=13)
ax = fig.add_subplot(111)
#fig,ax = fn.subplots(label='Tritium concentration')

font = {'family': 'serif',
    'color':  'darkgreen',
    'weight': 'normal',
    'size': 12,
    }  

ax.set_title(f'{profs.transpcid} Tritium concentration')

if profs.signal('NT'):
    ax.plot(profs.t+40.,
            profs.tttd*100, 
            color='k',
            linewidth=2,
            label='nt/(nt+nd+nh)',
            )
    ax.plot(kt5p.t,
            kt5p.tttd*100,
            color='olive',
            linewidth=2,
            linestyle=(0,(1,1)),
            label=f'{Ion_fraction_dda}/TTTD',
            )
if profs.signal('ND'):
    ax.plot(profs.t+40.,profs.dthd*100, color='b',linewidth=2,label='nd/(nt+nd+nh)')
    ax.plot(kt5p.t,kt5p.dthd*100, color='darkorange',linestyle=(0,(1,1)), label=f'{Ion_fraction_dda}/DTHD')

ax.plot(profs.t+40.,profs.ntne*100, color='r',linewidth=2,label='nt/ne',linestyle='dashed')
ax.set_xlim(profs.t[0]+40,profs.t[-1]+40)
ax.set_ylim(0., 1.4*np.amax(profs.tttd*100))
xleft,xright = ax.get_xlim()
ymin,ymax=ax.get_ylim()
ax.set_xlabel('Time [s]')
ax.set_ylabel('%')
ax.text(xleft+0.04*(xright-xleft),(ymax+ymin)/2+0.3*(ymax-ymin),text,fontdict=font,bbox=props)                                                                                   
#cornernote(f'{profs.transpcid}', '', ax=ax)                                                  
ax.legend() 
win.addPlot('T conc',fig)                                                                     
                                                                                              
#fig,ax = fn.subplots(label='Neutron rates vs T density')                                     
fig=plt.figure()
fig.suptitle(f'Neutron rates vs T density', fontsize=13)                                       
ax = fig.add_subplot(111)                                                                      
                                                                                               
text = (                                                                                       
        r'$T_{conc}=\left<\frac{n_T}{n_H+n_D+n_T}\right>=$'+f'{profs.tttd_av*100:.3f}\n'
        r'$D_{frac}=\left<\frac{n_D}{n_D+n_T}\right>=$'+f'{1/(1+profs.ntnd_av):.6f}\n'
        #r'$ \left<n_T/n_D\right>=$'+f'{profs.ntnd_av:.6f}\n'                                           
        #r'$ \left<R_{DT/DD}\right>=\left<\frac{R_{nDT}}{R_{nDD}}\right>=$'+f'{profs.Rdtdd_av:.4f}\n'   
        #r'$ \left<(n_T/n_D)/R_{DT/DD}\right> = $'+f'{ratio:.4f}'+'\n'                                  
        #r'$ \left< Wiesen\,\, coefficient\right>:$'+f' {ratio_W:.4f}\n'                                
        #r'$n_D,n_T,R_{nDT},R_{nTT}$ volume integrated.'                                                
       )                                                                                              
                                                                                               
ax.set_title(f'{profs.transpcid} Neutron rates vs T density')                                  
if any(profs.Rdtdd>0):                                                                                               
    ax.plot(profs.t+40.,
            profs.ntnd/profs.Rdtdd, 
            color='k',
            linewidth=2,
            label=r'$(nt/nd)/R_{DT/DD}$'
            )
if ratio>0: 
    ax.plot([profs.t[0]+40.,profs.t[-1]+40.],
            [ratio,ratio],
            linestyle='dashed',
            color='darkred',
            linewidth=2,
            label= 'average'
            )                                                                                      
                                                                                                   
ax.set_xlim(profs.t[0]+40,profs.t[-1]+40)
if any(profs.Rdtdd>0):
    ax.set_ylim(0.5*np.amin(profs.ntnd/profs.Rdtdd), 1.5*np.amax(profs.ntnd/profs.Rdtdd))          
xleft,xright = ax.get_xlim()                                                                   
ymin,ymax=ax.get_ylim()                                                                        
ax.set_xlabel('Time [s]')                                                                      
#ax.set_ylabel('%')                                                                            
                                                                                               
ax.text(xleft+0.04*(xright-xleft),(ymax+ymin)/2+0.12*(ymax-ymin),text,fontdict=font,bbox=props)
                                                                                               
#cornernote(f'{profs.transpcid}', '', ax=ax)                                                   
ax.legend()                                                                                    
win.addPlot('Neutron vs T density',fig)                                                        
                                                                                               
#fig,ax = fn.subplots(label=r'$R_{DT/DD}$')                                                    
fig = plt.figure()                                                                             
fig.suptitle(f'R_DT/DD', fontsize=13)                                                          
ax = fig.add_subplot(111)                                                                      
                                                                                               
                                                                                               
text = (                                                                                       
        r'$T_{conc}=\left<\frac{n_T}{n_H+n_D+n_T}\right>=$'+f'{profs.tttd_av*100:.3f}%\n'                  
        r'$ \left<n_T/n_D\right>=$'+f'{profs.ntnd_av:.6f}\n'                                           
        r'$ \left<R_{DT/DD}\right>=\left<\frac{R_{nDT}}{R_{nDD}}\right>=$'+f'{profs.Rdtdd_av:.4f}'     
       )                                                                                              
                                                                                               
ax.set_title(f'{profs.transpcid} Neutron rates vs T density')                                  
                                                                                               
ax.plot(profs.t+40, 
        profs.Rdtdd,
        color='b',
        linewidth=2,
        label=r'$R_{DT/DD}$'
        )                                                                                       
ax.plot([profs.t[0]+40.,profs.t[-1]+40.],
        [profs.Rdtdd_av,profs.Rdtdd_av],
        linestyle='dashed',
        color='darkred',
        linewidth=2,
        label= 'average'
        )
                                                            
ax.set_xlim(profs.t[0]+40,profs.t[-1]+40)
if all(profs.Rdtdd>0):
    ax.set_ylim(0.5*np.amin(profs.Rdtdd), 1.5*np.amax(profs.Rdtdd))
xleft,xright = ax.get_xlim()
ymin,ymax=ax.get_ylim()
ax.set_xlabel('Time [s]')
#ax.set_ylabel('%')

ax.text(xleft+0.04*(xright-xleft),(ymax+ymin)/2+0.26*(ymax-ymin),text,fontdict=font,bbox=props)
#cornernote(f'{profs.transpcid}', '', ax=ax)
ax.legend()
                                                                 
win.addPlot('R_DTDD',fig)
win.show()

#for key in profs._dat.variables.keys():
#    if 'MINORITY DENSITY' in profs._dat.variables[key].long_name: 
#        print(profs._dat.variables[key])
#        #print(profs._dat.variables[key].long_name)

