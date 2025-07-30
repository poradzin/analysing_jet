#!/usr/bin/env python                                                                               
import plotWindow as pw                                                                             
import numpy as np                                                                                  
import matplotlib.pyplot as plt                                                                     
import profiles as ps                                                                               
import sys   
import f90nml
import ppf

# Function to sum densities for an impurity across charge states
def sum_charge_state_densities(profs, impurity):
    total_density = np.zeros_like(profs.transp('NE'))  # Initialize with zeros
    for charge_state in profs.imp_charges[impurity]:  # Loop through charge states
        signal_name = f'NIMP_{impurity}_{charge_state}'  # Example: NIMP_NI_1, NIMP_NI_2, ...
        print(f'NIMP_{impurity}_{charge_state}')
        density = profs.transp(signal_name)  # Get the density for that charge state
        total_density += charge_state * density  # Multiply by charge state to get electron count
    return total_density

def sum_nimp_with_charges(profs):
    total_density = np.zeros_like(profs.transp('NE'))
    imp_density = {}
    imp_el_density = {}
    imp_zeff = {}
    for impurity in profs.impurities:   
        imp_el_density[impurity] = np.zeros_like(profs.transp('NE'))
        imp_density[impurity] = np.zeros_like(profs.transp('NE'))
        imp_zeff[impurity] = np.zeros_like(profs.transp('NE'))
        for charge_state in profs.imp_charges[impurity]:  # Loop through charge states
            signal_name = f'NIMP_{impurity}_{charge_state}'  # Example: NIMP_NI_1, NIMP_NI_2, ...
            print(f'NIMP_{impurity}_{charge_state}')
            density = profs.transp(signal_name)  # Get the density for that charge state
            imp_el_density[impurity] += charge_state * density
            imp_zeff[impurity] += density*charge_state**2
            imp_density[impurity] += density
            total_density += charge_state * density  # Multiply by charge state to get electron count
        imp_zeff[impurity] = imp_zeff[impurity]/profs.transp('NE')
    if np.all(total_density == 0): 
        total_density = profs.xzimp*profs.transp('NIMP')
        imp_el_density = total_density
        imp_density = profs.transp('NIMP')
    return (total_density, imp_el_density, imp_density,imp_zeff)

def sum_nimps_with_zimps(profs):
    total_density = np.zeros_like(profs.transp('NE'))
    imp_el = {}
    for impurity in profs.impurities:
        imp_el[impurity] = profs.transp(f'ZIMPS_{impurity}')*profs.transp(f'NIMPS_{impurity}')
        total_density +=  imp_el[impurity]
    return (total_density, imp_el)

# List of impurities and their maximum charge state
imp_charge = {
    'BE': 4,  # Beryllium has a maximum charge state of 4
    'NI': 28, # and so on..
    'NE': 10,
    'W' : 74
}

# NADVSIM: The NADVSIM(I) namelist variable controls this for 
#each impurity element.
#  NADVSIM(I) = 1     ! redistribute the charge states to satisfy the
#                     ! coronal equilibrium
#             = 0     ! do not modify the charge state distribution (default)
nadvsim = {'BE': 0, 'NE': 0, 'NI': 0}

################################################################

pulse=int(sys.argv[1])                                                                              
runid=str(sys.argv[2]) 
time = float(sys.argv[3])

profs = ps.Transp(pulse, runid)
tind = ps.gi(profs.transp_time,time)

signals_start = ['NE','NI','NMINI','BDENS']
signals = []
#ion_signals = ['NI','NMINI','BDENS']
for sig in signals_start: 
    if profs.check_signal(sig):
        signals.append(sig)

ion_signals_start = ['NI','BDENS','NMINI']
ion_signals = []
for sig in ion_signals_start:    
    if profs.check_signal(sig):
        ion_signals.append(sig)
print(f'Impurities: {profs.impurities}')

profs.add_data(*signals) 
profs.read_imp_charge_densities()

# Sum electron densities from all charge states

#profs._imp_charges = {'BE':[4], 'NI':[28],'NE':[10]}
    
#total_el_contribution, imp_el_density = sum_nimps_with_zimps(profs)
total_imp_el_density, imps_el_densities, _ , imp_zeff= sum_nimp_with_charges(profs)


print(f'np.all(total_imp_el_density == 0): {np.all(total_imp_el_density == 0)}')
if profs.xzimps is None:
    print(f'Single impurity Z: {profs.xzimp}')
#print(f'profs._namelist: {profs._namelist}')

nsum = np.sum([profs.transp(sig) for sig in ion_signals], axis=0) + total_imp_el_density# -0.5*profs.transp('NMINI')

#slice at time

############
win=pw.plotWindow()                                                                                 
                                                                                                        
fig = plt.figure()                                                                                  
fig.suptitle(f'{profs.runcid} densities at {profs.t[tind]+40:.3f}s', fontsize=13)                                                      
ax = fig.add_subplot(111)   
ax.plot(
        profs.x[ps.gi(profs.transp_time,time)],
        profs.transp('NE')[ps.gi(profs.transp_time,time)],
        color='red',
        linewidth=2,
        label=r'$n_e$'
        )


for sig in ion_signals:
    ax.plot(
        profs.x[ps.gi(profs.transp_time,time)],
        profs.transp(sig)[ps.gi(profs.transp_time,time)],
        linewidth=2,
        label=sig
        )


ax.plot(
        profs.x[ps.gi(profs.transp_time,time)],
        nsum[ps.gi(profs.transp_time,time)],
        color='black',
        linewidth=2,
        label=' + '.join(ion_signals)+r'+ $\sum_{i,imp} Z_i n_{i,imp}$'
        )

for imp in profs.impurities: 
    ax.plot(
        profs.x[ps.gi(profs.transp_time,time)],
        imps_el_densities[imp][ps.gi(profs.transp_time,time)],
        label=imp+':'+r'$\sum Z_i n_i$ ',
        linestyle='--',
        )
if not profs.xzimps:
    ax.plot(
        profs.x[ps.gi(profs.transp_time,time)],
        imps_el_densities[ps.gi(profs.transp_time,time)],
        label=r'$Z_i n_i$ '+': '+r'$Z_i$='+str(profs.xzimp),
        linestyle='--',
        ) 
ax.set_xlabel(r'$\rho_{tor}^{norm}$')                                                                           
ax.set_ylabel(r'$cm^{-3}$')                                                                          
leg=ax.legend()                                                                                     
leg.set_draggable(True) 
win.addPlot('Density profiles',fig)  

if profs.xzimps is not None:
    fig = plt.figure()                                                                                  
    fig.suptitle(f'{profs.runcid} imp. avg. effective charge at {profs.t[tind]+40:.3f}s', fontsize=13)                                                      
    ax = fig.add_subplot(111)  
    for imp in profs.impurities: 
        ax.plot(
            profs.x[ps.gi(profs.transp_time,time)],
            profs.transp(f'ZIMPS_{imp}')[ps.gi(profs.transp_time,time)],
            label=f'ZIMPS_{imp}',
            linestyle='-',
            )         
    ax.set_xlabel(r'$\rho_{tor}^{norm}$')                                                                           
    leg=ax.legend()                                                                                     
    leg.set_draggable(True) 
    win.addPlot('Imp. avg. effective charge profiles',fig)  
   
if profs.zeffp is not None: 
    fig = plt.figure()                                                                                  
    fig.suptitle(f'{profs.runcid} 2D Effective charge at {profs.t[tind]+40:.3f}s', fontsize=13)                                                      
    ax = fig.add_subplot(111)  
    ax.plot(
        profs.x[ps.gi(profs.transp_time,time)],
        profs.zeffp[ps.gi(profs.transp_time,time)],
        label='ZEFFP',
        color='r',
        linestyle='-',
        )    
    if profs.zeffpro is not None: 
        ax.plot(
            profs.x[ps.gi(profs.transp_time,time)],     
            profs.zeffpro[ps.gi(profs.transp_time,time)],
            label='ZEFFPRO', 
            color='k',
            linestyle='--'
            )
    for imp in profs.impurities:                                                                    
        ax.plot(                                                                                    
            profs.x[ps.gi(profs.transp_time,time)],                                                 
            imp_zeff[f'{imp}'][ps.gi(profs.transp_time,time)],                            
            label=f'ZEFF_{imp}',                                                                   
            linestyle='-',                                                                          
            )
    ax.set_xlabel(r'$\rho_{tor}^{norm}$')                                                                           
    leg=ax.legend()                                                                                     
    leg.set_draggable(True) 
    win.addPlot(' Effective charge profiles',fig)                                            
                                                                                    
win.show()     

