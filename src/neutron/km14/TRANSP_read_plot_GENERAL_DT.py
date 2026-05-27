# Script - reading CDF files for the analysis of TRANSP output
"""
Created on Tue Dec 06 15:42:45 2016

@author: k3nob1
"""

import os
import numpy as np
import scipy as sc
import scipy.interpolate as scinterpolate
import matplotlib as mpl
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import matplotlib.tri as tri
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.patches as patches
#import ppf


mpl.rcParams['font.family'] = 'Gulliver'
mpl.rcParams.update({'font.size': 16})

# Read the NetCDF file with the Dataset module

RunID = '104614M30'
pulse = int(RunID[:-3])

print(RunID)
# Read n


# Read the NetCDF file with the Dataset module
TRANSP_dir = f"/home/{os.environ['USER']}/jet/data/"
pulse_dir = TRANSP_dir+f'{RunID[:-3]}/'
run_dir = pulse_dir+f'{RunID[-3:]}/'

BasicDat = Dataset(run_dir+f'{RunID}.CDF', 'r')
BasicDat_list = []
for item in BasicDat.variables.keys():
    BasicDat_list.append(item)

TRANSPDat = Dataset(run_dir+f'{RunID}_fi_1.cdf', 'r')
TRANSPDat_list = []
for item in TRANSPDat.variables.keys():
    TRANSPDat_list.append(item)
print('FBM interval time: ', np.array(TRANSPDat.variables['TIME']))

#LimiterDat = Dataset('JET_Limiter_Eq_99972.cdf', 'r')
LimiterDat_list = []
#for item in LimiterDat.variables.keys():
#    LimiterDat_list.append(item)

NeutronDat = Dataset(run_dir+f'{RunID}_neut_1.cdf', 'r')
NeutronDat_list = []
for item in NeutronDat.variables.keys():
    NeutronDat_list.append(item)

# JET Shot information

BasicTime=np.array(BasicDat.variables['TIME'])
Time=np.array(TRANSPDat.variables['TIME'])

# Neutron emission time-trace data

TotalMeasNeutrons=np.array(BasicDat.variables['MNEUT'])
TotalNeutrons=np.array(BasicDat.variables['NEUTT'])
TotalTHNeutrons=np.array(BasicDat.variables['NEUTX'])
TotalBTNeutrons=np.array(BasicDat.variables['BTNTS'])
TotalBBNeutrons=np.array(BasicDat.variables['BBNTS'])

TotNeutSum=TotalTHNeutrons+TotalBTNeutrons+TotalBBNeutrons
TotNeutDiff=(TotalNeutrons*1.15/TotalMeasNeutrons-1)*100

# Neutron emissivity density profile data

NeutDensBB_DD=np.array(NeutronDat.variables['BBN4'])
NeutDensBB = NeutDensBB_DD
if 'BBN5' in NeutronDat_list:
    print('Neutron output BB TT present!')
    NeutDensBB_TT=np.array(NeutronDat.variables['BBN5'])
    NeutDensBB += NeutDensBB_TT
if 'BBN1' in NeutronDat_list:
    print('Neutron output BB DT present!')
    NeutDensBB_DT=np.array(NeutronDat.variables['BBN1'])
    NeutDensBB += NeutDensBB_DT

NeutDensBT_DD=np.array(NeutronDat.variables['BTN4'])
NeutDensBT = NeutDensBT_DD
if 'BTN5' in NeutronDat_list:
    print('Neutron output BT TT present!')
    NeutDensBT_TT=np.array(NeutronDat.variables['BTN5'])
    NeutDensBT += NeutDensBT_TT
if 'BTN1' in NeutronDat_list:
    print('Neutron output BT DT present!')
    NeutDensBT_DT=np.array(NeutronDat.variables['BTN1'])
    NeutDensBT += NeutDensBT_DT
if 'BTN7' in NeutronDat_list:
    print('Neutron output BT TD present!')
    NeutDensBT_TD=np.array(NeutronDat.variables['BTN7'])
    NeutDensBT += NeutDensBT_TD

NeutDensTH=np.array(NeutronDat.variables['THNTNT2d'])

NeutTotF2D=np.array(NeutronDat.variables['TOTNTNF2d'])

VolumeGrid=np.array(TRANSPDat.variables['BMVOL'])

# Read the R (major radius), Z (tokamak height) grid, limiter, plasma boundary and equilibrium into arrays

R2D=np.array(TRANSPDat.variables['R2D'])
Z2D=np.array(TRANSPDat.variables['Z2D'])
R2DN=np.array(NeutronDat.variables['R'])
Z2DN=np.array(NeutronDat.variables['Z'])
#RLIM=np.array(LimiterDat.variables['RLIM_CDF'])
#ZLIM=np.array(LimiterDat.variables['ZLIM_CDF'])
#RBND=np.array(LimiterDat.variables['RBND_Eq_CDF'])
#ZBND=np.array(LimiterDat.variables['ZBND_Eq_CDF'])
#RMAG=np.array(LimiterDat.variables['RMAG_Eq_CDF'])
#ZMAG=np.array(LimiterDat.variables['ZMAG_Eq_CDF'])

#########################################################################################
# Read the fast ion distribution function of D-NBI into arrays if present in output files

if 'F_D_NBI' in TRANSPDat_list:
    IonDist_DNBI=np.array(TRANSPDat.variables['F_D_NBI'])
    
    # Read the energy and pitch data - bin boundaries and centers - into arrays
    Ebound_DNBI=np.array(TRANSPDat.variables['EB_D_NBI'])
    Ecenter_DNBI=np.array(TRANSPDat.variables['E_D_NBI'])
    pbound_DNBI=np.array(TRANSPDat.variables['AB_D_NBI'])
    pcenter_DNBI=np.array(TRANSPDat.variables['A_D_NBI'])
    print(Ebound_DNBI[0], Ebound_DNBI[75], pbound_DNBI[0], pbound_DNBI[50])
    
    # Construct arrays containing the widths of energy and pitch bins
    # and calculate the fast ion densities for different ion species

    dE_DNBI=np.empty(len(Ebound_DNBI)-1) # DNBI component
    dp_DNBI=np.empty(len(pbound_DNBI)-1)
    for i in range (0,len(Ebound_DNBI)-1):
        dE_DNBI[i]=Ebound_DNBI[i+1]-Ebound_DNBI[i]
    for i in range (0,len(pbound_DNBI)-1):
        dp_DNBI[i]=pbound_DNBI[i+1]-pbound_DNBI[i]
    dpdE_DNBI=dE_DNBI*dp_DNBI[:,np.newaxis] # dp*dE matrix for integration
    n_DNBI=np.empty(len(R2D)) # Fast ion density array
    for i in range (0,len(R2D)): # Integration of the distribution function
        Int_dpdE_Temp=float(np.sum(dpdE_DNBI*0.5*IonDist_DNBI[i]))
        n_DNBI[i]=Int_dpdE_Temp*1000000 # Fast ion density distribution in FI/m^3

    # Construct arrays containing the centers of energy and pitch bins
    # and calculate the average fast ion energy and pitch distributions
    
    dpdEcenterE_DNBI=(dE_DNBI*Ecenter_DNBI)*dp_DNBI[:,np.newaxis] # dpdE*E matrix
    dpdEcenterp_DNBI=dE_DNBI*(dp_DNBI*pcenter_DNBI)[:,np.newaxis] # dpdE*p matrix
    
    AvgE_DNBI=np.empty(len(R2D)) # Average fast ion energy array
    for i in range (0,len(R2D)): # Integration of the distribution function
        Int_dpdE_Temp=float(np.sum(dpdEcenterE_DNBI*0.5*IonDist_DNBI[i])/np.sum(dpdE_DNBI*0.5*IonDist_DNBI[i]))
        AvgE_DNBI[i]=Int_dpdE_Temp # Average fast ion energy distribution in FI/m^3
    
    Avgp_DNBI=np.empty(len(R2D)) # Average fast ion energy array
    for i in range (0,len(R2D)): # Integration of the distribution function
        Int_dpdE_Temp=float(np.sum(dpdEcenterp_DNBI*0.5*IonDist_DNBI[i])/np.sum(dpdE_DNBI*0.5*IonDist_DNBI[i]))
        Avgp_DNBI[i]=Int_dpdE_Temp # Average fast ion energy distribution in FI/m^3
        
    # Construct arrays and calculate the fast ion energy and pitch distribution
    # function for specific grid location and averaged over whole plasma
    
    EmptyE=np.ones(len(Ecenter_DNBI)) # Dummy matrix containing ones with E dimension
    Emptyp=np.ones(len(pcenter_DNBI)) # Dummy matrix containing ones with p dimension
    N_FIR_E_DNBI={}
    N_FIR_EdV_DNBI={}
    N_FIR_p_DNBI={}
    N_FIR_pdV_DNBI={}
    N_FI_E_DNBI=np.zeros(len(Ecenter_DNBI))
    N_FI_p_DNBI=np.zeros(len(pcenter_DNBI))
    dp_DNBI=EmptyE*dp_DNBI[:,np.newaxis] # dp matrix
    dE_DNBI=dE_DNBI*Emptyp[:,np.newaxis] # dE matrix
    for i in range(0,len(R2D)):
        N_FIR_E_DNBI[i]=np.sum(dp_DNBI*0.5*IonDist_DNBI[i],axis=0) # Integration of f_i over dp
        N_FIR_EdV_DNBI[i]=N_FIR_E_DNBI[i]*VolumeGrid[i] # Multiplying with grid voxel volume
        N_FIR_p_DNBI[i]=np.sum(dE_DNBI*0.5*IonDist_DNBI[i],axis=1) # Integration of f_i over dE
        N_FIR_pdV_DNBI[i]=N_FIR_p_DNBI[i]*VolumeGrid[i] # Multiplying with grid voxel volume
        for j in range(0,len(Ecenter_DNBI)):
            N_FI_E_DNBI[j] += N_FIR_EdV_DNBI[i][j] # Fast ion energy distribution function
        for j in range(0,len(pcenter_DNBI)):
            N_FI_p_DNBI[j] += N_FIR_pdV_DNBI[i][j] # Fast ion pitch distribution function
    
    plt.figure()
    plt.plot(Ecenter_DNBI, N_FI_E_DNBI, 'r--')
    plt.show()
    plt.clf
    
    plt.figure()
    plt.plot(pcenter_DNBI, N_FI_p_DNBI, 'r--')
    plt.show()
    plt.clf

#########################################################################################
# Read the fast ion distribution function of T-NBI into arrays if present in output files

if 'F_T_NBI' in TRANSPDat_list:
    IonDist_TNBI=np.array(TRANSPDat.variables['F_T_NBI'])
    
    # Read the energy and pitch data - bin boundaries and centers - into arrays
    Ebound_TNBI=np.array(TRANSPDat.variables['EB_T_NBI'])
    Ecenter_TNBI=np.array(TRANSPDat.variables['E_T_NBI'])
    pbound_TNBI=np.array(TRANSPDat.variables['AB_T_NBI'])
    pcenter_TNBI=np.array(TRANSPDat.variables['A_T_NBI'])
    print(Ebound_TNBI[0], Ebound_TNBI[75], pbound_TNBI[0], pbound_TNBI[50])
    
    # Construct arrays containing the widths of energy and pitch bins
    # and calculate the fast ion densities for different ion species

    dE_TNBI=np.empty(len(Ebound_TNBI)-1) # TNBI component
    dp_TNBI=np.empty(len(pbound_TNBI)-1)
    for i in range (0,len(Ebound_TNBI)-1):
        dE_TNBI[i]=Ebound_TNBI[i+1]-Ebound_TNBI[i]
    for i in range (0,len(pbound_TNBI)-1):
        dp_TNBI[i]=pbound_TNBI[i+1]-pbound_TNBI[i]
    dpdE_TNBI=dE_TNBI*dp_TNBI[:,np.newaxis] # dp*dE matrix for integration
    n_TNBI=np.empty(len(R2D)) # Fast ion density array
    for i in range (0,len(R2D)): # Integration of the distribution function
        Int_dpdE_Temp=float(np.sum(dpdE_TNBI*0.5*IonDist_TNBI[i]))
        n_TNBI[i]=Int_dpdE_Temp*1000000 # Fast ion density distribution in FI/m^3

    # Construct arrays containing the centers of energy and pitch bins
    # and calculate the average fast ion energy and pitch distributions
    
    dpdEcenterE_TNBI=(dE_TNBI*Ecenter_TNBI)*dp_TNBI[:,np.newaxis] # dpdE*E matrix
    dpdEcenterp_TNBI=dE_TNBI*(dp_TNBI*pcenter_TNBI)[:,np.newaxis] # dpdE*p matrix
    
    AvgE_TNBI=np.empty(len(R2D)) # Average fast ion energy array
    for i in range (0,len(R2D)): # Integration of the distribution function
        Int_dpdE_Temp=float(np.sum(dpdEcenterE_TNBI*0.5*IonDist_TNBI[i])/np.sum(dpdE_TNBI*0.5*IonDist_TNBI[i]))
        AvgE_TNBI[i]=Int_dpdE_Temp # Average fast ion energy distribution in FI/m^3
    
    Avgp_TNBI=np.empty(len(R2D)) # Average fast ion energy array
    for i in range (0,len(R2D)): # Integration of the distribution function
        Int_dpdE_Temp=float(np.sum(dpdEcenterp_TNBI*0.5*IonDist_TNBI[i])/np.sum(dpdE_TNBI*0.5*IonDist_TNBI[i]))
        Avgp_TNBI[i]=Int_dpdE_Temp # Average fast ion energy distribution in FI/m^3
        
    # Construct arrays and calculate the fast ion energy and pitch distribution
    # function for specific grid location and averaged over whole plasma
    
    EmptyE=np.ones(len(Ecenter_TNBI)) # Dummy matrix containing ones with E dimension
    Emptyp=np.ones(len(pcenter_TNBI)) # Dummy matrix containing ones with p dimension
    N_FIR_E_TNBI={}
    N_FIR_EdV_TNBI={}
    N_FIR_p_TNBI={}
    N_FIR_pdV_TNBI={}
    N_FI_E_TNBI=np.zeros(len(Ecenter_TNBI))
    N_FI_p_TNBI=np.zeros(len(pcenter_TNBI))
    dp_TNBI=EmptyE*dp_TNBI[:,np.newaxis] # dp matrix
    dE_TNBI=dE_TNBI*Emptyp[:,np.newaxis] # dE matrix
    for i in range(0,len(R2D)):
        N_FIR_E_TNBI[i]=np.sum(dp_TNBI*0.5*IonDist_TNBI[i],axis=0) # Integration of f_i over dp
        N_FIR_EdV_TNBI[i]=N_FIR_E_TNBI[i]*VolumeGrid[i] # Multiplying with grid voxel volume
        N_FIR_p_TNBI[i]=np.sum(dE_TNBI*0.5*IonDist_TNBI[i],axis=1) # Integration of f_i over dE
        N_FIR_pdV_TNBI[i]=N_FIR_p_TNBI[i]*VolumeGrid[i] # Multiplying with grid voxel volume
        for j in range(0,len(Ecenter_TNBI)):
            N_FI_E_TNBI[j] += N_FIR_EdV_TNBI[i][j] # Fast ion energy distribution function
        for j in range(0,len(pcenter_TNBI)):
            N_FI_p_TNBI[j] += N_FIR_pdV_TNBI[i][j] # Fast ion pitch distribution function
    
    plt.figure()
    plt.plot(Ecenter_TNBI, N_FI_E_TNBI, 'r--')
    plt.show()
    plt.clf
    
    plt.figure()
    plt.plot(pcenter_TNBI, N_FI_p_TNBI, 'r--')
    plt.show()
    plt.clf

#########################################################################################
# Read the fast ion distribution function of HE4 into arrays if present in output files

if 'F_He4_FUSN' in TRANSPDat_list:
    IonDist_HE4FUSN=np.array(TRANSPDat.variables['F_He4_FUSN'])
    
    # Read the energy and pitch data - bin boundaries and centers - into arrays
    Ebound_HE4FUSN=np.array(TRANSPDat.variables['EB_He4_FUSN'])
    Ecenter_HE4FUSN=np.array(TRANSPDat.variables['E_He4_FUSN'])
    pbound_HE4FUSN=np.array(TRANSPDat.variables['AB_He4_FUSN'])
    pcenter_HE4FUSN=np.array(TRANSPDat.variables['A_He4_FUSN'])
    print(Ebound_HE4FUSN[0], Ebound_HE4FUSN[75], pbound_HE4FUSN[0], pbound_HE4FUSN[50])
    
    # Construct arrays containing the widths of energy and pitch bins
    # and calculate the fast ion densities for different ion species

    dE_HE4FUSN=np.empty(len(Ebound_HE4FUSN)-1) # HE4FUSN component
    dp_HE4FUSN=np.empty(len(pbound_HE4FUSN)-1)
    for i in range (0,len(Ebound_HE4FUSN)-1):
        dE_HE4FUSN[i]=Ebound_HE4FUSN[i+1]-Ebound_HE4FUSN[i]
    for i in range (0,len(pbound_HE4FUSN)-1):
        dp_HE4FUSN[i]=pbound_HE4FUSN[i+1]-pbound_HE4FUSN[i]
    dpdE_HE4FUSN=dE_HE4FUSN*dp_HE4FUSN[:,np.newaxis] # dp*dE matrix for integration
    n_HE4FUSN=np.empty(len(R2D)) # Fast ion density array
    for i in range (0,len(R2D)): # Integration of the distribution function
        Int_dpdE_Temp=float(np.sum(dpdE_HE4FUSN*0.5*IonDist_HE4FUSN[i]))
        n_HE4FUSN[i]=Int_dpdE_Temp*1000000 # Fast ion density distribution in FI/m^3

    # Construct arrays containing the centers of energy and pitch bins
    # and calculate the average fast ion energy and pitch distributions
    
    dpdEcenterE_HE4FUSN=(dE_HE4FUSN*Ecenter_HE4FUSN)*dp_HE4FUSN[:,np.newaxis] # dpdE*E matrix
    dpdEcenterp_HE4FUSN=dE_HE4FUSN*(dp_HE4FUSN*pcenter_HE4FUSN)[:,np.newaxis] # dpdE*p matrix
    
    AvgE_HE4FUSN=np.empty(len(R2D)) # Average fast ion energy array
    for i in range (0,len(R2D)): # Integration of the distribution function
        Int_dpdE_Temp=float(np.sum(dpdEcenterE_HE4FUSN*0.5*IonDist_HE4FUSN[i])/np.sum(dpdE_HE4FUSN*0.5*IonDist_HE4FUSN[i]))
        AvgE_HE4FUSN[i]=Int_dpdE_Temp # Average fast ion energy distribution in FI/m^3
    
    Avgp_HE4FUSN=np.empty(len(R2D)) # Average fast ion energy array
    for i in range (0,len(R2D)): # Integration of the distribution function
        Int_dpdE_Temp=float(np.sum(dpdEcenterp_HE4FUSN*0.5*IonDist_HE4FUSN[i])/np.sum(dpdE_HE4FUSN*0.5*IonDist_HE4FUSN[i]))
        Avgp_HE4FUSN[i]=Int_dpdE_Temp # Average fast ion energy distribution in FI/m^3
        
    # Construct arrays and calculate the fast ion energy and pitch distribution
    # function for specific grid location and averaged over whole plasma
    
    EmptyE=np.ones(len(Ecenter_HE4FUSN)) # Dummy matrix containing ones with E dimension
    Emptyp=np.ones(len(pcenter_HE4FUSN)) # Dummy matrix containing ones with p dimension
    N_FIR_E_HE4FUSN={}
    N_FIR_EdV_HE4FUSN={}
    N_FIR_p_HE4FUSN={}
    N_FIR_pdV_HE4FUSN={}
    N_FI_E_HE4FUSN=np.zeros(len(Ecenter_HE4FUSN))
    N_FI_p_HE4FUSN=np.zeros(len(pcenter_HE4FUSN))
    dp_HE4FUSN=EmptyE*dp_HE4FUSN[:,np.newaxis] # dp matrix
    dE_HE4FUSN=dE_HE4FUSN*Emptyp[:,np.newaxis] # dE matrix
    for i in range(0,len(R2D)):
        N_FIR_E_HE4FUSN[i]=np.sum(dp_HE4FUSN*0.5*IonDist_HE4FUSN[i],axis=0) # Integration of f_i over dp
        N_FIR_EdV_HE4FUSN[i]=N_FIR_E_HE4FUSN[i]*VolumeGrid[i] # Multiplying with grid voxel volume
        N_FIR_p_HE4FUSN[i]=np.sum(dE_HE4FUSN*0.5*IonDist_HE4FUSN[i],axis=1) # Integration of f_i over dE
        N_FIR_pdV_HE4FUSN[i]=N_FIR_p_HE4FUSN[i]*VolumeGrid[i] # Multiplying with grid voxel volume
        for j in range(0,len(Ecenter_HE4FUSN)):
            N_FI_E_HE4FUSN[j] += N_FIR_EdV_HE4FUSN[i][j] # Fast ion energy distribution function
        for j in range(0,len(pcenter_HE4FUSN)):
            N_FI_p_HE4FUSN[j] += N_FIR_pdV_HE4FUSN[i][j] # Fast ion pitch distribution function
    
    plt.figure()
    plt.plot(Ecenter_HE4FUSN, N_FI_E_HE4FUSN, 'r--')
    plt.show()
    plt.clf
    
    plt.figure()
    plt.plot(pcenter_HE4FUSN, N_FI_p_HE4FUSN, 'r--')
    plt.show()
    plt.clf

n_FI = n_DNBI + n_HE4FUSN # Combined fast ion density distribution in FI/m^3
      
# Construct arrays with neutron emission rate profiles - volume multiplication

NeutBT=NeutDensBT*VolumeGrid
NeutBB=NeutDensBB*VolumeGrid
NeutDensBTBB=NeutDensBT+NeutDensBB
NeutDensTH_Comp=NeutTotF2D-NeutDensBTBB

NeutDensTH_Mod=np.empty(len(NeutDensTH_Comp))
for i in range(0, len(NeutDensTH_Comp)):
    if NeutDensTH_Comp[i] <= 0:
        NeutDensTH_Mod[i]=0
    else:
        NeutDensTH_Mod[i]=NeutDensTH_Comp[i]

NeutDensTot_Comp=NeutDensBT+NeutDensBB+NeutDensTH
NeutDensTot_CompNorm=NeutDensTot_Comp/np.max(NeutDensTot_Comp)
NeutDensTHNorm=NeutDensTH/np.max(NeutDensTH)
NeutDensBT_RelDif=(NeutDensTot_CompNorm/NeutDensTHNorm-1)*100

# *******************  PLOTTING  *******************************

# Create an unstructured triangular grid for plotting with R and Z coordinates
R2DM=R2D/100.0
Z2DM=Z2D/100.0
print(np.min(R2DM))
print(np.min(Z2DM))
print(np.max(R2DM))
print(np.max(Z2DM))
triang = tri.Triangulation(R2DM, Z2DM)

# Fast ion GRID
plt.figure()
plt.gca().set_aspect('equal')
plt.scatter(R2DM, Z2DM, s=4, color='black')
#plt.plot(RLIM, ZLIM, linestyle='solid', color='black')
#plt.plot(RBND, ZBND, linestyle='solid', color='red')
#plt.plot(RMAG, ZMAG, marker='x', mew=1.5, ms=6, color='red')
#plt.scatter(R2DM[1], Z2DM[1], color='green')
plt.xticks(np.arange(2.0, 4.0, 0.5))
plt.ylabel('JET Z-axis [m]')
plt.xlabel('JET R-major radius [m]')
plt.title('Computational grid ({0})\n 99972' .format(len(R2DM)), fontsize='12')
#plt.savefig('JET_LimEq_FastIonDens_Shot{0}_T3.png' .format(RunID), dpi=600)
plt.show()

# Fast ion density ALL ION SPECIES
R=np.linspace(np.min(R2DM), np.max(R2DM), 51)
Z=np.linspace(np.min(Z2DM), np.max(Z2DM), 51)
n_FI_Int2=scinterpolate.griddata((R2DM, Z2DM), n_FI, (R[None,:], Z[:,None]), method='cubic')
plt.figure()
plt.gca().set_aspect('equal')
plt.tricontourf(triang, n_FI, 100, linewidths=3, cmap=plt.cm.YlOrRd, antialiased=False)
plt.colorbar()
#plt.contour(R,Z,n_FI_Int2, 20, linewidths=0.5, colors='black')
#plt.tricontour(triang, n_FI, 15, linewidths=0.5, colors='k')
#plt.triplot(triang, '-k')
#plt.plot(RLIM, ZLIM, linestyle='solid', color='black')
#plt.plot(RBND, ZBND, linestyle='solid', color='red')
#plt.plot(RMAG, ZMAG, marker='x', mew=1.5, ms=6, color='black')
plt.xticks(np.arange(2.0, 4.0, 0.5))
plt.ylabel('JET Z-axis [m]')
plt.xlabel('JET R-major radius [m]')
plt.title('Fast Ion Density [m$^{-3}$]'+'\n 99972 - ICRH', fontsize='16')
#plt.savefig('JET_LimEq_FastIonDens_Shot{0}_T3.png' .format(RunID), dpi=600)
plt.show()

# Average ion energy distribution
AvgE_DNBI_Int=scinterpolate.griddata((R2DM, Z2DM), AvgE_DNBI, (R[None,:], Z[:,None]), method='cubic')
plt.figure()
plt.gca().set_aspect('equal')
plt.tricontourf(triang, AvgE_DNBI, 100, linewidths=3, cmap=plt.cm.plasma)
plt.colorbar()
#plt.contour(R,Z,AvgE_DNBI_Int, 15, linewidths=0.5, colors='black')
#plt.plot(RLIM, ZLIM, linestyle='solid', color='black')
#plt.plot(RBND, ZBND, linestyle='solid', color='red')
#plt.plot(RMAG, ZMAG, marker='x', mew=1.5, ms=6, color='red')
plt.ylabel('JET Z-axis [m]')
plt.xlabel('JET R-major radius [m]')
plt.title('Average Fast Ion Energy [keV/m^3]\n Shot {0}' .format(RunID), fontsize='12')
#plt.savefig('JET_LimEq_FastIonEnergy_Shot{0}_T3.png' .format(RunID), dpi=600)
plt.show()

# Average ion pitch distribution
Avgp_DNBI_Int=scinterpolate.griddata((R2DM, Z2DM), Avgp_DNBI, (R[None,:], Z[:,None]), method='cubic')
plt.figure()
plt.gca().set_aspect('equal')
plt.tricontourf(triang, Avgp_DNBI, 100, linewidths=3, cmap=plt.cm.plasma)
plt.colorbar()
#plt.contour(R,Z,Avgp_DNBI_Int, 10, linewidths=0.5, colors='black')
#plt.plot(RLIM, ZLIM, linestyle='solid', color='black')
#plt.plot(RBND, ZBND, linestyle='solid', color='red')
#plt.plot(RMAG, ZMAG, marker='x', mew=1.5, ms=6, color='red')
plt.ylabel('JET Z-axis [m]')
plt.xlabel('JET R-major radius [m]')
plt.title('Average Fast Ion Pitch Angle [m^-3]\n Shot {0}' .format(RunID), fontsize='12')
#plt.savefig('JET_LimEq_FastIonPitch_Shot{0}_T3.png' .format(RunID), dpi=600)
plt.show()

# BT neutron emission density profile
NeutDensBT_Int=scinterpolate.griddata((R2DM, Z2DM), NeutDensBT, (R[None,:], Z[:,None]), method='cubic')
plt.figure()
plt.gca().set_aspect('equal')
plt.tricontourf(triang, NeutDensBT, 100, linewidths=3, cmap=plt.cm.plasma)
plt.colorbar()
plt.contour(R,Z,NeutDensBT_Int, 10, linewidths=0.5, colors='black')
#plt.scatter(R2DM, Z2DM, s=4, color='black')
#plt.plot(RLIM, ZLIM, linestyle='solid', color='black')
#plt.plot(RBND, ZBND, linestyle='solid', color='red')
#plt.plot(RMAG, ZMAG, marker='x', mew=1.5, ms=6, color='red')
plt.ylabel('JET Z-axis [m]')
plt.xlabel('JET R-major radius [m]')
plt.title('BT neutron emission rate density [cm^-3 s^-1]\n Shot {0}' .format(RunID), fontsize='12')
#plt.savefig('JET_LimEq_NeutDensBT_Shot{0}_T3.png' .format(RunID), dpi=600)
plt.show()

# BB neutron emission density profile
NeutDensBB_Int=scinterpolate.griddata((R2DM, Z2DM), NeutDensBB, (R[None,:], Z[:,None]), method='cubic')
plt.figure()
plt.gca().set_aspect('equal')
plt.tricontourf(triang, NeutDensBB, 100, linewidths=3, cmap=plt.cm.plasma)
plt.colorbar()
plt.contour(R,Z,NeutDensBB_Int, 10, linewidths=0.5, colors='black')
#plt.plot(RLIM, ZLIM, linestyle='solid', color='black')
#plt.plot(RBND, ZBND, linestyle='solid', color='red')
#plt.plot(RMAG, ZMAG, marker='x', mew=1.5, ms=6, color='red')
plt.ylabel('JET Z-axis [m]')
plt.xlabel('JET R-major radius [m]')
plt.title('BB neutron emission rate density [cm^-3 s^-1]\n Shot {0}' .format(RunID), fontsize='12')
#plt.savefig('JET_LimEq_NeutDensBB_Shot{0}_T3.png' .format(RunID), dpi=600)
plt.show()

# Thermal neutron emission density profile
NeutDensTH_Int=scinterpolate.griddata((R2DM, Z2DM), NeutDensTH, (R[None,:], Z[:,None]), method='cubic')
plt.figure()
plt.gca().set_aspect('equal')
plt.tricontourf(triang, NeutDensTH, 100, linewidths=3, cmap=plt.cm.plasma)
plt.colorbar()
plt.contour(R,Z,NeutDensTH_Int, 10, linewidths=0.5, colors='black')
#plt.plot(RLIM, ZLIM, linestyle='solid', color='black')
#plt.plot(RBND, ZBND, linestyle='solid', color='red')
#plt.plot(RMAG, ZMAG, marker='x', mew=1.5, ms=6, color='red')
plt.ylabel('JET Z-axis [m]')
plt.xlabel('JET R-major radius [m]')
plt.title('TH neutron emission rate density [cm^-3 s^-1]\n Shot {0}' .format(RunID), fontsize='12')
#plt.savefig('JET_LimEq_NeutDensTH_Shot{0}_T3.png' .format(RunID), dpi=600)
plt.show()

# Thermal neutron emission density profile normalized
NeutDensTHNorm_Int=scinterpolate.griddata((R2DM, Z2DM), NeutDensTHNorm, (R[None,:], Z[:,None]), method='cubic')
plt.figure()
plt.gca().set_aspect('equal')
plt.tricontourf(triang, NeutDensTHNorm, 100, linewidths=3, cmap=plt.cm.plasma)
plt.colorbar()
plt.contour(R,Z,NeutDensTHNorm_Int, 10, linewidths=0.5, colors='black')
#plt.plot(RLIM, ZLIM, linestyle='solid', color='black')
#plt.plot(RBND, ZBND, linestyle='solid', color='red')
#plt.plot(RMAG, ZMAG, marker='x', mew=1.5, ms=6, color='red')
plt.ylabel('JET Z-axis [m]')
plt.xlabel('JET R-major radius [m]')
plt.title('TH neutron emission rate density [cm^-3 s^-1]\n Shot {0}' .format(RunID), fontsize='12')
#plt.savefig('JET_LimEq_NeutDensTH_Shot{0}_T3.png' .format(RunID), dpi=600)
plt.show()

# Computed Total neutron emission density profile
NeutDensTot_Comp_Int=scinterpolate.griddata((R2DM, Z2DM), NeutDensTot_Comp, (R[None,:], Z[:,None]), method='cubic')
plt.figure()
plt.gca().set_aspect('equal')
plt.tricontourf(triang, NeutDensTot_Comp, 100, linewidths=3, cmap=plt.cm.plasma)
plt.colorbar()
plt.contour(R,Z,NeutDensTot_Comp_Int, 10, linewidths=0.5, colors='black')
#plt.scatter(R2DM, Z2DM, s=4, color='black')
#plt.contour(R,Z,NeutDensTot_Comp_Int, levels= [1.3E+9], linewidths=1.5, colors='white')
#plt.plot(RLIM, ZLIM, linestyle='solid', color='black')
#plt.plot(RBND, ZBND, linestyle='solid', color='red')
plt.plot(3.033,0.273, marker='x', mew=1.5, ms=6, color='white')
#plt.plot((2.965, 2.965), (np.min(ZBND)+0.4, np.max(ZBND)-0.2), 'w--', linewidth=1.5)
#plt.text(2.61, -0.9, r'$2\omega_{D}$', {'color': 'w', 'fontsize': 14})
plt.ylabel('JET Z-axis [m]')
plt.xlabel('JET R-major radius [m]')
plt.xticks(np.arange(2.0, 4.0, 0.5))
plt.title('Total neutron emissivity\n density profile [cm$^{-3}$ s$^{-1}$]', fontsize='12')
#plt.savefig('JET_LimEq_NeutDensTot_Comp_Shot{0}_T3.png' .format(RunID), dpi=600)
plt.show()

# Computed Total neutron emission density profile normalized
NeutDensTot_CompNorm_Int=scinterpolate.griddata((R2DM, Z2DM), NeutDensTot_CompNorm, (R[None,:], Z[:,None]), method='cubic')
plt.figure()
plt.gca().set_aspect('equal')
plt.tricontourf(triang, NeutDensTot_CompNorm, 100, linewidths=3, cmap=plt.cm.plasma)
#plt.colorbar()
#plt.contour(R,Z,NeutDensTot_CompNorm_Int, 8, linewidths=0.5, colors='black')
#plt.plot(RLIM, ZLIM, linestyle='solid', color='black')
#plt.plot(RBND, ZBND, linestyle='solid', color='red')
#plt.plot(RMAG, ZMAG, marker='x', mew=1.5, ms=6, color='red')
plt.ylabel('JET Z-axis [m]')
plt.xlabel('JET R-major radius [m]')
plt.title('TRANSP normalized emissivity\n profile [a.u.]' .format(RunID), fontsize='16')
#plt.savefig('JET_LimEq_NeutDensTot_Comp_Shot{0}_T3.png' .format(RunID), dpi=600)
plt.show()


# Computed Total neutron emission density profile for Zamir
NeutTotF2D_Int=scinterpolate.griddata((R2DM, Z2DM), NeutTotF2D, (R[None,:], Z[:,None]), method='cubic')

#NeutDensTot_CompNorm_Int=scinterpolate.griddata((R2DM, Z2DM), NeutDensTot_CompNorm, (R[None,:], Z[:,None]), method='cubic')
plt.figure()                                                                                        
plt.gca().set_aspect('equal')                                                                       
#plt.tricontourf(triang, NeutTotF2D_Int, 100, linewidths=3, cmap=plt.cm.plasma)
plt.tricontourf(triang, NeutTotF2D, 100, linewidths=3, cmap=plt.cm.plasma)                          
plt.colorbar()                                                                                      
plt.contour(R,Z,NeutTotF2D_Int, 10, linewidths=0.5, colors='black')  
#plt.colorbar()                                                                                     
#plt.contour(R,Z,NeutDensTot_CompNorm_Int, 8, linewidths=0.5, colors='black')                       
#plt.plot(RLIM, ZLIM, linestyle='solid', color='black')                                             
#plt.plot(RBND, ZBND, linestyle='solid', color='red')                                               
#plt.plot(RMAG, ZMAG, marker='x', mew=1.5, ms=6, color='red')                                       
plt.ylabel('JET Z-axis [m]')                                                                        
plt.xlabel('JET R-major radius [m]')                                                                
plt.title(f'{RunID} TOTNTNF2D [#/SEC/CM**3] emissivity\n profile ', fontsize='16')            
#plt.savefig('JET_LimEq_NeutDensTot_Comp_Shot{0}_T3.png' .format(RunID), dpi=600)                   
plt.show()   

# BT neutron emission density relative difference
NeutDensBT_RelDif_RZInt=scinterpolate.griddata((R2DM, Z2DM), NeutDensBT_RelDif, (R[None,:], Z[:,None]), method='cubic')
plt.figure()
plt.gca().set_aspect('equal')
plt.imshow(NeutDensBT_RelDif_RZInt, vmin=-40, vmax=40, origin='lower', extent=[R.min(), R.max(), Z.min(), Z.max()], cmap=plt.cm.bwr)
plt.colorbar()
plt.contour(R,Z,NeutDensTot_Comp_Int, 5, linewidths=1.5, colors='black')
#plt.plot(RBND, ZBND, linestyle='solid', color='red')
plt.ylabel('JET Z-axis [m]')
plt.xlabel('JET R-major radius [m]')
plt.title('Total vs. thermal fusion neutron emssion\n Shot {0}' .format(RunID), fontsize='12')
#plt.savefig('JET_LimEq_NeutDensBT_RelDif_Shot{0}_T3_Volume.png' .format(RunID), dpi=600)
plt.show()

# Plotting the neutron emission data - comparison between measured and TRANSP computed values

f, (ax1, ax2) = plt.subplots(2, sharex=True,gridspec_kw={'height_ratios':[2,1]})

ax1.plot(BasicTime+40.0, TotalMeasNeutrons, 'k-', label='Measured NR')
ax1.plot(BasicTime+40.0, TotalNeutrons*1.15, 'r-', label='Computed NR')
ax1.plot(BasicTime+40.0, TotalTHNeutrons*1.15, 'g:',linewidth=2.0, label='Thermal NR')
ax1.plot(BasicTime+40.0, TotalBTNeutrons*1.15, 'b:',linewidth=2.0, label='Beam-target NR')
ax1.plot(BasicTime+40.0, TotalBBNeutrons*1.15, 'm:',linewidth=2.0, label='Beam-beam NR')
ax1.set_ylabel('Neutron rate [$s^{-1}$]')
leg=ax1.legend(fancybox=True, loc='upper left', prop={'size':14})
leg.get_frame().set_alpha(0.5)
ax1.add_patch(patches.Rectangle((48.75, 0.0), 0.5, 4.0E+16,facecolor="gray",alpha=0.6))

ax2.plot(BasicTime+40.0, TotNeutDiff, 'k-', label='NR relative difference')
ax2.axhline(y=0.0, color='k', linestyle='--')
#ax3.set_ylim([0.0, 3.2E+16])
ax2.set_ylabel('Relative difference [%]')
ax2.set_xlabel('Time [s]')
ax2.set_ylim([-22, 45])
ax2.yaxis.set_ticks(np.arange(-20,60,20))
#ax3.legend(fancybox=True,prop={'size':14})
f.subplots_adjust(hspace=0.1)
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
axes = plt.gca()
axes.set_xlim([47.5,51])
plt.show()


f, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True)

# THERMAL neutron emission density profile
NeutDensTH_Int=scinterpolate.griddata((R2DM, Z2DM), NeutDensTH, (R[None,:], Z[:,None]), method='cubic')
#ax1.gca().set_aspect('equal')
con1=ax1.tricontourf(triang, NeutDensTH, 100, linewidths=0, cmap=plt.cm.plasma, antialiased=True)
f.colorbar(con1, ax=ax1, pad=0.01)
ax1.contour(R,Z,NeutDensTH_Int, 10, linewidths=0.5, colors='black')
#ax1.plot(RLIM, ZLIM, linestyle='solid', color='black')
#ax1.plot(RBND, ZBND, linestyle='solid', color='red')
ax1.plot(3.066,0.277, marker='x', mew=1.5, ms=6, color='white')
ax1.set_aspect('equal')
#ax1.set_axis_off()
ax1.set_ylabel('JET Z-axis [m]')
ax1.set_xlabel('JET R-major radius [m]')
ax1.set_title('Thermal neutron emission\n density [cm$^{-3}$ s$^{-1}$]', fontsize='12')

# BT neutron emission density profile
NeutDensBT_Int=scinterpolate.griddata((R2DM, Z2DM), NeutDensBT, (R[None,:], Z[:,None]), method='cubic')
#ax1.gca().set_aspect('equal')
con2=ax2.tricontourf(triang, NeutDensBT, 100, linewidths=0, cmap=plt.cm.plasma)
f.colorbar(con2, ax=ax2, pad=0.01)
ax2.contour(R,Z,NeutDensBT_Int, 10, linewidths=0.5, colors='black')
#ax2.plot(RLIM, ZLIM, linestyle='solid', color='black')
#ax2.plot(RBND, ZBND, linestyle='solid', color='red')
ax2.plot(3.066,0.277, marker='x', mew=1.5, ms=6, color='white')
ax2.set_aspect('equal')
#ax1.set_axis_off()
#ax2.set_ylabel('JET Z-axis [m]')
ax2.set_xlabel('JET R-major radius [m]')
ax2.set_title('BT neutron emission\n density [cm$^{-3}$ s$^{-1}$]', fontsize='12')

# Total neutron emission density relative difference
NeutDensBB_Int=scinterpolate.griddata((R2DM, Z2DM), NeutDensBB, (R[None,:], Z[:,None]), method='cubic')
#ax3.gca().set_aspect('equal')
con3=ax3.tricontourf(triang, NeutDensBB, 100, linewidths=0, cmap=plt.cm.plasma)
f.colorbar(con3, ax=ax3, pad=0.01)
ax3.contour(R,Z,NeutDensBB_Int, 10, linewidths=0.5, colors='black')
#ax3.plot(RLIM, ZLIM, linestyle='solid', color='black')
#ax3.plot(RBND, ZBND, linestyle='solid', color='red')
ax3.plot(3.066,0.277, marker='x', mew=1.5, ms=6, color='white')
ax3.set_aspect('equal')
#ax3.set_axis_off()
#ax3.set_ylabel('JET Z-axis [m]')
ax3.set_xlabel('JET R-major radius [m]')
ax3.set_title('BB neutron emission\n density [cm$^{-3}$ s$^{-1}$]', fontsize='12')

axes = plt.gca()
axes.set_ylim([-1.75,2.0])
#axes.set_xlim([1.75,4.0])

# AVERAGED Total neutron emission density profile
NeutTotF2D_Int=scinterpolate.griddata((R2DM, Z2DM), NeutTotF2D, (R[None,:], Z[:,None]), method='cubic')
plt.figure()
plt.gca().set_aspect('equal')
plt.tricontourf(triang, NeutTotF2D, 100, linewidths=3, cmap=plt.cm.plasma)
plt.colorbar()
plt.contour(R,Z,NeutTotF2D_Int, 10, linewidths=0.5, colors='black')
#plt.plot(RLIM, ZLIM, linestyle='solid', color='black')
#plt.plot(RBND, ZBND, linestyle='solid', color='red')
plt.plot(RMAG, ZMAG, marker='x', mew=1.5, ms=6, color='red')
plt.ylabel('JET Z-axis [m]')
plt.xlabel('JET R-major radius [m]')
plt.xticks(np.arange(2.0, 4.0, 0.5))
plt.title('Neutron emissivity profile\n averaged [cm$^{-3}$ s$^{-1}$]', fontsize='16')
#plt.savefig('JET_LimEq_NeutDensTot_Comp_Shot{0}_T3.png' .format(RunID), dpi=600)
plt.show()

NeutDensTot_RelDif=(NeutDensTot_Comp/NeutTotF2D-1)*100

# RELATIVE DIFFERENCE SAMPLING AVERAGE OR NUBEAM
NeutDensTot_RelDif_RZInt=scinterpolate.griddata((R2DM, Z2DM), NeutDensTot_RelDif, (R[None,:], Z[:,None]), method='cubic')
plt.figure()
plt.gca().set_aspect('equal')
plt.imshow(NeutDensTot_RelDif_RZInt, vmin=NeutDensTot_RelDif.min(), vmax=NeutDensTot_RelDif.max(), origin='lower', extent=[R.min(), R.max(), Z.min(), Z.max()], cmap=plt.cm.rainbow)
plt.colorbar
#plt.contour(R,Z,NeutTotF2D_Int, 10, linewidths=0.5, colors='black')
#plt.plot(RLIM, ZLIM, linestyle='solid', color='black')
#plt.plot(RBND, ZBND, linestyle='solid', color='red')
#plt.plot(RMAG, ZMAG, marker='x', mew=1.5, ms=6, color='red')
plt.ylabel('JET Z-axis [m]')
plt.xlabel('JET R-major radius [m]')
plt.xticks(np.arange(2.0, 4.0, 0.5))
plt.title('Relative difference [%]', fontsize='16')
#plt.savefig('JET_LimEq_NeutDensTot_Comp_Shot{0}_T3.png' .format(RunID), dpi=600)
plt.show()

TRANSPDat.close()
LimiterDat.close()
NeutronDat.close()
