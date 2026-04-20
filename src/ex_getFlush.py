import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
import py_flush as Flush
import ppf
import pdb
import warnings

class EquilibriumNotFoundError(RuntimeError):
    """Raised when EFIT equilibrium is not found in the database."""
    pass

class SignalNotFoundWarning(RuntimeWarning):
    """Warning raised when optional signal is missing."""
    pass



def ReadSignal(pulse, seq, user, dda, dtyp, fatal=False):
    # Initialise PPF routines
    ier = ppf.ppfgo(pulse, seq)
    if ier != 0:
        raise RuntimeError(
            f"Error initialising PPF routines for pulse {pulse}/{seq}"
        )

    # Set User ID for reading
    ier = ppf.ppfuid(user, rw="R")

    ihdat, iwdat, data, x, t, ier = ppf.ppfget(pulse, dda, dtyp)
    if ier != 0:
        msg = (
            f"Signal {pulse}/{dda}/{seq}/{dtyp} not found "
            f"(PPF error code {ier})"
        )
        if fatal:
            raise RuntimeError(msg)
        else:
            warnings.warn(msg, SignalNotFoundWarning)
    units = ihdat[0:7]
    xunits = ihdat[8:15]
    tunits = ihdat[16:23]

    return x, t, data, units, xunits, tunits, ier




def create_efit_psi(pulse,time):
    dda = 'efit'
    dty = 'jetppf'
    seq = 0
    #grid size in R direction
    #nw = 33
    nw = 65
    #grid size in z direction
    #nh = 33
    nh=65
    # Flush works in cm
    psir = 100*np.linspace(1.65,4.05,nw)
    psiz = 100*np.linspace(-1.9,2.15,nh)
    rv, zv = np.meshgrid(psir,psiz)
    
    npts = nw*nh
    time, ier = Flush.flushquickinit(pulse,time)
    
    f, ier = Flush.Flush_getFlux(npts,rv,zv)
    #reshape
    #f = np.reshape(f,(33,33)) #check
    return f, psir, psiz, ier

def normalize_psi_to_axis(psi, R, Z, Rmag, Zmag):
    """
    Subtract psi(Rmag, Zmag) so that psi=0 on magnetic axis.
    psi shape: (nZ, nR)
    """

    interp = RegularGridInterpolator(
        (Z, R), psi,
        bounds_error=False,
        fill_value=None
    )

    psi_axis = interp((Zmag, Rmag))
    return psi - psi_axis

def plot_psi(pulse, psi_flush, psir_flush, psiz_flush,nl=30):                                                                        
    """                                                                                             
    Plots equilibrium contour plot for a given time slice.                                          
    nl: number of contours                                                                          
    """                                                                                             
    R_grid = psir_flush                                                                               
    Z_grid = psiz_flush                                                                               
    ratio = (Z_grid[-1]-Z_grid[0])/(R_grid[-1]-R_grid[0])                                           
    # figure size                                                                                   
    fs=(6,6*ratio)                                                                                  
                                                                                                    
    nR = len(psir_flush)                                                                              
    nZ = len(psiz_flush)                                                                              
    psi_grid = np.reshape(psi_flush, (nR,nZ))                                                   
    R,Z = np.meshgrid(R_grid,Z_grid)                                                                
    #creating figure                                                                                
    fig, ax = plt.subplots(figsize=fs)                                                              
    levels = np.linspace(np.min(psi_grid),np.max(psi_grid),nl)                                      
    CS = ax.contour(R,Z,psi_grid,levels)                                                            
    ax.clabel(CS, fontsize=10)                                                                      
    #ax.plot(eq.Rmag[tind],eq.Zmag[tind],'ro',ms=2)                                                  
    #ax.plot(R,Z, 'ko',ms=1) #plots grid points on top                                              
    ax.set_xlabel('R[m]')                                                                           
    ax.set_ylabel('z[m]')                                                                           
    ax.set_title(f'EFIT pulse {pulse}')                                                                                                                                                                         
    plt.show()
def plot_psi_compare(
    psi_flush, psir_flush, psiz_flush,
    psi_efit,  psir_efit,  psiz_efit,
    nl=30
    ):
    """
    Side-by-side contour plots of FLUSH and EFIT psi.
    """

    # ---- reshape ----
    nR_f = len(psir_flush)
    nZ_f = len(psiz_flush)
    psi_flush = psi_flush.reshape(nZ_f, nR_f)

    nR_e = len(psir_efit)
    nZ_e = len(psiz_efit)
    psi_efit = psi_efit.reshape(nZ_e, nR_e)

    # ---- grids ----
    Rf, Zf = np.meshgrid(psir_flush, psiz_flush)
    Re, Ze = np.meshgrid(psir_efit, psiz_efit)

    # ---- common contour levels (important for comparison) ----
    psi_min = max(psi_flush.min(), psi_efit.min())
    psi_max = min(psi_flush.max(), psi_efit.max())
    levels = np.linspace(psi_min, psi_max, nl)

    # ---- figure ----
    fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

    # FLUSH
    cs0 = axs[0].contour(Rf, Zf, psi_flush, levels=levels)
    axs[0].clabel(cs0, fontsize=8)
    axs[0].set_title("FLUSH ψ")
    axs[0].set_xlabel("R [m]")
    axs[0].set_ylabel("Z [m]")

    # EFIT
    cs1 = axs[1].contour(Re, Ze, psi_efit, levels=levels)
    axs[1].clabel(cs1, fontsize=8)
    axs[1].set_title("EFIT ψ")
    axs[1].set_xlabel("R [m]")

    plt.tight_layout()
    plt.show()

def main():
    pulse = 105630
    seq = 0
    user = 'jetppf'
    time = 46.0
    dda = 'efit'

    # break if equilibrium not found
    _, efit_times, *__, ier = ReadSignal(pulse, seq, user, dda, 'Q',fatal=True)
    tind = np.abs(efit_times-time).argmin()                                                         
    time_efit = efit_times[tind]
    # if reading Q passes the Rmag and Zmag should be available as well
    _, _, Rmag, _, _, ier = ppf.ppfget(pulse, dda, 'RMAG',reshape=True)                         
    _, _, Zmag, _, _, ier = ppf.ppfget(pulse, dda, 'ZMAG',reshape=True)
    Rmag_t = Rmag[tind]
    Zmag_t = Zmag[tind]

    # PSI may not be present in older efit runs
    ihdat, iwdat, efit_psi, x, efit_times, ier = ppf.ppfget(pulse, dda, 'PSI',reshape=True)
    if ier==0:
        ifpsi=True
        _, _, efit_psir, _, _, ier = ppf.ppfget(pulse, dda, 'PSIR',reshape=True)
        _, _, efit_psiz, _, _, ier = ppf.ppfget(pulse, dda, 'PSIZ',reshape=True)
        
        print("EFIT ier status  :", ier)
        print("EFIT PSI units   :", ihdat[0:7])
        print("EFIT R units     :", ihdat[8:15])
        print("EFIT time units  :", ihdat[16:23])
        #if ier != 0:
        #    raise Exception('Error reading signal ' + DDA + '/PSIRZ')
    else:
        print(f'{pulse}/{dda}/{seq}/PSI not found. ')
        ifpsi = False
    
    #time, ier = Flush.flushquickinit(pulse,time_efit)

    #psi = []
    #for time in efit_times[tind-1:tind+1]:
    #    print(f'Time: {time:.3f}')
    #    add_psi, _ = create_efit_psi(pulse,time)
    #    psi.append(add_psi)
    #psi = np.column_stack(psi)
    
    # create example equilibrium
    psi_flush, psir, psiz, _ = create_efit_psi(pulse,time_efit)
    # convert FLUSH grid from cm → m
    psir = psir / 100.0
    psiz = psiz / 100.0

    # Normalize FLUX

    psi_flush_grid = psi_flush.reshape(len(psiz), len(psir))

    psi_flush_norm = normalize_psi_to_axis(
        psi_flush_grid,
        psir,
        psiz,
        Rmag_t,
        Zmag_t
    )


    # EFIT PSI if available
    if ifpsi:
        psi_efit = efit_psi[tind]   
   
        # Normalize EFIT
        psi_efit_grid = psi_efit.reshape(len(efit_psiz), len(efit_psir))
   
        psi_efit_norm = normalize_psi_to_axis(
            psi_efit_grid,
            efit_psir,
            efit_psiz,
            Rmag_t,
            Zmag_t
        )

    #print(f'ier: {ier}')
    #print(f'flush PSI shape {np.shape(psi_flush)}')
    #print(f'efit PSI shape {np.shape(psi_efit)}') 
    
    #print("FLUSH R range:", psir.min(), psir.max())
    #print("EFIT  R range:", efit_psir.min(), efit_psir.max())

    #print("FLUSH Z range:", psiz.min(), psiz.max())
    #print("EFIT  Z range:", efit_psiz.min(), efit_psiz.max())
    if ifpsi:
        print("ψ_axis EFIT :", psi_efit_grid.min())
    print("ψ_axis FLUSH:", psi_flush_grid.min())

    print("After normalization:")
    if ifpsi:
        print("ψ_axis EFIT :", psi_efit_norm.min())
    print("ψ_axis FLUSH:", psi_flush_norm.min())
    if ifpsi:
        plot_psi_compare(
            psi_flush_norm, psir, psiz,
            psi_efit_norm, efit_psir, efit_psiz,
            nl=30
        )
    else:
        plot_psi(pulse, psi_flush, psir, psiz,nl=30)
main()
