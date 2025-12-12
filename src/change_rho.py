import numpy as np
from scipy.interpolate import RegularGridInterpolator, griddata
from scipy.interpolate import PchipInterpolator
from scipy.interpolate import CubicSpline

import profiles as ps                                                                               
import sys                                                                                          
import matplotlib.pyplot as plt                                                                                                    
import argparse                                                                                     


def parse_arguments():                                                                                                    
    parser = argparse.ArgumentParser()                                                                  
    parser.add_argument('pulse', type=int)                                                              
    #parser.add_argument('runid', type=str) 
    parser.add_argument('-uid', '--uid', help='PPF user ID', type=str,default='jetppf')
    parser.add_argument('-dda', '--dda', type=str, default='eftp', help='default= eftp , PPF dda (Diagnostic Data Area) i.e. EFIT,EFTP,EFTM,...') 
    #parser.add_argument('-dty', '--datatype', help='PPF datatype') 
    parser.add_argument('-seq', '--sequence', type=int, default=0, help='PPF sequence, default=0') 
    parser.add_argument('-runid', '--runid', type=str, default=None, help="TRANSP runID")
    parser.add_argument('-prof', '--profile', type=str, default='TI', help="TRANSP output (t,x) signal")
    parser.add_argument("--save", action="store_true", help="Save profile data to text file")
    return parser.parse_args()

def plot_psi(eq,time,nl=30):
    """
    Plots equilibrium contour plot for a given time slice. 
    nl: number of contours
    """
    R_grid = eq._psir
    Z_grid = eq._psiz
    ratio = (Z_grid[-1]-Z_grid[0])/(R_grid[-1]-R_grid[0])
    # figure size
    fs=(6,6*ratio)

    nR = len(eq._psir)
    nZ = len(eq._psiz)
    tind = np.abs(eq.t-time).argmin()
    psi_grid = np.reshape(eq._psi[tind], (nR,nZ))
    R,Z = np.meshgrid(R_grid,Z_grid)
    #creating figure
    fig, ax = plt.subplots(figsize=fs)
    levels = np.linspace(np.min(psi_grid),np.max(psi_grid),nl)
    CS = ax.contour(R,Z,psi_grid,levels)
    ax.clabel(CS, fontsize=10)
    ax.plot(eq.Rmag[tind],eq.Zmag[tind],'ro',ms=2)
    #ax.plot(R,Z, 'ko',ms=1) #plots grid points on top
    ax.set_xlabel('R[m]')
    ax.set_ylabel('z[m]')
    ax.set_title(f'EFIT pulse {eq._pulse}')

    plt.show()

def RZ_to_psi(R_pts, Z_pts, eq, time,
    method='linear' # interpolation method: 'linear' or 'nearest' or 'cubic' (griddata)
    ):
    """
    Map points (R_pts, Z_pts) onto poloidal flux.
    
    Inputs (choose either grid version OR scattered version):
      - Grid version (preferred for EFIT output):
        psi_grid : 2D array shape (nZ, nR) or (nR, nZ) depending on your mesh
                   given by eq._psi
        R_grid : eq._psir 1D array of R coordinates (length nR)
        Z_grid : eq._psiz 1D array of Z coordinates (length nZ)

      - Scattered version:
        psi_scattered : 1D array of psi values at scattered (R,Z) points
                        then you must pass a matching points array below
      - psi_axis : scalar, flux on magnetic axis
      - psi_bnd  : scalar, flux on LCFS / boundary

      - R_pts, Z_pts : 1D arrays of points to map (same length)

      - method : interpolation method for RegularGridInterpolator or griddata
    
    Returns:
      s : 1D array of sqrt(normalized_flux) in [0,1] (or NaN for failed)
      norm : normalized flux before sqrt (clipped if clip=True)
      status : dict with masks e.g. 'inside' boolean mask where norm<=1 and >=0,
               'outside' boolean mask, 'nan' boolean mask.
    """
    
    R_pts = np.asarray(R_pts).ravel()
    Z_pts = np.asarray(Z_pts).ravel()
    assert R_pts.shape == Z_pts.shape, "R_pts and Z_pts must have same shape"
    n = R_pts.size
    
    # psi_grid dimensions
    
    R_grid = eq._psir
    Z_grid = eq._psiz

    nR = len(eq._psir)
    nZ = len(eq._psiz)
    tind = np.abs(eq.t-time).argmin()

    # reshaping, original eq._psi shape is (len(time), len(psir)*len(psiz))  
    # rechaping is in (R,Z) order
    psi_grid = np.reshape(eq._psi[tind], (nR,nZ))

    psi_scattered=None

    psi_axis = eq._faxs[tind]
    psi_bnd = eq._fbnd[tind]
    #print(f'PSI on AXIS = {psi_axis}')
    #print(f'PSI at BND = {psi_bnd}') 

    # interpolate psi at points
    if psi_grid is not None and R_grid is not None and Z_grid is not None:
        # Ensure R_grid and Z_grid are ascending
        R_grid = np.asarray(R_grid)
        Z_grid = np.asarray(Z_grid)
        psi_grid = np.asarray(psi_grid)
        #psi_gird = psi_grid.T

        # Detect psi_grid orientation: common EFIT provides psi as psi(z_index, r_index)
        # We will assume psi_grid.shape == (len(Z_grid), len(R_grid)).
        if psi_grid.shape != (len(R_grid), len(Z_grid)):
            # try transposing as fallback
            if psi_grid.T.shape == (len(R_grid), len(Z_grid)):
                psi_grid = psi_grid.T
            else:
                # do not raise, but try best effort: reshape if possible
                pass

        # Make sure grid axes are strictly ascending
        if (np.any(np.diff(R_grid) <= 0) or np.any(np.diff(Z_grid) <= 0)):
            print('Sorting (R,Z) points')
            # sort them
            R_sort = np.argsort(R_grid)
            Z_sort = np.argsort(Z_grid)
            R_grid = R_grid[R_sort]
            Z_grid = Z_grid[Z_sort]
            psi_grid = psi_grid[np.ix_(R_sort, Z_sort)]

        interp = RegularGridInterpolator((R_grid, Z_grid), psi_grid.T,
                                         method=method, bounds_error=False, fill_value=np.nan)
        #interp = Rbf((R_grid, Z_grid), psi_grid.T, epsilon=3)
        query = np.column_stack((R_pts, Z_pts))  # note the ordering (Z,R) because grid uses Z first
        print('Interpolating on (rpts,zpts)') 
        psi_at_pts = interp(query)

    elif psi_scattered is not None:
        # psi_scattered must come with scattered_points (R_scat, Z_scat)
        # We'll assume user provides psi_scattered as tuple (points, values) where points = (R_scat, Z_scat)
        # For generality allow psi_scattered be dict-like
        try:
            points, values = psi_scattered
            pts_array = np.column_stack(points)
        except Exception:
            raise ValueError("For scattered interpolation pass psi_scattered=(points, values) where points is (R_scat, Z_scat) arrays.")
        query = np.column_stack((R_pts, Z_pts))
        psi_at_pts = griddata(pts_array, values, query, method=method, fill_value=np.nan)
    else:
        raise ValueError("Provide either (psi_grid, R_grid, Z_grid) or psi_scattered=(points,values).")

    # include psi value on axis
    Rmag = eq._Rmag[tind]
    Zmag = eq._Zmag[tind]
    psi_axis = eq._faxs[tind]

    # distance threshold (you can adjust 1e-3 or 5e-3)
    axis_tol = 5e-3  

    d_axis = np.sqrt((R_pts - Rmag)**2 + (Z_pts - Zmag)**2)

    # Replace psi values close to the magnetic axis with exact psi_axis
    axis_mask = d_axis < axis_tol
    psi_at_pts[axis_mask] = psi_axis
    
    return psi_at_pts

def RZ_to_psin(R_pts, Z_pts, eq, time,
    method='linear' # interpolation method: 'linear' or 'nearest' or 'cubic' (griddata)
    ):
    """
    Map points (R_pts, Z_pts) onto poloidal flux.
    
    Inputs (choose either grid version OR scattered version):
      - Grid version (preferred for EFIT output):
        psi_grid : 2D array shape (nZ, nR) or (nR, nZ) depending on your mesh
                   given by eq._psi
        R_grid : eq._psir 1D array of R coordinates (length nR)
        Z_grid : eq._psiz 1D array of Z coordinates (length nZ)

      - Scattered version:
        psi_scattered : 1D array of psi values at scattered (R,Z) points
                        then you must pass a matching points array below
      - psi_axis : scalar, flux on magnetic axis
      - psi_bnd  : scalar, flux on LCFS / boundary

      - R_pts, Z_pts : 1D arrays of points to map (same length)

      - method : interpolation method for RegularGridInterpolator or griddata
    
    Returns:
      s : 1D array of sqrt(normalized_flux) in [0,1] (or NaN for failed)
      norm : normalized flux before sqrt (clipped if clip=True)
      status : dict with masks e.g. 'inside' boolean mask where norm<=1 and >=0,
               'outside' boolean mask, 'nan' boolean mask.
    """
    
    R_pts = np.asarray(R_pts).ravel()
    Z_pts = np.asarray(Z_pts).ravel()
    assert R_pts.shape == Z_pts.shape, "R_pts and Z_pts must have same shape"
    n = R_pts.size
 
    # psi_grid dimensions
    
    R_grid = eq._psir
    Z_grid = eq._psiz
    nR = len(R_grid)
    nZ = len(Z_grid)

    tind = np.abs(eq.t-time).argmin()


    

    # reshaping, original eq._psi_norm shape is (len(time), len(psir)*len(psiz))  
    # rechaping is in (R,Z) order

    psin_grid = None
    if method=='cubic':
        RZ = np.array([R_pts,Z_pts])
        # shape check
        RZ = RZ.T if RZ.shape[0] == 2 else RZ
        # EFIT grids
        Rg = eq._psirzmg[0]
        Zg = eq._psirzmg[1]
        psin = eq._psi_norm[tind]

        Rf = Rg.ravel()
        Zf = Zg.ravel()
        psinf = psin.ravel()

        # --- Add magnetic axis point ----------------------------------------
        Rmag = eq.Rmag[tind]
        Zmag = eq.Zmag[tind]
        rho_axis = 0.0  # by definition sqrt(normalized ψ)
    
        Rf2 = np.concatenate([Rf, [Rmag]])
        Zf2 = np.concatenate([Zf, [Zmag]])
        psinf2 = np.concatenate([psinf, [rho_axis]])
        # ---------------------------------------------------------------------

        # Primary interpolation
        psin_at_pts = griddata(
            (Rf2, Zf2),
            psinf2,
            RZ,
            method='cubic'
        )

        # Fallback for NaNs
        if np.any(np.isnan(psin_at_pts)):
            psin_at_pts = griddata(
                (Rf2, Zf2),
                psinf2,
                RZ,
                method='linear'
            )
            psin_at_pts = np.where(np.isnan(rho_interp), rho_interp2, rho_interp)
    else:
        psin_grid = np.reshape(eq._psi_norm[tind], (nR,nZ))


    # interpolate psi at points
    if psin_grid is not None and R_grid is not None and Z_grid is not None:
        # Ensure R_grid and Z_grid are ascending
        R_grid = np.asarray(R_grid)
        Z_grid = np.asarray(Z_grid)
        psin_grid = np.asarray(psin_grid)
        #psi_gird = psi_grid.T

        # Detect psi_grid orientation: common EFIT provides psi as psi(z_index, r_index)
        # We will assume psi_grid.shape == (len(Z_grid), len(R_grid)).
        if psin_grid.shape != (len(R_grid), len(Z_grid)):
            # try transposing as fallback
            if psin_grid.T.shape == (len(R_grid), len(Z_grid)):
                psin_grid = psin_grid.T
            else:
                # do not raise, but try best effort: reshape if possible
                pass

        # Make sure grid axes are strictly ascending
        if (np.any(np.diff(R_grid) <= 0) or np.any(np.diff(Z_grid) <= 0)):
            print('Sorting (R,Z) points')
            # sort them
            R_sort = np.argsort(R_grid)
            Z_sort = np.argsort(Z_grid)
            R_grid = R_grid[R_sort]
            Z_grid = Z_grid[Z_sort]
            psin_grid = psin_grid[np.ix_(R_sort, Z_sort)]

        interp = RegularGridInterpolator((R_grid, Z_grid), psin_grid.T,
                                         method=method, bounds_error=False, fill_value=np.nan)
        #interp = Rbf((R_grid, Z_grid), psi_grid.T, epsilon=3)
        query = np.column_stack((R_pts, Z_pts))  # note the ordering (Z,R) because grid uses Z first
        print('Interpolating on (rpts,zpts)') 
        psin_at_pts = interp(query)


        # include psi value on axis
        Rmag = eq._Rmag[tind]
        Zmag = eq._Zmag[tind]
        psin_axis = 0.0

        # distance threshold (you can adjust 1e-3 or 5e-3)
        axis_tol = 5e-3  

        d_axis = np.sqrt((R_pts - Rmag)**2 + (Z_pts - Zmag)**2)

        # Replace psi values close to the magnetic axis with exact psi_axis
        axis_mask = d_axis < axis_tol
        psin_at_pts[axis_mask] = psin_axis
    
    return psin_at_pts


def RZ_to_psin_2(RZ,eq,t1):
    # shape check
    RZ = RZ.T if RZ.shape[0] == 2 else RZ

    ix_t1 = (np.abs(eq.t - t1)).argmin()

    # EFIT grids
    Rg = eq._psirzmg[0]
    Zg = eq._psirzmg[1]
    rho = eq._psi_norm[ix_t1]

    Rf = Rg.ravel()
    Zf = Zg.ravel()
    rhof = rho.ravel()

    # --- Add magnetic axis point ----------------------------------------
    Rmag = eq.Rmag[ix_t1]
    Zmag = eq.Zmag[ix_t1]
    rho_axis = 0.0  # by definition sqrt(normalized ψ)
    
    Rf2 = np.concatenate([Rf, [Rmag]])
    Zf2 = np.concatenate([Zf, [Zmag]])
    rhof2 = np.concatenate([rhof, [rho_axis]])
    # ---------------------------------------------------------------------

    # Primary interpolation
    rho_interp = griddata(
        (Rf2, Zf2),
        rhof2,
        RZ,
        method='cubic'
    )

    # Fallback for NaNs
    if np.any(np.isnan(rho_interp)):
        rho_interp2 = griddata(
            (Rf2, Zf2),
            rhof2,
            RZ,
            method='linear'
        )
        rho_interp = np.where(np.isnan(rho_interp), rho_interp2, rho_interp)

    return rho_interp


def psi_to_psin(psi_at_pts, eq, time, clip=True, outside=1.3, sqrt=False):
    """
    Maps poloidal flux onto s = sqrt(normalized_poloidal_flux) at a chosen time slice 'time'

    From eq (Eq class instance) uses:
      eq._faxs - 'FAXS' psi at magn. axis (from efit ppf)
      eq._fbnd - 'FBND' psi at plasma boundary (from efit ppf)
      
    Returns:
      s : 1D array of sqrt(normalized_flux) in [0,1.3] (or NaN for failed)
      norm : normalized flux before sqrt (clipped if clip=True)
      status : dict with masks e.g. 'inside' boolean mask where norm<=1 and >=0,
               'outside' boolean mask, 'nan' boolean mask.
    """
    tind = np.abs(eq.t-time).argmin()
    psi_axis = eq._faxs[tind]
    psi_bnd = eq._fbnd[tind]
    psi_axis = float(psi_axis)
    psi_bnd = float(psi_bnd)
    denom = (psi_bnd - psi_axis)

    # Check psi_axis and psi_bnd
    if psi_axis is None or psi_bnd is None:
        raise ValueError("psi_axis and psi_bnd (flux on axis and on boundary) must be provided.")
    if denom == 0:
        raise ValueError("psi_bnd equals psi_axis; cannot normalize (denominator zero).")

    # normalized flux (rho^2 like)
    norm_raw = (psi_at_pts - psi_axis) / denom


    # masked status (0,1.3), range as in ProfileMaker
    nan_mask = np.isnan(norm_raw)
    outside_mask = norm_raw > outside
    negative_mask = norm_raw < 0.0
    inside_mask = (~nan_mask) & (~outside_mask) & (~negative_mask)

    # optionally clip to [0,1.3] before sqrt
    if clip:
        norm = np.clip(norm_raw, 0.0, outside)
    else:
        norm = norm_raw.copy()

    # sqrt=true  is default
    with np.errstate(invalid='ignore'):
        if sqrt:
            s = np.sqrt(norm)
        else: 
            s = norm

    status = {'inside': inside_mask, 'outside': outside_mask, 'negative': negative_mask, 'nan': nan_mask}
    
    #return psi_at_pts, norm_raw, status, psi_grid
    return s, norm_raw, status

def RZ_to_psin_3(Rpts, Zpts, eq, time, method='linear'):
    """
    (R,Z) points mapps onto rhop
    rhop: psi normalized
    """
    psi_at_pts= RZ_to_psi(Rpts, Zpts, eq, time, method=method)
    s, norm_raw, status = psi_to_psin(psi_at_pts, eq, time, sqrt=False)
    return s, norm_raw, status, psi_at_pts

def RZ_to_sqrt_psin(Rpts, Zpts, eq, time, method='linear'):
    """
    (R,Z) points mapped onto sqrt of psi normalized
    Wrapper mapping Rpts, Zpts (equal length) to sqrt of psi_norm
    rhop = sqrt_psin_norm
    """
    psi_at_pts= RZ_to_psi(Rpts, Zpts, eq, time, method=method) 

    s, norm_raw, status = psi_to_psin(psi_at_pts, eq, time, sqrt=True)

    return s, norm_raw, status, psi_at_pts    


def psin_to_ftor(pts_psin, eq, time):
    """
    ftor: 'FTOR' - Toroidal Flux(Psin) [Wb] (not normalized) 
    ftor_x: Psin
    """ 
    tind = np.abs(eq.t-time).argmin()
    ftor = eq._ftor
    ftor_x = eq._ftor_x
    f_psin_to_ftor = PchipInterpolator(ftor_x[tind], ftor[tind], extrapolate=True)
    return f_psin_to_ftor(pts_psin)
    
def psin_to_ftor_norm(pts_psin, eq, time):
    """
    ftor: 'FTOR' - Toroidal Flux(Psin) [Wb] (not normalized) 
    ftor_norm = ftor[time,:]/ftor[time,-1]
    ftor_x: Psin
    """ 
    tind = np.abs(eq.t-time).argmin()
    ftor_norm = eq._ftor_norm
    psin = eq._ftor_x
    f_psin_to_ftor = PchipInterpolator(spsin[tind], ftor_norm[tind], extrapolate=True)
    return f_psin_to_ftor(pts_psin)

def psin_to_sqrt_ftor_norm(pts_psin, eq, time):
    """
    ftor: 'FTOR' - Toroidal Flux(Psin) [Wb] (not normalized) 
    sqrt_ftor_norm = np.sqrt(ftor[time,:]/ftor[time,-1])
    ftor_x: Psin
    """ 
    tind = np.abs(eq.t-time).argmin()
    sqrt_ftor_norm = eq._sqrt_ftor_norm
    psin = eq._ftor_x
    # PchipInterpolater(x,y)
    f = PchipInterpolator(psin[tind], sqrt_ftor_norm[tind], extrapolate=True)
    return f(pts_psin)

def sqrt_ftor_norm_to_psin(pts_rhot, eq, time):
    """
    Inverse function to psin_to_sqrt_ftor_norm
    Maps from sqrt_ftor_norm to ftor_x= psin
    ftor: 'FTOR' - Toroidal Flux(Psin) [Wb] (not normalized) 
    sqrt_ftor_norm = np.sqrt(ftor[time,:]/ftor[time,-1]) i.e. rhot
    ftor_x: Psin
    """
    tind = np.abs(eq.t-time).argmin()
    sqrt_ftor_norm = eq._sqrt_ftor_norm[tind]
    psin = eq._ftor_x 
    f = PchipInterpolator(sqrt_ftor_norm, psin, extrapolate=True)
    return f(pts_rhot) 


def Rbnd_at_midplane(eq, time):
    """
    find Rbnd at midplane on the LFS
    """
    tind = np.abs(eq.t - time).argmin()

    # Midplane line (Z = Zmag)
    Z_mid = eq.Zmag[tind]
    Z_aux = np.full(200, Z_mid)
    R_aux = np.linspace(eq._Rmag[tind], eq._psir[-1], 200)

    # Compute psin along the midplane
    psin,*_ = RZ_to_psin_3(R_aux, Z_aux, eq, time)

    # Ensure psin is strictly increasing: keep only unique increasing values
    # otherwise inverse interpolation is not possible
    psin = np.asarray(psin)
    R_aux = np.asarray(R_aux)

    # Build monotonic mask
    # True if psin[i] > psin[i-1]
    monotonic_mask = np.concatenate(([True], np.diff(psin) > 0))

    psin_m = psin[monotonic_mask]
    R_m = R_aux[monotonic_mask]

    # Remove possible trailing plateau at 1.3 (EFIT's outside value)
    # Keep only psin <= 1.01 + small epsilon
    eps = 1e-5
    inside_mask = psin_m <= (1.01 + eps)
    psin_m = psin_m[inside_mask]
    R_m = R_m[inside_mask]

    # Check if we still have valid points
    if len(psin_m) < 3:
        raise ValueError("Not enough monotonic psin points to interpolate R(psin).")

    # Interpolate R(psin)
    f = PchipInterpolator(psin_m, R_m)

    psin_out = np.linspace(psin_m[0],1.0,100)
    # Return R at psin = 1
    return float(f(1.0))

def psin_to_RZ_midplane(psin_in, eq, time, method='cubic'):
    """
    Input: psin_in
    Maps input psin_in onto RZ midplane (LFS) 
    Maps RZ midplane on the LFS to psin from equilibirum.  
    May find R at the midplane where psin = 1 (i.e., outer midplane boundary).
    """
    tind = np.abs(eq.t - time).argmin()

    # Midplane line (Z = Zmag)
    Z_mid = eq.Zmag[tind]
    Z_aux = np.full(200, Z_mid)
    R_aux = np.linspace(eq._Rmag[tind], eq._psir[-1], 200)

    # Compute psin along the midplane
    psin = RZ_to_psin(R_aux, Z_aux, eq, time, method=method)

    # Ensure psin is strictly increasing: keep only unique increasing values
    # otherwise inverse interpolation is not possible
    psin = np.asarray(psin)
    R_aux = np.asarray(R_aux)

    # Build monotonic mask
    # True if psin[i] > psin[i-1]
    monotonic_mask = np.concatenate(([True], np.diff(psin) > 0))

    psin_m = psin[monotonic_mask]
    R_m = R_aux[monotonic_mask]

    # Remove possible trailing plateau at 1.3 (EFIT's outside value)
    # Keep only psin <= 1.01 + small epsilon
    eps = 1e-5
    inside_mask = psin_m <= (1.01 + eps)
    psin_m = psin_m[inside_mask]
    R_m = R_m[inside_mask]

    # Check if we still have valid points
    if len(psin_m) < 3:
        raise ValueError("Not enough monotonic psin points to interpolate R(psin).")

    # Interpolate R(psin)
    f = PchipInterpolator(psin_m, R_m)
    print(f'psin_m[:2]: {psin_m[:2]}')
    print(f'psin_in[:4]: {psin_in[:4]}')
    psin_in_mask = psin_in >= psin_m[0]
    #psin_out = np.linspace(psin_m[0],1.0,100)
    psin_out = psin_in[psin_in_mask]

    return f(psin_out)

def rhot_to_RZ_midplane(rhot_in, eq, time, norm=True,method='cubic'):
    """
    Input: sqrt_ftor_norm, eq, time slice
    Maps rhot onto RZ midplane (LFS) 
 
    """
    psin_in =  sqrt_ftor_norm_to_psin(rhot_in,eq, time)
    print(f'psin_in.shape: {psin_in.shape}') 
    Rpts = psin_to_RZ_midplane(psin_in, eq, time,method=method)
    Rmag = eq.Rmag[tind]
    if norm:
        Rbnd = Rbnd_at_midplane(eq, time)
        denom = Rbnd-Rmag
        Rpts_out = (Rpts-Rmag)/denom
    else:
        Rpts_out = Rpts

    return Rpts_out

args = parse_arguments()  

pulse = args.pulse 
time=50.9
eq=ps.Eq(pulse,dda = args.dda, uid=args.uid,seq=args.sequence)
dat = ps.Data(pulse, time)
#ex = ps.Exp(dat, eq)



tind = np.abs(eq.t-time).argmin()
print(f'tind = {tind}') 
rpts = np.linspace(2.988,3.82,100)
zpts = eq.Zmag[tind]* np.ones(100)


#X point of 104448 at 50.9s
rpts = [2.62543]
zpts=[-1.457]

# axis
rpts = [eq.Rmag[tind]+2e-3]
zpts = [eq.Zmag[tind]]
# axis to separatrix 
rpts = np.linspace(eq.Rmag[tind],Rbnd_at_midplane(eq, time),200)
zpts = eq.Zmag[tind]*np.ones(rpts.shape)
# sample point
#rpts = [3.888]
#zpts = [2.0618]


#s, norm_raw, status, psi_at_pts = RZ_to_rhop(rpts, zpts, eq, time)

RZ = np.array([rpts,zpts])
psin_at_pts = RZ_to_psin(rpts, zpts, eq, time)
psin_at_pts_2 = RZ_to_psin(rpts, zpts, eq, time, method='cubic')
psin_at_pts_3 = RZ_to_psin_2(RZ,eq,time)

#plt.plot(rpts,psin_at_pts)
#plt.plot(rpts,psin_at_pts_2,'-')
#plt.plot(rpts,psin_at_pts_3,'--')
#plt.show()

#print(f'RZ: {RZ}') 
#print(f'psin_at_pts = {psin_at_pts}')
#print(f'psin_at_pts_2 = {psin_at_pts_2}')
#plot_psi(eq,time,nl=30)


#sys.exit()
runID=args.runid
if runID is not None:
    tr=ps.Transp(pulse,runID)
# read Ti profiles

sig=args.profile
(t,x,vals) = (tr.t,tr.x,tr.profile(sig))
trind = np.abs(tr.t+40-time).argmin()
print(f'x.shape: {x.shape}')
print(f't.shape: {t.shape}')
plt.title(f'{pulse} TRANSP {runID} {sig} on {args.dda}/{args.uid}/{args.sequence} at {time}s')
Rmid = rhot_to_RZ_midplane(tr.x[trind], eq, time)

plt.plot(tr.x[trind], vals[trind], label=f'on rhot at {t[trind]+40.:.3f}s')
#plt.plot(tr.x[trind]**2, ti[trind], 'g', label=f'on ftor_norm at {t[trind]+40.:.3f}s')
plt.plot(Rmid, vals[trind],'r', label=f'on R midplane at {eq.t[tind]:.3f}s')
#plt.xlabel('eV')
plt.legend()
plt.show()

ne = tr.profile('NE') 

if args.save: 
    outfile = f"tmp/{pulse}_{runID}_{sig}_at_{40+tr.t[trind]:.2f}s_profile.txt"

    data = np.column_stack([tr.x[trind], Rmid, vals[trind]])

    np.savetxt(outfile, data,
           header=f'Remapped using {eq.pulse}/{args.dda}/{args.uid}/{args.sequence}\n'
                  f'At t={eq.t[tind]:.3f}s Rmag={eq.Rmag[tind]:.3f}m, Zmag={eq.Zmag[tind]:.3f}m\n'
                  f'R of separatrix at the midplane is Rbnd={Rbnd_at_midplane(eq, time):.3f}m\n'
                  f'    rhot          Rmid           {sig}',
           fmt="%.8e")


