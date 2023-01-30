from netCDF4 import Dataset
import ppf
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d 
import sys
import plotWindow as pw


class Transp():
    def __init__(self,pulse,runid):
        self._pulse = pulse
        self._runid = runid
        self._transpcid = None
        self._cdfloc = f'/common/transp_shared/Data/result/JET/{pulse}/{runid}/{pulse}{runid}.CDF'
        self._dat = Dataset(self._cdfloc)
        self._dat_list=[item for item in self._dat.variables.keys()]
        get = lambda key: np.array(self._dat.variables[key])
        self._x = get("X")
        self._t = get("TIME3")
        self._dt = np.concatenate(([(self._t[1] - self._t[0])/2], (self._t[2:] - self._t[:-2])/2, [(self._t[-1] - self._t[-2])]))
        self._dvol = get('DVOL')
        self._transp = {}
        self._transp_names = []
        self._transp_time = []
        
    @property
    def runcid(self):
        return self._transpcid
    @property
    def transpcid(self):
        return self._transpcid
    @property
    def runid(self):
        return self._runid
    @property
    def x(self):
        return self._x
    @property
    def t(self):
        return self._t
    @property
    def transp_time(self):
        return self._t
    @property
    def dt(self):
        return self._dt

    def transp(self,signal):
        if signal in self._transp.keys():
            return self._transp[signal]
        else:
            print(f'Signal {signal} not present')
            return None
    
    def get_ind(value, data):
        return np.abs(data - value).argmin()

    def profile(self,signal):
        return get(signal) 

    def taver(tvec):
        return np.average(tvec,weights=self._dt)

    def add_data(self, signal):
        get=lambda key:np.array(self._dat.variables[key])
        if type(signal)==list:
            for sig in signal:
                try:
                    data = get(sig)
                    self._transp[sig]=data
                except KeyError:
                    print(f'No {sig} found in TRANSP output')
                    pass
        elif type(signal)==str:            
            try:
                data = get(signal)
                self._transp[signal]=data
            except KeyError:
                print(f'No {signal} found in TRANSP output')
                pass
    

class Densities(Transp):
    def __init__(self,pulse,run):
        Transp.__init__(self,pulse,run)
        self._ntnd = None
        self._ntne = None
        self._tttd = None
        self._meff = None
        self._ntne_av = None
        self._ntnd_av = None
        self._tttd_av = None

    @property
    def ntnd(self):
        return self._ntnd
    @property
    def ntnd_av(self):
        return self._ntnd_av

    @property
    def ntne(self):
        return self._ntne
    @property
    def ntne_av(self):
        return self._ntne_av
    @property
    def tttd(self):
        if 'NT' in self._dat_list:
            return self._tttd
        else:
            return None
    @property
    def tttd_av(self):
        return self._tttd_av
    @property
    def meff(self):
        return self._meff

    def calculate_concentrations(self):
        signals = ['NH','ND','NT','NE']
        ions = ['NT','ND','NH']
        mass={'NT':3,'ND':2,'NH':1}

        self.add_data(signals)
	
        #integrate
        prof_int = {}
        for sig in signals:
            if sig in self._dat_list:
                prof_int[sig] = np.sum(self._dvol*self._transp[sig],axis=1)

        th_ion_density = sum(prof_int[ion] for ion in ions if ion in self._dat_list)

        # effective mass
        weighted_sum = sum(mass[ion]*prof_int[ion] for ion in ions if ion in self._dat_list)
        self._meff = weighted_sum/th_ion_density

        # time vectors
        if 'NT' in self._dat_list:
            self._tttd = prof_int['NT']/th_ion_density 
            self._ntne = prof_int['NT']/prof_int['NE']
        if 'NT' and 'ND' in self._dat_list:
            self._ntnd = prof_int['NT']/prof_int['ND']
        # average
        self._ntne_av = np.average(self._ntne,weights=self._dt)
        self._ntnd_av = np.average(self._ntnd,weights=self._dt)
        self._tttd_av = np.average(self._tttd,weights=self._dt)

        return None


class Neutrons(Densities):
    def __init__(self,pulse,run):
        Densities.__init__(self,pulse,run)
        self._Rdtdd = None
        self._Rdtdd_av = None
        self._RNT = None
        self._RNT_time  = None           
        self._R14 = None
        self._R14_time = None   
        self._wiesen = None
    
    @property
    def Rdtdd(self):
        return self._Rdtdd
    @property
    def Rdtdd_av(self):
        return self._Rdtdd_av
    @property
    def rnt(self):
        return (self._RNT_time,self._RNT)
    def rnt_on_transp(self):
        return self._rnt_on_transp
    @property
    def r14(self):
        return (self._R14_time, self._R14)
    @property
    def wiesen(self):
        return self._wiesen


    def get_transp_neutrons(self):
        self.calculate_concentrations()
        dd_signals = ['NEUTX_DD','BTNTS_DD','BBNTS_DD']
        dt_signals = ['NEUTX_DT', 'NEUTX_TT',
                      'BTNTS_DT','BTNTS_TT','BTNTS_TD',
                      'BBNTS_DT','BBNTS_TT'
                        ]
        self.add_data(dd_signals)
        self.add_data(dt_signals)
        self.add_data(['NEUTT','NEUTX','BBNTS','BTNTS'])

        dd_rate = sum(self.transp(sig) for sig in dd_signals if sig in self._transp.keys())
        dt_rate = sum(self.transp(sig) for sig in dt_signals if sig in self._transp.keys())
        dd_rate_wiesen = self.transp('NEUTX_DD')+ self.transp('BTNTS_DD')+self.transp('BBNTS_DD')
        dt_rate_wiesen = self.transp('NEUTX_DT')+self.transp('BTNTS_DT')        

        self._Rdtdd = dt_rate/dd_rate
        self._Rdtdd_av = np.average(self._Rdtdd,weights=self._dt)

        self._wiesen = self._ntnd/(dt_rate_wiesen/dd_rate_wiesen)

        return None
    
    def get_exp(self):
        PPFseq = 0
        user = 'JETPPF'
        #print(f'Pulse: {self._pulse}, PPF sequence = {PPFseq}, user: {user}')
        ier = ppf.ppfgo(self._pulse, PPFseq)
        #print(f'ier = {ier},user: {user}')
        if ier != 0:
            raise OMFITexception(f'ier = {ier}. Error initialising PPF routines. Aborting.')
        # Set User ID for reading
        ier = ppf.ppfuid(user, rw="R")
        self._RNT = ppf.ppfdata(self._pulse, 'TIN', 'RNT')[0]
        self._RNT_time = ppf.ppfdata(self._pulse, 'TIN', 'RNT')[2]        
        self._R14 = ppf.ppfdata(self._pulse, 'TIN', 'R14')[0]
        self._R14_time = ppf.ppfdata(self._pulse, 'TIN', 'R14')[2]
        # interpolate R14 onto RNT
        fr14 = interp1d(self._R14_time,self._R14)
        frnt = interp1d(self._RNT_time,self._RNT)

        self._RNT_on_transp = frnt(self._t+40.)
        self._RNT_on_R14 = frnt(self._R14_time)
        return None
        #return (self._NEUT_Time,self._NEUT, self._R14_time, self._R14 )
        

if __name__=='__main__': 

    pulse=int(sys.argv[1])
    runid=str(sys.argv[2])
    profs = Neutrons(pulse, runid)
    profs.get_transp_neutrons()
    #pulse averaged Wiesen coefficient 
    ratio = np.average(profs.ntnd/profs.Rdtdd,weights=profs.dt)
    ratio_W = np.average(profs.wiesen,weights=profs.dt)

    Meff = np.average(profs.meff,weights=profs.dt)
    text = (
    r'$T_{conc}=\left<\frac{n_T}{n_D+n_T}\right>=$'+f'{profs.tttd_av*100:.4f}%\n'
    r'$ \left<n_T/n_D\right>=$'+f'{profs.ntnd_av*100:.4f}%\n'
    r'$ \left<n_T/n_e\right>=$'+f'{profs.ntne_av*100:.4f}%'
    )

    print(f'Meff = {Meff}')

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
    ax.plot(profs.t+40.,profs.tttd*100, color='k',linewidth=2,label='nt/(nt+nd)')
    ax.plot(profs.t+40.,profs.ntne*100, color='r',linewidth=2,label='nt/ne',linestyle='dashed')
    ax.set_xlim(profs.t[0]+40,profs.t[-1]+40)
    ax.set_ylim(0., 2*np.amax(profs.tttd*100))
    xleft,xright = ax.get_xlim()
    ymin,ymax=ax.get_ylim()
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('%')

    ax.text(xleft+0.04*(xright-xleft),(ymax+ymin)/2+0.3*(ymax-ymin),text,fontdict=font,bbox=props)

    #cornernote(f'{profs.transpcid}', '', ax=ax)
    ax.legend()

    win.addPlot('T conc',fig)

    #fig,ax = fn.subplots(label='Neutron rates vs T density')
    fig = plt.figure()
    fig.suptitle(f'Neutron rates vs T density', fontsize=13)
    ax = fig.add_subplot(111)

    text = (
    r'$T_{conc}=\left<\frac{n_T}{n_D+n_T}\right>=$'+f'{profs.tttd_av*100:.3f}%\n'
    r'$ \left<n_T/n_D\right>=$'+f'{profs.ntnd_av:.6f}\n'
    r'$ \left<R_{DT/DD}\right>=\left<\frac{R_{nDT}}{R_{nDD}}\right>=$'+f'{profs.Rdtdd_av:.4f}\n'
    r'$ \left<(n_T/n_D)/R_{DT/DD}\right> = $'+f'{ratio:.4f}'+'\n'
    r'$ \left< Wiesen\,\, coefficient\right>:$'+f' {ratio_W:.4f}\n'
    r'$n_D,n_T,R_{nDT},R_{nTT}$ volume integrated.'
    )

    ax.set_title(f'{profs.transpcid} Neutron rates vs T density')

    ax.plot(profs.t+40.,profs.ntnd/profs.Rdtdd, color='k',linewidth=2,label=r'$(nt/nd)/R_{DT/DD}$')
    ax.plot([profs.t[0]+40.,profs.t[-1]+40.],[ratio,ratio],linestyle='dashed',color='darkred',linewidth=2,label= 'average')

    ax.set_xlim(profs.t[0]+40,profs.t[-1]+40)
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
    r'$T_{conc}=\left<\frac{n_T}{n_D+n_T}\right>=$'+f'{profs.tttd_av*100:.3f}%\n'
    r'$ \left<n_T/n_D\right>=$'+f'{profs.ntnd_av:.6f}\n'
    r'$ \left<R_{DT/DD}\right>=\left<\frac{R_{nDT}}{R_{nDD}}\right>=$'+f'{profs.Rdtdd_av:.4f}'
    )

    ax.set_title(f'{profs.transpcid} Neutron rates vs T density')


    ax.plot(profs.t+40.,profs.Rdtdd, color='b',linewidth=2,label=r'$R_{DT/DD}$')
    ax.plot([profs.t[0]+40.,profs.t[-1]+40.],[profs.Rdtdd_av,profs.Rdtdd_av],linestyle='dashed',color='darkred',linewidth=2,label= 'average')

    ax.set_xlim(profs.t[0]+40,profs.t[-1]+40)
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


