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
        self._transpcid = str(pulse)+str(runid)
        self._cdfloc = f'/common/transp_shared/Data/result/JET/{pulse}/{runid}/{pulse}{runid}.CDF'
        self._dat = Dataset(self._cdfloc)
        self._dat_list=[item for item in self._dat.variables.keys()]
        get = lambda key: np.array(self._dat.variables[key])
        self._x = get("X")
        self._t = get("TIME3")
        self._dt = np.concatenate(([(self._t[1] - self._t[0])/2], (self._t[2:] - self._t[:-2])/2, [(self._t[-1] - self._t[-2])]))
        self._dvol = get('DVOL')
        self._transp = {}
        self._units = {}
        self._long_name = {}
        self._transp_names = []
        self._transp_time = []
    @property
    def pulse(self):
        return self._pulse
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
    def all_keys(self):
        return self._dat_list
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
    def units(self,signal):
        if signal in self._transp.keys():                                                           
            return self._units[signal]
    def long_name(self,signal):
        if signal in self._transp.keys():                                                           
            return self._long_name[signal]
    def signal(self,signal):
        if signal in self._transp.keys():          
            return True
        else:
            print(f'Signal {signal} not present')
            return False
    def get_ind(value, data):
        return np.abs(data - value).argmin()

    def profile(self,signal):
        return get(signal) 

    def taver(tvec):
        return np.average(tvec,weights=self._dt)

    def add_data(self, *args):
        get=lambda key:np.array(self._dat.variables[key])
        for sig in args:
            try:
                #print(f'Signal  = {sig}, type(sig): {type(sig)}')
                self._transp[sig] = get(sig)
                self._units[sig] = self._dat.variables[sig].units.strip()
                self._long_name[sig] = self._dat.variables[sig].long_name.strip()
            except KeyError:
                print(f'No {sig} found in TRANSP output')
                pass
        #if type(signal)==list:
        #    for sig in signal:
        #        try:
        #            #data = get(sig)
        #            self._transp[sig]=get(sig)
        #        except KeyError:
        #            print(f'No {sig} found in TRANSP output')
        #            pass
        #elif type(signal)==str:            
        #    try:
        #        #data = get(signal)
        #        self._transp[signal]=get(signal)
        #    except KeyError:
        #        print(f'No {signal} found in TRANSP output')
        #        pass
    

class Densities(Transp):
    def __init__(self,pulse,run):
        Transp.__init__(self,pulse,run)
        self._ntnd = None
        self._ntne = None
        self._tttd = None
        self._dthd = None
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
            return False 
    @property              
    def dthd(self):
        if 'ND' in self._dat_list:
            return self._dthd
        else:
            return False

    @property
    def tttd_av(self):
        return self._tttd_av
    @property
    def meff(self):
        return self._meff

    def calculate_concentrations(self):
        signals = ('NH','ND','NT','NE','NMINI','NMINI_H')
        ions = ['NT','ND','NH','NMINI_H']
        mass={'NT':3,'ND':2,'NH':1,'NMINI_H':1}

        self.add_data(*signals)
	
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
            self._tttd_av = np.average(self._tttd,weights=self._dt)
            self._ntne = prof_int['NT']/prof_int['NE'] 
            self._ntne_av = np.average(self._ntne,weights=self._dt)
        
        if 'NT' and 'ND' in self._dat_list:
            self._ntnd = prof_int['NT']/prof_int['ND']
            self._ntnd_av = np.average(self._ntnd,weights=self._dt)
            self._dthd = prof_int['ND']/th_ion_density 

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
    def beam_dt(self):
        signals = ('BTNTS_DT','BBNTS_DT')
        return sum(self.transp(sig) for sig in signals if sig in self._transp.keys())
    @property                                                                           
    def beam_dd(self):                                                                  
        signals = ('BTNTS_DD','BBNTS_DD')                                               
        return sum(self.transp(sig) for sig in signals if sig in self._transp.keys())

    @property
    def beam_td(self):
        signals = ('BTNTS_TD','BBNTS_TD')
        return sum(self.transp(sig) for sig in signals if sig in self._transp.keys())
    @property
    def beam_tt(self):                                                                
        signals = ('BTNTS_TT','BBNTS_TT')                                    
        return sum(self.transp(sig) for sig in signals if sig in self._transp.keys())
    def thermal_dt(self):                                                                   
        signals = ('NEUTX_TT','NEUTX_DT')                                     
        return sum(self.transp(sig) for sig in signals if sig in self._transp.keys()) 

    @property
    def transp_dd(self):
        signals = ('NEUTX_DD','BTNTS_DD','BBNTS_DD')
        return sum(self.transp(sig) for sig in signals if sig in self._transp.keys())
    @property
    def transp_dt(self):
        signals = ('NEUTX_DT','BTNTS_DT','BTNTS_TD','BBNTS_DT')
        return sum(self.transp(sig) for sig in signals if sig in self._transp.keys())
    @property
    def transp_tt(self):
        signals = ('NEUTX_TT','BTNTS_TT','BBNTS_TT')
        return sum(self.transp(sig) for sig in signals if sig in self._transp.keys())
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
        dt_signals = ('NEUTX_DT', 'NEUTX_TT',
                      'BTNTS_DT','BTNTS_TT','BTNTS_TD',
                      'BBNTS_DT','BBNTS_TT'
                    )
        self.add_data(*dd_signals)
        self.add_data(*dt_signals)
        self.add_data('NEUTT','NEUTX','BBNTS','BTNTS')

        dd_rate = sum(self.transp(sig) for sig in dd_signals if sig in self._transp.keys())
        dt_rate = sum(self.transp(sig) for sig in dt_signals if sig in self._transp.keys())

        if 'NEUTX_DD' and 'BTNTS_DD' in self._transp.keys():
            dd_rate_wiesen = self.transp('NEUTX_DD')+ self.transp('BTNTS_DD')+self.transp('BBNTS_DD')
        if 'NEUTX_DT' and 'BTNTS_DT' in self._transp.keys(): 
            dt_rate_wiesen = self.transp('NEUTX_DT')+self.transp('BTNTS_DT')        

        self._Rdtdd = dt_rate/dd_rate
        self._Rdtdd_av = np.average(self._Rdtdd,weights=self._dt)
        if 'NEUTX_DD' and 'NEUTX_DT' in self._transp.keys(): 
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
        

class EXP():
    def __init__(self,pulse,runid):
        self._pulse = pulse
        self._runid = runid
        self._dda = 'KT5P'
        self._seq = 0
        self._ierr=ppf.ppfgo(pulse=self._pulse,seq=self._seq)
        ppf.ppfsetdevice("JET")
        ppf.ppfuid("JETPPF","r")
        get=lambda dtyp:ppf.ppfget(self._pulse,self._dda,dtyp)
        [_,_,self._tttd,_,self._time,_]=get('TTTD')
        [_,_,self._dthd,*_]=get('DTHD')
        [_,_,self._hthd,*_]=get('HTHD')

    @property
    def t(self):
        return self._time
    @property
    def tttd(self):
        return self._tttd
    @property              
    def dthd(self):
        return self._dthd
    @property              
    def hthd(self):
        return self._hthd

    




