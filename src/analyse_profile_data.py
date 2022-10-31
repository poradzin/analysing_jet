import ppf
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage, misc, fftpack
import plotWindow as pw
import sys

###### Reading EFIT++ equilibrium data
def get_data(pulse, dda, uid, seq, output=None):
    '''
    '''
    # Initialise PPF routines
    ier = ppf.ppfgo(pulse, seq)
    if ier != 0:
        raise exception('Error initialising PPF routines. Aborting.')

    # Set User ID for reading
    ier = ppf.ppfuid(uid, rw="R")
    
    if output=='TE' and dda=='KK3':
        KK3ind = np.arange(96)
        KK3={}
        ind=[]
        for ii in KK3ind:
            if ii<9:
                tmp_ind='0'+str(ii+1)
            else:
                tmp_ind=str(ii+1)
            tmpT = ppf.ppfdata(pulse,dda,'TE'+tmp_ind, seq=seq, uid=uid)
            tmpR = ppf.ppfdata(pulse,dda,'RC'+tmp_ind, seq=seq, uid=uid)
            if tmpT[-2]>0 and tmpR[-2]>0:
                KK3['TE'+tmp_ind]= tmpT
                KK3['RC'+tmp_ind]= tmpR
                ind.append(tmp_ind)
        tmpG=ppf.ppfdata(pulse,dda,'GEN', seq=seq, uid=uid)
        return (KK3,ind,tmpG)
    if output:
        dtyp = output
        ihdat, iwdat, data, x, t, ier = ppf.ppfget(pulse, dda, dtyp)
        if ier != 0:
            print('Error reading signal ' + str(dda) + '/' + str(dtyp) + '. PPF Error Code ' + str(ier))
            print()
        print(f'Fetching {dda}/{output}/{uid}/{seq}')
            
        value, x_axis, time, nd, nx, nt, dunits, xunits, tunits, desc, comm, seq, ier = ppf.ppfdata(pulse, dda, output, uid=uid, seq=seq)
    else:
        return None
    
    return [time,x_axis, value]



#q, efit_x, efit_t, nd, nx, nt, dunits, xunits, tunits, desc, comm, seq, ier = ppf.ppfdata(pulse, eqdda, 'Q', uid=equid, seq=eqseq)
#efit_Rmag, dump1, dump2, nd, nx, nt, dunits, xunits, tunits, desc, comm, seq, ier = ppf.ppfdata(pulse, eqdda, 'RMAG', uid=equid, seq=eqseq)


def statistics(data,dda='',units='',out=False, plot=False,verbose=0):
    '''
    input: 1D ordered numpy array i.e. time slices [times,x,ne]
    name: 
    '''
    

    #[times,x,ne] = get_data(pulse, dda, uid,seq, output=dty)
    #[times,x,ne]  = data
    #data=times
    
    interval_data = data[1:] - data[:-1]

    if out:
        return interval_data
    else:
        ldata = len(data)
        av_int_data = np.sum(interval_data) / (ldata-1)

 
        print('Average interval {0:2.3f}{1:s}'.format(av_int_data,units))
        print(f'{dda} average frequency: {1/av_int_data:2.2f} Hz.')
        if verbose ==2:
            sd = np.sqrt(np.sum((interval_data - av_int_data) ** 2) / (ldata - 1))
            skewness = np.sum((interval_data - av_int_data) ** 3) / (ldata - 1) / np.power(
                        np.sum((interval_data - av_int_data) ** 2) / (ldata - 2), 3. / 2.)
            print('{0:s} range: [{1:2.3f},{2:2.3f}]s. Length: {3}'.format(dda,data[0],data[-1],ldata))
            print('The shortest interval: {0:2.4}{1:s}'.format(np.min(interval_data),units))
            print('The longest interval:  {0:2.4}{1:s}'.format(np.max(interval_data),units) )
            print('Interval standard deviation {0:3.8f}{1}'.format(sd,units))
            print('Skewness b1: {:2.3f}\n'.format(skewness))
        
    if plot and verbose>0:
        plt.hist(statistics(data,out=True) , density=False, bins=20)  # density=False would make counts
        plt.show()
        
    return None



def get_ind(value, data):
    return np.abs(data-value).argmin()

def plot_time_slice(times,x,data,time):

    ind = get_ind(time,times)
    plt.plot(x,data[len(x)*(ind):len(x)*(ind+1)])
    plt.title("Time slice %6.4f s" % times[ind])
    plt.show()

def plot_time_trace(times,x,data,x_point):
    ind = get_ind(x_point,x)
    lx=len(x)
    lt=len(times)
    indices = np.arange(lt)*lx+ind
    plt.plot(times,data[indices])
    plt.title("x point {:5.3f} m".format(x_point))
    plt.show()

def get_no_slices(data,tstart,tend):
    ind1 = get_ind(tstart,data)
    ind2 = get_ind(tend,data)
    return np.abs(ind2-ind1)



def plot_KK3_profile(data,ind,time):

    x_axis=[]
    y_axis=[]
    points = np.array([])
    times = []
    for ii in ind:
        #print(ii)
        index = get_ind(time,data['RC'+ii][2])
        #print('First index: ',index)
        x_axis.append(data['RC'+ii][0][index])
        index = get_ind(time,data['TE'+ii][2])
        #print('Second index: ',index)
        y_axis.append(data['TE'+ii][0][index])
        times.append(data['TE'+ii][2][index])
        np.append(points,[x_axis[-1],y_axis[-1]],axis=0)
    print(x_axis)
    #print(y_axis)
    #print(times)
    #plt.scatter(np.arange(len(x_axis))+1,x_axis)
    plt.scatter(x_axis,y_axis)
    plt.title("Time slice %6.4f s" % times[0])
    plt.show()


def check_when_heating_active(pulse, plot=False,verbose=0,window=None):
    '''
    '''

    data=[ ['NBI','PTOT','blue',False],['ICRH','PTOT','red',False] ]
    times={}
    values={}
    dft={}
    for  ind, dat in enumerate(data):
        [dda,dty,col,active]= dat
        value, x_axis, time, nd, nx, nt, dunits, xunits, tunits, desc, comm, seq, ier = ppf.ppfdata(pulse, dda, dty, uid='jetppf', seq=0)

        if len(value)>0:
            data[ind][3]=True
            power_cut_off = 0.0
            above_cut_off = value>power_cut_off
            nonzerotimes = time[above_cut_off]
            times[dda] = nonzerotimes
            values[dda]=value[above_cut_off]
            if dda=='ICRH':
                dft[dda]={}
                dft[dda]['fft'] = fftpack.fft(values[dda])
                L=len(values[dda])
                dft[dda]['P2']=abs(dft[dda]['fft']/L)
                dft[dda]['P1']=abs(dft[dda]['P2'][:L//2])
                dft[dda]['Re']=np.real(dft[dda]['fft']/L)[1:]
                dft[dda]['Im']=np.imag(dft[dda]['fft']/L)[1:]
                print(f'L={L}')

                
                
            if verbose:
                print(f'Times {dda} active: {nonzerotimes[0]:2.4f}s to {nonzerotimes[-1]:2.4f}s.')
                if dda=='NBI':
                    #check species
                    octants = ['NBI4','NBI8']
                    species = {}
                    for octant in octants:
                        print(f'Checking {octant}...')
                        ppfdata = ppf.ppfdata(pulse,octant,'GAS', seq=0, uid='jetppf')
                        if len(ppfdata[0])==0:
                            print(f'Octant {octant[-1]} NOT active')
                        else:
                            species[octant] = ppfdata[9].split()[-1]
                            print(f'Octant {octant[-1]}, species: {species[octant]}')
                print('\n')

        else:
            if verbose: print('No {0} found in pulse {1}\n'.format(dda,pulse))

    if plot:
        fig = plt.figure()
        fig.suptitle(f'Heating', fontsize=13)
        ax1 = fig.add_subplot(111)
        for dat in data:
            if dat[3]:
                ax1.plot(times[dat[0]], values[dat[0]], color = dat[2])
        window.addPlot('heating',fig)
        
    if plot and data[1][3]:
        signal = [['Re','green'],['Im','red'],['P2','black']]
        fig = plt.figure()
        
        fig.suptitle(f'{data[1][0]} FFT', fontsize=13)
        ax1 = fig.add_subplot(111)
        for sig in signal:        
            ax1.plot(range(len(dft[data[1][0]][sig[0]])),dft[data[1][0]][sig[0]], color = sig[1],label=sig[0])
        ax1.legend()
        #ax.set_xlim(left=0, right=50)
        window.addPlot('ICRH FFT',fig)
    
    #if plot:
    #    plt.show()
    
    return (times,window)


def check_when_data_meaningful(pulse,data_in, eps=1.,plot=False,verbose=0,har=None,output=None,windows=None):
    '''
    input: [times,x,exp_data], where len(exp_data)=len(x)*len(times)
    return: time vector with times when data is gratear than eps
    '''
    #print('Namespace: ', __name__)
    (dda, dty, uid, seq, unit)= data_in
    
    data = get_data(pulse, dda, uid,seq, output=dty)
    [times,x,exp_data] = data
    if dda is 'KK3':
        print(f'x = {x}')
        sys.exit()
    
    lt=len(data[0])
    lx=len(data[1])
    ldata = len(data[2])
    print('{0:s} range: [{1:2.3f},{2:2.3f}]s.\nLength/no. of slices: {3}'.format(dda,data[0][0],data[0][-1],lt))
    statistics(times,dda=dda,units='s',out=False, plot=plot,verbose=verbose)
    
    #times may be of the length x*ne
    if ldata==lx*lt:
        if verbose:
            print(f'Converging data list of length {ldata} into 2D array ({lx},{lt})')
        data2D=[]
        for ind, time in enumerate(data[0]):
            data2D.append(data[2][ind*lx:(ind+1)*lx])
        data2D=np.array(data2D)

        
    dataT=data2D.transpose()
    if verbose:
        print(f'np.shape(dataT): {np.shape(dataT)}')
    #filtering data
    filter_size = {'HRTS':10, 'CXD6':3,'CXG6':3,'CX7C':5,'CX7D':5,'default':10}
    if dda in filter_size:
        size=filter_size[dda]
    else:
        size = filter_size['default']
    size = np.min([size,lt,lx])
    result = ndimage.median_filter(dataT, size=size)
    
    sumax= np.sum(result, axis=0)

    condition = sumax>eps

    first_time = np.min(times[condition])

    last_time = np.max(times[condition])
    print(f'{dda} range after filtering: ({first_time:2.3f},{last_time:2.3f})s.')
    print(r'Maximum value {:3.4g}{} at {:2.3f}s.'.format(np.max(result[0]),unit,times[np.argmax(result[0])] ))
    #print('Earliest filtered {0:s} time slice: {1:2.3f}s and last time slice {2:2.3f}s '.format(dda, first_time, last_time))
    ind1 = get_ind(first_time, times)
    ind2 = get_ind(last_time, times)
    no_of_time_slices=ind2-ind1
    print('Number of time slices in the filtered range: {}'.format(no_of_time_slices))
    
    # har: heating active range
    if har:
        ind1 = get_ind(har[0], times)
        ind2 = get_ind(har[1], times)
        no_of_time_slices=ind2-ind1
        print(f'Number of {dda} filtered time slices and with active heating: {no_of_time_slices}')
        print(f'in the time range ({times[ind1]:2.3f},{times[ind2]:2.3f})s\n')
    print('\n')
    
    if plot:
                 
        fig = plt.figure()
        fig.suptitle(f'{dda}/{dty}/{uid}/{seq}', fontsize=13)
        ax1 = fig.add_subplot(221)  # left side
        ax2 = fig.add_subplot(222)# right side
        ax3 = fig.add_subplot(223) # down left side
        ax4 = fig.add_subplot(224)

        
        xax=range(len(times))
        ax1.contourf(times,x,dataT)
        ax2.contourf(xax,x,result)
        ax3.plot(times,dataT[0])
        ax4.plot(times,result[0], color='r')        

        ax1.set_title('Raw data')
        ax2.set_title(f'filtered: size={size}')

        print('############################################')
        #check of object present
        #if windows:
        windows.addPlot(dda,fig)
        #plt.show()
    #out = {'time':data[0],'x':data[1],'values':dataT, 'data': np.array([time,x,dataT])}
    out = {'time':data[0],'x':data[1],'values':dataT}   
    
    if output in out:
        return (out[output],windows)
    else:
        return (None,windows)
    

def check_CX(pulse, plot=False,verbose=0,checkHCD=False):
    '''
    '''
    #print('Namespace: ', __name__)
    data=[['HRTS','NE','jetppf',0,'m^-3'],['HRTS','TE','jetppf',0,'eV'],#['KK3','TE02','jetppf',0,'eV'],
          #['CXG6','TI','jetppf',0,'eV'],['CXD6','TI','jetppf',0,'eV'],
          ['CX7C','TI','jetppf',0,'eV'],['CX7D','TI','jetppf',0,'eV']]
    
    no_times=np.array([])
    times_minmax = np.array([])
    CX_active = np.array([])
    heat_active_range = []

    if plot:
        print('We\'re plotting')
        win=pw.plotWindow()
    else:
        win=None
    
    if checkHCD:
        nbi = False
        icrh = False
        
        # check nbi and icrh times
        heat_active = np.array([])
        (heat_times,win) = check_when_heating_active(pulse,plot=plot, verbose=1,window=win)

        if 'NBI' in heat_times:
            nbi=True
            (nbi_t_min,nbi_t_max) = (heat_times['NBI'][0],heat_times['NBI'][-1])
            heat_active = np.append(heat_active, [nbi_t_min,nbi_t_max])
            
        if 'ICRH' in heat_times:
            icrh=True
            (ic_t_min,ic_t_max) = (heat_times['ICRH'][0],heat_times['ICRH'][-1])
            heat_active = np.append(heat_active, [ic_t_min,ic_t_max])
            
            
        

        heat_active_range = [np.min(heat_active), np.max(heat_active)]
        
        
    for dat in data:
        (time,win)=check_when_data_meaningful(pulse,dat,plot=plot,verbose=verbose,har = heat_active_range,output='time',windows=win)
        if dat[0] in ['CXD6','CXG6']:
            no_times= np.append(no_times, len(time))
            times_minmax= np.append(times_minmax, [time[0],time[-1]])
    if plot:
        print('Showing')
        win.show()
    print(f'CXD6 or CXG7 available from {np.min(times_minmax):2.3f}s to {np.max(times_minmax):2.3f}s. for {int(np.max(no_times))} slices.')  

    #return None
def get_no_slices_and_plot(pulse,dda,dty,uid,seq,t_start,t_end,plot=True):
    [t,x,dat] = get_data(pulse,dda,uid,seq,output=dty)
    data2D = np.reshape(dat, (len(t), len(x)))
    n0=data2D[:,0][ np.logical_and(t>=t_start, t<=t_end) ]
    t0=t[np.logical_and(t>=ne_tstart, t<=ne_tend)]
    if plot:
        print(f'{dda}/{dty} in {t0[0]:2.3f}s,{t0[-1]:2.3f}s, {len(t0)} slices')
        plt.plot(t0,n0)
        plt.title(f'{dda}/{dty} in {t0[0]:2.3f}s,{t0[-1]:2.3f}s, {len(t0)} slices')
        plt.show()
    else:
        print(f'{dda}/{dty} in {t0[0]:2.3f}s,{t0[-1]:2.3f}s, {len(t0)} slices')
    return (t0,x,n0)

if __name__=='__main__':
    #import ppf
    #import numpy as np
    #import matplotlib.pyplot as plt
    #from scipy import ndimage, misc

    #import plotWindow as pw
    
    pulse = 100854

    plot=True

    verbose=0

    print(f'PULSE {pulse}\n')
        
    #check_CX(pulse,checkHCD=True, plot=plot,verbose=verbose)
    
    ne_tstart= 47.008
    ne_tend = 54.5230

    (t,x,data) = get_no_slices_and_plot(pulse,'HRTS','NE','JETPPF',0,ne_tstart,ne_tend,plot=True)
    #(t,x,data) = get_no_slices_and_plot(pulse,'KK3','TE01','JETPPF',0,ne_tstart,ne_tend,plot=True)
    #(t,x,data) = get_no_slices_and_plot(pulse,'CXD6','TIFS','JETPPF',0,ne_tstart,ne_tend,plot=True)
# tests of class plotWindow
    def fun(f,x):
        return f(x)

    def add_p(x,y,name,windows):
            f = plt.figure()
            #plt.plot(x, y, '--')
            size=10
            ax1 = f.add_subplot(221)
            ax2 = f.add_subplot(222)
            ax3 = f.add_subplot(223)
            ax4 = f.add_subplot(224)
            ax1.plot(x, y, '--')
            ax2.plot(x,y/2.,'-')
            ax3.plot(x,np.sqrt(np.abs(y)),'-')
            ax4.plot(x,y**2.,'-')
            ax1.set_title('Raw data')
            ax2.set_title(f'filtered: size={size}')        
            
            windows.addPlot(str(name), f)

            return None
    def plot_funs(functions):
        win = pw.plotWindow()
        for function in functions:

            x = np.arange(0, 10, 0.001)
            y = fun(function,x)
            add_p(x,y,function.__name__,win)
        win.show()
        return None

        
    #funs = [np.sin,np.cos,np.exp]
    #plot_funs(funs)


    #shot=99811

    sequence=0

    owner='JETPPF'

    dda='KK3'

    #data,ind,gen = get_data(pulse, dda, uid,seq, output='TE')


    #plot_KK3_profile(data,ind,49.9)

    #plt.plot(ppfdata[2],ppfdata[0])
    #plt.show()
     
