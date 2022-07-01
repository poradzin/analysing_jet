import ppf
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage, misc


pulse = 99813
dda = 'HRTS'
uid='jetppf'
seq = 0


plot=True

plot_stats = False

print(f'PULSE {pulse}\n')
###### Reading EFIT++ equilibrium data
def get_data(pulse, dda, uid, seq, output=None):
    '''
    '''
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
        print(f'Fetching {dda}/{output}/{uid}/{seq}')
        value, x_axis, time, nd, nx, nt, dunits, xunits, tunits, desc, comm, seq, ier = ppf.ppfdata(pulse, dda, output, uid=uid, seq=seq)
    else:
        return None
    
    return [time,x_axis, value]

#q, efit_x, efit_t, nd, nx, nt, dunits, xunits, tunits, desc, comm, seq, ier = ppf.ppfdata(pulse, eqdda, 'Q', uid=equid, seq=eqseq)
#efit_Rmag, dump1, dump2, nd, nx, nt, dunits, xunits, tunits, desc, comm, seq, ier = ppf.ppfdata(pulse, eqdda, 'RMAG', uid=equid, seq=eqseq)


def statistics(pulse,dda = 'hrts',dty='NE',uid='jetppf',seq=0, units='',out=False, plot_hist=False):
    '''
    input: 1D ordered numpy array i.e. time slices
    name: 
    '''
    
    #print( 'Statistics: Data type: {}'.format(type(data)))
    #print( 'Statistics: len(data) = {0}'.format( len(data) ) )
    [times,x,ne] = get_data(pulse, dda, uid,seq, output=dty)
    data=times
    
    interval_data = data[1:] - data[:-1]

    if out:
        return interval_data
    else:
        ldata = len(data)
        av_int_data = np.sum(interval_data) / (ldata-1)
        sd = np.sqrt(np.sum((interval_data - av_int_data) ** 2) / (ldata - 1))
        skewness = np.sum((interval_data - av_int_data) ** 3) / (ldata - 1) / np.power(
        np.sum((interval_data - av_int_data) ** 2) / (ldata - 2), 3. / 2.)
        #print('{0:s} range: [{1:2.3f},{2:2.3f}]s. Length: {3}'.format(name,data[0],data[-1],ldata)) 
        print('Average interval {0:4.5f}{1:s}'.format(av_int_data,units))
        print('HRTS average frequency: {:2.2f} Hz.'.format(1/av_int_data))
        print('The shortest interval: {0:2.4}{1:s}'.format(np.min(interval_data),units))
        print('The longest interval:  {0:2.4}{1:s}'.format(np.max(interval_data),units) )
        print('Interval standard deviation {0:3.8f}{1}'.format(sd,units))
        print('Skewness b1: {:2.3f}\n'.format(skewness))
        
    if plot_hist:
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



def check_when_heating_active(pulse, plot=False,verbose=True):
    '''
    '''

    data=[['NBI','PTOT','blue'],['ICRH','PTOT','red']]
    times={}
    for  ind, data in enumerate(data):
        [dda,dty,col]= data
        value, x_axis, time, nd, nx, nt, dunits, xunits, tunits, desc, comm, seq, ier = ppf.ppfdata(pulse, dda, dty, uid='jetppf', seq=0)

        if len(value)>0:
            nonzerotimes = time[value>0]
            times[dda]=( nonzerotimes[0], nonzerotimes[-1] )
            if verbose:
                print('{0} start time: {1:2.3f}s'.format(dda,nonzerotimes[0]))
                print('{0} stop  time: {1:2.3f}s'.format(dda,nonzerotimes[-1]))
            if plot:
                plt.plot(nonzerotimes, value[value>0], color = col)
                
        else:
            if verbose: print('No {0} found in pulse {1}\n'.format(dda,pulse))
    if plot:
        plt.show()
    
    return times


def check_when_data_meaningful(pulse,dda='HRTS',dty ='NE',uid='jetppf',seq=0,eps=1.,plot=False,verbose=False,checkHCD=False,output=None):
    '''
    input: [times,x,exp_data], where len(exp_data)=len(x)*len(times)
    return: time vector with times when data is gratear than eps
    '''
    data = get_data(pulse, dda, uid,seq, output=dty)
    [times,x,exp_data] = data

    if checkHCD:
        nbi = False
        icrh = False
        
        # check nbi and icrh times
        heat_active = np.array([])
        heat_times = check_when_heating_active(pulse,plot=False, verbose=False)
        if 'NBI' in heat_times:
            nbi=True
            (nbi_t_min,nbi_t_max) = heat_times['NBI']
            heat_active = np.append(heat_active, [nbi_t_min,nbi_t_max])
            print('Times NBI active: {0:2.4f}s to {1:2.4f}s.'.format(nbi_t_min,nbi_t_max))
        if 'ICRH' in heat_times:
            icrh=True
            (ic_t_min,ic_t_max) = heat_times['ICRH']
            heat_active = np.append(heat_active, [ic_t_min,ic_t_max])
            print('Times ICRH active: {0:2.4f}s to {1:2.4f}s.'.format(ic_t_min,ic_t_max))
            
        

        heat_active_range = [np.min(heat_active), np.max(heat_active)]

    
    
    lt=len(data[0])
    lx=len(data[1])
    ldata = len(data[2])
    print('{0:s} range: [{1:2.3f},{2:2.3f}]s.\nLength/no. of slices: {3}'.format(dda,data[0][0],data[0][-1],lt))
    #times may of the length x*ne
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
    filter_size = {'HRTS':10, 'CXD6':3,'CXG6':3,'default':10}
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
    print(f'{dda} range after filtering: [{first_time:2.3f},{last_time:2.3f}]s.') 
    #print('Earliest filtered {0:s} time slice: {1:2.3f}s and last time slice {2:2.3f}s '.format(dda, first_time, last_time))
    ind1 = get_ind(first_time, times)
    ind2 = get_ind(last_time, times)
    no_of_time_slices=ind2-ind1
    print('Number of time slices in the filtered range: {}\n'.format(no_of_time_slices))
    if checkHCD:
        print(f'Heating active from {heat_active_range[0]:2.3f}s to {heat_active_range[1]:2.3f}s')
        ind1 = get_ind(heat_active_range[0], times)
        ind2 = get_ind(heat_active_range[1], times)
        no_of_time_slices=ind2-ind1
        print(f'Number of {dda} time slices in the filtered range and with active heating: {no_of_time_slices}\n')

    
    if plot:
         
        fig = plt.figure()
        ax1 = fig.add_subplot(221)  # left side
        ax2 = fig.add_subplot(222)# right side
        ax3 = fig.add_subplot(223) # down left side
        ax4 = fig.add_subplot(224)
        
        xax=range(len(times))
        ax1.contourf(times,x,dataT)
        ax2.contourf(xax,x,result)
        ax3.plot(times,dataT[7])
        ax4.plot(times,result[7], color='r')        
        
        plt.show()
    #out = {'time':data[0],'x':data[1],'values':dataT, 'data': np.array([time,x,dataT])}
    out = {'time':data[0],'x':data[1],'values':dataT}   
    
    if output in out:
        return out[output]
    else:
        return None
    

# check_CX_data is obsolete and replaced by check_CX
def check_CX_data(pulse, plot=False,verbose=True):
    '''
    '''

    data=[['CXG6','TI','blue'],['CXD6','TI','red'],['CX7C','TI','green'],['CX7D','TI','black']]
    no_times=np.array([])
    times_minmax = np.array([])
    CX_active = np.array([])
    print('\n')
    for  ind, data in enumerate(data):
        [dda,dty,col]= data
        print(f'Fetching {dda} TI measurements...')
        value, x_axis, time, nd, nx, nt, dunits, xunits, tunits, desc, comm, seq, ier = ppf.ppfdata(pulse, dda, dty, uid='jetppf', seq=0)
        CX_active = np.append(CX_active, [time[0],time[-1]])
        #print(f'{dda} collected from {time[0]:2.3f}s to {time[-1]:2.3f}s')
        print(f'No of times slices: {len(time)}\n')
        if dda in ['CXD6','CXG6']:
            no_times= np.append(no_times, len(time))
            times_minmax= np.append(times_minmax, [time[0],time[-1]])
    CX_active_min,CX_active_max = (np.min(CX_active), np.max(CX_active))
    if verbose:
        #print('CX measurements available from {:2.3f}s to {:2.3f}s'.format(CX_active_min,CX_active_max))
        print(f'CXD6 or CXG7 available from {np.min(times_minmax):2.3f}s to {np.max(times_minmax):2.3f}s. for {int(np.max(no_times))} slices.')
    return None

def check_CX(pulse, plot=False,verbose=False):
    '''
    '''
    data=[['CXG6','TI','jetppf',0,'blue'],['CXD6','TI','jetppf',0,'red'],['CX7C','TI','jetppf',0,'green'],['CX7D','TI','jetppf',0,'black']]
    
    no_times=np.array([])
    times_minmax = np.array([])
    CX_active = np.array([])
    
    for dat in data:
        time=check_when_data_meaningful(pulse,dda=dat[0],dty =dat[1],uid=dat[2],seq=dat[3],plot=plot,verbose=verbose,output='time')
        if dat[0] in ['CXD6','CXG6']:
            no_times= np.append(no_times, len(time))
            times_minmax= np.append(times_minmax, [time[0],time[-1]])

    print(f'CXD6 or CXG7 available from {np.min(times_minmax):2.3f}s to {np.max(times_minmax):2.3f}s. for {int(np.max(no_times))} slices.')  

    return None


# check HRTS data
check_when_data_meaningful(pulse,dda='HRTS',dty ='NE',uid='jetppf',seq=0,plot=plot,checkHCD=True)
statistics(pulse, dda='hrts',dty='ne', units='s', plot_hist=plot_stats)

#plot_time_slice(times,x, ne,time)
#plot_time_trace(times,x,te1, 3.4)

check_when_heating_active(pulse,plot=plot)

check_CX(pulse)

    
#shot=99811

sequence=0

owner='JETPPF'

dda='KK3'

#data,ind,gen = get_data(pulse, dda, uid,seq, output='TE')

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


#plot_KK3_profile(data,ind,49.9)

#plt.plot(ppfdata[2],ppfdata[0])
#plt.show()
 
