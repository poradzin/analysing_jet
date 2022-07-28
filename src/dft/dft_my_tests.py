# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 16:37:01 2022

@author: mporad
"""
from __future__ import division
import numpy as np
from matplotlib import pyplot as plt
from scipy import fftpack



# follow the example from 
# https://www.allaboutcircuits.com/technical-articles/an-introduction-to-the-discrete-fourier-transform/


def z(x):
    '''
    '''
    return np.exp(x*1j*2*np.pi)

def F(n,x):
    '''
    n in (0,len(x)-1)
    '''
    N=len(x)
    W = np.zeros((N,N),dtype=complex)
    
    for i in range(N):
        for j in range(N):
            W[i,j] = z(-i*j/N)

    args = [x[i]*W[i,n] for i in range(N)]
    return np.sum(args)


def dft(f):
    '''
    f is sampled function over a whole period
    x axis is not specified
    return: fourier transformed vector 
    '''
    N=len(f)

    result = np.zeros(N,dtype=complex)

    W = np.zeros((N,N),dtype=complex)
    for i in range(N):
        for j in range(N):
            W[i,j] = z(-i*j/N)  
    for n in range(N):
        args = [f[i]*W[i,n] for i in range(N)]
        result[n] = np.sum(args)
    return result


def idft(f):
    '''
    f is sampled function over a whole period
    x axis is not specified
    return: fourier transformed vector 
    '''
    N=len(f)
    result = np.zeros(N,dtype=complex)
    W = np.zeros((N,N),dtype=complex)
    
    for i in range(N):
        for j in range(N):
            W[i,j] = z(i*j/N)  
    for n in range(N):
        args = [f[i]*W[i,n] for i in range(N)]
        result[n] = np.sum(args)/N
    return result




#x_n = np.array([0.2165, 0.8321, 0.7835, 0.5821, 0.2165, -0.5821, -1.2165, -0.8321])

def f(t):
    '''
    Sample function
    '''
    #return 5.+2*np.cos(2*np.pi*t-np.pi/2) +3*np.cos(4*np.pi*t)
    return 5.+1*np.cos(3*np.pi*t-np.pi/2) +3*np.cos(4*np.pi*t) + 2*np.sin(0.5*np.pi*t)

def ft(t,dft,T,t0=0,nmax=-1):
    '''   
    Parameters
    ----------
    n : non negative integer Fourier mode, n=< len(dft)//2
    T : scalar float, sampling period 
    dft : 1D numpy complex array, Discrete Fourier transform 

    Returns
    interpolation mode derived from Fourrier transform
    -------
    None.

    '''
    
    N=len(dft)
    #sum_terms = np.ones(N//2+1)
    result = np.zeros(N//2+1,dtype=float)
    result[0]= abs(dft[0])/N
    for n in np.arange(1,N//2+1):
        # basis frequency
        omg = 2*np.pi/T
        phi = np.angle(dft[n])
        result[n]=2./N* abs(dft[n])*np.cos(omg*n*(t-t0) + phi)
    if N-1 & 1:
        result[-1]=result[-1]/2.
    if nmax in range(N//2+1):
        for i in range(N//2+1):
            if i>nmax:
                result[i]=0.0

    return np.sum(result)
    #return np.dot(sum_terms,result)
 


# period 
T=1.5
t0=0. # beginning of sampling period
sample_times=np.arange(0.,T,0.05)

#sampling
f_n = f(sample_times)

N = len(sample_times)
#fourier transform and inverse
df_n = dft(f_n)
idf_n = idft(df_n)

# fft and inverse for comparison
sdf_n = fftpack.fft(f_n)
isdf_n = fftpack.ifft(sdf_n)   


plot_t = T    
x_list = np.linspace(0.,plot_t, 200)  
y_list = [ft(x, df_n,T) for x in x_list]   
y_check = f(x_list)
    
if True:
    fig = plt.figure()
    
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    
    
    ax1.stem(range(N), np.absolute(df_n),"k")
    ax1.set_title('FFT')
    #ax1.set(aspect=0.1)
    
    ax2.scatter(sample_times, np.absolute(idf_n),color='red')
    ax2.plot(x_list,y_check,color='black')
    ax2.plot(x_list,y_list,'--',color='blue')
    ax2.set_title('iFFT')
    ax2.set_xlim(right=x_list[-1])
    
    plt.subplots_adjust(hspace=0.5)
    plt.show()

