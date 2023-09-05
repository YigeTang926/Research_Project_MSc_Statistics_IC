import os
import numpy as np
import pandas as pd
import scipy

# Get the path
args = {}
args['seed'] = 1234
args['cwd'] = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
args['derived_dir'] = os.path.join(args['cwd'], 'data', 'derived', 'czech')

# Read the data
wind = pd.read_csv(os.path.join(args['derived_dir'], 'wind_time_series.csv'),index_col=0)
solar = pd.read_csv(os.path.join(args['derived_dir'], 'solar_time_series.csv'),index_col=0)

# Periodogram
def periodogram(X,shift):
    '''function to compute the periodogram of time series X, outputting both the periodogram and the Fourier frequencies'''
    '''shift = 1 will do an fftshit onto [-1/2,1/2), shift = 0 will do an fft onto [0,1)'''
    N = len(X) # length of time series X
    if shift:
        Sp = (1/N)*abs(np.fft.fftshift(np.fft.fft(X)))**2
        f = -1/2 + np.arange(0,N)/N
    else:
        Sp = (1/N)*abs(np.fft.fft(X))**2
        f = np.arange(0,N)/N
    return f, Sp

# Direct spectral estimator
def direct(X,p,shift):
    '''The function that compute the direct spectrum estimator using the periodogram function
with p*100% cosine taper for a time series X'''
    N = len(X) # length of X
    h = np.ones(N)
    floor = int(p*N) # floor of pN
    pN2 = int(int(p*N)/2) # calculate the floor of pN and divide by 2
    # calculate the h value by definition
    for i in range(pN2):
        h[i] = 0.5*(1-np.cos((2*np.pi*(i+1))/(floor+1)))
    for i in range(N-pN2,N):
        h[i] = 0.5*(1-np.cos(2*np.pi*(N-i)/(floor+1)))
    C = 1/np.sum(h**2) # calculate the normalizing constant
    h = np.sqrt(C)*h
    f, S_pf = periodogram(h*X,shift) # using the periodigram()
    S_pf = S_pf*N # output spectrum estimate 
    return f, S_pf

# Find the peaks of frequency
def find_peaks(f,sp,n,threshold,distance):
    '''function to find the peaks of frequency'''
    peaks, _ = scipy.signal.find_peaks(x=sp,threshold=threshold,distance=distance)
    peak_f = f[peaks]
    peak_sp = sp[peaks]
    # order it according to the value of peak_sp
    peak_f = peak_f[np.argsort(peak_sp)[::-1]]
    peak_sp = peak_sp[np.argsort(peak_sp)[::-1]]
    peak_T = np.round(1/peak_f,2)
    if n <= len(peak_T):
        return peak_T[:n], peak_f[:n], peak_sp[:n]
    else:
        return peak_T, peak_f, peak_sp

# Get the periods for each node
wind_periods = pd.DataFrame()
solar_periods = pd.DataFrame()
for i in range(len(wind.columns)):
    wind_center = wind.iloc[:,i].values-wind.iloc[:,i].values.mean()
    f_wind, sp_wind = direct(wind_center,0.5,shift=1)
    peak_T_wind, _, _= find_peaks(f_wind[f_wind>0],sp_wind[np.where(f_wind>0)],8,5e5,20)
    peak_T_wind = peak_T_wind[peak_T_wind>3]
    wind_periods[wind.columns[i]] = peak_T_wind[0:5]
    solar_center = solar.iloc[:,i].values-solar.iloc[:,i].values.mean()
    f_solar, sp_solar = direct(solar_center,0.5,shift=1)
    peak_T_solar, _, _ = find_peaks(f_solar[f_solar>0],sp_solar[np.where(f_solar>0)],8,5e5,20)
    peak_T_solar = peak_T_solar[peak_T_solar>3]
    solar_periods[solar.columns[i]] = peak_T_solar[0:5]

# save the results to csv
wind_periods.to_csv(os.path.join(args['derived_dir'], 'wind_periods.csv'))
solar_periods.to_csv(os.path.join(args['derived_dir'], 'solar_periods.csv'))