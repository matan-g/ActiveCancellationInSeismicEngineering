# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 12:59:18 2017
@author: Matan Grossmann


***** STRUCTURAL DYNAMICS TERM PROJECT *****

This file is a numerical model of vibrations in a SDOF system:
    
The approach is to provide active damping by using destructive
interference from a massless vibration mechanism attached to our 
SDOF system.

*Convolution Integral translated from Matlab code from Paco
"""

import os
os.chdir('C:\Anaconda\StructuralDynamics_Term');

import numpy as np
import pandas as pd
from pandas import ExcelWriter
import scipy 
from scipy import signal
import heapq
import math
import matplotlib.pyplot as plt
import sklearn
from sklearn.neighbors import NearestNeighbors 


# ** Upload earthquake **
#eq_name = raw_input('Type name of Earthquake \n ::')
eq_name = 'Chalfant Valley' #fill in info from downloaded file
dt = .02 #time step & sampling rate (1/sec) from downloaded file
earthquake = pd.read_excel('ChalfantValley_EQ.xlsx', header = None)
earthquake_df = pd.DataFrame(earthquake)
earthquake = earthquake_df.values.reshape(earthquake.size,1)

# ** Initial Condiitons and Parameters **
k = 300
m = 200
wn = math.sqrt(k/m)
pi = math.pi
fn = wn /(2*pi)
Tn = 1/fn    
ic = 0.0, 0.0



# ***** NO DAMPING: BENCHMARK/COMPARITIVE CASE *****

#convolution integral
t_end = len(earthquake) * dt
tspan = np.arange(0,t_end,dt)
p = earthquake*m
tau = 0;
unitresp = np.zeros((len(tspan),len(tspan)))
u = np.zeros((len(tspan),1))
zeta = 0.05

for i in range(len(tspan)):
    tau = tspan[i]*.5
    u_t = (1/(m*wn))*np.exp(-zeta *wn*(tspan-tau))*np.sin((wn*(tspan-tau)))
    for n in range(len(u_t)):
        if tspan[n] < tau:
            u_t[n] = 0
    unitresp[i,:] = u_t

resp = p*unitresp*dt/m

for i in range(len(tspan)):
    u[i] = sum(resp[:,i])

# PLOTS
plt.figure(1)
plt.plot(tspan, u) #displacement plot over time
plt.figure(2)
plt.plot(tspan, earthquake)  #earthquake plot over time



# ***** ACTIVE DAMPING: DESTRUCTIVE INTERFERENCE ******
#Sampling Earthquake
upto = 3/dt #first 3 seconds
EQsamp = earthquake[:int(upto)]#sample

                    
# **Collecting Indicator Data**
                    
# Finding Dominant Frequency via Furier Transform                    
zer = np.zeros([1,int(upto)])#zero array for padding 
EQsamp_n = np.append([EQsamp],[zer]) #appending zeros 
fourier_transf = np.fft.fft(EQsamp_n)   #fast fourier transform
fourier_freq = np.fft.fftfreq(int(len(fourier_transf))) #frequency domain

plt.figure(3)
plt.plot(fourier_freq, fourier_transf) #magnitude plot -- look for phase plot

max_ind = np.where(fourier_transf==fourier_transf.max())
freq_est =float(fourier_freq[max_ind])

"""
raw_psd = abs(fourier_transf**2) #raw power spectrum density 
raw_psd = raw_psd[:(len(raw_psd)/2)] #half spectrum
max_ind = heapq.nlargest(1, range(len(raw_psd)), raw_psd.take) #finding index of max value
f_scale = np.arange(0,len(raw_psd/2))* dt/len(raw_psd)  #frequency scale
freq_est = float(f_scale[max_ind]) #dominant frequency estimate
"""

# Finding Peak Acceleration 
max_acc = float(max(abs((EQsamp))))

# Distance b/w peaks
peaks = scipy.signal.find_peaks_cwt(EQsamp.flatten(), np.arange(0.2,10,0.1))
peakOI_ind = heapq.nlargest(3, range(len(peaks)), EQsamp[peaks].take)
peak_dist = abs((peaks[peakOI_ind[0]] -peaks[peakOI_ind[1]])) * dt

#plotting peaks
plt.figure(4)
plt.plot(EQsamp)
plt.plot(peaks[peakOI_ind[0]],EQsamp[peaks[peakOI_ind[0]]],linestyle="",marker="o")
plt.plot(peaks[peakOI_ind[1]],EQsamp[peaks[peakOI_ind[1]]],linestyle="",marker="o")

#Max Disp of Structure
max_disp = float((max(abs(u[:150]))))

#Write all collected data into Data Frame (row)
new_row = pd.DataFrame({'Eq_name': [eq_name], \
                        'Eqsamp_Dom_Freq': [freq_est], \
                        'Eqsamp_Peak_Accel': [max_acc], \
                        'Peak_Dist': [peak_dist], \
                        'Max_Disp': [max_disp]}) 

#Load back end data --> append row to new data frame
BackEnd = pd.read_excel('BackEnd.xlsx')
BackEnd = pd.DataFrame(BackEnd)


"""
#Dataframe with comparable data
Compare = BackEnd[['Eq_name','Eqsamp_Dom_Freq', 'Eqsamp_Peak_Accel', \
                   'Peak_Dist', 'Max_Disp']]

                   
#Using Nearest Neighbor to find most similar earthquakes
neighb = NearestNeighbors(n_neighbors = 2, algorithm = 'brute', metric = 'euclidean')
neighb = neighb.fit(FullSet) # fits NN Alg
sim_EQs = neighb.kneighbors(FullSet,n_neighbors=2, return_distance = False) # runs NNalg
sim_EQs = pd.DataFrame(sim_EQs)
"""


#adjusted forcing function to account for active damping
#forcing function based on dominant frequency estimate
dest_func = - (max_acc/2)*np.sin(2*pi*freq_est*tspan)
p_n = p + dest_func
"""
p_n = np.add(p[:(len(p)/2)],dest_func[:(len(p)/2)])
p_n = np.append(p_n,p[len(p)/2:])
"""
"""
plt.figure(5)
plt.plot(tspan, p)
plt.plot(tspan, p_n)
"""

#convolution integral 
t_end = len(earthquake) * 0.02
dt = 0.02; #time step 
tspan = np.arange(0,t_end,dt)
tau = 0;


unitresp = np.zeros((len(tspan),len(tspan)))
u_n = np.zeros((len(tspan),1))

for i in range(len(tspan)):
    tau = tspan[i]*0.5
    u_t = (1/(m*wn))*np.exp(-zeta *wn*(tspan-tau))*np.sin((wn*(tspan-tau)))
    for n in range(len(u_t)):
        if tspan[n] < tau:
            u_t[n] = 0
    unitresp[i,:] = u_t

resp = p_n*unitresp*dt/m

for i in range(len(tspan)):
    u_n[i] = sum(resp[:,i])

plt.figure(6)
plt.plot(tspan, u_n)



# ***** POST PROCESSING ******
#**collecting data for storage** --> will mostly be used for creating clusters (if that is the method i decide on)

# Finding Dominant Frequency of entire EQ                    
zer_T = np.zeros([1,int(len(earthquake))])#zero array for padding 
EQ_n = np.append([earthquake],[zer_T]) #appending zeros 
fourier_transf_T = np.fft.fft(EQ_n)   #fast fourier transform
fourier_freq_T = np.fft.fftfreq(int(len(fourier_transf_T))) #frequency domain

plt.figure(7)
plt.plot(fourier_freq_T, fourier_transf_T) #magnitude plot -- look for phase plot

max_ind_T = np.where(fourier_transf_T==fourier_transf_T.max())
Tot_dom_freq =float(fourier_freq_T[max_ind_T])


#max accel of entire eq 
tot_max_accel = max(earthquake)

#attempted cancelling freq 
dest_freq = freq_est

# success or failure ? 1 or 0 
b_mark = max(u) + min(u) + np.average(u)
test = max(u_n) + min(u_n) + np.average(u_n)

if b_mark > test:
    success = 1
else:
    success = 0
 
extra_damp = b_mark - test

#optimal cancelling Freq --> this could be an optimization problem
#forcing function based on dominant frequency estimate
dest_func_opt = - (max_acc/2)*np.sin(2*pi*Tot_dom_freq*tspan)
p_opt = p + dest_func_opt
tau = 0;
unitresp = np.zeros((len(tspan),len(tspan)))
u_opt = np.zeros((len(tspan),1))

for i in range(len(tspan)):
    tau = tspan[i]*0.5
    u_t = (1/(m*wn))*np.sin((wn*(tspan-tau)))
    for n in range(len(u_t)):
        if tspan[n] < tau:
            u_t[n] = 0
    unitresp[i,:] = u_t

resp = p_n*unitresp*dt/m

for i in range(len(tspan)):
    u_opt[i] = sum(resp[:,i])

plt.figure(8)
plt.plot(tspan, u_opt)

b_mark2 = max(u_n) + min(u_n) + np.average(u_n)
test2 = max(u_opt) + min(u_opt) + np.average(u_opt)


if b_mark2 > test2:
    success2 = 1
else:
    success2 = 0
    
if success2 == 1:
    opt_freq = Tot_dom_freq
else:
    opt_freq = freq_est
  
    

#writing post processing data to file
new_row['dom_freq'] = Tot_dom_freq
new_row['max_accel'] = tot_max_accel
new_row['dest_freq'] = dest_freq
new_row['canc_success'] = success
new_row['extra_damp'] = extra_damp
new_row['opt_freq'] = opt_freq

#writing data to excel file 
BackEnd = BackEnd.append(new_row).reset_index(drop =1)
writer = ExcelWriter('BackEnd.xlsx')
BackEnd.to_excel(writer,'Sheet1') 
writer.save()
