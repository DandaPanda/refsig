# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 09:31:00 2017

@author: DandaPanda
"""

"""
for import use: *from refsig import ref*
"""

import numpy as np
from sklearn.decomposition import FastICA


def avg(unipol):
    
    """ creates average reference signal from a particular set of data
    
    unipol(  channels,samples  )      set of unipolar data obtained from iEEG
          
    """
    
    avg = np.zeros(unipol.shape[1])

    for i in range(unipol.shape[1]):
        avg[i] = np.mean(unipol[:,i])

    return avg


def m1(unipol, N_iterations = 20, p = 0.25):
    
    """ creates the reference signal of an iEEG set of data
    using a comparative method based on correlation of independent components
    
    x = ref.m1(unipol, N_iterations = 20, p = 0.25)
    
    unipol(  channels,samples  )      set (matrix) of unipolar data obtained from iEEG
         
          
    returns referential signal ref1...linear vector (the 1 indicates that method I was used)
    and a coefficient of accuracy R, which shows how accurate the calculated reference is

        if R (minimal correlation coefficient) is smaller than p (R<p), then the calculated reference is fairly accurate 
        and can be used in further equations
    """
    
    #calculate average reference
    
    avg = np.zeros(unipol.shape[1])

    for i in range(unipol.shape[1]):
        avg[i] = np.mean(unipol[:,i])
    
    #create bipolar montage
        
    bipol = np.zeros([unipol.shape[0]-1,unipol.shape[1]])
    for i in range(unipol.shape[0]-1):
        bipol[i] = unipol[i+1]-unipol[i]
    
    
    #ICA of unipols and bipols
    
    unipol = np.transpose(unipol) #->(columns,rows)
    bipol = np.transpose(bipol)   #->(columns,rows)
    
    
    Rs = []
    Refs = []
    for i in range(N_iterations):
        Rs.append(0)
        Refs.append(0)

    #cycling through N number of iterations and obtaining the best possible result
    for k in range(N_iterations):        
    
        ica = FastICA(max_iter=1000)
        S_uni = ica.fit_transform(unipol)  # Reconstruct signals
        Q = ica.mixing_  # Get estimated mixing matrix

        ica = FastICA(max_iter=1000)
        S_bi = ica.fit_transform(bipol)  # Reconstruct signals

        S_uni = np.transpose(S_uni) #rows,columns
        S_bi = np.transpose(S_bi)   #rows,columns
    
    
        #method I calculation

    
        R = np.zeros((S_uni.shape[0],S_bi.shape[0]))   
    
    
        for i in range(S_uni.shape[0]):
            for j in range(S_bi.shape[0]):
                r = np.corrcoef(S_uni[i],S_bi[j])
                R[i,j] = abs(r[0,1])
            
        R1 = np.zeros(R.shape[0])
    
    
        for i in range(R.shape[0]):
            R1[i] = max(R[i])
      
        
        Rs[k] = min(R1) #coefficient of accuracy
        Ri = np.argmin(R1)
        
        P = np.mean(Q[Ri])
        Refs[k] = np.dot(P,S_uni[Ri])  
    
    
    #pinpointing the best possible reference based on the lowest R
    ref1 = Refs[np.argmin(Rs)]
    R = min(Rs)   
    
    if R < p:
        print('The method has found a good estimation of the referential signal with the minimal correlation coefficient',R,'being smaller than the given condition p =',p)
    else:
        print('Unfortunately the method could not meet the condition of p =',p,'and has reconstructed the referential signal with the minimal correlation coefficient being',R)
    
    #check if the phase is fine, if not -> turn it upside down    
    C = np.corrcoef(ref1,avg)
    if C[1,0] < 0:
        ref1 = ref1*(-1)
        C = np.corrcoef(ref1,avg)
    
    return ref1             

    
def m2(unipol):
        
    """ creates the reference signal of an iEEG set of data
    using method II described in http://doi.org/10.1109/TBME.2007.892929
    The difference from method I is that this method is solely based on
    calculating the reference instead of locating it as it was in method I
    
    x = ref.m2(unipol)
    
    unipol(  rows,columns  )      set of unipolar data obtained from iEEG
          (channels,samples)
          
    returns only referential signal ref2 (the 2 indicates method II was used) 
    """
    
    #calculate average reference
    
    avg = np.zeros(unipol.shape[1])

    for i in range(unipol.shape[1]):
        avg[i] = np.mean(unipol[:,i])
    
    #create bipolar montage
        
    bipol = np.zeros([unipol.shape[0]-1,unipol.shape[1]])
    for i in range(unipol.shape[0]-1):
        bipol[i] = unipol[i+1]-unipol[i]
    
    
    #ICA calc only for bipols
    
    bipol = np.transpose(bipol)   #->(columns,rows)
    
    ica = FastICA(max_iter=1000, algorithm = 'parallel')
    S_bi = ica.fit_transform(bipol)  # Reconstruct signals
    
    S_bi = np.transpose(S_bi) #rows,columns
    
    
    #method II calc
    
    suma = 0
    R = np.zeros([unipol.shape[0],unipol.shape[1]])
    
    for i in range(unipol.shape[0]):
        for j in range(S_bi.shape[0]):
            citatel = np.multiply(unipol[i,:],S_bi[j,:])
            citatel = np.mean(citatel)
            jmenovatel = np.power(S_bi[j],2)
            jmenovatel = np.mean(jmenovatel)
            suma = suma + ((citatel/jmenovatel)*S_bi[j])
        R[i] = unipol[i] - suma
        suma = 0   
        
    ref2 = np.zeros(R.shape[1])  
    
    for i in range(R.shape[1]):
        ref2[i] = np.mean(R[:,i])
      
    #check if the phase is fine, if not -> turn it upside down
    
    print('The method has succesfully reconstructed the referential signal')
    
    C = np.corrcoef(ref2,avg)
    if C[1,0] < 0:
        ref2 = ref2*(-1)
        C = np.corrcoef(ref2,avg)
    
    return ref2    