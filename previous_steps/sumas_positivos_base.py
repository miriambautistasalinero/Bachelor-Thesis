
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import seaborn as sns
import astropy.stats as sts
from scipy.stats import kurtosis, skew, mode
import sys

#Insert path for load mat modules
sys.path.insert(0, 'C:/Users/miria/OneDrive/Documentos/Articulos UVA/data/Resultados senales base')

import loadmats as lm

pathbase = "C:/Users/miria/OneDrive/Documentos/Articulos UVA/data/Resultados senales base/mats"
def sumas_positivos_base(k):
    '''
    

    Parameters
    ----------
    k : str
        Frequency band.

    Returns
    -------
    suma_positivos: dict
                    Dictionary with each group as key. Matrix of the sum 
                    of positive values in each sample. 

    '''
    
    pacientesbase = os.listdir("C:/Users/miria/OneDrive/Documentos/Articulos UVA/data/Resultados senales base/mats")

    grupo = {'control':[], 'DCL':[], 'EA':[] }

    correlation = {'control':[], 'DCL':[], 'EA':[] }
    
    for p in pacientesbase:
        dsi= lm.loadmat(os.path.join(pathbase,p))['dynamicStateInfo']  #Load the mat and extract info
        group =dsi['group'] #Find group
        correlationp=dsi['correlation'][k]  #Extract correlation
        
        #Assing name and corr to each group
        if group == 'DCL1' or group =='DCL2' or group == 'DCL':
            grupo["DCL"].append(p)
            correlation['DCL'].append(correlationp)
            
        elif group == 'EA1' or group == 'EA2':
            grupo["EA"].append(p)
            correlation['EA'].append(correlationp)
            
        else:
            grupo["control"].append(p)
            correlation['control'].append(correlationp)

    #dict for positive values
    suma_positivos= {'control':  np.zeros((len(grupo['control']),correlationp.shape[1])), 'DCL': np.zeros((len(grupo['DCL']),correlationp.shape[1])), 'EA':np.zeros((len(grupo['EA']),correlationp.shape[1]))}
    
    for g in grupo.keys(): #control, dcl, ea
        for l in range(len(grupo[g])): #Number of patients in each group
            for i in range(correlationp.shape[1]): #Total number of samples
                c = correlation[g][l][:,i]  #Takes columns
                pos = []     #list where positive values are going to be saved
                for e in c:  
                    if e >= 0:  #if positive add to pos and to dict
                        pos.append(e)
                        suma_positivos[g][l,i] = np.sum(pos)
    
    #encontrar el binwidth               
    listapos=[]
    for kk in suma_positivos:
        listapos.append(suma_positivos[kk])
    todas_pos=np.concatenate(listapos).flatten()
                        
    
    ##Hallar el binwidth para todos los valores
    #positivo
    opt_width_scott_pos, bins_scott_pos = sts.scott_bin_width(todas_pos, return_bins=True)
    bpos=np.arange(np.min(todas_pos), np.max(todas_pos) + opt_width_scott_pos, opt_width_scott_pos)
   
    return suma_positivos, opt_width_scott_pos

postivos_base_delta = sumas_positivos_base('correlationDelta')