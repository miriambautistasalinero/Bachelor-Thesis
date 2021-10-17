#SUMA DE CORRELACIONES POSITIVAS SUBROGADAS
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import seaborn as sns
import astropy.stats as sts
import sys
import io
import sumas_positivos_base

#Insert path for load mat modules
sys.path.insert(0, 'C:/Users/miria/OneDrive/Documentos/Articulos UVA/data/Resultados senales base')

import loadmats as lm

#path subrogate signals
path = "C:/Users/miria/OneDrive/Documentos/Articulos UVA/data/Subrogadas"
pacientessub = os.listdir(path)

#k = 'correlationDelta'
grupo = {'control':[], 'DCL':[], 'EA':[] }
correlation = {'control':[], 'DCL':[], 'EA':[] }
sumaspacientes_totales = {'control':[], 'DCL':[], 'EA':[] } #ir guardando las sumas de las subrogsdas de todos los pacientes
media = {'control':[], 'DCL':[], 'EA':[] } #medias

def grups(pacientessub):
    '''
    

    Parameters
    ----------
    pacientessub : list
        list of patients

    Returns
    -------
    grupo: dict
            dictionary with the name of the patient corresponding to each group

    '''
    for p in pacientessub:
    
        ssi= lm.loadmat(os.path.join(path, p))['surrogateStateInfo'] 
        group =ssi['group']
        
        
        
        if group == 'DCL1' or group =='DCL2' or group == 'DCL':
             grupo["DCL"].append(p)
        
        
        elif group == 'EA1' or group == 'EA2':
             grupo["EA"].append(p)
        
        else:
            grupo["control"].append(p)
        
        return grupo
    
    
def subrogadas_por_grupo(patients, g, k):
    '''
    
    #NECESITO QUE RUPO SEA UNICAMANETE DE UN GRUPO
    Parameters
    ----------
    patients : dict
            dictionary with the name of the patient corresponding to each group
    g : string 'DCL', 'EA', 'control'
    
    k: Frequency band for sum of señales base
        
    Returns
    -------
    None.

    '''
    
    for p in grupo[g]:
        
        correlation = lm.loadmat(os.path.join(path, p))['surrogateStateInfo']['correlation']
        for i in range(len(correlation)):  #for 100 subrogates 
            ci =correlation[i].correlationDelta 
            sumasub_positivos= np.zeros((correlation.shape[0],ci.shape[1])) #se va creando para cada paciente
                    
            for l in range(ci.shape[1]): 
                        #print(l)
                    pos = []
                    e = ci[:,l]
                        #print(e)
                    for elem in e :
                        if elem >=0:
                            pos.append(elem)
                            sumasub_positivos[i,l] = np.sum(pos)
                            
        
            sumaspacientes_totales[g].append(sumasub_positivos)
            media[g].append(np.mean(sumasub_positivos))
            
        señales_base, binwidth = sumas_positivos_base(k)
        
            
    return media, señales_base, binwidth
       
#Para hacer los histogramas, voy a coger el valor optimo de binwidth de los positivos


