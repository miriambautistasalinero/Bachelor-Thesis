#FACTOR DE ANTAGONISMO 

import numpy as np
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import names_per_groups as nm #de señales base
import seaborn as sns
import time
from datetime import datetime
from scipy.stats import kruskal
from scipy.stats import mannwhitneyu
# df.rename(columns={"Gender":"MCI"})



#Insert path for load mat modules
sys.path.insert(0, 'C:/Users/miria/OneDrive/Documentos/Articulos UVA/data/Resultados senales base')

import loadmats as lm

path = "C:/Users/miria/OneDrive/Documentos/Articulos UVA/data/Resultados senales base/mats"
grupos = nm.name_per_groups() #Separation per groups

#b='correlationAlpha'

#NORMALES
def factor_antagonismo_base(b):
    n =  {'control':[], 'MCI':[], 'AD':[] }
    
    for key in grupos.keys():

        correlaciones_temporales = np.zeros((len(grupos[key]),12000))
        
        for l in range(len(grupos[key])): #number of patients in each group
            correlation = lm.loadmat(os.path.join(path, grupos[key][l]))['dynamicStateInfo']['correlation'][b]
            
            for i in range(correlation.shape[1]):
                correlaciones_total=0
                correlaciones_negativas = 0
                for c in correlation[:,i]:
                    correlaciones_total += abs(c)
                    
                    if c < 0:
                        correlaciones_negativas += c
                #En correlaciones_temporales acabo teniendo el resultado de todos los pacientes en cada instante
                correlaciones_temporales[l,i]= abs(correlaciones_negativas)/correlaciones_total
        
        #To have one value per patient. suma los 1200 valores de cada paciente y lo divide entre 12000
        #p = np.apply_along_axis(sum, 1, correlaciones_temporales)/correlaciones_temporales.shape[1] 
        #print(p)
        n[key]= np.apply_along_axis(sum, 1, correlaciones_temporales)/correlaciones_temporales.shape[1] 
       
    return n
    

#SUBROOGADAS
pathsub = "C:/Users/miria/OneDrive/Documentos/Articulos UVA/data/Subrogadas"


start = time.time()
def factor_antagonismo_sub():
    
    sub_medias_factorant = {'control':[], 'MCI':[], 'AD':[] }
    for key in grupos.keys(): #control, DCL, EA
   
        for l in range(len(grupos[key])): #para la longitud de cada grupo 
        
            correlation = lm.loadmat(os.path.join(pathsub, grupos[key][l]))['surrogateStateInfo']['correlation']
            correlaciones_temporales_100 = np.zeros((100,12000))
            
            for s in range(len(correlation)): #me meto en las 100 subrogadas
                
                ci = correlation[s].correlationAlpha #cargo las 100 subrogadas de cada banda
                
                for i in range(ci.shape[1]): #cada columna de la matrix 12000
                    
                    correlaciones_total=0 
                    correlaciones_negativas = 0
                    
                    for c in ci[:,i]:#dentro de todos los elementos de una columna
                        
                        
                        correlaciones_total += abs(c) #suma cada elemento en valor absoluto
                        
                        if c < 0:
                            correlaciones_negativas += c
                    
                    correlaciones_temporales_100[s,i]= abs(correlaciones_negativas)/correlaciones_total
            
            #Es lo mismo hacer el np.mean de toda una matriz que hacer la suma por instantes y luego hacer la media
            sub_medias_factorant[key].append(np.mean(correlaciones_temporales_100))
       
    return sub_medias_factorant
                
#fac_base_delta= factor_antagonismo_base(b='correlationDelta')
#fac_base_theta= factor_antagonismo_base(b='correlationTheta')
#fac_base_alpha= factor_antagonismo_base(b='correlationAlpha')
#fac_base_beta1= factor_antagonismo_base(b='correlationBeta1')
#fac_base_beta2= factor_antagonismo_base(b='correlationBeta2')

#antafact_sub_delta= factor_antagonismo_sub()
#antafact_sub_theta= factor_antagonismo_sub()
antafact_sub_alpha= factor_antagonismo_sub()

'''
negativos_total_bases = factor_antagonismo_base(b)

antafact_sub= factor_antagonismo_sub()

#Normalización

norm = {'control':[], 'MCI':[], 'AD':[] }
for hd in norm.keys():
    norm[hd] = np.divide(negativos_total_bases[hd],antafact_sub[hd])
    
krusk = kruskal(norm['control'],norm['DCL'],norm['EA'])

con_dcl= mannwhitneyu(norm['control'],norm['DCL'])
con_ea= mannwhitneyu(norm['control'],norm['EA'])
dcl_ea= mannwhitneyu(norm['DCL'],norm['EA'])


'''
'''
#PLOTEAR CON VIOLIN PLOTS BASES
df_bases = pd.DataFrame(dict([(h, pd.Series(j)) for h,j in negativos_total_bases.items()]))
plt.figure()
plt.title("Factor de Antagonismo Bases - Delta")
sns.set_theme(style="whitegrid")
ax = sns.violinplot(x="bandas", y="values", hue="grupos", data=df_bueno3)
ax1.axhline(1, ls='--, color='r'')
#ax = sns.violinplot(data=df_bases)
#plt.ylim(-0.05,0.4)
#sns.boxplot(data=df_medias)
plt.show()

#PLOTEAR CON VIOLIN PLOTS SUBROAGADAS
df_norm = pd.DataFrame(dict([(h, pd.Series(j)) for h,j in norm.items()]))
plt.figure()
plt.title("Factor de Antagonismo Normalizados - Delta")
sns.set_theme(style="whitegrid")
ax = sns.violinplot(data=df_norm)
plt.ylim((0.3,2))
#sns.boxplot(data=df_medias)
plt.show()

#PLOTEAR TODO JUNTO
plt.figure(figsize=[15,5])
plt.title("Degree of Antagonism")
sns.set_theme(style="whitegrid")
ax = sns.violinplot(x="Bands", y="Valores", hue="Grupos", data=df_norm)
#ax = sns.violinplot(data=df_bases)
ax.axhline(1, ls='--', color="red")
#plt.ylim(-0.05,0.4)
#sns.boxplot(data=df_medias)
plt.show()

'''
