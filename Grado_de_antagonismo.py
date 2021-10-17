#GRADO DE ANTAGONISMO
#CORR.POSITIVA MAYO - CORR. NEGATIVA MAYOR DE VALOR ABSOLUTO.

import numpy as np
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import names_per_groups as nm #de señales base
import seaborn as sns
from scipy.stats import kruskal
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import fdrcorrection

#Insert path for load mat modules
sys.path.insert(0, 'C:/Users/miria/OneDrive/Documentos/Articulos UVA/data/Resultados senales base')


import loadmats as lm

#b = 'correlationAlpha' #Banda de frequencia

def antagonismo_bases(b):
    
    path = "C:/Users/miria/OneDrive/Documentos/Articulos UVA/data/Resultados senales base/mats"
    grupos = nm.name_per_groups() #Separation per groups
    
    
    grado_total = {'control':[], 'MCI':[], 'AD':[] }
    prob_ocurrencia= {'control':[], 'MCI':[], 'AD':[] }
    
    
    for key in grupos.keys():
    #key = 'control'
        grados = np.zeros((len(grupos[key]),12000))
        
        for l in range(len(grupos[key])):
            count=0
            correlation = lm.loadmat(os.path.join(path, grupos[key][l]))['dynamicStateInfo']['correlation'][b]
            
            #print(grados.shape)
            for i in range(correlation.shape[1]): #1200
               
                maxi_positivo = max(correlation[:,i])
                min_negativo = min(correlation[:,i])
                
                if min_negativo >= 0:
                    min_negativo = 0
                    
                grados[l,i] = maxi_positivo - min_negativo
                #print(maxi_positivo)
                #print(min_negativo)
                if abs(maxi_positivo) < abs(min_negativo):
                    count += 1
            #indices.append(count/12000)
            
            prob_ocurrencia[key].append(count/12000)
            
            grado_total[key].append(np.sum(grados[l,:]) / grados.shape[1])  
                               
        #indices[key] = vector_indices
    return  grado_total, prob_ocurrencia

#las subrogadas hay que calcularlas por grupos para que no se sature            
#def antagomismo_subrogadas()

pathsub = "C:/Users/miria/OneDrive/Documentos/Articulos UVA/data/Subrogadas"
grupos = nm.name_per_groups() #Separation per groups

def ant_sub():
    sub_medias_grados = {'control':[], 'MCI':[], 'AD':[] }
    prob_ocurrencia_sub= {'control':[], 'MCI':[], 'AD':[] }
    #np.zeros((len(grupos['control']),1))
    #for key in grupos:
    #key = 'control'
    for key in grupos.keys():
        for l in range(len(grupos[key])): #range(1)
             #Load correlation
            correlation = lm.loadmat(os.path.join(pathsub, grupos[key][l]))['surrogateStateInfo']['correlation']
            
            #Final matrix for grados ant
            grados_sub= np.zeros((correlation.shape[0],12000))
            
            #Matrix for indices 100x12000
            #vector_indices= np.zeros((correlation.shape[0],12000))
            indices = []
            
            #SUBROGADAS
            for i in range(len(correlation)): #100 subrogadas
                count=0
                #Load each of the subrrogates for a specific band
                ci = correlation[i].correlationAlpha
                
                #COLUMNS
                for j in range(ci.shape[1]):
                    
                    #Find the max and the min of that column
                    maxi_positivo_sub = max(ci[:,j])
                    min_negativo_sub = min(ci[:,j])
                    
                    #If no negativo correlation it is set to 0 
                    if min_negativo_sub >= 0:
                            min_negativo_sub = 0
                            
                    #Formula of grados sub
                    grados_sub[i,j] = maxi_positivo_sub - min_negativo_sub
                    
                    #Set the indices matrix 
                    if abs(maxi_positivo_sub) < abs(min_negativo_sub):
                        count += 1
             
                indices.append(count/12000)
                
            prob_ocurrencia_sub[key].append(np.mean(indices)) 
            sub_medias_grados[key].append(np.mean(grados_sub))
        
    return sub_medias_grados, prob_ocurrencia_sub

#grado_total_base, vindices = antagonismo_bases(b)
#anta_subr = ant_sub()

#grad_delta_base, prob_delta_base = antagonismo_bases('correlationDelta')
#grad_theta_base, prob_theta_base = antagonismo_bases('correlationTheta')
#grad_alpha_base, prob_alpha_base = antagonismo_bases('correlationAlpha')
#grad_beta1_base, prob_beta1_base = antagonismo_bases('correlationBeta1')
#grad_beta2_base, prob_beta2_base = antagonismo_bases('correlationBeta2')

#sub_delta, prob_delta_sub =  ant_sub()
#sub_theta, prob_theta_sub =  ant_sub()
sub_apha, prob_alpha_sub =  ant_sub()
#sub_beta1, prob_beta1_sub =  ant_sub()
#sub_beta2, prob_beta2_sub =  ant_sub()




#Normalización
'''
norm = {'control':[], 'MCI':[], 'AD':[] }
for hd in norm.keys():
    norm[hd] = np.divide(grado_total_base[hd],anta_subr[hd])

krusk = kruskal(norm['control'],norm['DCL'],norm['EA'])


con_dcl_beta1= mannwhitneyu(prob_norm_beta1['control'],prob_norm_beta1['DCL'])
con_ea_beta1= mannwhitneyu(prob_norm_beta1['control'],prob_norm_beta1['EA'])
dcl_ea_beta1= mannwhitneyu(prob_norm_beta1['DCL'],prob_norm_beta1['EA'])

pval_mann_beta1=[con_dcl_beta1[1], con_ea_beta1[1], dcl_ea_beta1[1]]

fdr_mann_beta1= fdrcorrection(pval_mann_beta1)



#LISTA DE GRUPOS
nombres_grupos = []
for k in grupos.keys():
    for i in range(len(grupos[k])):
        nombres_grupos.append(k)

#LISTA DE BANDAS
bandas = []
for k in grupos.keys():
    for ii in range(len(grupos[k])):
        bandas.append('Alpha')
for k in grupos.keys():
    for ii in range(len(grupos[k])):
        bandas.append('Beta1')

#VALORES DE LA MÉTRICA
valores_norm = []

for k in grupos.keys():
    for ii in range(len(grupos[k])):
        valores_norm.append(prob_norm[k][ii])
        
for k in grupos.keys():
    for ii in range(len(grupos[k])):
        valores_norm.append(prob_norm_beta1[k][ii])

dict_norm_alphabeta1= {'Bands': bandas, 'Probability_of_occurrence': valores_norm, 'Groups': nombres_grupos }
df_norm_prob = pd.DataFrame(dict([(h, pd.Series(j)) for h,j in dict_norm_alphabeta1.items()]))

'''
'''
#Plotear las bases
df_bases = pd.DataFrame(dict([(h, pd.Series(j)) for h,j in grado_total_base.items()]))
plt.figure()
plt.title("Grado de Antagonismo Bases - Delta")
sns.set_theme(style="whitegrid")
plt.plot(sns.violinplot(data=df_gradant_delta))
plt.plot(sns.violinplot(data=df_gradant_theta))
plt.ylim(0.3,1.5)
#sns.boxplot(data=df_medias)
plt.show()


plt.figure()
sns.heatmap(vindices['control'], yticklabels='auto', xticklabels='auto', cmap="viridis")
plt.show()

plt.figure()
sns.heatmap(vindices['DCL'], yticklabels='auto', xticklabels='auto', cmap="viridis")
plt.show()

plt.figure()
sns.heatmap(vindices['EA'], yticklabels='auto', xticklabels='auto', cmap="viridis")
plt.show()
'''

