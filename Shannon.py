#ENTROPIA DE SHANNON CON NUEVO ALFABETO
#1. Cambiar las nuevas matrices con alfabeto
#Debería elimar las rows de correlación que no son ni 1,2,3
#2. Calcular la probabilidad de ocurrencia
#3. Entropia de Shannon

import numpy as np
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import names_per_groups as nm #de señales base
import seaborn as sns
from scipy.stats import kruskal, pearsonr
from scipy.stats import mannwhitneyu
import statsmodels.api as sm
import astropy.stats as sts
from natsort import natsorted
from collections import Counter
import matplotlib.ticker as ticker
import pickle
import itertools
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import skew, kurtosis

sys.path.insert(0, 'C:/Users/miria/OneDrive/Documentos/Articulos UVA/data/Resultados senales base')
import loadmats as lm
import scipy.io as sio
path = "C:/Users/miria/OneDrive/Documentos/Articulos UVA/data/Resultados senales base/mats"
pathsub = "C:/Users/miria/OneDrive/Documentos/Articulos UVA/data/Subrogadas"

grupos = nm.name_per_groups() #Separation per groups

##PARA ABRIR UN PICKLE
#pickle_in = open("delta_sub_prob","rb")
#example_dict = pickle.load(pickle_in)

##PARA GUARDAR UN PICKLE
#pickle_file = open('alphabet_alpha', 'wb')
#pickle.dump(alphabet_alpha, pickle_file)

def alphabet_change(index_band, b):

    grupos = nm.name_per_groups()
    new_alphabet = {'control':[], 'MCI':[], 'AD':[] }
    prob_ocurrencia = {'control':[], 'MCI':[], 'AD':[] }
    prob_sinnormalizar =  {'control':[], 'MCI':[], 'AD':[] }
    tas_numbers= [1,2,3]
    numeros = [10,11,12,13,14,15,20,21,22,23,24,25,30,31,32,33,34,35]

    
    for key in grupos.keys():
        
        alphabet_mat= np.zeros((len(grupos[key]), 12000))
        p = np.zeros((len(grupos[key]), len(numeros)))
        
        for l in range(len(grupos[key])):
                       
            correlation_pat = lm.loadmat(os.path.join(path, grupos[key][l]))['dynamicStateInfo']['correlation'][b]
            tas= lm.loadmat(os.path.join(path, grupos[key][l]))['dynamicStateInfo']['temporalActivation'][index_band]
            #p = np.zeros((1, len(numeros)))
            position_corr = []
           
            for n in tas_numbers:
                #Find the column where the metastates are dominant
                position_tas= np.where(n == tas)[0][0]
                #Find the row where that metastate is
                position_corr.append(np.where(np.max(correlation_pat[:,position_tas])==correlation_pat[:,position_tas])[0][0])
    
            #If present, delete meta-states 4 and 5
            #Create the total list of meta-states present
            total_states = list(range(0, correlation_pat.shape[0]))
            
            #Substract the actual meta-states lists
            missing = list(Counter(total_states)-Counter(position_corr).keys())
            
            #Cheeck if missing is not empty:
            if len(missing) > 0:
                #delete the row of the correlation matrix that are not going to be used
                correlation = np.delete(correlation_pat, missing, 0)
            else:
                correlation = correlation_pat
            #If NOT 0, please continue
            
            #print(position_corr)
            #For each column:
            for c in range(correlation.shape[1]):
               
                #1.Obtain the columns
                column = natsorted(correlation[:,c])[::-1]
                #2. Find the position of the second and third values
                row1_c = np.where(column[0] == correlation_pat[:,c])[0][0] #Para evitar problemas con el tas
                row2_c = np.where(column[1] == correlation_pat[:,c])[0][0]
                row3_c = np.where(column[2] == correlation_pat[:,c])[0][0]
                #print(row2_c)
                #Find the tas position where it is.The meta-state is +1 the index
                metastado_1 = position_corr.index(row1_c) + 1 
                metastado_2 = position_corr.index(row2_c) + 1 
                metastado_3 = position_corr.index(row3_c) + 1 
                
                combinacion = [metastado_1, metastado_2, metastado_3]
                #print(row2_c)
                #print(combinacion)
                
                #FOR TAS 1
                if metastado_1 == 1:
                    
                    #Second and Third States are positive
                    if column[1]>0 and column [2]>0:
                        #If second position is state 2
                        if combinacion[1] == 2:
                            alphabet = 10 #123
                                                    
                        #If thrid position is state 3
                        elif combinacion[1] == 3:
                            alphabet = 11 #132
                    
                    #If Second and third value are negative       
                    elif column[1]<0 and column [2]<0:
                        if combinacion[1] == 2:
                            alphabet = 14 # 1-2-3
                        
                        elif combinacion[1] == 3:
                            alphabet = 15# 1-3-2
                            
                    #Third position is negative
                    elif column[2]<0:
                        if combinacion[2] == 2:
                            alphabet = 13 #13-2
                            
                        else: 
                            alphabet = 12 #12-3
                            
                    
                #FOR TAS 2    
                elif metastado_1 == 2:
                    
                    #Second and Third States are positive
                    if column[1]>0 and column [2]>0:
                        #If second position is state 2
                        if combinacion[1] == 1:
                            alphabet = 20 #213
                            
                        #If thrid position is state 3
                        elif combinacion[1] == 3:
                            alphabet = 21 #231
                    
                    #If Second and third value are negative       
                    elif column[1]<0 and column [2]<0:
                        if combinacion[1] == 1:
                            alphabet = 24 # 2-1-3
                            
                        elif combinacion[1] == 3:
                            alphabet = 25# 2-3-1
                    
                    #Third position is negative
                    elif column[2]<0:
                        if combinacion[2] == 1:
                            alphabet = 23 #23-1
                            
                        else: 
                            alphabet = 22 #21-3
                            
                            
                #TAS 3    
                else:
                    #Second and Third States are positive
                    if column[1]>0 and column [2]>0:
                        #If second position is state 1
                        if combinacion[1] == 1:
                            alphabet = 30 #312
                            
                        #If thrid position is state 2
                        elif combinacion[1] == 2:
                            alphabet = 31 #321
                                                
                    #If Second and third value are negative       
                    elif column[1]<0 and column [2]<0:
                        if combinacion[1] == 1:
                            alphabet = 34 #3-1-2
                            
                        elif combinacion[1] == 2:
                            alphabet = 35 #3-2-1
                        
                    #Third position is negative
                    elif column[2]<0:
                        if combinacion[2] == 1:
                            alphabet = 33 #32-1
                            
                        else: 
                            alphabet = 32 #31-2
                            
                
                alphabet_mat[l,c] = alphabet

                             
                for e in range(len(numeros)):
                    if numeros[e] == alphabet_mat[l,c]:
                        p[l,e] +=1
               
        new_alphabet[key].append(alphabet_mat)
        prob_ocurrencia[key].append(p/correlation.shape[1])
        prob_sinnormalizar[key].append(p)
        #print(key)
        #print(len(grupos[key]))

    return new_alphabet, prob_ocurrencia #, prob_sinnormalizar

#PROB DE OCURRENCIA
#Valores: 10, 11, 12, 13, 14, 15, 16, 17, 20-27, 20-27

#alphabet_alpha, prob_alpha = alphabet_change(index_band=2, b='correlationAlpha')

#ENTROPÍA DE SHANNON
# Mientras menos probable, más sorpresa y más información contiene.

#H(X) = - sum(n=(i=1)) p(x1)log2p(x1)
#def entropia_shannon
#PARA UN SOLO SUJETO DE UN ÚNICO GRUPO, HABRÍA QUE HACER LA MEDIA
def shannon(prob):
    entropy= {'control':[], 'MCI':[], 'AD':[] }
    for k in grupos.keys():
        numeros = [10,11,12,13,14,15,20,21,22,23,24,25,30,31,32,33,34,35]
        for l in range(len(prob[k][0])):
          
            h=0
            for i in range(len(numeros)):
               
                if prob[k][0][l][i] != 0:
                    multiplication = prob[k][0][l][i]*np.log(prob[k][0][l][i])
                    #print(prob[k][0][l][i])
                    h += multiplication
                else:
                    h += 0
                #print(h)
            entropy[k].append(-h/np.log(len(numeros)))
        
    return entropy

#PROBABILIDADES

alphabet_delta, prob_delta = alphabet_change(index_band=0, b='correlationDelta')
#alphabet_theta, prob_theta = alphabet_change(index_band=1, b='correlationTheta')
#alphabet_alpha, prob_alpha = alphabet_change(index_band=2, b='correlationAlpha')
#alphabet_beta1, prob_beta1 = alphabet_change(index_band=3, b='correlationBeta1')
#alphabet_beta2, prob_beta2 = alphabet_change(index_band=4, b='correlationBeta2')



#SHANNON

#shannon_delta = shannon(prob = prob_delta)
#shannon_theta = shannon(prob = prob_theta)
#shannon_alpha = shannon(prob = prob_alpha)
#shannon_beta1 = shannon(prob = prob_beta1)
#shannon_beta2 = shannon(prob = prob_beta2)


#MAPAS DE CALOR SIN HIST: Hacer una lista con todos los valores de la banda
'''
dc = list(np.concatenate(prob_delta['control'][0]))
dm= list(np.concatenate(prob_delta['MCI'][0]))
dad= list(np.concatenate(prob_delta['AD'][0]))
lista_delta = list(itertools.chain(dad, dc, dm))
dmax= np.max(lista_delta)
dmin=np.min(lista_delta)

#MAPA DE CALOR 
plt.figure()
cmap= sns.color_palette("viridis", as_cmap=True)
ax= sns.heatmap(prob_alfa['MCI'][0], yticklabels=False, xticklabels='auto', cmap=cmap, vmin=dmin, vmax=dmax, cbar_kws={'label': 'Probabilidad de ocurrencia'})
#ax.figure.axes[-1].yaxis.label.set_size(16)
plt.xticks(ticks= [i for i in range(18)], labels= combi_sin,rotation=75,fontsize=14) 
plt.title("DCL",fontsize=21)
plt.xlabel("Patrones Simbólicos", fontsize=16)
plt.ylabel("Sujetos", fontsize=16)
plt.show()


##DISTRIBUTION GRAFICA 
yticks= np.arange(0, 0.16, step=0.02)

plt.figure(figsize=(15, 6))
#sns.set_theme()
plt.plot(x,distri_ac, label='Control', marker="o", linestyle="-")
plt.plot(x, distri_am, label='DCL', marker="o", linestyle="-")
plt.plot(x, distri_aad, label='EA', marker="o", linestyle="-")
plt.xlabel("Patrones Simbólicos", fontsize=18)
plt.ylabel("Probabilidad de ocurrencia", fontsize=16)
plt.xticks(ticks= [i for i in range(18)], labels= combi_sin, fontsize=15) 
plt.yticks(ticks=yticks, fontsize=15)
plt.legend()
plt.show()

#Mapas de Calor
plt.figure()
cmap= sns.color_palette("viridis", as_cmap=True)
sns.heatmap(prob_delta['AD'][0], yticklabels=False, xticklabels='auto', cmap=cmap, vmin=dmin, vmax=dmax)
plt.title("AD Delta")
plt.xticks(rotation = 75) 
plt.xlabel("Alphabet combinations")
plt.ylabel("Subjects")
plt.show()
'''
'''
###FUNCION DE DISTRIBUCIÓN 
distri_tc = np.apply_along_axis(sum, 0, prob_theta['control'][0])/prob_theta['control'][0].shape[0]
distri_tm = np.apply_along_axis(sum, 0, prob_theta['MCI'][0])/prob_theta['MCI'][0].shape[0]
distri_tad = np.apply_along_axis(sum, 0, prob_theta['AD'][0])/prob_theta['AD'][0].shape[0]


combi_sin =['[123]','[132]','[12-3]','[13-2]','[1-2-3]','[1-3-2]','[213]','[231]','[21-3]','[23-1]','[2-1-3]','[2-3-1]','[312]','[321]','[31-2]','[3-21]','[3-1-2]','[3-2-1]']

#GRUPAL
plt.figure(figsize=(15, 6))
plt.plot(distri_dc, label='control')
plt.plot(distri_dm, label='MCI')
plt.plot(distri_dad, label='AD')
plt.title('Distribution Functions - Delta')
plt.xlabel("Alphabet combination")
plt.ylabel("Ocurrence probability")
plt.xticks([i for i in range(18)]) 
plt.xlim([0, 17])
plt.legend()
plt.show()

#INDIVIDUAL
plt.figure(figsize=(13, 3))
plt.plot(distri_dc)
plt.title('Distribution Function Control - Delta')
plt.xlabel("Alphabet combination")
plt.ylabel("Ocurrence probability")
plt.xticks([i for i in range(18)]) 
plt.xlim([0, 17])
plt.ylim([0, 0.15])
plt.show()

plt.figure(figsize=(13, 3))
plt.plot(distri_tm, color='orange')
plt.title('Distribution Function MCI - Theta')
plt.xlabel("Alphabet combination")
plt.ylabel("Ocurrence probability")
plt.xticks([i for i in range(18)]) 
plt.xlim([0, 17])
plt.ylim([0, 0.15])
plt.show()


plt.figure(figsize=(13, 3))
plt.plot(distri_tad, color='green')
plt.title('Distribution Function AD - Theta')
plt.xlabel("Alphabet combination")
plt.ylabel("Ocurrence probability")
plt.xticks([i for i in range(18)]) 
plt.xlim([0, 17])
plt.ylim([0, 0.15])
plt.show()


norm_shannon_delta= {'control':[], 'MCI':[], 'AD':[] }
for k in norm_shannon_delta.keys():
    norm_shannon_delta[k] = np.divide(shannon_delta['k'][0],sub_shannon_delta['k'][0])

'''

def new_alphabet_sub(index_band):
    
    #print('Hola estoy en Delta')
    grupos = nm.name_per_groups()
    numeros = [10,11,12,13,14,15,20,21,22,23,24,25,30,31,32,33,34,35]
    tas_numbers= [1,2,3]
    prob_ocurrencia_sub = {'control':[], 'MCI':[], 'AD':[] }
    
    for k in grupos.keys():
        
        #alphabet_mat= np.zeros((len(grupos[k]), 12000))
        p = np.zeros((len(grupos[k]), len(numeros)))
        
        for l in range(len(grupos[k])):
            
            correlation_sub = lm.loadmat(os.path.join(pathsub, grupos[k][l]))['surrogateStateInfo']['correlation']
            
            #SUBROGADAS
            for i in range(len(correlation_sub)): #100 subrogadas 
            
                
                #Load tas for each surrogate
                tas_sub= lm.loadmat(os.path.join(pathsub, grupos[k][l]))['surrogateStateInfo']['temporalActivation'][i]
                #Load each surrogate 
                surrogate_corr_band = correlation_sub[i].correlationDelta
                
                #Matrix for 100 surrogates to add the alphabet and from this calculate the prob 
                alphabet_mat_temp= np.zeros((len(correlation_sub), surrogate_corr_band.shape[1]))
                
                #Matrix for the prob of ocurrence temporal
                p_temp = np.zeros((len(correlation_sub), len(numeros)))
                
                #p = np.zeros((1, len(numeros)))
                position_corr = []
               
                for n in tas_numbers:
                    #Find the column where the metastates are dominant
                    position_tas= np.where(n == tas_sub[index_band])[0][0]
                    #Find the row where that metastate is
                    position_corr.append(np.where(np.max(surrogate_corr_band[:,position_tas])==surrogate_corr_band[:,position_tas])[0][0])
        
                #If present, delete meta-states 4 and 5
                #Create the total list of meta-states present
                total_states = list(range(0, surrogate_corr_band.shape[0]))
                
                #Substract the actual meta-states lists
                missing = list(Counter(total_states)-Counter(position_corr).keys())
                
                #Cheeck if missing is not empty:
                if len(missing) > 0:
                    #delete the row of the correlation matrix that are not going to be used
                    correlation = np.delete(surrogate_corr_band, missing, 0)
                else:
                    correlation = surrogate_corr_band
                
                #For each column:
                for c in range(correlation.shape[1]):
               
                    #1.Obtain the columns
                    column = natsorted(correlation[:,c])[::-1]
                    #2. Find the position of the second and third values
                    row1_c = np.where(column[0] == surrogate_corr_band[:,c])[0][0] #Para evitar problemas con el tas
                    row2_c = np.where(column[1] == surrogate_corr_band[:,c])[0][0]
                    row3_c = np.where(column[2] == surrogate_corr_band[:,c])[0][0]
                    #print(row2_c)
                    #Find the tas position where it is.The meta-state is +1 the index
                    metastado_1 = position_corr.index(row1_c) + 1 
                    metastado_2 = position_corr.index(row2_c) + 1 
                    metastado_3 = position_corr.index(row3_c) + 1 
                
                    combinacion = [metastado_1, metastado_2, metastado_3]
                
                    #FOR TAS 1
                    if metastado_1 == 1:
                        
                        
                        #Second and Third States are positive
                        if column[1]>0 and column [2]>0:
                            #If second position is state 2
                            if combinacion[1] == 2:
                                alphabet = 10 #123
                                                        
                            #If thrid position is state 3
                            elif combinacion[1] == 3:
                                alphabet = 11 #132
                        
                        #If Second and third value are negative       
                        elif column[1]<0 and column [2]<0:
                            if combinacion[1] == 2:
                                alphabet = 14 # 1-2-3
                            
                            elif combinacion[1] == 3:
                                alphabet = 15# 1-3-2
                                
                        #Third position is negative
                        elif column[2]<0:
                            if combinacion[2] == 2:
                                alphabet = 13 #13-2
                                
                            else: 
                                alphabet = 12 #12-3
                                
                    #FOR TAS 2    
                    elif metastado_1 == 2:
                        
                        
                        #Second and Third States are positive
                        if column[1]>0 and column [2]>0:
                            #If second position is state 2
                            if combinacion[1] == 1:
                                alphabet = 20 #213
                                
                            #If thrid position is state 3
                            elif combinacion[1] == 3:
                                alphabet = 21 #231
                        
                        #If Second and third value are negative       
                        elif column[1]<0 and column [2]<0:
                            if combinacion[1] == 1:
                                alphabet = 24 # 2-1-3
                                
                            elif combinacion[1] == 3:
                                alphabet = 25# 2-3-1
                        
                        #Third position is negative
                        elif column[2]<0:
                            if combinacion[2] == 1:
                                alphabet = 23 #23-1
                                
                            else: 
                                alphabet = 22 #21-3
                                
                    #FOR TAS 3    
                    else:
                        
                        #Second and Third States are positive
                        if column[1]>0 and column [2]>0:
                            #If second position is state 1
                            if combinacion[1] == 1:
                                alphabet = 30 #312
                                
                            #If thrid position is state 2
                            elif combinacion[1] == 2:
                                alphabet = 31 #321
                                                    
                        #If Second and third value are negative       
                        elif column[1]<0 and column [2]<0:
                            if combinacion[1] == 1:
                                alphabet = 34 #3-1-2
                                
                            elif combinacion[1] == 2:
                                alphabet = 35 #3-2-1
                            
                        #Third position is negative
                        elif column[2]<0:
                            if combinacion[2] == 1:
                                alphabet = 33 #32-1
                                
                            else: 
                                alphabet = 32 #31-2
                    
                    #Add the alphabet for each 100 surrogate of each sujeto
                    alphabet_mat_temp[i,c] = alphabet
            #print(alphabet_mat_temp)
                    
                    #Calculate the prob of ocurrence for the 100 surrogates
                    for n in range(len(numeros)):
                        if numeros[n] == alphabet_mat_temp[i,c]:
                            p_temp[i,n] +=1
                
                #Obtain the probabilities of each combination for each subrogate
                p_norm = p_temp/correlation.shape[1]
                #Do the mean of the probs 
                p[l]=np.mean(p_norm, axis=0)
                
                
        #Add to the dict of that group
        prob_ocurrencia_sub[k].append(p)
        

    return prob_ocurrencia_sub



#REMEMBER TO CHANGE THE BAND INSIDE DE CODE!!!
#sub_prob_delta = new_alphabet_sub(index_band = 0)
#sub_prob_theta = new_alphabet_sub(index_band = 1)
#sub_prob_alpha = new_alphabet_sub(index_band = 2)
#sub_prob_beta1 = new_alphabet_sub(index_band = 3)
#sub_prob_beta2 = new_alphabet_sub(index_band = 4)

