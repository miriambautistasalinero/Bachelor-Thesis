#Mapas de calor para correlaciones NO dominantes por bandas
import loadmats as lm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import seaborn as sns
import astropy.stats as sts
from scipy.stats import kurtosis, skew
import pandas as pd
from scipy.stats import kruskal
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import fdrcorrection

#SI QUIERES CAMBIAR DE BANDA
#1.CAMBIA k
#2.CAMBIA tasp
#3.CAMBIA LOS MULTIPLOS: Delta, Alpha, Beta1, =12;Beta2 ,Theta 11


k='correlationAlpha'

pacientesbase = os.listdir("C:/Users/miria/OneDrive/Documentos/Articulos UVA/data/Resultados senales base/mats")
grupo = {'control':[], 'MCI':[], 'AD':[] }
tas = {'control':[], 'MCI':[], 'AD':[] }
correlation = {'control':[], 'MCI':[], 'AD':[] }
corr_ordenadas = {'control':[], 'MCI':[], 'AD':[] }
suma_nodom= {'control':[], 'MCI':[], 'AD':[] }
for p in pacientesbase:
    
    dsi= lm.loadmat(p)['dynamicStateInfo']
    group =dsi['group']
    tasp = dsi['temporalActivation'][2,:] #El 0 es porque estoy en banda DELTA
    correlationp=dsi['correlation'][k]
    #corr_ordenadas= np.zeros(correlation.shape)
    elemintas = list(set(tasp))

    if group == 'DCL1' or group =='DCL2' or group == 'DCL':
            grupo["MCI"].append(p)
            tas['MCI'].append(tasp)
            correlation['MCI'].append(correlationp)
            nodominantes = np.delete(np.sort(correlationp, axis=0), elemintas[-2], axis=0) 
            suma_nodom['MCI'].append( np.apply_along_axis(sum, 0, abs(nodominantes)))
            
    
    elif group == 'EA1' or group == 'EA2':
        grupo["AD"].append(p)
        tas['AD'].append(tasp)
        correlation['AD'].append(correlationp)
        nodominantes = np.delete(np.sort(correlationp, axis=0), elemintas[-2], axis=0) 
        suma_nodom['AD'].append( np.apply_along_axis(sum, 0, abs(nodominantes)))
        
    else:
        grupo["control"].append(p)
        tas['control'].append(tasp)
        correlation['control'].append(correlationp)
        nodominantes = np.delete(np.sort(correlationp, axis=0), elemintas[-2], axis=0) 
        suma_nodom['control'].append( np.apply_along_axis(sum, 0, abs(nodominantes)))


#Propiedades
skw ={'control':[], 'MCI':[], 'AD':[] }
kurt = {'control':[], 'MCI':[], 'AD':[] }
std ={'control':[], 'MCI':[], 'AD':[] }
mean={'control':[], 'MCI':[], 'AD':[] }
for k in suma_nodom.keys():
    for ii in range(len(suma_nodom[k])):
        skw[k].append(skew(suma_nodom[k][ii]))
        #skw[k] = np.mean(skew(suma_nodom[k]))
        kurt[k].append(kurtosis(suma_nodom[k][ii]))
        #kurt[k] = np.mean(kurtosis(suma_nodom[k]))
        mean[k].append(np.mean(suma_nodom[k][ii]))
        std[k].append(np.std((suma_nodom[k][ii])))
    
    
         
#Hacer lista con todos los valores de las correlaciones
lista1=[]
for kk in suma_nodom:
    lista1.append(suma_nodom[kk])
todas_nodom=np.concatenate(lista1).flatten()
        
##Hallar el binwidth para todos los valores
opt_width_scott_corre, bins_scott_corre = sts.scott_bin_width(todas_nodom, return_bins=True)
b=np.arange(np.min(todas_nodom), np.max(todas_nodom) + opt_width_scott_corre, opt_width_scott_corre)  

#Hacer los histogramas para cada grupo   
m=[]
hisporgrupos = {'control':np.zeros((len(grupo["control"]),len(b)-1)), 'MCI': np.zeros((len(grupo["MCI"]),len(b)-1)), 'AD': np.zeros((len(grupo["AD"]),len(b)-1)) }     
for cx in suma_nodom.keys():
    for pt in range(len(grupo[cx])):
        n,bgrupo=np.histogram(suma_nodom[cx][pt],bins=b)
        #print(np.sum(n))
        hisporgrupos[cx][pt]=n #/np.sum(n)
    #hisporgrupos[cx] = hisporgrupos[cx]/np.sum(hisporgrupos[cx])
    m.append(hisporgrupos[cx].max())
maxi = min(m) #cojo el más restrictivo
#Normalización del histograma dividiendo entre el numero total de elementos

#Disposión de los números de correlaciones. 
#ESTE VALOR DEPENDE DEL TAMAÑO DE LOS PLOTS!!!!!!!!
#PARA COMPROBARLO HAY QUE HACER UNA GRÁFICA SIN XTICKS
multiplos=[valor for valor in range(len(b)) if valor%20==0]
g=[]
for elem in multiplos:
    g.append(np.round(b[elem], decimals=3))

xticks= ticker.FixedLocator(multiplos)    
#Mapas de Calor

#Ensure same color scale

#plotear las propiedades
df_medias= pd.DataFrame(dict([(h, pd.Series(j)) for h,j in mean.items()]))
df_kurt= pd.DataFrame(dict([(h, pd.Series(j)) for h,j in kurt.items()]))
df_std= pd.DataFrame(dict([(h, pd.Series(j)) for h,j in std.items()]))
df_skw= pd.DataFrame(dict([(h, pd.Series(j)) for h,j in skw.items()]))

plt.figure(figsize=(13,8))
plt.suptitle("Properties of sum of non-dominant correlations - Alpha", fontsize=13)
sns.set_theme(style="whitegrid")

plt.subplot(2,2,1)
sns.violinplot(data=df_medias)
plt.title('Mean')

plt.subplot(2,2,2)
sns.violinplot(data=df_std)
plt.title('Standard Deviation')

plt.subplot(2,2,3)
sns.violinplot(data=df_kurt)
plt.title('Kurtosis')

plt.subplot(2,2,4)
sns.violinplot(data=df_skw)
plt.title('Skewness')

plt.show()

#Mapas de Calor
plt.figure(figsize=(12,4))
cmap= sns.color_palette("viridis", as_cmap=True)
xticks= ticker.FixedLocator(multiplos)
formater = ticker.FixedFormatter(g)
#plt.suptitle("Sum of non-dominant correlations - Alpha", fontsize=20)
plt.xlabel('Correlations')
plt.xticks(rotation = 75) 

ax1= plt.subplot(1,3,1)
sns.heatmap(hisporgrupos['control'], yticklabels=False, xticklabels='auto', cmap=cmap, vmin=0, vmax=maxi)
plt.title('Control')
ax1.xaxis.set_major_locator(xticks)
ax1.xaxis.set_major_formatter(formater)
plt.ylabel('Subjects')
plt.xticks(rotation = 75) 

ax2= plt.subplot(1,3,2)
sns.heatmap(hisporgrupos['MCI'], yticklabels=False, xticklabels='auto', cmap=cmap, vmin=0, vmax=maxi)
plt.title('MCI')
ax2.xaxis.set_major_locator(xticks)
ax2.xaxis.set_major_formatter(formater)
plt.xticks(rotation = 75) 
plt.xlabel('Correlations')

ax3=plt.subplot(1,3,3)
sns.heatmap(hisporgrupos['AD'], yticklabels=False, xticklabels='auto', cmap=cmap, vmin=0, vmax=maxi)
plt.title('AD')
ax3.xaxis.set_major_locator(xticks)
ax3.xaxis.set_major_formatter(formater)
plt.xticks(rotation = 75) 

plt.tight_layout()
plt.show()
