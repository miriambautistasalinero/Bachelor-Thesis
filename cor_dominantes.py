import loadmats as lm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import seaborn as sns
import astropy.stats as sts

#k='correlationBeta2'
grupo = {'control':[], 'MCI':[], 'AD':[] }
tas={'control':[], 'MCI':[], 'AD':[] }
correlation={'control':[], 'MCI':[], 'AD':[] }
cxbanda={'control':[], 'MCI':[], 'AD':[] }

#p='091_B_EA_Norm.mat'
pacientesbase = os.listdir("C:/Users/miria/OneDrive/Documentos/Articulos UVA/data/Resultados senales base/mats")
todas_correlaciones=[]

for p in pacientesbase:

    dsi= lm.loadmat(p)['dynamicStateInfo']
    
    #tas = dsi['temporalActivation'][0,:]
    #correlation=dsi['correlation'][k]
    group =dsi['group']
    #todas_correlaciones.append(dsi['correlation'][k])
    dsicorre = dsi['correlation']
    #for k in dsicorre.keys():
     #   print(k)

    if group == 'DCL1' or group =='DCL2' or group == 'DCL':
            grupo["MCI"].append(p)
            #tas['DCL'].append(dsi['temporalActivation'][0,:]) #EL 0 ES POR ESTAR EN LA BANDA DELTA
            correlation['MCI'].append(dsi['correlation'])
            
    elif group == 'EA1' or group == 'EA2':
        grupo["AD"].append(p)
        #tas['EA'].append(dsi['temporalActivation'][0,:]) #EL 0 ES POR ESTAR EN LA BANDA DELTA
        correlation['AD'].append(dsi['correlation'])
    else:
        grupo["control"].append(p)
        #tas['control'].append(dsi['temporalActivation'][0,:]) #EL 0 ES POR ESTAR EN LA BANDA DELTA
        correlation['control'].append(dsi['correlation'])
        
        
#maximos={'control':np.zeros((len(grupo['control']),12000)), 'DCL':np.zeros((len(grupo['DCL']),12000)), 'EA':np.zeros((len(grupo['EA']),12000)) }

#To save all the dominant correaltions
pruebadm= {'control':[], 'MCI':[], 'AD':[] } 

k='correlationAlpha'
for g in grupo.keys(): 
    for ii in range(len(grupo[g])):
        #pruebadm= {'control':[], 'MCI':[], 'AD':[] }
        #for k in dsicorre.keys():
           
            
        m = np.max(correlation[g][ii][k], axis=0)
        pruebadm[g].append(abs(m))
        #Va sumando la correlacion al anterior pero no sé para qué
        cxbanda[g].append((np.apply_along_axis(sum, 0, pruebadm[g])))
           #for c in range(12000):
               #   correlation[g].append(max(dsi['correlation'][k][:,c]))



#SUBROGADAS


#MAPAS DE CALOR               
lista1=[]            
#Hacer lista con todos los valores de las correlaciones
for kk in cxbanda:
    lista1.append(pruebadm[kk]) #pruebadm[kk]
lista=np.concatenate(lista1).flatten()

#Hallar el binwidth para todos los valores
opt_width_scott_corre, bins_scott_corre = sts.scott_bin_width(lista, return_bins=True)
b=np.arange(np.min(lista), np.max(lista) + opt_width_scott_corre, opt_width_scott_corre)

#Hacer los histogramas para cada grupo    
hisporgrupos = {'control':np.zeros((len(grupo["control"]),len(b)-1)), 'MCI': np.zeros((len(grupo["MCI"]),len(b)-1)), 'AD': np.zeros((len(grupo["AD"]),len(b)-1)) }     
for cx in cxbanda.keys():
    print(cx)
    for pt in range(len(grupo[cx])):
        n,bgrupo=np.histogram(pruebadm[cx][pt],bins=b) 
        hisporgrupos[cx][pt]=n      
        
#Encontrar el valor máximo y mínimo de todos los grupos para poner límites
#a las barras de colores
lista_hist=[]
for h in hisporgrupos:
    lista_hist.append(hisporgrupos[kk]) #pruebadm[kk]
lista_final=np.concatenate(lista_hist).flatten()
valor_max= np.max(lista_final)
valor_min= np.min(lista_final)
        
#MAPAS DE CALOR
#Disposión de los números de correlaciones. 
#ESTE VALOR DEPENDE DEL TAMAÑO DE LOS PLOTS!!!!!!!!
#PARA COMPROBARLO HAY QUE HACER UNA GRÁFICA SIN XTICKS
multiplos=[valor for valor in range(len(b)) if valor%20==0]
g=[]
for elem in multiplos:
    g.append(np.round(b[elem], decimals=3))
    
    

    
plt.figure(figsize=(12,4))
cmap= sns.color_palette("viridis", as_cmap=True)
xticks= ticker.FixedLocator(multiplos)
formater = ticker.FixedFormatter(g)
#plt.suptitle("Distribution of Dominant Correlations - Alpha", fontsize=20)
plt.xlabel('Correlations')
plt.xticks(rotation = 75) 


ax1= plt.subplot(1,3,1)
sns.heatmap(hisporgrupos['control'], yticklabels=False, xticklabels='auto', cmap=cmap, vmin=valor_min, vmax=valor_max)
plt.title('Control')
ax1.xaxis.set_major_locator(xticks)
ax1.xaxis.set_major_formatter(formater)
plt.ylabel('Subjects')
plt.xticks(rotation = 75) 

ax2= plt.subplot(1,3,2)
sns.heatmap(hisporgrupos['MCI'], yticklabels=False, xticklabels='auto', cmap=cmap, vmin=valor_min, vmax=valor_max)
plt.title('MCI')
ax2.xaxis.set_major_locator(xticks)
ax2.xaxis.set_major_formatter(formater)
plt.xticks(rotation = 75) 
plt.xlabel('Correlations')

ax3=plt.subplot(1,3,3)
sns.heatmap(hisporgrupos['AD'], yticklabels=False, xticklabels='auto', cmap=cmap, vmin=valor_min, vmax=valor_max)
plt.title('AD')
ax3.xaxis.set_major_locator(xticks)
ax3.xaxis.set_major_formatter(formater)
plt.xticks(rotation = 75) 

plt.tight_layout()
plt.show()

        
    