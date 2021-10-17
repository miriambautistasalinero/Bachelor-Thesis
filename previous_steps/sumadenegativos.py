import loadmats as lm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import seaborn as sns
import astropy.stats as sts
from scipy.stats import kurtosis, skew, mode

pacientesbase = os.listdir("C:/Users/miria/OneDrive/Documentos/Articulos UVA/data/Resultados senales base/mats")

#k='correlationBeta2'
grupo = {'control':[], 'DCL':[], 'EA':[] }
numeros_neg = {'control':[], 'DCL':[], 'EA':[] }
numeros_pos = {'control':[], 'DCL':[], 'EA':[] }
correlation = {'control':[], 'DCL':[], 'EA':[] }
corr_ordenadas = {'control':[], 'DCL':[], 'EA':[] }
suma_negativos= {'control':[], 'DCL':[], 'EA':[] }
suma_positivos= {'control':[], 'DCL':[], 'EA':[] }
#negativos= {'control':[], 'DCL':[], 'EA':[] }
cxgrupo={'control':[], 'DCL':[], 'EA':[] }

for p in pacientesbase:
    
    dsi= lm.loadmat(p)['dynamicStateInfo']
    group =dsi['group']
    #tasp = dsi['temporalActivation'][0,:] #El 0 es porque estoy en banda DELTA
    correlationp=dsi['correlation']
    #corr_ordenadas= np.zeros(correlation.shape)
    #elemintas = list(set(tasp))
    
    if group == 'DCL1' or group =='DCL2' or group == 'DCL':
            grupo["DCL"].append(p)
            correlation['DCL'].append(correlationp)
            
    elif group == 'EA1' or group == 'EA2':
        grupo["EA"].append(p)
        correlation['EA'].append(correlationp)
        
    else:
        grupo["control"].append(p)
        correlation['control'].append(correlationp)

for b in grupo.keys():
    
    for k in correlationp.keys():
        negativos= {'control':[], 'DCL':[], 'EA':[] }
        positivos= {'control':[], 'DCL':[], 'EA':[] }
        
        
        elemneg=0
        elempos=0
        for i in range(len(grupo[b])):
             c = correlation[b][i][k].flatten()
             for e in range(len(c)):
                 if c[e] < 0:
                     negativos[b].append(c[e])
                     elemneg+=1
                 else:
                     positivos[b].append(c[e])
                     elempos+=1
                    
                     
        #Normalizando
        numeros_neg[b].append(elemneg/len(grupo[b]))
        suma_negativos[b].append(sum(negativos[b])/len(grupo[b]))
        
        numeros_pos[b].append(elempos/len(grupo[b]))
        suma_positivos[b].append(sum(positivos[b])/len(grupo[b]))
        
   
plt.figure()
plt.title('Sum of negative values')                
plt.plot([1,2,3,4,5],suma_negativos['control'], label='Control')
plt.plot([1,2,3,4,5],suma_negativos['DCL'],label='DCL')
plt.plot([1,2,3,4,5],suma_negativos['EA'], label='EA')
labels= ['Delta','Theta', 'Alpha', 'Beta1', 'Beta2']
ticks=[1,2,3,4,5]
plt.xticks(ticks, labels)
plt.legend()
plt.xlabel('Frequency bands')

plt.figure()
plt.title('Sum of positive values')                
plt.plot([1,2,3,4,5],suma_positivos['control'], label='Control')
plt.plot([1,2,3,4,5],suma_positivos['DCL'],label='DCL')
plt.plot([1,2,3,4,5],suma_positivos['EA'], label='EA')
labels= ['Delta','Theta', 'Alpha', 'Beta1', 'Beta2']
ticks=[1,2,3,4,5]
plt.xticks(ticks, labels)
plt.legend()
plt.xlabel('Frequency bands')


#Boxplots
#POSITIVE VALUES

# set width of bar
barWidth = 0.25

# Set position of bar on X axis
r1 = np.arange(len(suma_positivos['control']))
r2 = [x + barWidth for x in r1]
r3 = [x + barWidth for x in r2]


plt.figure()
plt.title('Sum of positive values')  
plt.bar(r1, suma_positivos['control'],  edgecolor='white', width = 0.25, label='control')
plt.bar(r2, suma_positivos['DCL'],  edgecolor='white',width = 0.25, label='DCL')
plt.bar(r3, suma_positivos['EA'],  edgecolor='white', width = 0.25, label='EA')

# Add xticks on the middle of the group bars
#plt.xlabel('group', fontweight='bold')
plt.xticks([r + barWidth for r in range(len(suma_positivos['control']))], ['Delta','Theta', 'Alpha', 'Beta1', 'Beta2'])
 
# Create legend & Show graphic
plt.legend()
plt.show()

#NEGATIVES VALUES
# Set position of bar on X axis

plt.figure()
plt.title('Sum of negative values')  
plt.bar(r1, suma_negativos['control'],  edgecolor='white', width = 0.25, label='control')
plt.bar(r2, suma_negativos['DCL'],  edgecolor='white',width = 0.25, label='DCL')
plt.bar(r3, suma_negativos['EA'],  edgecolor='white', width = 0.25, label='EA')

# Add xticks on the middle of the group bars
#plt.xlabel('group', fontweight='bold')
plt.xticks([r + barWidth for r in range(len(suma_negativos['control']))], ['Delta','Theta', 'Alpha', 'Beta1', 'Beta2'])
 
# Create legend & Show graphic
plt.legend()
plt.show()
