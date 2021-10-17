import os
import loadmats as lm
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#for dirpath, dirnames, filenames in os.walk("C:/Users/miria/Documents/Artículos UVA/data/Resultados señales base"):
#    print("Ruta actual:", dirpath)
#    print("Carpetas:", ", ".join(dirnames))
#    print("Archivos:", ", ".join(filenames))

def grupos():
    
     pacientes = os.listdir("C:/Users/miria/Documents/Artículos UVA/data/Resultados señales base/mats")
     DCL =[]
     EA= []
     CONTROL = []
     b=0
     d=[]
     ddic ={}
     path = []
     for i in pacientes:
        # print(i)
        #os.path.join((pacientes i))
         path.append(os.path.join((pacientes + i)))
         mat_pac = lm.loadmat(i)
         rawdata =mat_pac['dynamicStateInfo']
         group = rawdata['group']
         d.append(group)
         #ddic={i:group}
         
         #df= pd.DataFrame(d, index=(i),columns=('grupo'))
         if group == 'DCL1' or group =='DCL2' or group == 'DCL':
             
             DCL.append(i)
         elif group == 'EA1' or group == 'EA2':
             EA.append(i)
         else:
             CONTROL.append(i)
         b += 1  
         
         
     return path, DCL,  EA, CONTROL, d
         
         
def info(p):

    pacientesbase = os.listdir("C:/Users/miria/Documents/Artículos UVA/data/Resultados señales base")
    

    mat_pac = lm.loadmat(pacientesbase[p])
    rawdata =mat_pac['dynamicStateInfo']
    groups = rawdata['group']
    print("This is patient " + pacientesbase[p] + "from group " + groups)
    return rawdata

def plots(rawdata):

#Tengo que limpiar las variables para que al volverlo a ejecutar no tenga problema    
    
    #Matriz de correlaciones
    Matcor = rawdata['correlation']
    #Matriz del TAS
    TAS = rawdata['temporalActivation']
    #Matriz de clusters
    Matclusters = rawdata['clusters']
    
    #Obtener la lista de los keys del dicccionario de correlacion
    lista =[]
    for n in Matcor.keys():
        lista.append(n)
    
    #Obtener las matrices de correlacion por banda y la suma 
    suma=[]
    
    for k in lista:    
        
        #axis= 0 me suma las columnas
        suma.append(np.apply_along_axis(sum, 0, abs( Matcor[k][:,0:2000])))
        Matcor[k] = Matcor[k][:,0:2000].T
        r,c = Matcor[k].shape
        print(np.mean(suma))
    
        
    #Matriz de clusters
    clusters =[]
    for n in Matclusters.keys():
        clusters.append(n)
    
    #Ploteando todo     
    for i in range(len(clusters)):
    
       
       plt.figure(figsize=(20,10))
       sample = np.arange(0,2000) 
       
       #SUBPLOT DEL TAS
       Tas= TAS[i, 0:2000]
       plt.subplot(3,1,1)
       plt.plot(sample, Tas)
       plt.xlim(0,2000)
       plt.xlabel('Samples');
       plt.ylabel('TAS in' + ' ' + str(clusters[i]))
   
   
       #SUBPLOT DEL ICT
      
           
       plt.subplot(3,1,2)
    
       plt.plot(sample, Matcor[lista[i]])
       plt.xlim(0,2000)
       plt.xlabel('Samples');
       plt.ylabel('Correlation')
        #for l in range(c):
            #print(l)
       #plt.legend()
   
    
       #SUBPLOT DEL VALOR ABSOLUTO DE LA CORRELACION
       plt.subplot(3,1,3)
       plt.plot(sample, suma[i])
       plt.xlim(0,2000)
       plt.xlabel('Samples');
       plt.ylabel('Abs Corr' + '' + str(clusters[i]))
    
       plt.show()
       
            