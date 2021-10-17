#SEPERACIÓN DE PACIENTES BASE POR GRUPO, UNICAMENTE SE SEPARA EL NOMBRE

import os
import sys

#Insert path for load mat modules
sys.path.insert(0, 'C:/Users/miria/OneDrive/Documentos/Articulos UVA/data/Resultados senales base')
sys.path.insert(0, 'C:/Users/miria/OneDrive/Documentos/Articulos UVA/data/Resultados senales base')
import loadmats as lm

def name_per_groups():
    
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
    path = "C:/Users/miria/OneDrive/Documentos/Articulos UVA/data/Resultados senales base/mats"
    pacientesbase = os.listdir(path)

    grupo = {'control':[], 'MCI':[], 'AD':[] }
    
    for p in pacientesbase:
    
        group= lm.loadmat(os.path.join(path, p))['dynamicStateInfo']['group']
        #group =dsi['group']
        
        
        
        if group == 'DCL1' or group =='DCL2' or group == 'DCL':
             grupo["MCI"].append(p)
        
        
        elif group == 'EA1' or group == 'EA2':
             grupo["AD"].append(p)
        
        else:
            grupo["control"].append(p)
        
    return grupo