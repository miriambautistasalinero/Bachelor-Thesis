import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pickle
import os
import names_per_groups as nm #de señales base
from sklearn.model_selection import train_test_split

import loadmats as lm

grupos = nm.name_per_groups() #Separation per groups

#LOAD DATA AND CONVERT IT INTO DATA FRAME
#Grado de Ant
path_grado = 'C:/Users/miria/OneDrive/Documentos/Articulos UVA/data/var Grado_Ant'
path_factor='C:/Users/miria/OneDrive/Documentos/Articulos UVA/data/var Factor_Ant'

def extract_base(path):
    with open(os.path.join(path, 'Negativos_ocurrencia_delta'), 'rb') as g:
       dict_delta = pickle.load(g)
    
    with open(os.path.join(path, 'Negativos_ocurrencia_theta'), 'rb') as g1:
       dict_theta = pickle.load(g1)
    
    with open(os.path.join(path, 'Negativos_ocurrencia_alpha'), 'rb') as g2:
       dict_alpha = pickle.load(g2)
    
    with open(os.path.join(path, 'Negativos_ocurrencia_beta1'), 'rb') as g3:
      dict_beta1 = pickle.load(g3)
       
    with open(os.path.join(path, 'Negativos_ocurrencia_beta2'), 'rb') as g4:
       dict_beta2 = pickle.load(g4)
    return dict_delta, dict_theta, dict_alpha, dict_beta1, dict_beta2

def extract_sub(path):
    with open(os.path.join(path, 'Sub_negativos_ocurrencia_delta'), 'rb') as g:
       dict_delta = pickle.load(g)
    
    with open(os.path.join(path, 'Sub_negativos_ocurrencia_theta'), 'rb') as g1:
       dict_theta = pickle.load(g1)
    
    with open(os.path.join(path, 'Sub_negativos_ocurrencia_alpha'), 'rb') as g2:
       dict_alpha = pickle.load(g2)
    
    with open(os.path.join(path, 'Sub_negativos_ocurrencia_beta1'), 'rb') as g3:
      dict_beta1 = pickle.load(g3)
       
    with open(os.path.join(path, 'Sub_negativos_ocurrencia_beta2'), 'rb') as g4:
       dict_beta2 = pickle.load(g4)
    return dict_delta, dict_theta, dict_alpha, dict_beta1, dict_beta2



dict_delta, dict_theta, dict_alpha, dict_beta1, dict_beta2= extract_base(path_grado)
sub_delta, sub_theta, sub_alpha, sub_beta1, sub_beta2= extract_sub(path_grado)

def normalize(base, sub):
    norm = {'control':[], 'MCI':[], 'AD':[] }
    for hd in norm.keys():
        norm[hd] = np.divide(base[hd],sub[hd])
    return norm

norm_delta = normalize(dict_delta, sub_delta)
norm_theta = normalize(dict_theta, sub_theta)
norm_alpha = normalize(dict_alpha, sub_alpha)
norm_beta1 = normalize(dict_beta1, sub_beta1)
norm_beta2 = normalize(dict_beta2, sub_beta2)