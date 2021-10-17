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

path='C:/Users/miria/OneDrive/Documentos/Articulos UVA/data/var prob_ocurrencia base'
path_sub= 'C:/Users/miria/OneDrive/Documentos/Articulos UVA/data/var sub prob_ocurrencia'
grupos = nm.name_per_groups() #Separation per groups
def extract_prob(path):
    with open(os.path.join(path, 'prob_delta'), 'rb') as g:
       dict_delta = pickle.load(g)
    
    with open(os.path.join(path, 'prob_theta'), 'rb') as g1:
       dict_theta = pickle.load(g1)
    
    with open(os.path.join(path, 'prob_alpha'), 'rb') as g2:
       dict_alpha = pickle.load(g2)
    
    with open(os.path.join(path, 'prob_beta1'), 'rb') as g3:
      dict_beta1 = pickle.load(g3)
       
    with open(os.path.join(path, 'prob_beta2'), 'rb') as g4:
       dict_beta2 = pickle.load(g4)
    return dict_delta, dict_theta, dict_alpha, dict_beta1, dict_beta2

prob_delta, prob_theta, prob_alpha, prob_beta1, prob_beta2 = extract_prob(path)

def extract_prob_sub(path_sub):
    with open(os.path.join(path_sub, 'delta_sub_prob'), 'rb') as g:
       dict_delta = pickle.load(g)
    
    with open(os.path.join(path_sub, 'theta_sub_prob'), 'rb') as g1:
       dict_theta = pickle.load(g1)
    
    with open(os.path.join(path_sub, 'alpha_sub_prob'), 'rb') as g2:
       dict_alpha = pickle.load(g2)
    
    with open(os.path.join(path_sub, 'beta1_sub_prob'), 'rb') as g3:
      dict_beta1 = pickle.load(g3)
       
    with open(os.path.join(path_sub, 'beta2_sub_prob'), 'rb') as g4:
       dict_beta2 = pickle.load(g4)
    return dict_delta, dict_theta, dict_alpha, dict_beta1, dict_beta2

sub_prob_delta, sub_prob_theta, sub_prob_alpha, sub_prob_beta1, sub_prob_beta2 = extract_prob_sub(path_sub)
#INDICE DE VARIACIÓN QUALITATIVA
def ivq(prob_bandas):
    ivq={'control':[], 'MCI':[], 'AD':[]}
    B = {'control':[], 'MCI':[], 'AD':[]}
    for k in prob_bandas.keys():
        for i in range(len(prob_bandas[k][0])):
            B_control= 1- (np.sum(np.square(prob_bandas[k][0][i])))
            indice= B_control*len(prob_bandas[k][0][i])/(len(prob_bandas[k][0][i])-1)
            #B_control*((len(prob_bandas[k][0][i])-1)/len(prob_bandas[k][0][i]))
            ivq[k].append(indice)
            B[k].append(B_control)
    return B, ivq

B_delta,ivq_delta=ivq(prob_delta)
B_theta,ivq_theta=ivq(prob_theta)
B_alpha,ivq_alpha=ivq(prob_alpha)
B_beta1,ivq_beta1=ivq(prob_beta1)
B_beta2,ivq_beta2=ivq(prob_beta2)

#SUBROGADAS
sub_B_delta,sub_ivq_delta=ivq(sub_prob_delta)
sub_B_theta,sub_ivq_theta=ivq(sub_prob_theta)
sub_B_alpha,sub_ivq_alpha=ivq(sub_prob_alpha)
sub_B_beta1,sub_ivq_beta1=ivq(sub_prob_beta1)
sub_B_beta2,sub_ivq_beta2=ivq(sub_prob_beta2)

def normalizar(base,sub):
    norm = {'control':[], 'MCI':[], 'AD':[] }
    for hd in norm.keys():
        norm[hd] = np.divide(base[hd],sub[hd])
    return norm

#Normalización
norm_delta= normalizar(ivq_delta,sub_ivq_delta)
norm_theta= normalizar(ivq_theta,sub_ivq_theta)
norm_alpha= normalizar(ivq_alpha,sub_ivq_alpha)
norm_beta1= normalizar(ivq_beta1,sub_ivq_beta1)
norm_beta2= normalizar(ivq_beta2,sub_ivq_beta2)



lista_ivq=[norm_delta, norm_theta, norm_alpha, norm_beta1, norm_beta2]
p_values_krusk=[]
for ii in lista_ivq:
    #print(ii)
    krusk = kruskal(ii['control'],ii['MCI'],ii['AD'])
    p_values_krusk.append(krusk[1])
fdr_bandas=fdrcorrection(p_values_krusk)

condcl_alpha= mannwhitneyu(norm_alpha['control'],norm_alpha['MCI'])
conea_alpha= mannwhitneyu(norm_alpha['control'],norm_alpha['AD'])
dclea_alpha= mannwhitneyu(norm_alpha['MCI'],norm_alpha['AD'])

p_man_alpha =[condcl_alpha[1],conea_alpha[1],dclea_alpha[1]]
fdr_grupos=fdrcorrection(p_man_alpha)

#ivq = (1 - sum(fdp(:).^2))*length(fdp(:))/(length(fdp(:))-1);