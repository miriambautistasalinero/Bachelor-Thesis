from sklearn.model_selection import train_test_split
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.model_selection import cross_val_score
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis 
from sklearn.model_selection import LeaveOneOut
from sklearn import datasets
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pickle
import os
import names_per_groups as nm #de señales base
from sklearn.model_selection import train_test_split
from sklearn.metrics import ConfusionMatrixDisplay
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
from sklearn.metrics import plot_confusion_matrix
from sklearn.naive_bayes import GaussianNB
import seaborn as sns

import loadmats as lm

#LOAD DATA AND CONVERT IT INTO DATA FRAME
#Grado de Ant
path_grado = 'C:/Users/miria/OneDrive/Documentos/Articulos UVA/data/var Grado_Ant'
path_factor='C:/Users/miria/OneDrive/Documentos/Articulos UVA/data/var Factor_Ant'
path_norm = 'C:/Users/miria/OneDrive/Documentos/Articulos UVA/data/Normalizados'

def extract_factor(path):
    with open(os.path.join(path, 'fac_delta'), 'rb') as g:
       dict_delta = pickle.load(g)
    
    with open(os.path.join(path, 'fac_theta'), 'rb') as g1:
       dict_theta = pickle.load(g1)
    
    with open(os.path.join(path, 'fac_alpha'), 'rb') as g2:
       dict_alpha = pickle.load(g2)
    
    with open(os.path.join(path, 'fac_beta1'), 'rb') as g3:
      dict_beta1 = pickle.load(g3)
       
    with open(os.path.join(path, 'fac_beta2'), 'rb') as g4:
       dict_beta2 = pickle.load(g4)
    return dict_delta, dict_theta, dict_alpha, dict_beta1, dict_beta2

def extract_grado(path):
    with open(os.path.join(path, 'grad_delta'), 'rb') as g:
       dict_delta = pickle.load(g)
    
    with open(os.path.join(path, 'grad_theta'), 'rb') as g1:
       dict_theta = pickle.load(g1)
    
    with open(os.path.join(path, 'grad_alpha'), 'rb') as g2:
       dict_alpha = pickle.load(g2)
    
    with open(os.path.join(path, 'grad_beta1'), 'rb') as g3:
      dict_beta1 = pickle.load(g3)
       
    with open(os.path.join(path, 'grad_beta2'), 'rb') as g4:
       dict_beta2 = pickle.load(g4)
    return dict_delta, dict_theta, dict_alpha, dict_beta1, dict_beta2

def extract_shannon(path):
    with open(os.path.join(path, 'shannon_delta'), 'rb') as g:
       dict_delta = pickle.load(g)
    
    with open(os.path.join(path, 'shannon_theta'), 'rb') as g1:
       dict_theta = pickle.load(g1)
    
    with open(os.path.join(path, 'shannon_alpha'), 'rb') as g2:
       dict_alpha = pickle.load(g2)
    
    with open(os.path.join(path, 'shannon_beta1'), 'rb') as g3:
      dict_beta1 = pickle.load(g3)
       
    with open(os.path.join(path, 'shannon_beta2'), 'rb') as g4:
       dict_beta2 = pickle.load(g4)
    return dict_delta, dict_theta, dict_alpha, dict_beta1, dict_beta2

def extract_iqv(path):
    with open(os.path.join(path, 'iqv_delta'), 'rb') as g:
       dict_delta = pickle.load(g)
    
    with open(os.path.join(path, 'iqv_theta'), 'rb') as g1:
       dict_theta = pickle.load(g1)
    
    with open(os.path.join(path, 'iqv_alpha'), 'rb') as g2:
       dict_alpha = pickle.load(g2)
    
    with open(os.path.join(path, 'iqv_beta1'), 'rb') as g3:
      dict_beta1 = pickle.load(g3)
       
    with open(os.path.join(path, 'iqv_beta2'), 'rb') as g4:
       dict_beta2 = pickle.load(g4)
    return dict_delta, dict_theta, dict_alpha, dict_beta1, dict_beta2

f_delta, f_theta, f_alpha, f_beta1, f_beta2 = extract_factor(path_norm)
g_delta, g_theta, g_alpha, g_beta1, g_beta2 = extract_grado(path_norm)
s_delta, s_theta, s_alpha, s_beta1, s_beta2 = extract_shannon(path_norm)
i_delta, i_theta, i_alpha, i_beta1, i_beta2 = extract_iqv(path_norm)

metrics=[g_delta, g_theta, g_alpha, g_beta1, g_beta2, f_delta, f_theta, f_alpha, f_beta1, f_beta2, s_delta, s_theta, s_alpha, s_beta1, s_beta2, i_delta, i_theta, i_alpha, i_beta1, i_beta2]
name_metrics=['g_delta', 'g_theta', 'g_alpha', 'g_beta1', 'g_beta2', 'f_delta', 'f_theta', 'f_alpha', 'f_beta1', 'f_beta2', 's_delta', 's_theta', 's_alpha', 's_beta1', 's_beta2', 'i_delta', 'i_theta', 'i_alpha', 'i_beta1', 'i_beta2']

#DO A LIST WITH THE NAME AND THE GROUP
def list_grupos():
    path = "C:/Users/miria/OneDrive/Documentos/Articulos UVA/data/Resultados senales base/mats"
    pacientesbase = os.listdir(path)
    dic_final={}
    target=[]

    #Take each metric 
    for item in range(len(metrics)):
        var_metrics=[] #to upload the value and update the data fram
        grupo=[] #save the group of the patient
        subject=[] #save the name of the patient
        target=[] #save the category
        print(name_metrics[item]) #follow up

        #Number of the patient for each group
        a=0 #AD
        c=0 #control
        m=0 #MCI
        
        for p in pacientesbase:
            
            group= lm.loadmat(os.path.join(path, p))['dynamicStateInfo']['group']
            subject.append(p)

            
            
            if group == 'DCL1' or group =='DCL2' or group == 'DCL':
                grupo.append('MCI')
                #Extracts the value of one metric for the patient m in group MCI
                elemento_de_var=metrics[item]['MCI'][m]
                #delta= dict_delta['MCI'][m]
                
                t=1
                m+=1
            
            elif group == 'EA1' or group == 'EA2':
                grupo.append('AD')
                #Extracts the value of one metric for the patient a in group control
                elemento_de_var=metrics[item]['AD'][a]
                
                t=2
                a+=1
            
            else:
                grupo.append('control')
                #Extracts the value of one metric for the patient c in group control
                elemento_de_var=metrics[item]['control'][c]
                
                t=0
                c+=1
                
            target.append(str(t))
            #The list is upladated for each subject
            var_metrics.append(elemento_de_var)
        #update the final dictionary
        dic_final[name_metrics[item]] = var_metrics
    
    dic_final['target'] = target
    dic_final['group'] = grupo
    dic_final['subject']=subject
    
    
    return dic_final



dic =list_grupos()

#CREATE DATA FRAME
df = pd.DataFrame(dict([(h, pd.Series(j)) for h,j in dic.items()]))


#Fit the LDA model
X = df[['g_delta', 'g_theta', 'g_alpha', 'g_beta1', 'g_beta2', 
        'f_delta', 'f_theta', 'f_alpha', 'f_beta1', 'f_beta2', 
        's_delta', 's_theta', 's_alpha', 's_beta1', 's_beta2', 
        'i_delta', 'i_theta', 'i_alpha', 'i_beta1', 'i_beta2']]
y = df['target']

#NAIVE BAYES
model = GaussianNB()
cv = LeaveOneOut()

scores_nb = cross_val_score(model, X, y, scoring='accuracy', cv=cv)
print("ACC total", np.mean(scores_nb))  

X_selec = df[['g_alpha', 'g_beta1', 'f_alpha', 'f_beta1',  's_alpha', 'i_alpha']]
y_selec = df['target']

scores_nb_selec = cross_val_score(model, X_selec, y_selec, scoring='accuracy', cv=cv)
print("ACC total", np.mean(scores_nb_selec)) 

def NB_fit(df):
    predicted=[]
    for ii in range(len(df)):
        #Create the new df extracting one patient
        new_df= df.drop([ii])
        single_df = df.iloc[[ii]]
        #From the new df create X and Y
        X = new_df[['g_delta', 'g_theta', 'g_alpha', 'g_beta1', 'g_beta2', 'f_delta', 'f_theta', 'f_alpha', 'f_beta1', 'f_beta2', 's_delta', 's_theta', 's_alpha', 's_beta1', 's_beta2', 'i_delta', 'i_theta', 'i_alpha', 'i_beta1', 'i_beta2']]
        y = new_df['target']
        #Eliminate the group category for the subject used to predict
        X_predict = single_df[['g_delta', 'g_theta', 'g_alpha', 'g_beta1', 'g_beta2', 'f_delta', 'f_theta', 'f_alpha', 'f_beta1', 'f_beta2', 's_delta', 's_theta', 's_alpha', 's_beta1', 's_beta2', 'i_delta', 'i_theta', 'i_alpha', 'i_beta1', 'i_beta2']]
        #Fit the model an create prediction
        model.fit(X, y)
        y_hat= model.predict(X_predict)
        predicted.append(str(y_hat[0]))
    return predicted

def NB_fit_selec(df):
    predicted=[]
    for ii in range(len(df)):
        #Create the new df extracting one patient
        new_df= df.drop([ii])
        single_df = df.iloc[[ii]]
        #From the new df create X and Y
        X = new_df[['g_alpha', 'g_beta1', 'f_alpha', 'f_beta1',  's_alpha', 'i_alpha']]
        y = new_df['target']
        #Eliminate the group category for the subject used to predict
        X_predict = single_df[['g_alpha', 'g_beta1', 'f_alpha', 'f_beta1',  's_alpha', 'i_alpha']]
        #Fit the model an create prediction
        model.fit(X, y)
        y_hat= model.predict(X_predict)
        predicted.append(str(y_hat[0]))
    return predicted

#Call the model fit and predict
yp =  NB_fit(df)
confusionmatrix = confusion_matrix(y, yp,  labels=["0", "1", "2"])
print("Normal", confusionmatrix)
print(classification_report(y, yp, labels=["0", "1", "2"]))

yp_selec  =  NB_fit_selec(df)
confusion_matrix_selec = confusion_matrix(y, yp_selec,  labels=["0", "1", "2"])
print("Seleccion de características", confusion_matrix_selec)
print(classification_report(y, yp_selec, labels=["0", "1", "2"]))

#Plot confusion matrix
#Create df
df_cm= pd.DataFrame(confusionmatrix, columns=['Control', 'MCI', 'AD'], index=['Control', 'MCI', 'AD'])
df_cm.index.name = 'Predicted'
df_cm.columns.name = 'Actual'

#To ensure equal color scale
maximos = max([np.max(confusionmatrix), np.max(confusion_matrix_selec)])

plt.figure()
sns.heatmap(df_cm, annot=True, vmin=0, vmax=maximos)
plt.show()

df_cm_selec= pd.DataFrame(confusion_matrix_selec, columns=['Control', 'MCI', 'AD'], index=['Control', 'MCI', 'AD'])
df_cm_selec.index.name = 'Predicted'
df_cm_selec.columns.name = 'Actual'

plt.figure()
sns.heatmap(df_cm_selec, annot=True, vmin=0, vmax=maximos)
plt.show()