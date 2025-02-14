## Study of the temporal architecture of neuronal activity in Alzheimer's Disease using dynamic functional connectivity techniques

### Description

This is the complete repository for my bachelor thesis, carried out in collaboration with the Biomedical Engineering Research Group from Universidad de Valladolid. The aim of the project was the investigation of neuronal alterations produced by dementia due to Alzheimer’s Disease through MEG and EEG signals. The final strategy on this study was to applied a predictive model able to differentia between the three groups under study; i.e, controls, patients with Mild Cognitive Impairment and patients with Alzheimer's Disease. For this, we implemented a Naïve Bayes classifier. 

Some of the results of this research where published in the Spanish Biomedical Engineering Society Conference 2021. ISBN: 978-84-09-36054-3 pg:22-25 [^1]

This project was the continuation of a previous research [^2], where a new method for detecting recurrent brain patterns (meta-states) was introduced, based on the calculation of instantaneous functional connectivity, recurrence plots, and community detection algorithms. The data for the current project can also be found in this reference. 

### Data

The data used in this bachelor thesis was the one provided from [^2]. For the current project we worked with two tensors of temporal representation at subject-level meta-states. 
* **Temporal Activation Sequence** (TAS): The TAS corresponds a discrete representation of the meta-state with the highest correlation value at each time instant. Therefore, the TAS is shaped as 1  T, where T represents the total temporal samples of the 60 second EEG.
* **Instantaneous Correlation Tensor** (ICT). The ICT shows the Spearman correlation of each meta-state with the IAC for each temporal sample. Its dimensions are K  T, where K represents the number of metastates and T the total temporal samples of the 60 second EEG.

### Repository Structure

* **Previous Steps folder**: This directory contains the code used to understand the data provided and make representations. It also includes simple metrics that lead to the final metrics. Dominant meta-states are those that presented the highest absolute value correlation in the ICT per instant. The remaining meta-states in that same instant are called non-dominant meta-states.

* *load_mats.py*: Load the mat_structs as python dictionaries. Checks for duplicates.
* *normalizar.py*: Normalizes each metric by the surrogate in each frequency band. 

* Metrics:
  * Degree of Antagonism: Reflects the force difference between the most attracting force from the ICT (positive correlation with highest absolute value), and the most repelling force (negative correlation with highest absolute value)
  * Antagonism Factor: Closely related to the degree of antagonism. Provides information about the amount of repulsive forces present against all correlations. 
  * IVG: Computes the index of qualitative variation to study the dispersion of the categorical distribution of meta-states. 
  * Shannon's Entropy: Computes the Shannon's entropy for each subject of the meta-states activation sequence. 
  * Naive Bayes: Implementation of a Naive Bayes Classifier using sklearn package.
    * Load metrics and select which ones to fit
    * Uses Leave One out cross validation
    * Returns the confussion matrix 


### References
[^1]: https://seib.org.es/wp-content/uploads/2024/05/Libro-de-actas-Caseib-2021.pdf
[^2]: P. Núñez et al., “Abnormal meta-state activation of dynamic brain networks across the Alzheimer spectrum,” Neuroimage, vol. 232, no. October 2020, p.117898, 2021.

