# EPID600_Final-Project
Data Science for Biomedical Informatics Final Project

## Seizure Detection in Canine Subjects
EEG from four canine subjects will be analyzed. EEG data will be clipped into 1 second segments and classified as 'Ictal' or Interictal'. Features derived from the literature will be extracted from the clips, analyzed for importance, and ran through a seizure detection algorithm. The performance of the algorithm will be evaluated by ROC and AUC metrics

## Sample code
Sample code and data is provided to evaluate one subject in the study. Please use the sample R script file along with the data stored in data/Dog1. This project must be your working directory path for the script to execute correctly.

## Method
1. EEG data is pulled from ieeg.org and clipped into 1-second segments using the ieeg.org toolbox in Matlab.  
2. .mat files are loaded into R and processed.
3. The following features are extracted:
  - electrode channel correlations
  - correlation eigenvalues
  - energy
  - log10 magnitudes of frequencies in the range of 1-47Hz
  - frequency correlations
  - frequency correlation eigenvalues
4. Subject feature matrices are processed by seizure detector (classification algorithm)
  - Random Forest Classifier
5. Algorithms are evaluated using ROC and AUC metrics.

## Feature Visualizations
1. Time Domain
![energy and correlations](https://raw.github.com/tyblevins/EPID600_Final-Project/master/pics/time_feats.png)
2. Frequency Domain
![eigenvalues and frequency correlations](https://raw.github.com/tyblevins/EPID600_Final-Project/master/pics/freq_feats.png)

## Results - Model Performance
![ROC curves](https://raw.github.com/tyblevins/EPID600_Final-Project/master/pics/results.png)
