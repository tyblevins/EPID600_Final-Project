# EPID600_Final-Project
Data Science for Biomedical Informatics Final Project

## Seizure Detection in Canine Subjects
EEG from four canine subjects will be analyzed. EEG data will be clipped into 1 second segments and classified as 'Ictal' or Interictal'. Features derived from the literature will be extracted from the clips, analyzed for importance, and ran through a seizure detection algorithm. The performance of the algorithm will be evaluated by ROC and AUC metrics

## Method
1. EEG data is pulled from ieeg.org and clipped into 1-second segments using the ieeg.org toolbox in Matlab.  
2. .mat files are loaded into R and processed.
3. The following features are extracted:
  - electrode channel correlations
  - correlation eigenvalues
  - log10 magnitudes of frequencies in the range of 1-47Hz
  - frequency correlations
  - frequency correlation eigenvalues
4. Subject feature matrices are processed by seizure detector.....TBD
