library(R.matlab)
library(stats)
library(gtools)

######################################## FUNCTIONS #######################################

## Function to load .mat files into R
loadMat = function(fileName)
{
  # Load in mat file and read the fields into variables
  rawData = readMat(fileName)
  out = list(data=rawData['data'],freq = rawData['freq'],latency = rawData['latency'])
}

## function to extract features from EEG data
extractFeats = function(data)
{
  
  # normalize channels
  norm_data = data.frame((data-rowMeans(data))/apply(data,1,sd))
  
  # calculate time frequency correlation matrix
  corr_matrix = cor(t(norm_data))
  
  # sorted eigen values
  eigen_values = eigen(corr_matrix)$values
  
  # upper right triangular to eliminate redundant features
  uppr_right = corr_matrix[upper.tri(corr_matrix,diag = FALSE)]
  
  # Frequency Domain
  # fft magnitudes of 1-47Hz
  fdata = abs(fft(data))
  fdata = fdata + 2e-13
  fdata = log10(t(fdata[,1:47]))
  
  # normalize by frequency bin
  ddata = data.frame((fdata-rowMeans(fdata))/apply(fdata,1,sd))
  fdata = t(ddata)
  logdata = as.vector(fdata)
  
  # calculate time frequency correlation matrix
  fcorr_matrix = cor(t(fdata))
  
  # sorted eigen values
  feigen_values = eigen(fcorr_matrix)$values
  
  # upper right triangular to eliminate redundant features
  fuppr_right = fcorr_matrix[upper.tri(fcorr_matrix,diag = FALSE)]
  
  # concatenate features into a vector
  feature_vec = t(c(eigen_values,feigen_values,logdata,uppr_right,fuppr_right))
  
}

# clip dimensions
#size = as.data.frame(lapply(clip$data, dim))


#################################### PROCESS SUBJECTS ####################################

### Subjects in study
subjects = c('Dog_1','Dog_2','Dog_3','Dog_4')

## Process subjects
for i in 1:length(subjects){
  
  # ictal clip paths
  clips.ictal = list.files(paste('Documents/Volumes/Seagate/seizure_detection/competition_data/clips/',subjects[i],sep = ''), pattern="*ictal*", full.names=TRUE)
  clips.ictal = clips.ictal[mixedorder(clips.ictal)]
  
  # interictal clip paths
  clips.interictal = list.files(paste('Documents/Volumes/Seagate/seizure_detection/competition_data/clips/',subjects[i],sep = ''), pattern="*interictal*", full.names=TRUE)
  clips.interictal = clips.interictal[mixedorder(clips.interictal)]
  
  # test clip paths
  clips.test = list.files(paste('Documents/Volumes/Seagate/seizure_detection/competition_data/clips/',subjects[i],sep = ''), pattern="*test*", full.names=TRUE)
  clips.test = clips.test[mixedorder(clips.test)]
  
  ## process ictal clips
  for ict in 1:length(clips.ictal){
    
    # load clip
    clip = loadMat(clips.ictal[ict])
  
    # convert to matrix
    data = matrix(unlist(clip$data),round(as.numeric(clip$freq)), byrow = TRUE)
    
    # extract features and create feature matrix
    if ict == 1{
      ictal_feats = data.frame(extract_feats(data))
    } else {
    ictal_feats = rbind(ictal_feats,extract_feats(data))
    }
   
  }
  
  ## process interictal clips
  for intict in 1:length(clips.interictal){
    
    # load clip
    clip = loadMat(clips.interictal[intict])
    
    # convert to matrix
    data = matrix(unlist(clip$data),round(as.numeric(clip$freq)), byrow = TRUE)
    
    # extract features and create feature matrix
    if intict == 1{
      interictal_feats = data.frame(extract_feats(data))
    } else {
      interictal_feats = rbind(interictal_feats,extract_feats(data))
    }
    
  }
  
  ## process test clips
  for tes in 1:length(clips.test){
    
    # load clip
    clip = loadMat(clips.test[tes])
    
    # convert to matrix
    data = matrix(unlist(clip$data),round(as.numeric(clip$freq)), byrow = TRUE)
    
    # extract features and create feature matrix
    if tes == 1{
      test_feats = data.frame(extract_feats(data))
    } else {
      test_feats = rbind(test_feats,extract_feats(data))
    }
    
  }
  
  ## store feature matrices in list of subjects
  if i == 1 {
    features.ictal = list(ictal_feats)
    features.interictal = list(interictal_feats)
    features.test = list(test_feats)
  } else {
    features.ictal = c(features.ictal,ictal_feats)
    features.interictal = c(features.interictal,interictal_feats)
    features.test = c(features.test,test_feats)
  }
  
}

########################## SEIZURE DETECTION ALGORITHMS/SUBJECT ##########################




