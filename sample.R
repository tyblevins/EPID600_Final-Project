library(R.matlab)
library(stats)
library(gtools)
library(randomForest)
library(e1071)
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
  norm_data = data.frame((data-colMeans(data))/apply(data,2,sd))
  
  # calculate time frequency correlation matrix
  corr_matrix = cor(norm_data)
  
  # sorted eigen values
  eigen_values = eigen(corr_matrix)$values
  
  # upper right triangular to eliminate redundant features
  uppr_right = corr_matrix[upper.tri(corr_matrix,diag = FALSE)]
  
  # Frequency Domain
  # fft magnitudes of 1-47Hz
  fdata = abs(fft(data))
  fdata = fdata + 2e-13
  fdata = log10(fdata[1:47,])
  logdata = as.vector(fdata)
  
  # normalize by frequency bin
  ddata = ((fdata-rowMeans(fdata))/apply(fdata,1,sd))
  
  # calculate time frequency correlation matrix
  fcorr_matrix = cor(fdata)
  
  # sorted eigen values
  feigen_values = eigen(fcorr_matrix)$values
  
  # upper right triangular to eliminate redundant features
  fuppr_right = fcorr_matrix[upper.tri(fcorr_matrix,diag = FALSE)]
  
  # concatenate features into a vector
  feature_vec = c(eigen_values,feigen_values,logdata,uppr_right,fuppr_right)
  
}

# clip dimensions
#size = as.data.frame(lapply(clip$data, dim))


#################################### PROCESS SUBJECTS ####################################

### Subjects in study
subjects = c('Dog_1')

## Process subjects
for (i in 1:length(subjects)){
  print(i)
  # ictal clip paths
  clips.ictal = list.files(paste('./data/',subjects[i],sep = ''), pattern="*ictal*", full.names=TRUE)
  clips.ictal = clips.ictal[mixedorder(clips.ictal)]
  
  # interictal clip paths
  clips.interictal = list.files(paste('./data/',subjects[i],sep = ''), pattern="*interictal*", full.names=TRUE)
  clips.interictal = clips.interictal[mixedorder(clips.interictal)]
  
  # test clip paths
  clips.test = list.files(paste('./data/',subjects[i],sep = ''), pattern="*test*", full.names=TRUE)
  clips.test = clips.test[mixedorder(clips.test)]
  
  ## process ictal clips
  for (ict in 1:length(clips.ictal)){
    
    # load clip
    clip = loadMat(clips.ictal[ict])
    
    # convert to matrix
    data = matrix(unlist(clip$data),round(as.numeric(clip$freq)), byrow = TRUE)
    
    # extract features and create feature matrix
    if (ict == 1){
      ictal_feats = extractFeats(data)
    } else {
      ictal_feats = rbind(ictal_feats,extractFeats(data))
    }
    
  }
  print('ictal clips done!')
  
  ## process interictal clips
  for (intict in 1:length(clips.interictal)){
    
    # load clip
    clip = loadMat(clips.interictal[intict])
    
    # convert to matrix
    data = matrix(unlist(clip$data),round(as.numeric(clip$freq)), byrow = TRUE)
    
    # extract features and create feature matrix
    if (intict == 1){
      interictal_feats = extractFeats(data)
    } else {
      interictal_feats = rbind(interictal_feats,extractFeats(data))
    }
    
  }
  
  print('interictal clips done!')
  ## process test clips
  for (tes in 1:length(clips.test)){
    # load clip
    clip = loadMat(clips.test[tes])
    
    # convert to matrix
    data = matrix(unlist(clip$data),round(as.numeric(clip$freq)), byrow = TRUE)
    
    # extract features and create feature matrix
    if (tes == 1){
      test_feats = extractFeats(data)
    } else {
      test_feats = rbind(test_feats,extractFeats(data))
    }
    
  }
  print('test clips done!')
  ## store feature matrices in list of subjects
  if (i == 1) {
    features.ictal = list(ictal_feats)
    features.interictal = list(interictal_feats)
    features.test = list(test_feats)
  } else {
    features.ictal[[i]] = ictal_feats
    features.interictal[[i]] = interictal_feats
    features.test[[i]] = test_feats
  }
  print('features stored!')
}



# import test labels


########################## SEIZURE DETECTION ALGORITHMS/SUBJECT ##########################
# random forest
pred.rf = list()
for (i in 1:length(subjects)){
  print('')
  
  # Make labels
  labels.ictal = rep(1,dim(features.ictal[[i]])[1])
  labels.interictal = rep(0,dim(features.interictal[[i]])[1])
  
  # random forest
  detector.rf = randomForest(rbind(features.ictal[[i]],features.interictal[[i]]),factor(c(labels.ictal,labels.interictal)),ntree=100,importance=TRUE) 
  
  # make predictions
  pred.rf[[i]] = predict(detector.rf,features.test[[i]], type="prob")
}

# svm
pred.svm = list()


# train
detector.svm = svm(rbind(features.ictal[[i]],features.interictal[[i]]),factor(c(labels.ictal,labels.interictal)), scale = TRUE, kernel = "radial")

# test
pred.svm = predict(detector.svm,features.test[[i]], type="prob")



