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
  
  # calculate energy
  energy = log10(colSums(data*data)/dim(data)[1])
  
  # normalize channels
  norm_data = data.frame(data-matrix(rep(colMeans(data),400),ncol=16,byrow=TRUE)/matrix(rep(apply(data,2,sd),400),ncol=16,byrow=TRUE))
  
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
  fdata = log10(fdata[2:48,])
  logdata = as.vector(fdata)
  
  # normalize by frequency bin
  ddata = ((fdata-matrix(rep(rowMeans(fdata),16),nrow=47,byrow=FALSE))/matrix(rep(apply(fdata,1,sd),16),nrow=47,byrow=FALSE))
  
  # calculate time frequency correlation matrix
  fcorr_matrix = cor(fdata)
  
  # sorted eigen values
  feigen_values = eigen(fcorr_matrix)$values
  
  # upper right triangular to eliminate redundant features
  fuppr_right = fcorr_matrix[upper.tri(fcorr_matrix,diag = FALSE)]
  
  # concatenate features into a vector
  feature_vec = c(eigen_values,feigen_values,logdata,uppr_right,fuppr_right,energy)
  
}

# clip dimensions
#size = as.data.frame(lapply(clip$data, dim))


#################################### PROCESS SUBJECTS ####################################

### Subjects in study
subjects = c('Dog_1','Dog_2','Dog_3','Dog_4')

## Process subjects
for (i in 1:length(subjects)){
  print(i)
  # ictal clip paths
  clips.ictal = list.files(paste('Documents/Volumes/Seagate/seizure_detection/competition_data/clips/',subjects[i],sep = ''), pattern="*_ictal_*", full.names=TRUE)
  clips.ictal = clips.ictal[mixedorder(clips.ictal)]
  
  # interictal clip paths
  clips.interictal = list.files(paste('Documents/Volumes/Seagate/seizure_detection/competition_data/clips/',subjects[i],sep = ''), pattern="*_interictal_*", full.names=TRUE)
  clips.interictal = clips.interictal[mixedorder(clips.interictal)]
  
  # test clip paths
  clips.test = list.files(paste('Documents/Volumes/Seagate/seizure_detection/competition_data/clips/',subjects[i],sep = ''), pattern="*_test_*", full.names=TRUE)
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
true_labels = read.csv("Documents/Github/revised_solution.csv")
dog1_tlabels = true_labels[1:3181,2]
dog1_tlabels[dog1_tlabels==-1] = 0
dog2_tlabels = true_labels[3182:6178,2]
dog2_tlabels[dog2_tlabels==-1] = 0
dog3_tlabels = true_labels[6179:10628,2]
dog3_tlabels[dog3_tlabels==-1] = 0
dog4_tlabels = true_labels[10629:13641,2]
dog4_tlabels[dog4_tlabels==-1] = 0

dog_labels = list(dog1_tlabels,dog2_tlabels,dog3_tlabels,dog4_tlabels)

########################## SEIZURE DETECTION ALGORITHMS/SUBJECT ##########################


# Random forest classifier to predict test set labels
pred.rf = list()
for (i in 1:length(subjects)){
  print(i)
  
  # Make labels
  labels.ictal = rep(1,dim(features.ictal[[i]])[1])
  labels.interictal = rep(0,dim(features.interictal[[i]])[1])
  
  # random forest
  detector.rf = randomForest(rbind(features.ictal[[i]],features.interictal[[i]]),factor(c(labels.ictal,labels.interictal)),ntree=1000,importance=TRUE,mtry=1) 
  
  # make predictions
  pred.rf[[i]] = predict(detector.rf,features.test[[i]], type="prob")
}

# sample for SVM classifier
# train
detector.svm = svm(rbind(features.ictal[[i]],features.interictal[[i]]),factor(c(labels.ictal,labels.interictal)), scale = TRUE, kernel = "radial")

# test
pred.svm = predict(detector.svm,features.test[[i]], type="prob")

# heatmap of correlations
heat.ictal = features.ictal[[1]][,785:904]
heat.interictal = features.interictal[[1]][,785:904]
k = rbind(heat.interictal,heat.ictal)

heat.if = features.ictal[[1]][,905:1024]
heat.iif = features.interictal[[1]][,905:1024]
kk = rbind(heat.iif,heat.if)

# feature index values (for visualizations)
ev_ref = 1:16

fev_ref = 17:32

log_ref = 33:784

corr_ref = 485:904

fcorr_ref = 905:1024

energy_ref = 1025:1040


ev.if = features.ictal[[1]][,ev_ref]
ev.iif = features.interictal[[1]][,ev_ref]

fev.if = features.ictal[[1]][,fev_ref]
fev.iif = features.interictal[[1]][,fev_ref]

log.if = features.ictal[[1]][,log_ref]
log.iif = features.interictal[[1]][,log_ref]

energy.if = features.ictal[[1]][,energy_ref]
energy.iif = features.interictal[[1]][,energy_ref]


a = apply(heat.if,2,mean)

b = matrix(1,16,16)
b[lower.tri(b, diag = FALSE)] = a
b = t(b)
b[lower.tri(b, diag = FALSE)] = a


aa = apply(heat.iif,2,mean)
bb = matrix(1,16,16)
bb[lower.tri(bb,diag = FALSE)] = aa
bb = t(bb)
bb[lower.tri(bb,diag = FALSE)] = aa



# Calculate AUC and ROC metrics
library(AUC)
roc_curve = list()
auc_val = list()
for (i in 1:length(subjects)){
  roc_curve[[i]] = roc(as.vector(pred.rf[[i]][,2]),factor(dog_labels[[i]]))
  auc_val[[i]] = auc(roc_curve[[i]],min=0,max=1) 
}


# Dog 1
pred <- prediction(as.vector(pred.rf[[1]][,2]), factor(dog_labels[[1]]))
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]

# Dog 2
pred.2 <- prediction(as.vector(pred.rf[[2]][,2]), factor(dog_labels[[2]]))
perf.2 <- performance(pred.2, measure = "tpr", x.measure = "fpr")
auc.2 <- performance(pred.2, measure = "auc")
auc.2 <- auc.2@y.values[[1]]

# Dog 3
pred.3 <- prediction(as.vector(pred.rf[[3]][,2]), factor(dog_labels[[3]]))
perf.3 <- performance(pred.3, measure = "tpr", x.measure = "fpr")
auc.3 <- performance(pred.3, measure = "auc")
auc.3 <- auc.3@y.values[[1]]

# Dog 4
pred.4 <- prediction(as.vector(pred.rf[[4]][,2]), factor(dog_labels[[4]]))
perf.4 <- performance(pred.4, measure = "tpr", x.measure = "fpr")
auc.4 <- performance(pred.4, measure = "auc")
auc.4 <- auc.4@y.values[[1]]

# ROC data frames
roc.data <- data.frame(fpr=unlist(perf@x.values),
                       tpr=unlist(perf@y.values),
                       model="GLM")

roc.data.2 <- data.frame(fpr=unlist(perf.2@x.values),
                       tpr=unlist(perf.2@y.values),
                       model="GLM")

roc.data.3 <- data.frame(fpr=unlist(perf.3@x.values),
                         tpr=unlist(perf.3@y.values),
                         model="GLM")

roc.data.4 <- data.frame(fpr=unlist(perf.4@x.values),
                         tpr=unlist(perf.4@y.values),
                         model="GLM")

# Plot ROCs
ggplot(roc.data,aes(fpr,tpr))+geom_line(aes(color=sprintf("Dog 1 (AUC = %.4f)",auc)))+
  geom_line(data=roc.data.2,aes(color=sprintf("Dog 2 (AUC = %.4f)",auc.2)))+
  geom_line(data=roc.data.3,aes(color=sprintf("Dog 3 (AUC = %.4f)",auc.3)))+
  geom_line(data=roc.data.4,aes(color=sprintf("Dog 4 (AUC = %.4f)",auc.4)))+
  labs(color="Subject")+
  ggtitle("ROC Curves")

# example of single ROC
ggplot(roc.data, aes(x=fpr, ymin=0, ymax=tpr)) +
  geom_ribbon(alpha=0.2) +
  geom_line(aes(y=tpr)) +
  geom_line(aes(y=))
  ggtitle(paste0("ROC Curve w/ AUC=", auc))

# feature visualizations
ggplot(data=melt(fev.iif), aes(as.factor(X2), value))+xlab('eigenvalue') + geom_boxplot() + ggtitle("Interictal Eigenvalues")+ ylim(0,15)

ggplot(data=melt(fev.if), aes(as.factor(X2), value))+xlab('eigenvalue') + geom_boxplot() + ggtitle("Ictal Eigenvalues")+ ylim(0,15)


