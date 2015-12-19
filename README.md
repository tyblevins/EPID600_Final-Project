# Seizure Detection in Canine Subjects
### Tyler Blevins

## Project Description
EEG from four canine subjects was analyzed. The EEG data is publicly available (and nicely packaged/organized) on Kaggle.com as part of the Seizure Detection challenge. The data was clipped into 1 second segments and classified as 'Ictal' (seizure) or 'Interictal' (between seizure). Features were extracted from the clips, analyzed for importance, and ran through a seizure detection algorithm. The performance of the algorithm was evaluated by ROC and AUC metrics.

## Directory Information
`data` contains the training and testing data for only Dog 1 in the study. This is for use with the sample code. The data for Dogs 2-4 can be found on [Kaggle](https://www.kaggle.com/c/seizure-detection/data). Download the clips.tar file and unzip. Please be aware that the clips.tar file is ~50GB uncompressed.

`pics` contains pictures of the figures used in the final report.

`report` contains the .html and .RMD files for the final report.

`Final_Project.R` is the raw code for processing and analyzing the data.

`project-info.txt` details the instructions for the final project in EPID 600.

`Sample.R` is the raw code for running this analysis on the sample data in `data`. This project must be your working directory path for the script to execute correctly, but otherwise should be ready to go out-of-the-box.



