# EEG-network-ML
Machine learning classification of schizophrenia using EEG brain networks


rawconv.py : convert mne formats to .npy or .mat [RESULT e.g. pat6140.npy]

bashrun6140.py : example matrix calculation (tigramite) for patient 6140 to be run in bash

bashrun6140.m : example matrix calculation (pmime) for patient 6140 to be run in bash

pbsscripts.py : example PBS scripts to iteratively calculate matrices using the two files above ('i' range must be continually updated as job time limited) [RESULT e.g. Pmat_pat6140_xxx000.txt]

datagen.m : calculates network metrics from all matrices (uses functions downloaded from https://sites.google.com/site/bctnet/ but Python versions exist too at https://pypi.org/project/bctpy/ ) [RESULT e.g. Pmetrics_pat6140.txt]

M4R_mlfinal.Rmd : performs all machine learning and visualisations (messy due to existence of backboned 'Tmats' batch which was added later)

M4R_nnonly.Rmd : only performs neural network and made to expect only one batch of matrices. Also compares two training/testing subsetting methods.

metricdata_upc.txt : metric data for unthresholded parcorr batch

metricdata_bpc.txt : metric data for backboned parcorr batch

metricdata_pmime.txt : metric data for pmime batch
