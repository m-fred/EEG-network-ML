# -*- coding: utf-8 -*-
"""
Created on Sat Oct 27 23:32:55 2018

@author: Matthew
"""
import numpy as np
import mne
import mat4py


#%%
#healthy
mne6140 = mne.io.read_raw_eeglab(r"C:\Users\Matthew\Documents\Matthew Imperial\Fourth Year\M4R\Stefan_Data\Healthy Control\mara_sasica_chanrej_filt_6140_020711_RestEyesClosed_MAPM_6.set")
mne6227 = mne.io.read_raw_eeglab(r"C:\Users\Matthew\Documents\Matthew Imperial\Fourth Year\M4R\Stefan_Data\Healthy Control\mara_sasica_chanrej_filt_6227_010810_RestEyesClosed_Unknown_2.set")
mne6232 = mne.io.read_raw_eeglab(r"C:\Users\Matthew\Documents\Matthew Imperial\Fourth Year\M4R\Stefan_Data\Healthy Control\mara_sasica_chanrej_filt_6232_010710_RestEyesClosed_Unknown_2.set")
mne6255 = mne.io.read_raw_eeglab(r"C:\Users\Matthew\Documents\Matthew Imperial\Fourth Year\M4R\Stefan_Data\Healthy Control\mara_sasica_chanrej_filt_6255_090910_RestEyesClosed_AC_4.set")
mne6383a = mne.io.read_raw_eeglab(r"C:\Users\Matthew\Documents\Matthew Imperial\Fourth Year\M4R\Stefan_Data\Healthy Control\mara_sasica_chanrej_filt_6383_012012_RestEyesClosed_PM_3.set")
mne6383b = mne.io.read_raw_eeglab(r"C:\Users\Matthew\Documents\Matthew Imperial\Fourth Year\M4R\Stefan_Data\Healthy Control\mara_sasica_chanrej_filt_6383_012012_RestEyesClosed_PM_4.set")
mne6395 = mne.io.read_raw_eeglab(r"C:\Users\Matthew\Documents\Matthew Imperial\Fourth Year\M4R\Stefan_Data\Healthy Control\mara_sasica_chanrej_filt_6395_10252010_RestEyesClosed_NT_4.set")
mne6396a = mne.io.read_raw_eeglab(r"C:\Users\Matthew\Documents\Matthew Imperial\Fourth Year\M4R\Stefan_Data\Healthy Control\mara_sasica_chanrej_filt_6396_061410_RestEyesClosed_NT_3.set")
mne6396b = mne.io.read_raw_eeglab(r"C:\Users\Matthew\Documents\Matthew Imperial\Fourth Year\M4R\Stefan_Data\Healthy Control\mara_sasica_chanrej_filt_6396_061410_RestEyesClosed_NT_7.set")
mne7577 = mne.io.read_raw_eeglab(r"C:\Users\Matthew\Documents\Matthew Imperial\Fourth Year\M4R\Stefan_Data\Healthy Control\mara_sasica_chanrej_filt_7577_10142010_RestEyesClosed_PM_7.set")
mne7890 = mne.io.read_raw_eeglab(r"C:\Users\Matthew\Documents\Matthew Imperial\Fourth Year\M4R\Stefan_Data\Healthy Control\mara_sasica_chanrej_filt_7890_02232012_RestEyesClosed_KS_5.set")


#schizophrenic
mne6527 = mne.io.read_raw_eeglab(r"C:\Users\Matthew\Documents\Matthew Imperial\Fourth Year\M4R\Stefan_Data\Proband with Schizophrenia\mara_sasica_chanrej_filt_6527_011210_RestEyesClosed_Unk_2.set")
mne7063 = mne.io.read_raw_eeglab(r"C:\Users\Matthew\Documents\Matthew Imperial\Fourth Year\M4R\Stefan_Data\Proband with Schizophrenia\mara_sasica_chanrej_filt_7063_050310_RestEyesClosed_NT_3.set")
mne7574 = mne.io.read_raw_eeglab(r"C:\Users\Matthew\Documents\Matthew Imperial\Fourth Year\M4R\Stefan_Data\Proband with Schizophrenia\mara_sasica_chanrej_filt_7574_081210_RestEyesClosed_NT_3.set")
mne7608 = mne.io.read_raw_eeglab(r"C:\Users\Matthew\Documents\Matthew Imperial\Fourth Year\M4R\Stefan_Data\Proband with Schizophrenia\mara_sasica_chanrej_filt_7608_071510_RestEyesClosed_NT_3.set")
mne7634 = mne.io.read_raw_eeglab(r"C:\Users\Matthew\Documents\Matthew Imperial\Fourth Year\M4R\Stefan_Data\Proband with Schizophrenia\mara_sasica_chanrej_filt_7634_01142011_RestEyesClosed_AC_4.set")
mne7771 = mne.io.read_raw_eeglab(r"C:\Users\Matthew\Documents\Matthew Imperial\Fourth Year\M4R\Stefan_Data\Proband with Schizophrenia\mara_sasica_chanrej_filt_7771_9152010_RestEyesClosed_ACNT_4.set")
mne7943 = mne.io.read_raw_eeglab(r"C:\Users\Matthew\Documents\Matthew Imperial\Fourth Year\M4R\Stefan_Data\Proband with Schizophrenia\mara_sasica_chanrej_filt_7943_111510_RestEyesClosed_PM_6.set")


#%%
#chunk to create files to export (takes ~15min)

pats = [6140, 6227, 6232, 6255, '6383a', '6383b', 6395, '6396a', '6396b', 7577, 7890,
            6527, 7063, 7574, 7608, 7634, 7771, 7943]
pats = [str(i) for i in pats]

npydict = {}

for i in range(len(pats)):
    pat = pats[i]
    mnename = 'mne'+pat
    mnedata = globals()[mnename]
    
    #export channel names as .npy
    np.save('ch_names'+pat, mnedata.ch_names[:-1]) #(exclude last row as this contains reference node)
    
    #export data as .npy
    np.save('pat'+pat, mnedata[:-1,:][0]) #(exclude last row as this contains reference node)
    
    #import the .npy data files in dictionary
    npydict['raw'+pat] = np.load('pat'+pat+'.npy')
    
    #export data as .mat
    n = len(mnedata.ch_names) - 1 #(to exclude reference node)
    rawdict = {}
    for i in range(n):
        rawdict[mnedata.ch_names[i].replace(" ", "")] = mnedata[i, :][0].tolist()
        
    mat4py.savemat('pat'+pat+'.mat', rawdict)
    
    