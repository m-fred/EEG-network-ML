# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 00:48:26 2019

@author: Matthew
"""

import numpy as np
import pandas as pd
#import tigramite
from tigramite import data_processing as pp
#from tigramite import plotting as tp
from tigramite.pcmci import PCMCI
from tigramite.independence_tests import ParCorr#, GPDC, CMIknn, CMIsymb
#from tigramite.models import LinearMediation, Prediction
#import cython
#import sklearn
import sys
import os
from scipy.signal import butter, lfilter
import backboning

#ids = [1461, 2117, 3464, 3757, 6527, 6568, 7063, 7574, 7608, 7634, 7771,
 #      7943, 2261, 2645, 6140, 6227, 6232, 6255, 6383, 6395, 6396, 7577, 7890,
  #     8368] #patient number

#Band Pass filtering functions
def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

def getrmat(path):
    m = np.loadtxt(path)
    m = m.T
    m = m - np.mean(m,axis=0)
    m = m/np.std(m,axis=0)
    e = np.shape(m)[1]
    print("Filtering...")
    for j in range(e):
        m[:,j] = butter_bandpass_filter(m[:,j],3,42,1000,order=5)
    return m
  
def parcorr_fn(dataframe,var_names,tau_max):
    parcorr = ParCorr(significance='analytic')
    pcmci = PCMCI(dataframe=dataframe, cond_ind_test=parcorr, var_names=var_names, verbosity=0)
    results = pcmci.run_pcmci(tau_max=tau_max, pc_alpha=None)
    val_matrix = results['val_matrix']
    q_matrix = pcmci.get_corrected_pvalues(p_matrix=results['p_matrix'], fdr_method='fdr_bh')
    link_matrix = pcmci._return_significant_parents(pq_matrix=q_matrix, val_matrix=val_matrix, alpha_level=0.01)['link_matrix']
    return(val_matrix, link_matrix)

def matconv_fn(val_matrix, link_matrix):
    dims = np.shape(val_matrix)[1]
    combined_matrix = np.multiply(val_matrix, link_matrix) #numeric*boolean matrix
    final_matrix = np.zeros((dims,dims)) #initialise
    for i in np.arange(dims):
        for j in np.arange(dims):
            final_matrix[i,j] = max(combined_matrix[i,j,:], key=abs) #selects biggest magnitude value across all lags
    return(final_matrix)

def mat2edge(m):
    N = np.shape(m)[0]
    elist = np.zeros(3)
    for i in range(N):
        for j in range(N):
            if(m[i][j]!=0):
                edg = np.array([i,j,m[i][j]])
                elist = np.vstack((elist,edg))
    elist = elist[1:]
    return elist

def edge2mat(m):
    N = 61#int(np.max(m[:,:2]))
    mat = np.zeros((N+1,N+1))
    for i in range(len(m)):
        mat[int(m[i,0])][int(m[i,1])] = m[i,2]
    return mat

def backbone(mat):
    elist = mat2edge(mat)
    table = pd.DataFrame(elist)
    #nnodes = np.shape(mat)[0]
    #nnedges = len(elist)
    threshold_value = 0.05
    table.columns = ["src","trg","nij"]
    nc_table = backboning.noise_corrected(table)
    nc_backbone = backboning.thresholding(nc_table, threshold_value)
    elist_backbone = np.array(nc_backbone.iloc[:,:3])
    mat_backbone = edge2mat(elist_backbone)
    return mat_backbone

pat = int(sys.argv[1]) - 1
datapath = os.path.join(os.getcwd(),'Data')
flist = os.listdir(datapath)
file = flist[pat]
print(file)
print('Getting Data')
raw = getrmat(os.path.join(datapath,file))
print('Got Data')
data = raw[:,:]
#data = data.swapaxes(0, 1) #flip shape of array
#var_names = np.load('ch_names'+str(pat)+'.npy').tolist()

n = np.shape(data)[0]
dims = np.shape(data)[1]
var_names = np.arange(dims)
tau_max = 5

i = int(sys.argv[2]) - 1

if(int(i*1000)>n-1000):
    exit(0)

start = int(i*1000)
end = start + 1000
data1 = data[start:end,]
dataframe1 = pp.DataFrame(data1)
print('Runing ParCorr')
val_matrix, link_matrix = parcorr_fn(dataframe1,var_names,tau_max)
final_matrix = matconv_fn(val_matrix, link_matrix)
print('Backboning')
final_matrix = backbone(final_matrix)
dest = os.path.join(os.getcwd(),'Networks',file[:-4]+'_'+str(i)+'_.txt')
print('Saving')
np.savetxt(dest, final_matrix)


'''    
def cmi_knn_fn(dataframe,var_names,tau_max):
    cmi_knn = CMIknn(significance='shuffle_test', knn=0.1, shuffle_neighbors=5)
    pcmci_cmi_knn = PCMCI(dataframe=dataframe, cond_ind_test=cmi_knn, var_names=var_names, verbosity=0)
    results = pcmci_cmi_knn.run_pcmci(tau_max=tau_max, pc_alpha=0.05)
    val_matrix = results['val_matrix']
    p_matrix = results['p_matrix']
    link_matrix = pcmci._return_significant_parents(pq_matrix=p_matrix, val_matrix=val_matrix, alpha_level=0.01)['link_matrix']
    return(val_matrix, link_matrix)
'''    

    

