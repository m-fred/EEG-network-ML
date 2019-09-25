# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 00:48:26 2019

@author: Matthew
"""

import numpy as np
import tigramite
from tigramite import data_processing as pp
from tigramite import plotting as tp
from tigramite.pcmci import PCMCI
from tigramite.independence_tests import ParCorr, GPDC, CMIknn, CMIsymb
from tigramite.models import LinearMediation, Prediction
import cython
import sklearn
import sys

pat = '6140' #patient number


raw = np.load('pat'+str(pat)+'.npy')
data = raw[:,:]
data = data.swapaxes(0, 1) #flip shape of array
var_names = np.load('ch_names'+str(pat)+'.npy').tolist()


def parcorr_fn(dataframe,var_names,tau_max):
    parcorr = ParCorr(significance='analytic')
    pcmci = PCMCI(dataframe=dataframe, cond_ind_test=parcorr, var_names=var_names, verbosity=0)
    results = pcmci.run_pcmci(tau_max=tau_max, pc_alpha=None)
    val_matrix = results['val_matrix']
    q_matrix = pcmci.get_corrected_pvalues(p_matrix=results['p_matrix'], fdr_method='fdr_bh')
    link_matrix = pcmci._return_significant_parents(pq_matrix=q_matrix, val_matrix=val_matrix, alpha_level=0.01)['link_matrix']
    return(val_matrix, link_matrix)
    
def cmi_knn_fn(dataframe,var_names,tau_max):
    cmi_knn = CMIknn(significance='shuffle_test', knn=0.1, shuffle_neighbors=5)
    pcmci_cmi_knn = PCMCI(dataframe=dataframe, cond_ind_test=cmi_knn, var_names=var_names, verbosity=0)
    results = pcmci_cmi_knn.run_pcmci(tau_max=tau_max, pc_alpha=0.05)
    val_matrix = results['val_matrix']
    p_matrix = results['p_matrix']
    link_matrix = pcmci._return_significant_parents(pq_matrix=p_matrix, val_matrix=val_matrix, alpha_level=0.01)['link_matrix']
    return(val_matrix, link_matrix)
    
def matconv_fn(val_matrix, link_matrix):
    dims = np.shape(val_matrix)[1]
    combined_matrix = np.multiply(val_matrix, link_matrix) #numeric*boolean matrix
    final_matrix = np.zeros((dims,dims)) #initialise
    for i in np.arange(dims):
        for j in np.arange(dims):
            final_matrix[i,j] = max(combined_matrix[i,j,:], key=abs) #selects biggest magnitude value across all lags
    return(final_matrix)
    
n = np.shape(data)[0]
dims = np.shape(data)[1]
tau_max = 5

i = int(sys.argv[1])

data1 = data[i:(i+1000),]
dataframe1 = pp.DataFrame(data1)
val_matrix, link_matrix = parcorr_fn(dataframe1,var_names,tau_max)
final_matrix = matconv_fn(val_matrix, link_matrix)
np.save('T_mat'+str(i)+'_pat'+str(pat), final_matrix)
