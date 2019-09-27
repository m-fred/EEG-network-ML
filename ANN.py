# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 13:55:00 2019

@author: hardi
"""

import tensorflow as tf
#import pickle
from sklearn import preprocessing
import pandas as pd
import numpy as np

n_nodes_hl1 = 300
n_nodes_hl2 = 300
n_nodes_hl3 = 300

n_features = 23
n_classes = 2
batch_size = 1000

x = tf.placeholder('float',[None,n_features])
y = tf.placeholder('float')

def neural_network_model(data):
    hidden_1_layer = {'weights':tf.Variable(tf.random_normal([n_features,n_nodes_hl1])),
                      'biases':tf.Variable(tf.random_normal([n_nodes_hl1]))}
    
    hidden_2_layer = {'weights':tf.Variable(tf.random_normal([n_nodes_hl1,n_nodes_hl2])),
                      'biases':tf.Variable(tf.random_normal([n_nodes_hl2]))}
    
    hidden_3_layer = {'weights':tf.Variable(tf.random_normal([n_nodes_hl2,n_nodes_hl3])),
                      'biases':tf.Variable(tf.random_normal([n_nodes_hl3]))}
    
    output_layer = {'weights':tf.Variable(tf.random_normal([n_nodes_hl3,n_classes])),
                      'biases':tf.Variable(tf.random_normal([n_classes]))}
    
    l1 = tf.add(tf.matmul(data,hidden_1_layer['weights']), hidden_1_layer['biases'])
    l1 = tf.nn.relu(l1)
    
    l2 = tf.add(tf.matmul(l1,hidden_2_layer['weights']), hidden_2_layer['biases'])
    l2 = tf.nn.relu(l2)
    
    l3 = tf.add(tf.matmul(l2,hidden_3_layer['weights']), hidden_3_layer['biases'])
    l3 = tf.nn.relu(l3)
    
    output = tf.matmul(l3,output_layer['weights']) + output_layer['biases']
    
    return output

def train_neural_network(x_train,y_train,x_test,y_test):
    prediction = neural_network_model(x)
    cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits_v2(logits=prediction,labels=y))
    optimizer = tf.train.AdamOptimizer().minimize(cost)
    hm_epochs = 20
    with tf.Session() as sess:
        sess.run(tf.initialize_all_variables())
        
        for epoch in range(hm_epochs):
            epoch_loss = 0
            for i in range(int(len(x_train)/batch_size)):
                start = int(i*batch_size)
                end = start + batch_size
                epoch_x = x_train[start:end,:]
                epoch_y = y_train[start:end]
                _, c = sess.run([optimizer,cost], feed_dict = {x: epoch_x, y: epoch_y})
                epoch_loss += c
            #print('Epoch',epoch,'completed out of',hm_epochs,'loss:',epoch_loss)
        
        correct = tf.equal(tf.argmax(prediction,1),tf.argmax(y,1))
        accuracy = tf.reduce_mean(tf.cast(correct,'float'))
        print('Accuracy:',accuracy.eval({x:x_test,y:y_test}))
        return accuracy.eval({x:x_test,y:y_test})

def get_class(a):
    c = np.zeros((len(a),2))
    for i in range(len(a)):
        if (a[i]==0):
            c[i] = np.array([1,0])
        else:
            c[i] = np.array([0,1])
    return c

def train_test_split(data,test_sub=1):
    ids = np.array(data['subject_id'])
    ids = np.unique(ids)
    #data.loc[:,:] = preprocessing.scale(np.array(data.loc[:,:])) 
    #N_sub = len(ids)
    train_sub = np.delete(ids,test_sub)
    #train_sub = np.random.choice(ids,N_sub-test_sub,replace=False)
    test_sub = np.array([ids[test_sub]])#
    #test_sub = np.array([i for i in ids if i not in train_sub])
    train_set = data.loc[data['subject_id'].isin(train_sub)]
    train_set = train_set.sample(frac=1)
    test_set = data.loc[data['subject_id'].isin(test_sub)]
    test_set = test_set.sample(frac=1)
    x_train = np.array(train_set.drop(columns=['subject_id','schizophrenia'],axis=0))
    #x_train = preprocessing.scale()
    x_test = np.array(test_set.drop(columns=['subject_id','schizophrenia'],axis=0))
    y_train = np.array(train_set['schizophrenia'])
    y_train = get_class(y_train)
    y_test = np.array(test_set['schizophrenia'])
    y_test = get_class(y_test)
    return x_train,y_train,x_test,y_test

    
if __name__ == "__main__":
    df = pd.read_csv('metricsdata_pc_recon.txt',delimiter='\t',header=None)
    df.columns = ['subject_id','state','in_deg_mean','in_deg_std','out_deg_std','deg_dif_std','in_str_mean','in_str_std','out_str_std','str_dif_std','EBC_mean','EBC_std','NBC_mean','NBC_std','dens','K','E_glob','Q','C_mean','C_std','R_oi','R_io','R_oo','R_ii','schizophrenia']
    ids = np.array(df['subject_id'])
    #ids = [i[:4] for i in ids]
    #df['subject_id'] = ids
    df.iloc[:,:-2] = preprocessing.scale(np.array(df.iloc[:,:-2])) 
    #x_train,y_train,x_test,y_test = train_test_split(df,test_sub=2)
    #train_neural_network(x_train,y_train,x_test,y_test)
    a = np.zeros(16)
    for i in range(16):
        print(i)
        x_train,y_train,x_test,y_test = train_test_split(df,test_sub=i)
        a[i] = train_neural_network(x_train,y_train,x_test,y_test)
    

