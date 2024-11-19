from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
import keras
import tensorflow as tf
from keras.models import Sequential
from keras.layers import Dense, LSTM, Dropout, Bidirectional, LeakyReLU, BatchNormalization
from keras.layers import Convolution1D, MaxPooling1D
from keras.layers.core import Flatten
from keras.layers.embeddings import Embedding
from keras.callbacks import ModelCheckpoint, EarlyStopping
from keras import regularizers
from keras.regularizers import L1L2
from keras import optimizers
from sklearn import metrics

from keras.constraints import max_norm,MinMaxNorm,NonNeg
#from keras.utils import multi_gpu_model
np.random.seed(12345)

# loading the input sequence vector and labels for the testing sets
import sys

bin = sys.argv[1]

acgt2vec = {'A': [1, 0, 0, 0],
            'C': [0, 1, 0, 0],
            'G': [0, 0, 1, 0],
            'T': [0, 0, 0, 1],
            'N': [0, 0, 0, 0]}

allWindows_permut = [] ### vector of sequence
#y = []

folder = "/path/to/DLinput/data/"

prefix = "1Kgenomes_phase3_commonAFR_or_commonEUR_snps.commonChrs.top5percentFst.slidingWindowOverlapSNP.DLinput"
file = folder + prefix + ".bin_" + bin
#file = folder + "test.input"
with open(file) as IN:
    for line in IN:
        seq1 = line.rstrip().split("\t")[0]   

        ntVecMat1 = []
        for nt in seq1:
            ntVec = acgt2vec[nt]
            ntVecMat1.append(ntVec)
        
            
        allWindows_permut.append(ntVecMat1)
        
        
        
X_permute = np.array(allWindows_permut)
print(X_permute.shape)

from keras import backend as K

def recall_m(y_true, y_pred):
        true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
        possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
        recall = true_positives / (possible_positives + K.epsilon())
        return recall

def precision_m(y_true, y_pred):
        true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
        predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
        precision = true_positives / (predicted_positives + K.epsilon())
        return precision

def f1_m(y_true, y_pred):
    precision = precision_m(y_true, y_pred)
    recall = recall_m(y_true, y_pred)
    return 2*((precision*recall)/(precision+recall+K.epsilon()))
    


model = Sequential()

model.add(Convolution1D(320, 8, activation = "relu", input_shape = (1001, 4),
                        kernel_regularizer=L1L2(l1=1e-8, l2=5e-7), kernel_constraint=max_norm(1)))

model.add(Dropout(0.2))
model.add(MaxPooling1D(4, 4))


model.add(Convolution1D(240, 4, activation = "relu",
                        kernel_regularizer=L1L2(l1=1e-8, l2=5e-7), kernel_constraint=max_norm(0.9)))

model.add(Dropout(0.2))
model.add(MaxPooling1D(2, 2))


model.add(Convolution1D(240, 4, activation = "relu", kernel_regularizer=L1L2(l1=1e-8, l2=5e-7),
                        kernel_constraint=max_norm(1)))

model.add(Dropout(0.2))
model.add(MaxPooling1D(2,2))


model.add(Convolution1D(320, 4, activation = "relu", kernel_regularizer=L1L2(l1=1e-8, l2=5e-7),
                        kernel_constraint=max_norm(1)))

model.add(Dropout(0.2))
model.add(MaxPooling1D(2,2))

model.add(Convolution1D(480, 4, activation = "relu", kernel_regularizer=L1L2(l1=1e-8, l2=5e-7),
                        kernel_constraint=max_norm(1)))

model.add(Dropout(0.2))



model.add(Flatten())

model.add(Dense(180, activation='relu',
                kernel_regularizer=L1L2(l1=1e-8, l2=5e-7),kernel_constraint=max_norm(1)))
model.add(Dense(1, activation='sigmoid',
                kernel_regularizer=L1L2(l1=1e-8, l2=5e-7),kernel_constraint=max_norm(1)))


# compile the model
#model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['acc',f1_m,precision_m, recall_m])
model.compile(optimizer='adam', loss='binary_crossentropy', metrics=[f1_m, precision_m, recall_m, 'acc'])
print(model.summary())




model_enh="./LNCap_after_DHT_all_DHS.overlapH3K27ac.center1KB_DHSctrl.bestmodel.considerRC.hdf5"
model.load_weights(model_enh)
y_pred = model.predict(X_permute, verbose=1)


#outPath = "/data/lis11/CS23/"
output = folder + prefix + ".bin_" + bin + ".slidingWindowOverlapSNP.DLscore"
with open(output, 'w') as OUT:
    
    for i in range(len(y_brain_pred)):
        OUT.write(str(y_brain_pred[i]) + "\n")
