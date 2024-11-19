from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
import sys
import keras
import tensorflow as tf
from keras.models import Sequential
from keras.layers import Dense, LSTM, Dropout, Bidirectional, LeakyReLU, BatchNormalization
from keras.layers import Convolution1D, MaxPooling1D
from keras.layers.core import Flatten
from keras.callbacks import ModelCheckpoint, EarlyStopping
from keras import regularizers
from keras.regularizers import L1L2
from keras import optimizers
from sklearn import metrics
from keras.constraints import max_norm,MinMaxNorm,NonNeg
np.random.seed(12345)



### training set
# loading the input sequence vector and labels

acgt2vec = {'A': [1, 0, 0, 0],
            'C': [0, 1, 0, 0],
            'G': [0, 0, 1, 0],
            'T': [0, 0, 0, 1],
            'N': [0, 0, 0, 0]}

allWindows = [] ### vector of sequence
trainLabels = []

allWindows0 = [] ### vector of sequence
validLabels = []



folder = sys.argv[1] #### /path/to/your/data/
train_input = sys.argv[2] #### "LNCap_after_DHT_all_DHS.overlapH3K27ac.center1KB.withDHSctrl.training.DLinput.shuffled"
valid_input = sys.argv[3] #### "LNCap_after_DHT_all_DHS.overlapH3K27ac.center1KB.withDHSctrl.validation.DLinput.shuffled"
### training set
file = folder +  train_input #### training input, you can replace it with your file name

with open(file) as IN1:
    for line in IN1:
        seq = line.rstrip().split("\t")[0]
        label = line.rstrip().split("\t")[2]
        trainLabels.append(label)
        
        ntVecMat = []
        for nt in seq:
            ntVec = acgt2vec[nt]
            ntVecMat.append(ntVec)
        allWindows.append(ntVecMat)
        
X_train = np.array(allWindows)
print(X_train.shape)

intTrainLabels = [int(x) for x in trainLabels]
y_train = np.array(intTrainLabels)
print(y_train.shape)

### validation set
file0 = folder + valid_input #### validation input, you can replace it with your file name

with open(file0) as IN0:
    for line in IN0:
        seq = line.rstrip().split("\t")[0]
        label = line.rstrip().split("\t")[2]
        validLabels.append(label)
        
        ntVecMat = []
        for nt in seq:
            ntVec = acgt2vec[nt]
            ntVecMat.append(ntVec)
        allWindows0.append(ntVecMat)
        
X_valid = np.array(allWindows0)
#print(X_valid.shape)

intValidLabels = [int(x) for x in validLabels]
y_valid = np.array(intValidLabels)
#print(y_valid.shape)



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

# building the model

model = Sequential()

model.add(Convolution1D(320, 8, activation = "relu", input_shape = (X_train.shape[1], 4),
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

model.compile(optimizer='adam', loss='binary_crossentropy', metrics=[f1_m, precision_m, recall_m, 'acc'])
print(model.summary())


import h5py
#from sklearn.utils import class_weight
#class_weights = class_weight.compute_class_weight('balanced',
#                                                 np.unique(intTrainLabels),
#                                                 intTrainLabels)
#print(class_weights)

outPath = "./LNCap_after_DHT_all_DHS.overlapH3K27ac.center1KB_DHSctrl.bestmodel.considerRC.hdf5"
checkpointer = ModelCheckpoint(filepath=outPath, verbose=1, save_best_only=True)
earlystopper = EarlyStopping(monitor='val_loss', patience=50, verbose=1)
print('Training model...')
history = model.fit(X_train, y_train, epochs=1000, batch_size=1000, shuffle=True, 
                    validation_data=(X_valid, y_valid), 
                    callbacks=[checkpointer,earlystopper], verbose=1)
                    





