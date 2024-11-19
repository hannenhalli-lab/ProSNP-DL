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
from keras.callbacks import ModelCheckpoint, EarlyStopping
from keras import regularizers
from keras.regularizers import L1L2
from keras import optimizers
from sklearn import metrics
from keras.constraints import max_norm,MinMaxNorm,NonNeg
np.random.seed(12345)


# loading the input sequence vector and labels

acgt2vec = {'A': [1, 0, 0, 0],
            'C': [0, 1, 0, 0],
            'G': [0, 0, 1, 0],
            'T': [0, 0, 0, 1],
            'N': [0, 0, 0, 0]}



allWindows1 = [] ### vector of sequence
testLabels = []

folder = "/path/to/your/data/"


## testing set
file1 = folder + "LNCap_after_DHT_all_DHS.overlapH3K27ac.center1KB.withDHSctrl.testing.DLinput.shuffled" ### testing input, you can replace it with your file name

with open(file1) as IN2:
    for line in IN2:
        seq = line.rstrip().split("\t")[0]
        label = line.rstrip().split("\t")[2]
        testLabels.append(label)
        
        ntVecMat = []
        for nt in seq:
            ntVec = acgt2vec[nt]
            ntVecMat.append(ntVec)
        allWindows1.append(ntVecMat)
        
X_test = np.array(allWindows1)
#print(X_test.shape)

intTestLabels = [int(x) for x in testLabels]
y_test = np.array(intTestLabels)
#print(y_test.shape)

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

model.add(Convolution1D(320, 8, activation = "relu", input_shape = (X_test.shape[1], 4),
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

                    

#test model
import h5py
from sklearn.metrics import roc_curve
from sklearn.metrics import auc
print('Testing model...')

model.load_weights(outPath)
tresults = model.evaluate(X_test, y_test, batch_size=1000)
print(tresults)
y_pred = model.predict(X_test, batch_size=1000, verbose=1)

print('Calculating AUC...')
auroc = metrics.roc_auc_score(y_test, y_pred)
#auprc = metrics.average_precision_score(y_test, y_pred)
print(auroc)
#print(auprc)


import matplotlib
matplotlib.get_backend()
from matplotlib import pyplot as plt
### plot ROC curve

    
fpr, tpr, threshold = metrics.roc_curve(y_test, y_pred)
plt.plot(fpr,tpr,lw=2,label="auROC="+str(auroc))
#print(len(fpr2), len(tpr2), len(threshold2))
    
legend = 'auROC = %0.2f' % auroc


plt.plot([0, 1], [0, 1], 'k--')
plt.xlabel('False positive rate')
plt.ylabel('True positive rate')
plt.title('ROC curve')
plt.legend([legend], loc='best')

plt.savefig("./ROC_curve_LNCap_after_DHT.png", dpi=300)
plt.show()


