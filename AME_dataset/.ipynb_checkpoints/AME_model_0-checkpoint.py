import tensorflow as tf
from tensorflow.keras.layers import Conv2D, Input
from sklearn.model_selection import train_test_split
#############################################
import numpy as np
import healpy as hp
from astropy.io import fits
import matplotlib.pyplot as plt
import sys
###############################
sys.path.insert(1, '/media/BINGODATA1/ComponentSeparation/beam_analyzes/scripts')
import handling_data as hdata
import beam_modelling as model

# Let's define the Peak to Signal Ratio metric
def psnr_metric(y_true, y_pred):
    return tf.reduce_mean(tf.image.psnr(y_true, y_pred, max_val=1.0))

# Load data
X = np.load('/media/BINGODATA1/BeamAnalysis/AME_dataset/g_smoothed_maps.npy', mmap_mode='r') # features (inputs)
y = np.load('/media/BINGODATA1/BeamAnalysis/AME_dataset/og_maps.npy', mmap_mode='r') # labels (outputs)
print('Data loaded')

print('Spliting data...')
# Split train and test data
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1)
print('Data split')

# Reshape data
X_train = np.expand_dims(np.squeeze(X_train, axis=1), axis=-1)
y_train = np.expand_dims(np.squeeze(y_train, axis=1), axis=-1)
X_test = np.expand_dims(np.squeeze(X_test, axis=1), axis=-1)
y_test = np.expand_dims(np.squeeze(y_test, axis=1), axis=-1)
print('Data reshaped')

# Create model
input_shape = (256, 256, 1)

model_0 = tf.keras.Sequential([
    Input(shape=input_shape, name='Input'),
    
    # Encoder
    Conv2D(64, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same', name='Conv2D_1'),
    Conv2D(64, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same', name='Conv2D_2'),
    Conv2D(64, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same', name='Conv2D_3'),
    Conv2D(64, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same', name='Conv2D_4'),
    Conv2D(64, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same', name='Conv2D_5'),

    # Output layer
    Conv2D(1, (3, 3), activation='sigmoid', padding='same', name='Output')
])

# Compile model
model_0.compile(
    optimizer='adam',
    loss='mean_squared_error',
    metrics=['accuracy', psnr_metric]
)

# Fit model
print('Fitting model...')
model_0_history = model_0.fit(X_train, y_train, epochs=200, validation_data=(X_test, y_test))

# Save model
print('Saving model...')
save_path = '/media/BINGODATA1/BeamAnalysis/AME_dataset/AME_model_0.ckpt'
model_0.save(save_path=save_path)

# Create dataframe of history
print('Saving history dataframe...')
model_0_dataframe = pd.DataFrame(model_0_history.history)
save_df = '/media/BINGODATA1/BeamAnalysis/AME_dataset/AME_model_0_history.csv'
model_0_dataframe.to_csv(save_df, sep='\t')

print('All done!')