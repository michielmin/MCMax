#!/home/lucia3/anaconda/bin/python

import numpy as np
import scipy
import matplotlib as mpl
import ml_read_data as mlr
import ml_plot_model_CD as mlp
#mpl.use('ps')
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LogNorm

# path to the folder containing all subfolders
common_path = '/home/klarmann/Desktop/MCMax_examples/'

# list of all subfolders, e.g. sh1, sh2, sh3 etc.
subfolders = ['contineous_evaporation_dust_composition']

# list of models per subfolder, e.g. rim1, rim2, rim3
model_names =[ 'a1e-1_la1e2_r0.0001_L40', 'f_1e0_r0.01_L8']

#it is assumed that each subfolder contains the same list of models
# e.g. models with 3 different scaheights and 3 different rim positions => 9 models
# subfolder sh1 contains rim1, rim2, rim3
# subfolder sh2 contains rim1, rim2, rim3
# subfolder sh3 contains rim1, rim2, rim3
# In this example, there is only one subfolder

for subfolder in subfolders:
    path = common_path + subfolder + '/'
    print subfolder
    print type(path), 'path'
    for model in model_names:
        spectrum_name = 'no_spectrum'# no raytraced SPECTRUM data has been calculated with MCMax, replace with name of SPECTRUM file. Like that, RT spectrum at 44.2 deg is plotted
        #basevis = '' # uncomment if there is no basevis file
        basevis = 'basevis45.0MATISSE_L.dat' # visibity vs baseline output has been created (basevis file)
        timestep = '' #shows the last iteration by default
        PAH_free = '' #allows and additional raytraced SPECTRUM to be plotted, provide name of spectrum file
        print model
        # if-else example if not all models are in all subfolders
        # if subfolder == 'mdust3.2e-10_1ym_flat':
        #     if (model == 'rin0.005' or model == 'rin0.008'):
        #         pass
        #     else:
        #         mlp.model_multiplot(path, model, spectrum_name)
        # else:
        #     mlp.model_multiplot(path, model, spectrum_name)

        mlp.model_multiplot(path, model, spectrum_name=spectrum_name, timestep = timestep, basevis=basevis, PAH_free=PAH_free)
