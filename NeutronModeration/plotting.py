# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 09:38:05 2022

@author: Jan Chojnacki

Neutron moderation plots 
"""

import os
import numpy as np
import pandas as pd
#import random as rand
from matplotlib import pyplot as plt

path = "Desktop/NeutronModeration"
os.chdir(path)

data = pd.read_csv("raw_data3.csv")
properties = ['Width','Length','Iterations']
data[properties] = data[properties].fillna('')
LENGTH = data.loc[0].Length
WIDTH = data.loc[0].Width
#%%
#histogram .that traces number of times neutron appeared in the bin
counts, yedges, xedges = np.histogram2d(np.array(data.y),np.array(data.x), bins=[40,40])
#total energy tracing
Total_energy, _ ,_ = np.histogram2d(np.array(data.y),np.array(data.x),weights=np.array(data.E), bins=[40,40])
#mean energy per bin
#nan_to_num replaces invalid values coming from E/0
plt.figure(figsize=(10,8))
plt.rc("axes",titlesize=20)
plt.rc("axes",labelsize=14)
plt.rc("xtick",labelsize=14)
plt.rc("ytick",labelsize=14)
plt.rc("legend",fontsize=15)
plt.pcolormesh(xedges,yedges, np.nan_to_num(Total_energy/counts, nan=0.0, posinf=0.0, neginf=0.0), cmap='coolwarm')
plt.colorbar().set_label('Energy (eV)', rotation=270,labelpad=20)
plt.xlim(0,LENGTH)
plt.ylim(0,WIDTH)
plt.title("Average neutron energy histogram \n undirected source")
plt.show()
#%%
#Histogram of neutron counts per bin 
import matplotlib.colors as colors
plt.figure(figsize=(10,8))

plt.rc("axes",titlesize=20)
plt.rc("axes",labelsize=14)
plt.rc("xtick",labelsize=14)
plt.rc("ytick",labelsize=14)
plt.rc("legend",fontsize=15)
plt.hist2d(np.array(data.x), np.array(data.y),norm=colors.LogNorm(), bins=[40,40],cmap='coolwarm')
plt.colorbar().set_label('Number of neutrons', rotation=270,labelpad=10)

plt.xlim(0,LENGTH)
plt.ylim(0,WIDTH)
plt.title("Neutron distribution histogram \n undirected source")
plt.show()     
#%%
#Mean free path plot
import matplotlib
mean_free_path_data = pd.read_csv("mean_free_path.csv")

plt.figure(figsize=(10,8))
plt.rc("axes",titlesize=20)
plt.rc("axes",labelsize=14)
plt.rc("xtick",labelsize=14)
plt.rc("ytick",labelsize=14)
plt.rc("legend",fontsize=15)
#plt.figure().set_xticks([i*10**6 for i in range(1,11)])
#plt.set_xticklabels([f"{i}" for i in range(1,11)])
#plt.ticklabel_format(style="scientific",axis="x",scilimits=(0,0))
plt.grid()
plt.plot(mean_free_path_data.Energy,mean_free_path_data['Mean free path'])
plt.title("Mean free path of $n^0$ in $H_2O$")
plt.xlabel("Energy (eV)")
plt.ylabel("Path (cm)")
