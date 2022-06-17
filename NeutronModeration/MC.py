# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 17:01:18 2022

@author: Jan Chojnacki

Mean free path approach taken from https://www.researchgate.net/publication/318413306_Monte_Carlo_studies_on_neutron_interactions_in_radiobiological_experiments

"""

import os
import numpy as np
import pandas as pd
import random as rand
from matplotlib import pyplot as plt



MASS_H = 931.49432e6
MASS_O = 14891.26551081036e6
MASS_n = 939.56542052e6 
E_init = 1.641063e+06
MCS = 2.108162e-24
sigma_O = 2.108162e-24
sigma_H = 3.186791e-24
E_min = 0.1*E_init

WIDTH = 10
LENGTH = 10
X_init = 0.01
Y_init = WIDTH/2

b = 10**(-24) #barn to cm^2
M_H20 = 18 #g/mol
rho_H20 = 1 #g/cm^3

class Element:
    N_Av = 6.023e23 #1/mol
    
    def __init__(self, MASS, CROSS_SECTION_DATA):
       self.Mass = MASS
       self.cross_section_data = CROSS_SECTION_DATA
       self.cross_section_data["Cross-section"] = b*np.array(self.cross_section_data["Cross-section"])
    def cross_section(self,E):
       return float(self.cross_section_data.iloc[(self.cross_section_data["Incident energy"]-E).abs().argsort()][:1]["Cross-section"])
       

#cd Desktop\NeutronModeration

path = "Desktop/NeutronModeration"
os.chdir(path)

Hydrogen = Element(1.00794, pd.read_csv("H-n.csv"))
Oxygen =  Element(15.9994, pd.read_csv("O-n.csv"))


#uniformly distributed post-collision angle of the neutron 
def random_angle(E):
    angle =np.random.choice([-1,1])* np.arccos(np.sqrt(np.abs(E)/E_init)) 
    
    return(angle + np.random.normal(loc = angle, scale = 0.1))
#average path lenght between collisions
def mean_free_path(MCS):
    return(- 1/MCS*np.log(rand.random()))
#needed in in mean free path
def macroscopic_cross_section(sigma_H,sigma_O):
    return(rho_H20*Element.N_Av/M_H20 *(2*sigma_H+sigma_O))
#Energy spectrum of Americium
def initial_distribution():
    return(E_init)
#picking target nuclei either H or O
def random_target(E):
    a = np.array([MASS_H, MASS_O])
    p_H = 2 *Hydrogen.cross_section(E)/(2*Hydrogen.cross_section(E)+ Oxygen.cross_section(E))
    p_O = Oxygen.cross_section(E)/(2*Hydrogen.cross_section(E)+ Oxygen.cross_section(E))
    return np.random.choice(a, p=[p_H,p_O])
#Energy of a neutron scattered at "angle" and initial energy E 
def energy(E,angle,M):
    m = MASS_n
    return(m/2*  (1/(1 + m/M) * (np.sqrt(2*m*E)/M *np.cos(angle) + np.abs(2* m*E/M**2 *np.cos(angle)**2 + 2*(1/m - 1/M)*E)**(1/2)))**2)

def inbox_condition(x,y,E):
    return(0 <= x <= LENGTH and 0 <= y <= WIDTH and E > E_min)

#Data Frame to which we will pass generation number, x coordinate, y coordinate, particle's energy
data = pd.DataFrame(columns = ["#", "x", "y", "E"])

#%%
#progress bar
iteration_number = 50000
import progressbar
bar = progressbar.ProgressBar(maxval=iteration_number, \
    widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
bar.start()

#main part of the program
for i in range(iteration_number):
    x = X_init
    y = Y_init
    E = initial_distribution()
    Theta = np.random.uniform(-np.pi/2,np.pi/2)
    bar.update(i)
    #data.loc[len(data.index)] = [i,x,y,E]
    while inbox_condition(x,y,E):
        data.loc[len(data.index)] = [i,x,y,E]
        Lambda = mean_free_path(macroscopic_cross_section(Hydrogen.cross_section(E), Oxygen.cross_section(E)))
        
        Theta_n = random_angle(E) + Theta
        x_n = x + Lambda*np.cos(Theta_n)
        y_n = y + Lambda*np.sin(Theta_n)
        E_n = energy(E,Theta_n,random_target(E))
        
        Theta = Theta_n
        x = x_n
        y = y_n
        E = E_n
bar.finish()        

#path plot
for i in range(10):
    plt.plot(data[data["#"]==i]["x"],data[data["#"]==i]["y"])
    plt.ylim((0,WIDTH))
    plt.xlim((0,LENGTH))        

#%%
properties = ['Width','Length','Iterations']
data.loc[0,properties] = [WIDTH, LENGTH, iteration_number]
data[properties] = data[properties].fillna('')

from os import path
data.to_csv(path.join("C:/Users/Jan Chojnacki/Desktop/NeutronModeration/",'raw_data3.csv'), index=False) 
#%%
# mean free path as a function of energy

energy_vector = np.arange(4e5,1e7, 1000)
mean_free_path_df = pd.DataFrame({"Energy":energy_vector, "Mean free path": 1/macroscopic_cross_section(np.vectorize(Hydrogen.cross_section)(energy_vector), np.vectorize(Oxygen.cross_section)(energy_vector)) })
mean_free_path_df.to_csv(path.join("C:/Users/Jan Chojnacki/Desktop/NeutronModeration/",'mean_free_path.csv'), index=False) 
#plt.plot(energy_vector, 1/macroscopic_cross_section(np.vectorize(Hydrogen.cross_section)(energy_vector), np.vectorize(Oxygen.cross_section)(energy_vector)))        
#%%
angle_data = [random_angle(E_init) for i in range(1000)]
plt.hist(angle_data)