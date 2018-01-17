#!/usr/bin/python
#-------------------------------------------------------------------------------
#Giving the path where gameModelPars5 and gamePlotUtils6 are located
#%cd ~/Documents/Leiden_University/3_Semester/MSc_thesis/Codes/Code_v15
'''
Arecibo:
Sources that satisfy Arecibo limits= 76 (Observations)
Index factor of the observed flux density function -0.471161157475 (Observation)
################################################################################
MMB:
Sources that satisfy MMB limits= 908 (Observations)
Index factor of the observed flux density function -0.581896115547 (Observation)
'''
#-------------------------------------------------------------------------------
#Importing Libraries
import matplotlib
import sys
import numpy
import csv
import ephem
global debug, interactive
debug = 0
interactive = 1
#if (interactive):
#    matplotlib.use('MacOSX')
#else:
#    matplotlib.use('pdf')
from pylab import *
from math import *
from random import *
from time import clock, time, localtime
from aux.gameModelPars31 import *
from aux.gamePlotUtils31 import *
from aux.gameFunctions31 import *
from aux.GaMe_v31 import *
from aux.export_txt_Reid31 import *
import os
from multiprocessing import Process
import subprocess
import shutil
from joblib import Parallel, delayed  
import multiprocessing
#-------------------------------------------------------------------------------
def calculate_index(inp):
    from_num=900
    to_num=1800
    from_ind=-1.1
    to_ind=-2.0
    total=100
    num_masers=numpy.linspace(from_num,to_num,total)
    ind=numpy.linspace(from_ind,to_ind,total)
    #for k in range(number_iterations):
    Matrix=numpy.zeros([total**2,10])
    w=0
    for i in range(len(num_masers)):
        for j in range(len(ind)):
            #print('Iteration:',inp+1,'of',number_iterations)
            print('Simulation:',w+1,'of',len(ind)*len(num_masers))
            print('Number:',int(num_masers[i]))
            print('Index:',ind[j])
            #gsample,par,file_masers,MMB,Arecibo=full_model(num_masers[i],ind[j])
            gsample,par,file_masers,fname_directory,MMB,Arecibo,prt1,prt2,prt3,prt7,prt8,prt9,prt10=full_model(file_name,int(num_masers[i]),ind[j])
            #gsample,par,file_masers,fname_directory,MMB,Arecibo,prt1,prt2,prt3,prt7,prt8,prt9,prt10=full_model(file_name,num_masers,ind)
            MMB_num=MMB[0]
            MMB_ind=MMB[1]
            MMB_chi2=MMB[2]
            MMB_KLD=MMB[3]            
            Arecibo_num=Arecibo[0]
            Arecibo_ind=Arecibo[1]
            Arecibo_chi2=Arecibo[2]
            Arecibo_KLD=Arecibo[3]
            Matrix[w,:]=[int(num_masers[i]),ind[j],MMB_num,MMB_ind,MMB_chi2,MMB_KLD,Arecibo_num,Arecibo_ind,Arecibo_chi2,Arecibo_KLD]
            w+=1    
    numpy.savetxt('../Matrix_ind_num/Calc_Num_Ind_'+str(inp)+'.txt',Matrix)
    return

           
def calculate_err(iterations):
    MMB_numbers=[]
    MMB_alpha=[]    
    Arecibo_numbers=[]    
    Arecibo_alpha=[]    
    w=0
    for i in range(iterations):
        print('Simulation %s of %s'%(w+1,iterations))
        #gsample,par,file_masers,MMB,Arecibo=full_model(num_masers[i],ind[j])
        gsample,par,file_masers,fname_directory,MMB,Arecibo,prt1,prt2,prt3,prt7,prt8,prt9,prt10=full_model(file_name)
        MMB_numbers.append(MMB[0])
        MMB_alpha.append(MMB[1])
        Arecibo_numbers.append(Arecibo[0])
        Arecibo_alpha.append(Arecibo[1])
        w+=1    
        matrix=[MMB_numbers,MMB_alpha,Arecibo_numbers,Arecibo_alpha]
    numpy.savetxt('../surface_bestnalpha/One_same_galaxy/Cube_One_Galaxy.txt',matrix)
    return                          

#Questions                                            
file_name = raw_input("Please enter the parameter's file: ")
number_iterations=int(raw_input("How many iterations? (height of the cube): "))
confirmation_function=raw_input("grid or one-same-galaxy?: ")
inputs=range(number_iterations)

#Calculation
if(confirmation_function=='grid'):
    #Number of cores
    num_cores_max = multiprocessing.cpu_count()
    ans = raw_input("How many cores do you want to use (max=%s): " % (num_cores_max)) 
    if (int(ans)>num_cores_max):
        print('Max cores reached!')
        sys.exit("Error message")
    else:    
        num_cores=int(ans)    
    print("numCores = " + str(num_cores)) 
    results = Parallel(n_jobs=num_cores)(delayed(calculate_index)(i) for i in inputs)
elif(confirmation_function=='one-same-galaxy'):
    results = calculate_err(number_iterations)
else:
    print('Bad parameters selected please check your inputs')
'''
#To test without parallel        
for i in inputs:
    calculate_index(i)
'''       