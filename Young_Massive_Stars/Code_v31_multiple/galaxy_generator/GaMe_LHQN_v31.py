# -*- coding: utf-8 -*-
################################################################################################################
#Importing Libraries
################################################################################################################
import matplotlib
import sys
import numpy
import csv
import ephem
global debug, interactive
debug = 0
interactive = 1
from pylab import *
from math import *
from random import *
import os
from time import clock, time, localtime
import shutil
import astropy
from astropy.io import ascii
from astropy.table import Table, Column
import scipy as scipy
from aux.gameModelPars31 import *
from aux.gamePlotUtils31 import *
from aux.gameFunctions31 import *
from aux.GaMe_v31 import *
from aux.export_txt_Reid31 import *
import re
################################################################################################################
#Defining the parameter txt file
################################################################################################################
file_name = ''
if (len(sys.argv) == 2):
    file_name = sys.argv[1]
else:
    print 'File name expected'
folder_=os.getcwd()+'/parameters/'
#file_name=os.getcwd()+'/parameters/'+file_name
try:
    control = parse_com(file_name,folder_)    
except:
    print 'Done nothing'
#os.chdir(path_in)
################################################################################################################
#Creating one galaxy
################################################################################################################
gsample,par,file_masers,fname_directory,MMB,Arecibo,prt1,prt2,prt3,prt7,prt8,prt9,prt10=full_model(file_name,folder_)
path=os.getcwd()
################################################################################################################
#Taking values of the control txt file
################################################################################################################
sam_selec_control                              =float(control['sam_selec'])
sam_complete_control                           =float(control['sam_complete'])
sam_bri_dec_control                            =float(control['sam_bri_dec'])
many_control                                   =int(control['many'])
Bri_faint_control                              =int(control['Bri_faint'])
include_BeSSeL_control                         =int(control['BeSSeL_include'])
samp_rand_control                              =float(control['samp_rand'])
many_rand_control                              =int(control['many_rand'])
plots_control                                  =float(control['plots'])
Specific_Unc_control                           =float(control['Unc_confirmation_spec'])
Fudge_Factor_control                           =float(control['Fudge_Factor'])
################################################################################################################
#Creating Plots directories
################################################################################################################
if(plots_control==1):
    os.chdir(path+'/output/'+fname_directory)
    os.makedirs('output_plots')
################################################################################################################
#Creating Sample Selection
################################################################################################################
if(sam_selec_control==1):
    os.chdir(path+'/output/'+fname_directory)
################################################################################################################
#By bright of sources
################################################################################################################
    if(sam_bri_dec_control==1):
        Ns =create_array_control(control,'N_sour')        
        dds=create_array_control(control,'down_limit')        
        dus=create_array_control(control,'up_limit')        
        brisamples={}          
        if(Specific_Unc_control==2):
                unc_mu_xy=create_array_control(control,'Values_Unc')        
        if(Fudge_Factor_control==1):
                fudge_values=create_array_control(control,'Values_Fudge')        
        for i in range(many_control):
            directory_nsources='brig_dec_sample'+file_masers[9:-4]+'_'+str(i)
            os.chdir(path+'/output/'+fname_directory)        
            os.makedirs(directory_nsources) 
            print('Brightest sources with certain declination sample saved in folder: %s' % directory_nsources)
            globals()['variable{}'.format(i)] = 0
            #Creating the subsamples (bright or faint)
            sample_small,sample_small_file=Sample_Selection(gsample,Ns[i],dds[i],dus[i],file_masers,i,Bri_faint_control,include_BeSSeL_control)
            #Printing txt files
            if(Specific_Unc_control==2):
                if(Fudge_Factor_control==1):                            
                    Create_txt_Reid(sample_small_file,directory_nsources,Fudge_Factor_control,err_mu_xy[i],fudge_values[i])
                elif(Fudge_Factor_control==0):
                    Create_txt_Reid(sample_small_file,directory_nsources,Fudge_Factor_control,unc_mu_xy[i])
            elif(Specific_Unc_control==0 or Specific_Unc_control==1):
                if(Fudge_Factor_control==1):
                    Create_txt_Reid(sample_small_file,directory_nsources,Fudge_Factor_control,0,fudge_values[i])                                    
                elif(Fudge_Factor_control==0):
                    Create_txt_Reid(sample_small_file,directory_nsources,Fudge_Factor_control)
            if(plots_control==1):
                create_dir_plots(path,fname_directory,'bri_sample',i)
                do_plots(sample_small,par)
    os.chdir(path+'/output/'+fname_directory)
################################################################################################################    
#Random selection
################################################################################################################
    if(samp_rand_control==1):
        Nsrand=create_array_control(control,'N_sour_rand')
        for i in range(many_rand_control):
            directory_random='random_sample'+file_masers[9:-4]+'_'+str(i)
            os.chdir(path+'/output/'+fname_directory)
            os.makedirs(directory_random)
            print('Random sample saved in folder: %s' % directory_random) 
            #Creating sample
            random_sampl,random_sampl_file=Sample_Selection_random(gsample,Nsrand[i],file_masers,i)
            #Printing txt files
            Create_txt_Reid(file_masers,directory_complete,Fudge_Factor_control)
            #PLots
            if(plots_control==1):
                create_dir_plots(path,fname_directory,'rand_sample',i)
                do_plots(random_sampl,par)
    os.chdir(path+'/output/'+fname_directory)
################################################################################################################
#Complete sample
################################################################################################################
    if (sam_complete_control==1):
        #Still creating the whole sample without fudge value or unc assigment
        directory_complete='complete_sample'+file_masers[9:-4]
        os.chdir(path+'/output/'+fname_directory)
        os.makedirs(directory_complete)
        print('Complete sample saved in folder: %s' % directory_complete)
        Create_txt_Reid(file_masers,directory_complete,Fudge_Factor_control)
        if(plots_control==1):
            create_dir_plots(path,fname_directory,'complete_sample','')            
            do_plots(gsample,par)
        os.chdir(path+'/output/'+fname_directory)
else:
    print('No Sample Selection, making complete sample')
    directory_complete='complete_sample'+file_masers[9:-4]
    os.chdir(path+'/output/'+fname_directory)
    os.makedirs(directory_complete)
    print('Complete sample saved in folder: %s' % directory_complete)    
    Create_txt_Reid(file_masers,directory_complete,Fudge_Factor_control)
    if(plots_control==1):
        create_dir_plots(path,fname_directory,'complete_sample','')
        do_plots(gsample,par)       
    os.chdir(path+'/output/'+fname_directory)
################################################################################################################
#Creating the output txt file with the input and output parameteres
################################################################################################################
os.chdir(path)
shutil.copy2(path+'/parameters/'+file_name,path+'/output/'+fname_directory+'/outpara_'+file_masers[9:-4]+'.txt')
os.chdir(path+'/output/'+fname_directory+'/')
with open("outpara_"+file_masers[9:-4]+".txt", "a") as myfile:
    myfile.write('\n')
    myfile.write('-------------------------------------------------------------------------------')
    myfile.write('\n Output results for simulation number: \n')
    myfile.write(file_masers[9:-4])    
    myfile.write('\n ------------------------------------------------------------------------------- \n')
    myfile.write(prt1)
    myfile.write(prt2)    
    myfile.write(prt3)
    myfile.write(prt7)    
    myfile.write(prt8)
    myfile.write(prt9)    
    myfile.write(prt10)
    myfile.write('\n ------------------------------------------------------------------------------- \n')
    myfile.write('Comparison with Surveys: \n')    
    myfile.write('Sources that satisfy Arecibo limits: 76 (Observed) vs Simulated=\n')
    myfile.write(str(Arecibo[0]))   
    myfile.write('\n Flux density slope in Arecibo: -0.471 (Observed) vs Simulated=\n')
    myfile.write(str(round(Arecibo[1],3)))       
    myfile.write('\n Sources that satisfy MMB limits: 666 (Observed) vs Simulated=\n')
    myfile.write(str(MMB[0]))       
    myfile.write('\n Flux density slope in MMB: -0.533 (Observed) vs Simulated=\n')
    myfile.write(str(round(MMB[1],3)))  
    myfile.write('\n ------------------------------------------------------------------------------- \n')
os.chdir(path)
################################################################################################################
print ('Methanol masers simulation %s appears to be end successfully') %file_masers[10:-4]
print ('Please check the folder output/%s where all the outputs of this simulations were saved.' %fname_directory)