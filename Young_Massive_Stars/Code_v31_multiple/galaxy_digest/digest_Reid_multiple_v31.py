#/Users/luishenryquiroganunez/Documents/Leiden_University/3_Semester/MSc_thesis/Codes/
#Libraries
'''
Digests the prtfiles and creates a text comparison,
as well a csv file for plotting
'''
#------------------------------------------------------------------------------------              
import os
import csv
#Libraries
import numpy as np
import astropy
from astropy.io import ascii
from astropy.table import Table, Column
import scipy as scipy
from scipy.stats import kde
from scipy.stats.stats import pearsonr
import scipy.stats as stats
from scipy.stats import norm as norm
from numpy import trapz
from pylab import *
import re
import sys
import itertools
import aux.HTML as HTML
import webbrowser
from collections import OrderedDict
from scipy.stats import norm as norm
from matplotlib.patches import Rectangle
from scipy.spatial import distance

 #------------------------------------------------------------------------------------     
def removeComments(string):
    '''
    For remove /*COMMENT */ in control txt file
    '''
    string = re.sub(re.compile("/\*.*?\*/",re.DOTALL ) ,"" ,string)
    string = re.sub(re.compile("//.*?\n" ) ,"" ,string) 
    return string
    
def parse_com(name):
    # who line regexps
    # keyword = value, value
    repair = re.compile(r'^(\s*(\S+)\s*=\s*(\S+(?:\s*\,\s*\S+)*)\s*)+$')
    # hash in first 
    recomment = re.compile(r'^\#')
    #control={}
    control = OrderedDict()
    file = open(name,'r')
    for line in file:
        line = line.rstrip('\n')
        line = removeComments(line)  
        yescomment  = recomment.match(line)      
        yespair = repair.match(line)
        if yescomment:
            pass
            #print 'Reading input properly'
            #print 'Comment',line
        elif yespair:
            input = yespair.groups()
            control[input[1]]=input[2]
        else:
            None
            #print 'crap, no match:', line
    return control
    
    #------------------------------------------------------------------------------------              
def appendfit(dofile,fit,cov,No_data_files):
    curf = open(dofile)
    content = curf.readlines()
    results = []
    #print "Opening ",fitfile
    matching = [s for s in content if "Best parameter values for trial" in s]
    if(len(matching) != 0):
        #assuming the next line has best values
        kk=content.index(matching[0])
        results = content[kk+1].split()
        for i in range(len(fitpars)):
            #They must appear in fitpars order
            fit[fitpars[i]].append(float(results[i]))
        #Find matrix with cov estimates
        matching = [s for s in content if "Correlation Matrix" in s]    
        kk=content.index(matching[0])
        results = []
        for i in range(len(fitpars)):
            ipar = fitpars[i]
            results = content[kk+2+i].split()
            for j in range(len(fitpars)):
                jpar = fitpars[j]
                cov[ipar,jpar].append(float(results[j]))
        names_prt_files.append(fitfile)
    else:
        No_data_files.append(fitfile)
    curf.close()
    return
def sqrt_array_norm(arr):
    summ=[]
    for i in range(len(arr)):
        summ.append(arr[i]*arr[i])
    k=np.sqrt(sum(summ))
    return(k)       
#------------------------------------------------------------------------------------    
#Reading the files                         
path_origin=os.getcwd()+'/'
varq = raw_input("Where are the prt folders (path)?: ")
if(varq[-1]!='/'):
    varq=varq+'/'
prt_folders = os.listdir(varq) #including control.inp and parameters.txt
prt_folders_com=[]
for i in range(len(prt_folders)):
    if (prt_folders[i][0:4]=='prt_'):
        prt_folders_com.append(prt_folders[i])
        
####################################################################################
#Getting data for each parameter
####################################################################################
simu=[]
for kk in range(len(prt_folders_com)): 
    if (prt_folders_com[kk]=='prt_files__312_std_v17'):
        vass=kk
        print kk
    var=varq+prt_folders_com[kk]
    if(var[-1]!='/'):
        var=var+'/'
    dummy5=var.index('prt_files')+10
    name_ouputfolder=var[dummy5:-1]
    #------------------------------------------------------------------------------------
    os.chdir(var)
    prt_files_incomplete = os.listdir(var) #including control.inp and parameters.txt
    prtfiles=[]
    control_file=[]
    parameters_file=[]
    for i in range(len(prt_files_incomplete)):
        dummy=prt_files_incomplete[i]
        if(dummy[-4:]=='.inp'):
            control_file.append(dummy)
        elif(dummy[-4:]=='.txt'):
            parameters_file.append(dummy)
        elif(dummy[:7]=='output_'):
            pass
        elif(dummy=='.DS_Store'):
            pass
        else:
            prtfiles.append(dummy)
    os.chdir(path_origin)
    directory='output'
    if not os.path.exists(directory):
        os.makedirs(directory)
    os.chdir(directory)
    path_results=os.getcwd()
    if not os.path.exists('output_'+name_ouputfolder):
            os.makedirs('output_'+name_ouputfolder)
    os.chdir(path_results+'/output_'+name_ouputfolder)
    path_outputs=os.getcwd()
    os.chdir(var)
    #------------------------------------------------------------------------------------    
    #Declaring initial values
    initial_values = parse_com(parameters_file[0])
    dummy=[]
    f=open(control_file[0],'rb')
    for line in f:
        dummy.append(line)
    f.close()
    fitpars=['R0','Th0','dTd','U0','V0','W0','VS','US']
    #Initial Values given by the model
    modval={'R0':float(initial_values['r0']), 'Th0':float(initial_values['v0t']), 'dTd':float(initial_values['Rotation_Curve']), 'U0':float(initial_values['usun']), 'V0':float(initial_values['vsun']), 'W0':float(initial_values['wsun']), 'US':float(initial_values['v0r_us']), 'VS':(float(initial_values['v0t_vs'])*-1)}
    #sigval={'R0':, 'Th0':, 'dTd':float(initial_values['Rotation_Curve']), 'U0':float(initial_values['usun']), 'V0':float(initial_values['vsun']), 'W0':float(initial_values['wsun']), 'US':float(initial_values['v0r_us']), 'VS':(float(initial_values['v0t_vs'])*-1)}
    if (prt_folders_com[kk]=='prt_files__312_std_v17'):
        print modval    
    #Apriori Values        
    first_value={'R0':float(dummy[11][0:5]),   'Th0':float(dummy[12][0:5]), 'dTd':float(dummy[13][0:5]), 'U0':float(dummy[14][0:5]), 'V0':float(dummy[15][0:5]), 'W0':float(dummy[16][0:5]), 'VS':float(dummy[17][0:5]), 'US':float(dummy[18][0:5])}
    aprior_err={'R0':float(dummy[11][12:17]),   'Th0':float(dummy[12][12:17]), 'dTd':float(dummy[13][12:17]), 'U0':float(dummy[14][12:17]), 'V0':float(dummy[15][12:17]), 'W0':float(dummy[16][12:17]), 'VS':float(dummy[17][12:17]), 'US':float(dummy[18][12:17])}
    final_unc={'R0':float(dummy[11][21:27]),   'Th0':float(dummy[12][21:27]), 'dTd':float(dummy[13][21:27]), 'U0':float(dummy[14][21:27]), 'V0':float(dummy[15][21:27]), 'W0':float(dummy[16][21:27]), 'VS':float(dummy[17][21:27]), 'US':float(dummy[18][21:27])}       
    names_prt_files=[]
    #------------------------------------------------------------------------------------              
    #Make dict of lists to store simulation results
    #fit={}
    fit=OrderedDict()
    mean_vals=OrderedDict()
    error_para=OrderedDict()
    for par in fitpars:
        fit[par]=[]
        mean_vals[par]=[]   
        error_para[par]=[] 
    #Make a dict to store all the cov estimates, make it double-keyed
    #cov={}
    cov=OrderedDict()
    for ipar in fitpars:
        for jpar in fitpars:
            cov[ipar,jpar]=[]
    #List of files that do not fit
    No_data_files=[]
    #Extract values
    for fitfile in prtfiles:
        dofile=var+fitfile
        appendfit(dofile,fit,cov,No_data_files)
        
    for par in fitpars:
        mean_vals[par]=np.mean(fit[par])
    
    modval['VS']=modval['VS']*-1
    dummy1=[]  
    comp_fitpars=['R0', 'Th0','V0','VS']
    for par in fitpars:
        dummy1.append((mean_vals[par]-modval[par])/modval[par])
    os.chdir(path_origin)
        
    #simu.append(distance.euclidean(dummy1,dummy2))
    simu.append(sqrt_array_norm(dummy1))

########################################################################################################
#Plot
########################################################################################################            
os.chdir(path_origin)
matplotlib.rcParams.update({'font.size': 23})
fig1 = plt.figure(figsize=(16,8))
ax1 = fig1.add_subplot(111)
ax1.plot(simu,'o',color='b')
#ax1.plot(simu, linestyle='--',color='k',linewidth=2.5)
ax1.axhline(y=simu[vass],linewidth=3,color='k',ls='dotted',label='Reid 2014 initial parameters')
ax1.axhline(y=np.median(simu),linewidth=3,color='b',ls='dashed',label='Median value')
#ax1.set_ylim([8.12,8.58]) 
#ax1.set_yticks([8.2,8.3,8.4,8.5,8.6])
ax1.set_xlabel('Set of Simulation',labelpad=-4)
ax1.set_ylabel('$\\kappa$',labelpad=-4)
ax1.legend(ncol=2)
#ax1.set_yscale('log')
plt.savefig('errr_simu.pdf',bbox_inches='tight')
