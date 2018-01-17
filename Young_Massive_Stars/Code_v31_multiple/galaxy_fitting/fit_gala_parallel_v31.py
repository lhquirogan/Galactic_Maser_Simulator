##############################################################################################################################################
#Importing library
##############################################################################################################################################
import sys
import numpy
from numpy import f2py
import os
from multiprocessing import Process
import subprocess
import shutil
from joblib import Parallel, delayed  
import multiprocessing
path_ini=os.getcwd()
##############################################################################################################################################
#Getting the paths
##############################################################################################################################################
zaar = raw_input("Do you want to make multiple file parameters (yes or no)?:")
if (zaar=='yes' or zaar=='y'):
    vaar = raw_input("Where are folder with multiple galaxies (path): ")
    if(vaar[-1]!='/'):
        vaar=vaar+'/'                            
    list_multiple=os.listdir(vaar)
    try:
        list_multiple.remove('.DS_Store')
    except ValueError:
        pass
        
    var3 = raw_input("Name Folder that contains the control_file_bayesian.inp to be used (inside ./control_files/)?:")
    var3=path_ini+'/control_files/'+var3
    if(var3[-1]!='/'):
        var3=var3+'/'
    var4 = raw_input("Name Folder that contains the fit_galaxy_bayesian_16.f and a.out to be used (inside path ./fortran_codes/)?:")
    var4=path_ini+'/fortran_codes/'+var4
    if(var4[-1]!='/'):
        var4=var4+'/'
    ##############################################################################################################################################
    #Number of cores
    ##############################################################################################################################################
    num_cores_max = multiprocessing.cpu_count()
    ans = raw_input("How many cores do you want to use (max=%s): " % (num_cores_max)) 
    if (int(ans)>num_cores_max):
        print('Max cores reached!')
        sys.exit("Error message")
    else:    
        num_cores=int(ans)        
    for ku in range(len(list_multiple)):           
        var = vaar+list_multiple[ku]+'/'        
        dum=os.listdir(var4)
        for i in range(len(dum)):
            if(dum[i][-2:]=='.f'):
                file_fortran=dum[i]
        ##############################################################################################################################################        
        #Creating folders and csv files
        ##############################################################################################################################################
        os.chdir(var)
        list_galaxies=os.listdir(var)
        try:
            list_galaxies.remove('.DS_Store')
        except ValueError:
            pass
        os.chdir(path_ini)
        directory='output'
        if not os.path.exists(directory):
            os.makedirs(directory)
        os.chdir(directory)
        dummy2=var.index('/si')+6
        dummy23=var3.index('/control_files/')+15
        dummy24=var4.index('/fortran_codes/')+15
        os.makedirs('prt_files_'+var[dummy2:-1]+'_'+var3[dummy23:-1]+'_'+var4[dummy24:-1])
        os.makedirs('csv_files_'+var[dummy2:-1]+'_'+var3[dummy23:-1]+'_'+var4[dummy24:-1])
        os.chdir('csv_files_'+var[dummy2:-1]+'_'+var3[dummy23:-1]+'_'+var4[dummy24:-1])
        path_outputfiles_csv=os.getcwd()
        os.chdir('..')
        os.chdir('prt_files_'+var[dummy2:-1]+'_'+var3[dummy23:-1]+'_'+var4[dummy24:-1])
        path_outputfiles=os.getcwd()
        shutil.copy2(var3+'control_file_bayesian.inp', path_outputfiles+'/control_file_bayesian.inp')
        try:
            shutil.copy2(var+list_galaxies[0]+'/outpara_@'+list_galaxies[0]+'.txt', path_outputfiles+'/outpara_@'+list_galaxies[0]+'.txt')
        except IOError:
            try:
                shutil.copy2(var+list_galaxies[0]+'/outpara_'+list_galaxies[0]+'.txt', path_outputfiles+'/outpara_@'+list_galaxies[0]+'.txt')
            except:
                pass        
        os.chdir(var)
        ##############################################################################################################################################
        #Declaring paths and output names
        ##############################################################################################################################################
        inputs=[]
        for i in range(len(list_galaxies)):
            os.chdir(var+list_galaxies[i])
            list_innerfiles=os.listdir(os.getcwd())
            csv_tmp=[]
            for kk in range(len(list_innerfiles)):
                if(list_innerfiles[kk][-4:]=='.csv'):
                    csv_tmp.append(list_innerfiles[kk])
            for ll in range(len(csv_tmp)):
                shutil.copy2(var+list_galaxies[i]+'/'+csv_tmp[ll], path_outputfiles_csv+'/'+csv_tmp[ll])
            list_samples=os.walk('.').next()[1]
            for plots in list_samples:
                if(plots=='output_plots'):
                    aa=list_samples.index("output_plots")
                    del list_samples[aa]
            for j in range(len(list_samples)):
                os.chdir(var+list_galaxies[i]+'/'+list_samples[j])
                tmp_sample=os.getcwd()
                shutil.copy2(var3+'control_file_bayesian.inp', tmp_sample)
                shutil.copy2(var4+file_fortran, tmp_sample)
                shutil.copy2(var4+'a.out', tmp_sample)
                dummy=list_samples[j]
                path_gala_sam=os.getcwd()
                output_file='fit_'+dummy[16:]+'_'+dummy[0:3]+'.prt'
                inputs.append([path_gala_sam,output_file])
        ##############################################################################################################################################
        #Function to be parallelize
        ##############################################################################################################################################
        def procesinput(inp):
            print('Fitting Galaxy simulation number %s' % (inp[0]))
            os.chdir(inp[0])
            with open(inp[1], "w+") as output:
                subprocess.call(['./a.out'],stdout=output)
            output.closed
            shutil.copy2(inp[0]+'/'+inp[1], path_outputfiles)
            return
        '''   
        #To test         
        for i in inputs:
            procesinput(i)
        '''       
        print("numCores = " + str(num_cores)) 
        results = Parallel(n_jobs=num_cores)(delayed(procesinput)(i) for i in inputs)  
        os.chdir(path_ini)
        
else: 
    var = raw_input("Where are folder with galaxies (path): ")
    if(var[-1]!='/'):
        var=var+'/'
        
    var3 = raw_input("Name Folder that contains the control_file_bayesian.inp to be used (inside ./control_files/)?:")
    var3=path_ini+'/control_files/'+var3
    if(var3[-1]!='/'):
        var3=var3+'/'
    
    var4 = raw_input("Name Folder that contains the fit_galaxy_bayesian_16.f and a.out to be used (inside path ./fortran_codes/)?:")
    var4=path_ini+'/fortran_codes/'+var4
    if(var4[-1]!='/'):
        var4=var4+'/'
        
    dum=os.listdir(var4)
    for i in range(len(dum)):
        if(dum[i][-2:]=='.f'):
            file_fortran=dum[i]
    ##############################################################################################################################################        
    #Creating folders and csv files
    ##############################################################################################################################################
    os.chdir(var)
    list_galaxies=os.listdir(var)
    os.chdir(path_ini)
    directory='output'
    if not os.path.exists(directory):
        os.makedirs(directory)
    os.chdir(directory)
    dummy2=var.index('/si')+6
    dummy23=var3.index('/control_files/')+15
    dummy24=var4.index('/fortran_codes/')+15
    os.makedirs('prt_files_'+var[dummy2:-1]+'_'+var3[dummy23:-1]+'_'+var4[dummy24:-1])
    os.makedirs('csv_files_'+var[dummy2:-1]+'_'+var3[dummy23:-1]+'_'+var4[dummy24:-1])
    os.chdir('csv_files_'+var[dummy2:-1]+'_'+var3[dummy23:-1]+'_'+var4[dummy24:-1])
    path_outputfiles_csv=os.getcwd()
    os.chdir('..')
    os.chdir('prt_files_'+var[dummy2:-1]+'_'+var3[dummy23:-1]+'_'+var4[dummy24:-1])
    path_outputfiles=os.getcwd()
    shutil.copy2(var3+'control_file_bayesian.inp', path_outputfiles+'/control_file_bayesian.inp')
    try:
        shutil.copy2(var+list_galaxies[0]+'/outpara_@'+list_galaxies[0]+'.txt', path_outputfiles+'/outpara_@'+list_galaxies[0]+'.txt')
    except IOError:
        try:
            shutil.copy2(var+list_galaxies[0]+'/outpara_'+list_galaxies[0]+'.txt', path_outputfiles+'/outpara_@'+list_galaxies[0]+'.txt')
        except:
            pass        
    os.chdir(var)
    ##############################################################################################################################################
    #Number of cores
    ##############################################################################################################################################
    num_cores_max = multiprocessing.cpu_count()
    ans = raw_input("How many cores do you want to use (max=%s): " % (num_cores_max)) 
    if (int(ans)>num_cores_max):
        print('Max cores reached!')
        sys.exit("Error message")
    else:    
        num_cores=int(ans)
    ##############################################################################################################################################
    #Declaring paths and output names
    ##############################################################################################################################################
    inputs=[]
    for i in range(len(list_galaxies)):
        os.chdir(var+list_galaxies[i])
        list_innerfiles=os.listdir(os.getcwd())
        csv_tmp=[]
        for kk in range(len(list_innerfiles)):
            if(list_innerfiles[kk][-4:]=='.csv'):
                csv_tmp.append(list_innerfiles[kk])
        for ll in range(len(csv_tmp)):
            shutil.copy2(var+list_galaxies[i]+'/'+csv_tmp[ll], path_outputfiles_csv+'/'+csv_tmp[ll])
        list_samples=os.walk('.').next()[1]
        for plots in list_samples:
            if(plots=='output_plots'):
                aa=list_samples.index("output_plots")
                del list_samples[aa]
        for j in range(len(list_samples)):
            os.chdir(var+list_galaxies[i]+'/'+list_samples[j])
            tmp_sample=os.getcwd()
            shutil.copy2(var3+'control_file_bayesian.inp', tmp_sample)
            shutil.copy2(var4+file_fortran, tmp_sample)
            shutil.copy2(var4+'a.out', tmp_sample)
            dummy=list_samples[j]
            path_gala_sam=os.getcwd()
            output_file='fit_'+dummy[16:]+'_'+dummy[0:3]+'.prt'
            inputs.append([path_gala_sam,output_file])
    ##############################################################################################################################################
    #Function to be parallelize
    ##############################################################################################################################################
    def procesinput(inp):
        print('Fitting Galaxy simulation number %s' % (inp[0]))
        os.chdir(inp[0])
        with open(inp[1], "w+") as output:
            subprocess.call(['./a.out'],stdout=output)
        output.closed
        shutil.copy2(inp[0]+'/'+inp[1], path_outputfiles)
        return
    '''   
    #To test         
    for i in inputs:
        procesinput(i)
    '''       
    print("numCores = " + str(num_cores)) 
    results = Parallel(n_jobs=num_cores)(delayed(procesinput)(i) for i in inputs)  
    os.chdir(path_ini)    