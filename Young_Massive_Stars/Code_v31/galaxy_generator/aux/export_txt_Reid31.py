#Importing Libraries
import numpy
from numpy import genfromtxt
from math import *
import csv
import heapq
import random
import os
import shutil
from shutil import copyfile
import copy
from random import *
#------------------------------------------------------------------------------------------------------------
path=os.getcwd()
def decdeg2dms(dd):
   is_positive = dd >= 0
   dd = abs(dd)
   minutes,seconds = divmod(dd*3600,60)
   degrees,minutes = divmod(minutes,60)
   degrees = degrees if is_positive else -degrees
   if (degrees<10 and degrees>=0):
       degrees='{num:02d}'.format(num=int(degrees))
   elif(degrees<0 and degrees>-10):
       degrees='{num:03d}'.format(num=int(degrees))
   elif(degrees):
       degrees=str(int(degrees))
   if (minutes<10):
       minutes='{num:02d}'.format(num=int(minutes))
   elif(minutes):
       minutes=str(int(minutes))
   if (seconds<10):
       seconds='{num:02d}'.format(num=int(seconds))+'.'+str(int((round(seconds,2)-int(seconds)+0.0000001)*100))
   elif(seconds):
       seconds=str(round(seconds,2))
   return (degrees,minutes,seconds)
#------------------------------------------------------------------------------------------------------------
#Reading data and writing files
def Create_txt_Reid(filename,folder_number,Fudge_Factor_control,error_value=0,fudge_value=1.0):
    '''
    If error_value =0, then the error is calucluated. Otehrwise, the value given,
    will be settle for mu_x and mu_y
    '''
    with open(filename) as csvfile:
        #Opening the file
        reader = csv.DictReader(csvfile)
        row_count = sum(1 for row in reader)
    with open(filename) as csvfile:
        #Opening the file
        reader = csv.DictReader(csvfile)        
        #Declaring variables
        err_mud=numpy.zeros(row_count)
        para=numpy.zeros(row_count)        
        err_para=numpy.zeros(row_count)        
        vlsr=numpy.zeros(row_count)
        mu_x=numpy.zeros(row_count)
        err_mux=numpy.zeros(row_count)                
        mu_dc=numpy.zeros(row_count)      
        RAs=numpy.zeros(row_count)
        DECs=numpy.zeros(row_count)
        arm_no=numpy.zeros(row_count)
        err_vlsr=numpy.zeros(row_count)
        dummy_para=numpy.zeros(row_count)
        dummy_mux=numpy.zeros(row_count)
        dummy_muy=numpy.zeros(row_count)
        dummy_vlsr=numpy.zeros(row_count)
        i=0
        for row in reader:
            #print row['err_mud(mas/yr)']
            if(error_value != 0.0):
                err_mud[i]=error_value
                err_mux[i]=error_value
                err_para[i]=error_value
            elif(error_value == 0.0):
                err_mud[i]=row['err_mud(mas/yr)']
                err_mux[i]=row['err_mu_x(mas/yr)']            
                err_para[i]=row['err_parallax(mas)']
            err_vlsr[i]=row['err_vls(km/s)']  
            RAs[i]=row['ra(rad)']
            DECs[i]=row['dc(rad)']
            arm_no[i]=row['ic(Arm)']              
            if(Fudge_Factor_control==1):
                dummy_para[i]=row['parallax(mas)']           
                para[i]=gauss(dummy_para[i],err_para[i]*fudge_value)             
                dummy_mux[i]=row['mux(mas/yr)']
                mu_x[i]=gauss(dummy_mux[i],err_mux[i]*fudge_value)
                dummy_muy[i]=row['mud(mas/yr)']
                mu_dc[i]=gauss(dummy_muy[i],err_mud[i]*fudge_value) 
                dummy_vlsr[i]=row['vl(km/s)']            
                vlsr[i]=gauss(dummy_vlsr[i],err_vlsr[i]*fudge_value) 
            elif(Fudge_Factor_control==0):                            
                mu_dc[i]=row['mud_obs(mas/yr)']
                mu_x[i]=row['mux_obs(mas/yr)']
                para[i]=row['parallax_obs(mas)']
                vlsr[i]=row['vl_obs(km/s)']                
            i+=1
    RA=['' for x in range(row_count)]
    for i in range(row_count):
        y=decdeg2dms((RAs[i]*180/pi)*24/360)
        RA[i]=''.join(y)
    DEC=['' for x in range(row_count)]
    for i in range(row_count):
        y=decdeg2dms(DECs[i]*180/pi)
        DEC[i]=''.join(y)
    arm=['' for x in range(row_count)]
    for i in range(row_count):
        if(arm_no[i] == 0.0):
            arm[i]='Nor'
        if(arm_no[i] == 1.0):
            arm[i]='Sgr'
        if(arm_no[i] == 2.0):
            arm[i]='Per'
        if(arm_no[i] == 3.0):
            arm[i]='Sct'
        if(arm_no[i] == 4.0):
            arm[i]='Loc'            
#    for i in range(data.shape[0]):
    for i in range(row_count):            
        output=numpy.empty((9,5), object)
        output[0]=['! parallax and proper motion results for','','','','']
        output[1]=['! Quiroga-Nunez & van Langevelde Simulation','','','','']
        output[2]=[str(i),str(1.0),arm[i],'Quiroga-Nunez&Langevelde','! Name/#, Near/far code, Arm, Refs']
        output[3]=[RA[i],'','','','! RA  (hhmmss.s)']
        output[4]=[DEC[i],'','','','! Dec (ddmmss.s)']
        output[5]=[str('{:f}'.format(round(para[i],7))),str('{:f}'.format(round(err_para[i],7))),'','','! Parallax (mas), unc']
        output[6]=[str('{:f}'.format(round(mu_x[i],7))),str('{:f}'.format(round(err_mux[i],7))),'','','! X proper motion (mas/yr), unc']
        output[7]=[str('{:f}'.format(round(mu_dc[i],7))),str('{:f}'.format(round(err_mud[i],7))),'','','! Y proper motion (mas/yr), unc']
        output[8]=[str('{:f}'.format(round(vlsr[i],7))),str('{:f}'.format(round(err_vlsr[i],7))),'','','! VLSR (km/s), unc']
        numpy.savetxt(str(folder_number)+'/'+filename[10:-4]+'ppm_'+str(i+1)+'.txt', output, delimiter='\t', fmt='%s')
        os.rename(str(folder_number)+'/'+filename[10:-4]+'ppm_'+str(i+1)+'.txt',str(folder_number)+'/'+filename[10:-4]+'ppm_'+str(i+1)+'.dat')
    list_files_dat=os.listdir(folder_number)
    numpy.savetxt(str(folder_number)+'/source_files.inp',list_files_dat, delimiter='\t', fmt='%s')
    return
#------------------------------------------------------------------------------------------------------------
def Sample_Selection(sample,brigthest_number,dec_min,dec_max,csv_file,version,control_bri_faint,include_BeSSeL):
    #First mimicking BeSSeL results
    if(control_bri_faint==0):    
        if(include_BeSSeL==1):            
            BeSSeL_brigthest_number=100
            BeSSeL_dec_min=-30
            BeSSeL_dec_max=70
            BeSSeL_drad= [i['dc(rad)'] for i in sample]
            BeSSeL_drad_degrees=numpy.array([BeSSeL_drad[i]*180/pi for i in range(len(BeSSeL_drad))])
            BeSSeL_indexs=numpy.where((BeSSeL_drad_degrees>BeSSeL_dec_min) * (BeSSeL_drad_degrees<BeSSeL_dec_max))    
            BeSSeL_fbright_limited= numpy.array([i['f_MMB(Jy)'] for i in sample])    
            BeSSeL_brightest=heapq.nlargest(BeSSeL_brigthest_number,BeSSeL_fbright_limited[BeSSeL_indexs])    
            BeSSeL_indexes=numpy.where((BeSSeL_drad_degrees>BeSSeL_dec_min)*(BeSSeL_drad_degrees<BeSSeL_dec_max)*(BeSSeL_fbright_limited>=BeSSeL_brightest[-1]))    
            BeSSeL_small_sample=[sample[i] for i in BeSSeL_indexes[0]]    
            a=BeSSeL_indexes[0]
            aa=list(a)
            aa.reverse()
            '''
            queda=set(range(len(sample))).difference(a.tolist())
            z=list(queda)
            '''
            new_sample=copy.deepcopy(sample)#no Reid sample included
            for i in aa:
                del new_sample[i]    
            drad= [i['dc(rad)'] for i in new_sample]
            drad_degrees=numpy.array([drad[i]*180/pi for i in range(len(drad))])
            indexs=numpy.where((drad_degrees>dec_min) * (drad_degrees<dec_max))
            fbright_limited= numpy.array([i['f_MMB(Jy)'] for i in new_sample])
            brightest=heapq.nlargest(brigthest_number,fbright_limited[indexs])
            indexes=numpy.where((drad_degrees>dec_min)*(drad_degrees<dec_max)*(fbright_limited>=brightest[-1]))
            small_sample=[new_sample[i] for i in indexes[0]]
            final_sample=small_sample+BeSSeL_small_sample
        elif(include_BeSSeL==0):
            new_sample=copy.deepcopy(sample)#no Reid sample included
            drad= [i['dc(rad)'] for i in new_sample]
            drad_degrees=numpy.array([drad[i]*180/pi for i in range(len(drad))])
            indexs=numpy.where((drad_degrees>dec_min) * (drad_degrees<dec_max))
            fbright_limited= numpy.array([i['f_MMB(Jy)'] for i in new_sample])
            brightest=heapq.nlargest(brigthest_number,fbright_limited[indexs])
            indexes=numpy.where((drad_degrees>dec_min)*(drad_degrees<dec_max)*(fbright_limited>=brightest[-1]))
            small_sample=[new_sample[i] for i in indexes[0]]
            final_sample=small_sample                        
    elif(control_bri_faint==1):    
        new_sample=copy.deepcopy(sample)#no Reid sample included
        drad= [i['dc(rad)'] for i in new_sample]
        drad_degrees=numpy.array([drad[i]*180/pi for i in range(len(drad))])
        indexs=numpy.where((drad_degrees>dec_min) * (drad_degrees<dec_max))
        fbright_limited= numpy.array([i['f_MMB(Jy)'] for i in new_sample])
        faintest=heapq.nsmallest(brigthest_number,fbright_limited[indexs])
        indexes=numpy.where((drad_degrees>dec_min)*(drad_degrees<dec_max)*(fbright_limited<=faintest[-1]))
        small_sample=[new_sample[i] for i in indexes[0]]
        final_sample=small_sample
    #directory='output'
    path_csv='./'+csv_file[:-4]+'_brisample_'+str(version)+'.csv'
    fcsv = open(path_csv, 'wb')
    wdict = csv.DictWriter(fcsv, final_sample[0].keys() )
    wdict.writerow(dict(zip(final_sample[0].keys(),final_sample[0].keys())))
    wdict.writerows(final_sample)
    fcsv.close()
    return(final_sample,csv_file[:-4]+'_brisample_'+str(version)+'.csv') 
#------------------------------------------------------------------------------------------------------------    
def Sample_Selection_random(sample,number_rand_sour,csv_file,version):
    my_randoms = random.sample(xrange(len(sample)), number_rand_sour)
    random_sample=[sample[i] for i in my_randoms]
    #directory='output'
    path_csv='./'+csv_file[:-4]+'_ransample_'+str(version)+'.csv'
    fcsv = open(path_csv, 'wb')
    wdict = csv.DictWriter(fcsv, random_sample[0].keys() )
    wdict.writerow(dict(zip(random_sample[0].keys(),random_sample[0].keys())))
    wdict.writerows(random_sample)
    fcsv.close()
    return(random_sample,csv_file[:-4]+'_ransample_'+str(version)+'.csv') 
#------------------------------------------------------------------------------------------------------------    