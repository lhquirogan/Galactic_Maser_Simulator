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
#Reading the files                         
path_origin=os.getcwd()
var = raw_input("Where are the .prt files or fitting files (path)?: ")
if(var[-1]!='/'):
    var=var+'/'
dummy5=var.index('prt_files')+10
name_ouputfolder=var[dummy5:-1]
#------------------------------------------------------------------------------------
#Question Plots
colork2=['k','g','r','b','sienna','m','c','orange','chocolate','firebrick','rosybrown','crimson',
               'cadetblue','tomato']
plots_confirmation = raw_input("Do you want plots for parameters distributions and PDFs for coupled parameters (only 100 galaxies with 100 sources)? (yes or no):")
plots_confirmation_nsources = raw_input("Do you want plots for evolution of parameters, Pearson and errors? (yes or no):?")
latex_tables = raw_input("Do you want tables for latex (yes or no-->recommened):")
if plots_confirmation_nsources=='yes':
    b_w_n = raw_input("Northern sources (b) or whole galaxy sources (w)? (b or w or none):")
    errors_conf = raw_input("Compare errors between Northern and Whole? (yes or no--->recommened):")
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
#Apriori Values
first_value={'R0':float(dummy[11][0:5]),   'Th0':float(dummy[12][0:5]), 'dTd':float(dummy[13][0:5]), 'U0':float(dummy[14][0:5]), 'V0':float(dummy[15][0:5]), 'W0':float(dummy[16][0:5]), 'VS':float(dummy[17][0:5]), 'US':float(dummy[18][0:5])}
aprior_err={'R0':float(dummy[11][12:17]),   'Th0':float(dummy[12][12:17]), 'dTd':float(dummy[13][12:17]), 'U0':float(dummy[14][12:17]), 'V0':float(dummy[15][12:17]), 'W0':float(dummy[16][12:17]), 'VS':float(dummy[17][12:17]), 'US':float(dummy[18][12:17])}
final_unc={'R0':float(dummy[11][21:27]),   'Th0':float(dummy[12][21:27]), 'dTd':float(dummy[13][21:27]), 'U0':float(dummy[14][21:27]), 'V0':float(dummy[15][21:27]), 'W0':float(dummy[16][21:27]), 'VS':float(dummy[17][21:27]), 'US':float(dummy[18][21:27])}
names_prt_files=[]
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
#------------------------------------------------------------------------------------              
#Make dict of lists to store simulation results
#fit={}
fit=OrderedDict()
for par in fitpars:
    fit[par]=[]
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
#New valid names
#names_prt_files=list(set(names_prt_files) - set(No_data_files))
#------------------------------------------------------------------------------------              
#Getting data 
if(initial_values['sam_selec']=='1'):
    if(initial_values['sam_bri_dec']=='1'):
        #message_bri_n={}
        message_bri_n=OrderedDict()
        sam_bri=int(initial_values['many'])
        N=initial_values['N_sour']
        Nr= N.replace (",", " ")
        N_bri_sources = [int(x) for x in Nr.split()]          
        for x in range(sam_bri):
            message_bri_n['bri_n_'+str(x)]="""
            <p>     Specifications:            </p>
            <ul>
                <li> Number of sources:  %s per galaxy </li>
            </ul>
            """ % (N_bri_sources[x]) 
        #message_bri_d={}
        message_bri_d=OrderedDict()
        ddl=initial_values['down_limit']
        ddlr= ddl.replace (",", " ")
        Down_dec_limit = [int(x) for x in ddlr.split()]
        for x in range(sam_bri):
            message_bri_d['bri_d_'+str(x)]="""
            <ul>
                <li>Lower declination limit:  %s degrees </li>        
            </ul>
            """ % (Down_dec_limit[x]) 
        #message_bri_u={}            
        message_bri_u=OrderedDict()
        udl=initial_values['up_limit']
        udlr= udl.replace (",", " ")
        Up_dec_limit = [int(x) for x in udlr.split()] 
        for x in range(sam_bri):
            message_bri_u['bri_u_'+str(x)]="""
            <ul>
                <li> Upper declination limit:  %s  degrees </li>        
            </ul>
            """ % (Up_dec_limit[x])         
        if(float(initial_values['Unc_confirmation_spec'])==2):   
                message_unc=OrderedDict()
                Values_Unc=initial_values['Values_Unc']
                uncmuxy=Values_Unc.replace (",", " ")
                unc_mu_xy = [float(x) for x in uncmuxy.split()]
                for x in range(sam_bri):
                    message_unc['unc_'+str(x)]="""
                    <ul>
                        <li> Uncertinty for parallax and proper motions:  %s  mas [/yr] </li>        
                    </ul>
                    """ % (unc_mu_xy[x])                                                                 
        if(float(initial_values['Fudge_Factor'])==1):
                message_fudge=OrderedDict()
                Values_fudge=initial_values['Values_Fudge']
                fudge=Values_fudge.replace (",", " ")
                fudge_fact = [float(x) for x in fudge.split()] 
                for x in range(sam_bri):
                    message_fudge['fudge_'+str(x)]="""
                    <ul>
                        <li> Fudge Factor for parallax and proper motions distributions (0= perfect data):  %s   </li>        
                    </ul>
                    """ % (fudge_fact[x])                                                                 
        
    if(initial_values['samp_rand']=='1'):
        #message_ran_n={}
        message_ran_n=OrderedDict()
        sam_ran=int(initial_values['many_rand'])
        nran=initial_values['N_sour_rand']
        nranr= nran.replace (",", " ")
        N_ran_sources = [int(x) for x in nranr.split()]  
        for x in range(sam_ran):
            message_ran_n['ran_n_'+str(x)]="""
            <p> Specifications: </p>
            <ul>
                <li> Number of sources:  %s per galaxy </li>
                <li> No limitaton in location </li>        
            </ul>
            """ % (N_ran_sources[x]) 
    else:
        sam_ran=0

if (float(initial_values['BeSSeL_include'])==1.0):
    message_note="""     
        <p> <b> Important Note </b>: These specifications DO NOT the first sample which include the current stage 
            of BeSSeL porject: 100 sources in -30 < dec <70, which were included in ALL simulations.
        </p>                  
        """
else:
    message_note="""     
        <p> <b> Notes: N/A </b>.
        </p>                  
        """
    
if(initial_values['sam_selec']=='1'):
    #zz = {}
    zz=OrderedDict()
    if(initial_values['sam_bri_dec']=='1'):
        for x in range(sam_bri):
            if(float(initial_values['Unc_confirmation_spec'])==2):
                if(float(initial_values['Fudge_Factor'])==1):
                    zz[str(x)]=message_bri_n['bri_n_'+str(x)]+message_bri_d['bri_d_'+str(x)]+message_bri_u['bri_u_'+str(x)]+message_unc['unc_'+str(x)]+message_fudge['fudge_'+str(x)]+message_note
                elif(float(initial_values['Fudge_Factor'])==0):
                    zz[str(x)]=message_bri_n['bri_n_'+str(x)]+message_bri_d['bri_d_'+str(x)]+message_bri_u['bri_u_'+str(x)]+message_unc['unc_'+str(x)]+message_note
            elif(float(initial_values['Unc_confirmation_spec'])==0 or float(initial_values['Unc_confirmation_spec'])==1):
                if(float(initial_values['Fudge_Factor'])==1):
                    zz[str(x)]=message_bri_n['bri_n_'+str(x)]+message_bri_d['bri_d_'+str(x)]+message_bri_u['bri_u_'+str(x)]+message_fudge['fudge_'+str(x)]+message_note
                elif(float(initial_values['Fudge_Factor'])==0):                
                    zz[str(x)]=message_bri_n['bri_n_'+str(x)]+message_bri_d['bri_d_'+str(x)]+message_bri_u['bri_u_'+str(x)]+message_note            
    if(initial_values['samp_rand']=='1'):
        for x in range(sam_ran):
            zz[str(x+sam_bri)]=message_ran_n['ran_n_'+str(x)]   
    zz[str(sam_bri+sam_ran)]="""
    """
#------------------------------------------------------------------------------------              
#Sorting by samples
com_samples=[]
bri_samples=[]
ran_samples=[]
for l in range(len(names_prt_files)):
    dummys2=names_prt_files[l]
    if(dummys2[-7:]=='bri.prt'):
        bri_samples.append(dummys2)
    elif(dummys2[-7:]=='ran.prt'):
        ran_samples.append(dummys2)
    else:
        com_samples.append(dummys2)     
#Creating dictionaries
maxs_bri=[]
if(initial_values['sam_bri_dec']=='1'):
    for l in range(len(bri_samples)):
        dummys3=bri_samples[l]
        maxs_bri.append(int(dummys3[-9:-8]))     
    #Number of bright samples    
    n_bri=max(maxs_bri)+1    

maxs_ran=[]
if(initial_values['samp_rand']=='1'):
    for l in range(len(ran_samples)):
        dummys4=ran_samples[l]
        maxs_ran.append(int(dummys4[-9:-8]))     
    #Number of random samples    
    n_ran=max(maxs_ran)+1 
#------------------------------------------------------------------------------------
#Creating a dictionary of dictionaries
samples_names=[]
if(initial_values['sam_bri_dec']=='1'):
    for n in range(n_bri):
        samples_names.append('bri_'+str(n))
if(initial_values['samp_rand']=='1'):
    for n in range(n_ran):
        samples_names.append('ran_'+str(n))
samples_names.append('com')    
dic={name:{ keypar:[] for keypar in fitpars} for name in samples_names}
#------------------------------------------------------------------------------------
#Adding all the results to a dictionary of dictionaries
for i in range(len(names_prt_files)):
    name_file=names_prt_files[i]
    if(name_file[-7:]=='bri.prt'):
        for j in range(n_bri):
            if(name_file[-9:-8]==str(j)):
                for keypar in fitpars:   
                    if keypar=='VS':
                        dic['bri_'+str(j)][keypar].append(fit[keypar][i]*-1)                    
                    else:                 
                        dic['bri_'+str(j)][keypar].append(fit[keypar][i])
    elif(name_file[-7:]=='ran.prt'):
        for j in range(n_ran):
            if(name_file[-9:-8]==str(j)):
                for keypar in fitpars:                    
                    dic['ran_'+str(j)][keypar].append(fit[keypar][i])
    else:
        for keypar in fitpars:
            dic['com'][keypar].append(fit[keypar][i])
#------------------------------------------------------------------------------------  
#reid matruices
def append_cov(dofile,fit,cov,No_data_files,sample):
    curf = open(dofile)
    content = curf.readlines()
    #print "Opening ",fitfile
    matching = [s for s in content if "Best parameter values for trial" in s]
    if(len(matching) != 0):
        #Find matrix with cov estimates
        matching = [s for s in content if "Correlation Matrix" in s]    
        kk=content.index(matching[0])
        results = []
        for i in range(len(fitpars)):
            ipar = fitpars[i]
            results = content[kk+2+i].split()
            for j in range(len(fitpars)):
                jpar = fitpars[j]
                cov_cov[sam_ple][ipar,jpar].append(float(results[j]))
        names_prt_files.append(fitfile)
    else:
        No_data_files.append(fitfile)
    curf.close()
    return

cov_cov=OrderedDict()
No_data_files_cov=[]
fit_cov=OrderedDict()
for par in fitpars:
    fit_cov[par]=[]
if(initial_values['sam_bri_dec']=='1'):
    for k in range(n_bri):
        cov_cov['bri_'+str(k)]=OrderedDict()    
if(initial_values['samp_rand']=='1'):    
    for k in range(n_ran):
        cov_cov['ran_'+str(k)]=OrderedDict()   
cov_cov['com']=OrderedDict()           
ssamples=cov_cov.keys()
for sample in ssamples:
    for ipar in fitpars:
        for jpar in fitpars:
            cov_cov[sample][ipar,jpar]=[]
        
for fitfile in prtfiles:
    #if (fitfile=='.DS_Store'):
    dofile=var+fitfile
    if(initial_values['sam_bri_dec']=='1'):
        for ii in range(n_bri):
            if(fitfile[-9:-4]==str(ii)+'_bri'):
                sam_ple='bri_'+str(ii)
    if(initial_values['samp_rand']=='1'):
        for ii in range(n_bri):
            if(fitfile[-9:-4]==str(ii)+'_ran'):
                sam_ple='ran_'+str(ii)    
    if(fitfile[-7:-4]=='com'):
        sam_ple='com'     
    append_cov(dofile,fit_cov,cov_cov,No_data_files_cov,sam_ple)
#------------------------------------------------------------------------------------                  
#Path and file
path=os.getcwd()
html_filename='Table_Results.html'
#------------------------------------------------------------------------------------
#starting html file
message_ini="""<html>
<body>
    <h1>
        <center>Simulations developed for the BeSSeL VLBA Key Science Project </center>        
    </h1>
    <a href="http://bessel.vlbi-astrometry.org/">   
        <center>
            <img src="https://science.nrao.edu/science/key-science-projects/images/bessel.jpg" alt="bessel.vlbi-astrometry.org">
        </center>
    </a>
    
    <p> These codes were developed by:</p>
    <ul>
        <li> <a href="https://www.cfa.harvard.edu/~reid/" target="_blank">Dr. Mark Reid (CfA)</a> </li>
        <li> <a href="http://jive.eu/~huib/" target="_blank"> Prof. Dr. Huib van Langevelde (JIVE/Leiden) </a> </li>
        <li> <a href="https://www.strw.leidenuniv.nl/~quiroganunez" target="_blank"> Luis Henry Quiroga-Nunez (Leiden/JIVE) </a> </li>
    </ul>
    <h1>
        <center> Results for simulations called : %s </center>
    </h1>
""" % (name_ouputfolder)
message_ini2="""
    <p>
        Prt files used for these results are saved in: %s
    </p>
""" % (var)
'''
if(plots_function_confirmation=='yes'):
    message_plots="""
    <p>  </p>
    <p>  </p>
    <p>  PLOTS:</p>

    <center> <img src="Pearson_%s.png" alt="Pearson N sources" style="width:1000px;height:500px;"> </center> <br>
    <center> <img src="pline_%s.png" alt="P value Soruces" style="width:1000px;height:500px;"> </center> <br>
    """ %(clue,clue)
else: 
    message_plots="""
        <p> NO PLOTS MADE! </p>
    """
'''
message_final="""
    <p> These codes were developed by: </p>
    <ul>
        <li> <a href="https://www.cfa.harvard.edu/~reid/" target="_blank">Dr. Mark Reid (CfA)</a> </li>
        <li> <a href="http://jive.eu/~huib/" target="_blank"> Prof. Dr. Huib van Langevelde (JIVE/Leiden) </a> </li>
        <li> <a href="https://www.strw.leidenuniv.nl/~quiroganunez" target="_blank"> Luis Henry Quiroga-Nunez (Leiden/JIVE) </a> </li>
    </ul>
    </body>
</html>"""
#------------------------------------------------------------------------------------
#Creating Table of results    
names_samples=[]
head_row=['Parameter','Result','No. Galaxies','Model Value','A priori','Expect Unc.']
if(initial_values['sam_bri_dec']=='1'):
    for z in range(n_bri):
        names_samples.append('Sample limited by Brightness & Location No. '+str(z))
if(initial_values['samp_rand']=='1'):   
    for z in range(n_ran):
        names_samples.append('Random Sample No. '+str(z))
names_samples.append('Complete Sample (no limitations)')
#dic_table={}
dic_table=OrderedDict()
dic_table_names=[]
if(initial_values['sam_bri_dec']=='1'):
    for x in range(n_bri):
            dic_table['t_bri_'+str(x)]=HTML.Table(header_row=head_row)
            dic_table_names.append('t_bri_'+str(x))
if(initial_values['samp_rand']=='1'):        
    for x in range(n_ran):
            dic_table['t_ran_'+str(x)]=HTML.Table(header_row=head_row)        
            dic_table_names.append('t_ran_'+str(x))
            
dic_table['t_com']=HTML.Table(header_row=head_row)
dic_table_names.append('t_com')
#------------------------------------------------------------------------------------
#Saving message results
message_results=[]
for k in range(len(dic)):        
    message_results.append('''
    <h1>Final results for: %s </h1>
    ''' % (names_samples[k]))
#------------------------------------------------------------------------------------
#Saving results in Table
for i in range(len(dic_table_names)):
    for keypar in fitpars:
        if keypar != 'void':
            if(aprior_err[keypar]<0.0):
                dic_table[dic_table_names[i]].rows.append(['{:>4s}:'.format(keypar),'{:6.2f} +/- {:5.2f}'.format(np.mean(dic[dic_table_names[i][2:]][keypar]),np.std(dic[dic_table_names[i][2:]][keypar])),'{:4}'.format(len(dic[dic_table_names[i][2:]][keypar])),'{:6.2f}'.format(modval[keypar]),'{:6.2f} +/- {:5.3s}'.format(first_value[keypar],'NaN'),'{:6.2f}'.format(final_unc[keypar])])
            else:
                dic_table[dic_table_names[i]].rows.append(['{:>4s}:'.format(keypar),'{:6.2f} +/- {:5.2f}'.format(np.mean(dic[dic_table_names[i][2:]][keypar]),np.std(dic[dic_table_names[i][2:]][keypar])),'{:4}'.format(len(dic[dic_table_names[i][2:]][keypar])),'{:6.2f}'.format(modval[keypar]),'{:6.2f} +/- {:5.2f}'.format(first_value[keypar],aprior_err[keypar]),'{:6.2f}'.format(final_unc[keypar])])
   
#------------------------------------------------------------------------------------
#Saving results of the tables
results_html_code=[]
for i in range(len(dic_table_names)):
    results_html_code.append(str(dic_table[dic_table_names[i]]))
#results_html_code=[str(t_100),str(t_200),str(t_complete),str(t_fit)]
#------------------------------------------------------------------------------------   
#Reid values to compare with
Reid_pearson=[[ 1.000, 0.465, 0.103, 0.452, 0.023,-0.003, -0.002, 0.517],
              [ 0.465, 1.000, 0.136, 0.243,-0.796,-0.009, -0.809, 0.171],
              [ 0.103, 0.136, 1.000,-0.124,-0.009, 0.025, -0.018,-0.094],
              [ 0.452, 0.243,-0.124, 1.000,-0.014,-0.017,0.025, 0.839],
              [ 0.023,-0.796,-0.009,-0.014, 1.000, 0.011,0.990,-0.006],
              [-0.003,-0.009, 0.025,-0.017, 0.011, 1.000,0.010,-0.002],
              [ -0.002, -0.809, -0.018,0.025,0.990,0.010, 1.000,0.028],
              [ 0.517, 0.171,-0.094, 0.839,-0.006,-0.002,0.028, 1.000]]
#------------------------------------------------------------------------------------ 
#Saving message results
message_pearson=[]
for k in range(len(dic)):        
    message_pearson.append('''
    <h1>Pearson coefficients for: %s</h1>
    ''' % (names_samples[k]))
#--------------------------------------------------------------------------------
#Creating Table of Pearson   
pearson_head_row=[]
pearson_head_row.append('Pearson')
for keypar in fitpars:
    pearson_head_row.append(keypar)

#dic_table_p={}
dic_table_p=OrderedDict()
dic_table_names_p=[]
if(initial_values['sam_bri_dec']=='1'):
    for x in range(n_bri):
            dic_table_p['p_bri_'+str(x)]=HTML.Table(header_row=pearson_head_row)
            dic_table_names_p.append('p_bri_'+str(x))
if(initial_values['samp_rand']=='1'): 
    for x in range(n_ran):
            dic_table_p['p_ran_'+str(x)]=HTML.Table(header_row=pearson_head_row)        
            dic_table_names_p.append('p_ran_'+str(x))
        
dic_table_p['p_com']=HTML.Table(header_row=pearson_head_row)
dic_table_names_p.append('p_com')
#--------------------------------------------------------------------------------
#Saving Table of Pearson  
#final_coeefs={}
final_coeefs=OrderedDict()
#final_coeefs2={}
final_coeefs2=OrderedDict()
for kk in range(len(dic_table_p)):
    matx= [[0 for x in range(len(fitpars))] for x in range(len(fitpars))] 
    mat2x = [[0 for x in range(len(fitpars))] for x in range(len(fitpars))] 
    for i in range(len(fitpars)):
        for j in range(len(fitpars)):
            ipar=fitpars[i]
            jpar = fitpars[j]        
            (matx[i][j],mat2x[i][j])=pearsonr(dic[dic_table_names_p[kk][2:]][ipar],dic[dic_table_names_p[kk][2:]][jpar])

    for i in range(len(fitpars)):    
        dummy_row=[]   
        dummy_row.append(fitpars[i]) 
        for j in range(len(fitpars)):
            if (abs(matx[i][j]-Reid_pearson[i][j]) > 0.5):
                colorname='RED'
            elif(0.3<abs(matx[i][j]-Reid_pearson[i][j])<0.5):
                colorname='YELLOW'   
            else:           
                colorname='transparent'
            dummy_color= HTML.TableCell('{:6.2f}'.format(matx[i][j])+'({:4.2f})'.format(Reid_pearson[i][j]), bgcolor=colorname)
            dummy_row.append(dummy_color)
        dic_table_p[dic_table_names_p[kk]].rows.append(dummy_row)
    final_coeefs[dic_table_names_p[kk]]=matx
    final_coeefs2[dic_table_names_p[kk]]=mat2x
#------------------------------------------------------------------------------------
#Saving results of the tables
results_pearson_html=[]
for i in range(len(dic_table_names_p)):
    results_pearson_html.append(str(dic_table_p[dic_table_names_p[i]]))

#results_pearson_html=[str(p_100),str(p_200),str(p_complete)]
#------------------------------------------------------------------------------------
#Saving html value
fhu=open(path_outputs+'/'+html_filename,'w') 
fhu.write(message_ini)
fhu.write(message_ini2)
for i in range(len(dic)):
    fhu.write(message_results[i])
    fhu.write(zz[str(i)])
    fhu.write(results_html_code[i])
    fhu.write(message_pearson[i])
    fhu.write(results_pearson_html[i])
#fhu.write(message_plots)
fhu.write(message_final)
fhu.close()
'''
import webbrowser
file_name=path+'/'+html_filename
webbrowser.open_new_tab(file_name)
'''

#-----------------------------------------------------------------------------------
num_sources_per_galaxy=OrderedDict()        
for i in range(len(N_bri_sources)):
    #if(==1):
    #add=100
    #elif(==0):
    #add=0    
    num_sources_per_galaxy['p_bri_'+str(i)]=N_bri_sources[i]+100
    N_bri_sources[i]=N_bri_sources[i]+100
dummy_7=num_sources_per_galaxy.keys()


###################################     Plots      ##########################################
##########################################################################################################################
###################################     Galac. Paramters Histograms BeSSeL      ##########################################
##########################################################################################################################
if (plots_confirmation=='yes'):
    # Eight Plots
    fig1 = plt.figure(figsize=(15,10))
    matplotlib.rcParams.update({'font.size': 16})
    zz=0 #which sample?
    label_height=35
    ax = fig1.add_subplot(111)
          
    ax1 = fig1.add_subplot(241)
    mean_R0,sigmas_R0_Reid = modval['R0'],final_unc['R0']  # initial values
    (mu_R0, sigmaR0) = norm.fit(dic[dic_table_names[zz][2:]]['R0'])
    normed_R0 = [x for x in dic[dic_table_names[zz][2:]]['R0']]    
    (mu_normed_R0, sigma_normed_R0) = norm.fit(normed_R0)
    a,bins_data,c=ax1.hist(normed_R0,bins=6,color=colork2[1], alpha=0.7) #     
    area = sum(np.diff(bins_data)*a)
    xmin,xmax=ax1.get_xlim()
    #xmin,xmax=ax5.get_xlim()
    x = np.linspace(xmin, 2*xmax, 1000)
    y=mlab.normpdf(x, mu_normed_R0, sigma_normed_R0)*area
    y1=mlab.normpdf(x, mean_R0, sigmas_R0_Reid)*area
    ax1.plot(x,y1, color='grey',linewidth=2.5)    
    ax1.plot(x,y, linestyle='--',color='k',linewidth=2.5)
    ax1.fill_between(x, 0, y1, facecolor='grey', alpha=0.4)
    ax1.set_yticks([0,10,20,30,40])
    #ax1.set_xticks([0.96,0.98,1.0,1.02,1.04])
    #ax1.set_xticks([0.7,0.85,1.0,1.15,1.3])
    ax1.set_xlim([8.12,8.58]) 
    ax1.set_xticks([8.2,8.3,8.4,8.5,8.6])
    ax1.grid(True)
    #ax1.set_xlabel('$R_0$ ($\\rm{kpc}$)',labelpad=-4)
    #plt.setp(ax1.get_xticklabels(), visible=False) 
    ax1.text( mu_R0,label_height , '$R_0$ $\\rm{(kpc)}$', verticalalignment='center', horizontalalignment='center', color=colork2[1], fontsize=25)
    #ax1.axvline(x=mean_R0,linewidth=3,color='k',ls='dotted')
    
    ax2 = fig1.add_subplot(242,sharey=ax1)
    mean_Th0,sigmas_Th0_Reid = modval['Th0'], final_unc['Th0']# initial values    
    (mu_Th0, sigmaTh0) = norm.fit(dic[dic_table_names[zz][2:]]['Th0'])
    normed_Th0 = [x for x in dic[dic_table_names[zz][2:]]['Th0']]
    (mu_normed_Th0, sigma_normed_Th0) = norm.fit(normed_Th0)
    a,bins_data,c=ax2.hist(normed_Th0,bins=7,color=colork2[2], alpha=0.7) # 
    area = sum(np.diff(bins_data)*a) 
    xmin,xmax=ax2.get_xlim()
    #xmin,xmax=ax6.get_xlim()   
    x = np.linspace(xmin/2, 2*xmax, 1000)
    y=mlab.normpdf(x, mu_normed_Th0, sigma_normed_Th0)*area
    y1=mlab.normpdf(x, mean_Th0, sigmas_Th0_Reid)*area
    ax2.plot(x,y1, color='grey',linewidth=2.5) 
    ax2.plot(x,y, linestyle='--',color='k',linewidth=2.0)     
    ax2.fill_between(x, 0, y1, facecolor='grey', alpha=0.4)    
    #ax2.set_xticks([0.6,0.8,1.0,1.2,1.4])
    ax2.set_xlim([228,252]) 
    ax2.set_xticks([230,235,240,245,250])
    ax2.grid(True)
    #ax2.set_xlabel('$\\Theta_0$ ($\\rm{km \, s^{-1}}$)',labelpad=-4)
    plt.setp(ax2.get_yticklabels(), visible=False)
    #plt.setp(ax2.get_xticklabels(), visible=False)
    ax2.text( mu_Th0,  label_height, '$\\Theta_0$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[2], fontsize=25)         
    #ax2.axvline(x=mean_Th0,linewidth=3,color='k',ls='dotted')
    
    ax3 = fig1.add_subplot(243,sharey=ax1)
    mean_dTd,sigmas_dTd_Reid = modval['dTd'],final_unc['dTd'] # initial values
    (mu_mean_dTd, sigma_mean_dTd) = norm.fit(dic[dic_table_names[zz][2:]]['dTd'])
    normed_dTd = [x for x in dic[dic_table_names[zz][2:]]['dTd']]    
    (mu_normed_dTd, sigma_normed_dTd) = norm.fit(normed_dTd)
    a,bins_data,c=ax3.hist(normed_dTd,bins=7,color=colork2[8], alpha=0.7) #     
    area = sum(np.diff(bins_data)*a)
    xmin,xmax=ax3.get_xlim()    
    x = np.linspace(-2.7, xmax, 100)
    y1=mlab.normpdf(x, mean_dTd, sigmas_dTd_Reid)*area
    ax3.plot(x,y1, color='grey',linewidth=2.5)     
    ax3.fill_between(x, 0, y1, facecolor='grey', alpha=0.4)
    y=mlab.normpdf(x, mu_normed_dTd, sigma_normed_dTd)*area
    ax3.plot(x,y, linestyle='--',color='k',linewidth=2.0)
    #ax3.set_xlim([-12,15]) 
    ax3.set_xlim([-2.7,xmax])     
    ax3.set_xticks([-2,-1,0,1,2])  
    ax3.grid(True)
    plt.setp(ax3.get_yticklabels(), visible=False)     
    ax3.text(mu_mean_dTd, label_height, '$\\frac{d \\Theta}{dR}$ $\\rm{(km/s \\cdot kpc)}$', verticalalignment='center', horizontalalignment='center', color=colork2[8], fontsize=25)
    #ax3.axvline(x=mean_dTd,linewidth=3,color='k',ls='dotted')
          
    ax4 = fig1.add_subplot(244,sharey=ax1)
    mean_U0,sigmas_U0_Reid = modval['U0'],final_unc['U0'] # initial values    
    (mu_U0, sigmaU0) = norm.fit(dic[dic_table_names[zz][2:]]['U0'])
    normed_U0 = [x for x in dic[dic_table_names[zz][2:]]['U0']]
    (mu_normed_U0, sigma_normed_U0) = norm.fit(normed_U0)
    a,bins_data,c=ax4.hist(normed_U0,bins=9,color=colork2[5], alpha=0.7) # 
    area = sum(np.diff(bins_data)*a)
    xmin,xmax=ax4.get_xlim()
    x = np.linspace(xmin/2, xmax*2, 1000)
    y=mlab.normpdf(x, mu_normed_U0, sigma_normed_U0)*area
    y1=mlab.normpdf(x, mean_U0, sigmas_U0_Reid)*area
    ax4.plot(x,y1, color='grey',linewidth=2.5)  
    ax4.plot(x,y, linestyle='--',color='k',linewidth=2.0)    
    ax4.fill_between(x, 0, y1, facecolor='grey', alpha=0.4)    
    ax4.set_yticks([0,10,20,30,40])
    ax4.set_ylim([0,39])        
    ax4.set_xlim([6.5,14.6])     
    ax4.set_xticks([7,9,11,13])
    #ax4.set_xlim([0.75,1.35])    
    #ax4.set_xlabel('$U_{\\odot}$ ($\\rm{km \, s^{-1}}$)',labelpad=-4)     
    ax4.grid(True)
    plt.setp(ax4.get_yticklabels(), visible=False) 
    ax4.text( mu_U0,  label_height, '$U_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[5], fontsize=25)         
    #ax4.axvline(x=mean_U0,linewidth=3,color='k',ls='dotted')
         
    ax5 = fig1.add_subplot(245)    
    mean_V0,sigmas_V0_Reid = modval['V0'],final_unc['V0'] # initial values    
    (mu_V0, sigmaV0) = norm.fit(dic[dic_table_names[zz][2:]]['V0'])
    normed_V0 = [x for x in dic[dic_table_names[zz][2:]]['V0']]
    (mu_normed_V0, sigma_normed_V0) = norm.fit(normed_V0)
    a,bins_data,c=ax5.hist(normed_V0,bins=7,color=colork2[6], alpha=0.7) #     
    area = sum(np.diff(bins_data)*a)
    xmin,xmax=ax5.get_xlim()    
    x = np.linspace(xmin/2, xmax*2, 1000)
    y=mlab.normpdf(x, mu_normed_V0, sigma_normed_V0)*area
    y1=mlab.normpdf(x, mean_V0, sigmas_V0_Reid)*area
    ax5.plot(x,y1, color='grey',linewidth=2.5)     
    ax5.plot(x,y, linestyle='--',color='k',linewidth=2.0)    
    ax5.fill_between(x, 0, y1, facecolor='grey', alpha=0.4)    
    ax5.set_xlim([8,22]) 
    ax5.set_xticks([9,12,15,18,21])  
    ax5.set_yticks([0,10,20,30,40])    
    ax5.grid(True)
    ax5.text( mu_V0, label_height, '$V_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[6], fontsize=25)        
    #ax5.axvline(x=mean_V0,linewidth=3,color='k',ls='dotted')

    ax6 = fig1.add_subplot(246,sharey=ax5)    
    mean_W0,sigmas_W0_Reid = modval['W0'],final_unc['W0'] # initial values
    (mu_W0, sigmaW0) = norm.fit(dic[dic_table_names[zz][2:]]['W0'])
    normed_W0 = [x for x in dic[dic_table_names[zz][2:]]['W0']]
    (mu_normed_W0, sigma_normed_W0) = norm.fit(normed_W0)
    a,bins_data,c=ax6.hist(normed_W0,bins=6,color=colork2[7], alpha=0.7) # 
    area = sum(np.diff(bins_data)*a)
    xmin,xmax=ax6.get_xlim()
    x = np.linspace(xmin, xmax, 100)
    y=mlab.normpdf(x, mu_normed_W0, sigma_normed_W0)*area
    y1=mlab.normpdf(x, mean_W0, sigmas_W0_Reid)*area
    ax6.plot(x,y1, color='grey',linewidth=2.5,label='Results Model A5')     
    ax6.plot(x,y, linestyle='--',color='k',linewidth=2.0)    
    ax6.fill_between(x, 0, y1, facecolor='grey', alpha=0.4)    
    ax6.set_xlim([xmin,xmax])     
    ax6.set_xticks([8,8.5,9,9.5,10])
    #ax6.set_xlim([-0.3,2.8])  
    ax6.grid(True)
    plt.setp(ax6.get_yticklabels(), visible=False) 
    ax6.text( mu_W0, label_height, '$W_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[7], fontsize=25)         
    #ax6.axvline(x=mean_W0,linewidth=3,color='k',ls='dotted')
    ax6.legend(bbox_to_anchor=(0.7, -0.05))        

    ax7 = fig1.add_subplot(247,sharey=ax5)    
    mean_US,sigmas_US_Reid = modval['US'],final_unc['US']# initial values    
    (mu_US, sigmaUS) = norm.fit(dic[dic_table_names[zz][2:]]['US'])
    normed_US = [x for x in dic[dic_table_names[zz][2:]]['US']]
    (mu_normed_US, sigma_normed_US) = norm.fit(normed_US)
    a,bins_data,c=ax7.hist(normed_US,bins=7,color=colork2[3], alpha=0.7) # 
    area = sum(np.diff(bins_data)*a)
    xmin,xmax=ax7.get_xlim()
    #xmin,xmax=ax7.get_xlim()
    x = np.linspace(xmin*2, xmax*2, 1000)
    y=mlab.normpdf(x, mu_normed_US, sigma_normed_US)*area
    y1=mlab.normpdf(x, mean_US, sigmas_US_Reid)*area
    ax7.plot(x,y1, color='grey',linewidth=2.5) 
    ax7.plot(x,y, linestyle='--',color='k',linewidth=2.0,label='Gaussian Fitting')     
    ax7.fill_between(x, 0, y1, facecolor='grey', alpha=0.4)    
    #ax7.set_xticks([-0.5,0.0,0.5,1.0,1.5,2.0,2.5])
    #ax7.set_xticks([-1.0,0.0,1.0,2.0,3.0])
    ax7.set_xlim([-1.2,7.9])     
    ax7.set_xticks([0,2,4,6])
    ax7.grid(True)
    #ax7.set_xlabel('$U_s$ ($\\rm{km \, s^{-1}}$)',labelpad=-4)
    plt.setp(ax7.get_yticklabels(), visible=False)    
    #plt.setp(ax7.get_xticklabels(), visible=False)    
    ax7.text( mu_US,  label_height, '$\\bar{U}_s$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[3], fontsize=25)               
    #ax7.axvline(x=mean_US,linewidth=3,color='k',ls='dotted') 
    ax7.legend(bbox_to_anchor=(1.4, -0.05))                
    
    ax8 = fig1.add_subplot(248,sharey=ax5)
    mean_VS,sigmas_VS_Reid = modval['VS'],final_unc['VS'] # initial values          
    (mu_VS, sigmaVS) = norm.fit(dic[dic_table_names[zz][2:]]['VS'])
    normed_VS = [x  for x in dic[dic_table_names[zz][2:]]['VS']]
    (mu_normed_VS, sigma_normed_VS) = norm.fit(normed_VS)
    a,bins_data,c=ax8.hist(normed_VS,bins=6,color=colork2[4], alpha=0.7) # 
    area = sum(np.diff(bins_data)*a)
    xmin,xmax=ax8.get_xlim()
    #xmin,xmax=ax8.get_xlim()
    x = np.linspace(xmin, xmax, 100)
    y=mlab.normpdf(x, mu_normed_VS, sigma_normed_VS)*area
    y1=mlab.normpdf(x, mean_VS, sigmas_VS_Reid)*area
    ax8.plot(x,y1, color='grey',linewidth=2.5)     
    ax8.plot(x,y, linestyle='--',color='k',linewidth=2.0)    
    ax8.fill_between(x, 0, y1, facecolor='grey', alpha=0.4)    
    #ax8.set_xticks([-14,-9,-4,1,6,11,16])
    ax8.set_ylim([0,39])  
    ax8.set_yticks([0,10,20,30,40])    
    ax8.set_xlim([-12,8])     
    ax8.set_xticks([-9,-6,-3,0,3,6])
    ax8.grid(True)
    #ax8.set_xlabel('$V_s$ ($\\rm{km \, s^{-1}}$)',labelpad=-4)    
    plt.setp(ax8.get_yticklabels(), visible=False)
    #plt.setp(ax8.get_xticklabels(), visible=False)
    ax8.text( mean_VS,  label_height, '$\\bar{V}_s$ $\\rm{(km/s)}$ ', verticalalignment='center', horizontalalignment='center', color=colork2[4], fontsize=25)            
    #ax8.axvline(x=mean_VS,linewidth=3,color='k',ls='dotted')
                                                            
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    #plt.setp(ax.get_yticklabels(), visible=False)  
    #plt.setp(ax.get_xticklabels(), visible=False)    
    #plt.axes(frameon=False)      
    #ax.axis('off')
    #fig1.axes.get_xaxis().set_visible(False)
    #fig1.axes.get_yaxis().set_visible(False)
    ax.set_ylabel('Number Counts',labelpad=32)
    #ax.set_xlabel('Normed Parameter',labelpad=22)
    ax.tick_params(
    axis='both',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off',right='off', left='off', labelleft='off')
    #fig1.patch.set_visible(False)
    fig1.tight_layout()
        
    plt.savefig(path_outputs+'/para_'+name_ouputfolder+'.pdf',bbox_inches='tight')
    
    modval2=modval.copy()
    parameters=modval2.keys()
    means=[mu_R0,mu_U0,mu_US,mu_mean_dTd,mu_V0,mu_VS,mu_W0,mu_Th0]
    sigmas=[sigmaR0,sigmaU0,sigmaUS,sigma_mean_dTd,sigmaV0,sigmaVS,sigmaW0,sigmaTh0]
    for i in range(len(parameters)):
        modval2[parameters[i]]=[]
        modval2[parameters[i]].append(modval[parameters[i]])        
        modval2[parameters[i]].append(means[i])
        modval2[parameters[i]].append(sigmas[i])
    para=['$R_0 \, \\rm{(kpc)}$','$U_{\sun} \, \\rm{(km \, s^{-1})}$','$U_s  \, \\rm{(km \, s^{-1})} $','d$\\Theta$/d$R  \, \\rm{(km \, s^{-1} \, kpc^{-1})}$','$V_{\sun}  \, \\rm{(km \, s^{-1})} $','$V_s  \, \\rm{(km \, s^{-1})} $','$W_{\sun}  \, \\rm{(km \, s^{-1})} $','$\\Theta_0  \, \\rm{(km \, s^{-1})} $']
    if (latex_tables=='yes'):
        print('---------------------------------------------------------------')
        print('------------Final Values---------------------------------------------')
        print('---------------------------------------------------------------')
        print('Parameter & Initial & BeSSeL & & & &&')
        print(' & Value & Sample & Sample & Sample & Sample &&')    
        print('\\hline')
        for j in range(len(parameters)):
            print('%s %s %s %s' % (para[j]+'&', str(modval2[parameters[j]][0])+'  &$',str(round(modval2[parameters[j]][1],2))+' \pm',str(round(modval2[parameters[j]][2],2))+'$ & & & &&'))    
    

    
##########################################################################################################################
###################################     Evolution. Paramters with nsources      ##########################################
##########################################################################################################################    
if (plots_confirmation_nsources=='yes'):     
    f = raw_input("Expected quality of data (bad, good), if you don't know then bad:")
    if f=='good':
        f='f01'
    else:
        f='none'
    means={}
    for keypar in fitpars:
        means[keypar]=[]
    for i in range(len(dic_table_names)-1):#putting out teh last = complete sample
        for keypar in fitpars:
            if keypar != 'void':
                means[keypar].append(np.mean(dic[dic_table_names[i][2:]][keypar]))

    mean_complete={}
    for keypar in fitpars:
        mean_complete[keypar]=[]
    for keypar in fitpars:            
        mean_complete[keypar].append(np.mean(dic[dic_table_names[-1][2:]][keypar]))
        
    errors={}
    for keypar in fitpars:
        errors[keypar]=[]
    for i in range(len(dic_table_names)-1):#putting out teh last = complete sample
        for keypar in fitpars:
            if keypar != 'void':
                errors[keypar].append(np.std(dic[dic_table_names[i][2:]][keypar]))        

    errors_complete={}
    for keypar in fitpars:
        errors_complete[keypar]=[] 
    for keypar in fitpars:
        errors_complete[keypar].append(np.std(dic[dic_table_names[-1][2:]][keypar]))                       
         
    if(b_w_n=='b'):
        x_tickts=[0,100,200,300,400,500]
        #y_tickts=[0.4,0.7,1.0,1.3,1.6]
        xpos_label=250
        #ylim=[0.5,1.5]
        xlim=[0,550]
        #ylim_4=[-3.5,5.5]        
        #ylim_8=[-7,8] 
        #y_tickts_4=[-3.0,-1.0,1.0,3.0,5.0]      
        #y_tickts_8=[7,4,1,-2,-5]
        N_BeSSeL=110
        N_BeSSeL_past=30
    elif(b_w_n=='w'):
        x_tickts=[0,250,500,750,1000]
        #y_tickts=[-1.8,-1.1,-0.4,0.3,1.0,1.7]
        xpos_label=500
        #ylim=[0.25,1.75]
        xlim=[0,1050]
        #ylim_4=[-43,7]        
        #ylim_8=[-150,11]
        #y_tickts_4=[1.0,-10,-20,-30,-40]
        #y_tickts_8=[-120,-80,-40,1.0]
        N_BeSSeL=120
        N_BeSSeL_past=30
    else:
        x_tickts=[]
        #y_tickts=[]
        xpos_label=440
        ylim=[]
        xlim=[]
        ylim_4=[]        
        ylim_8=[]
        y_tickts_4=[]
        y_tickts_8=[]

    normed_R0=means['R0']
    normed_errR0=errors['R0']
    normed_Th0=means['Th0']
    normed_errTh0=errors['Th0']    
    normed_US=means['US']
    normed_errUS=errors['US']
    normed_VS=means['VS']
    normed_errVS=errors['VS']
    normed_Usun=means['U0']
    normed_errUsun=errors['U0']
    normed_Vsun=means['V0']
    normed_errVsun=errors['V0']
    normed_Wsun=means['W0']
    normed_errWsun=errors['W0']
    normed_dTd=means['dTd']
    normed_errdTd=errors['dTd']

    pconf= raw_input("Do you want to include 16 sources results? (yes or no):")                                                
    if (pconf == 'yes'):
        #Run only once if 
        N_bri_sources.append(16)
        #normed_R0[0]= 8.34 #1000 galaxies results
        #normed_errR0[0]=  0.08  #1000 galaxies results
        normed_R0.append(8.34)   #16 sources 1000 gaalxies 
        normed_errR0.append(0.24) #16 sources 1000 gaalxies 
    
        #normed_Th0[0]= 239.8 #1000 galaxies results
        #normed_errTh0[0]=  4.2  #1000 galaxies results
        normed_Th0.append(240.32)   #16 sources 1000 gaalxies 
        normed_errTh0.append(9.11) #16 sources 1000 gaalxies     
        
        #normed_US[0]= 2.9 #1000 galaxies results
        #normed_errUS[0]=  1.6  #1000 galaxies results    
        normed_US.append(2.99)  #16 sources 1000 gaalxies 
        normed_errUS.append(3.66) #16 sources 1000 gaalxies     
    
        #normed_VS[0]= 1.0 #1000 galaxies results
        #normed_errVS[0]=  4.9  #1000 galaxies results    
        normed_VS.append(1.21)   #16 sources 1000 gaalxies 
        normed_errVS.append(9.78) #16 sources 1000 gaalxies         
    
        #normed_Usun[0]= 10.5 #1000 galaxies results
        #normed_errUsun[0]=  1.2  #1000 galaxies results    
        normed_Usun.append(10.79)   #16 sources 1000 gaalxies 
        normed_errUsun.append(1.58) #16 sources 1000 gaalxies             
        
        #normed_Vsun[0]= 15.5 #1000 galaxies results
        #normed_errVsun[0]=  3.1  #1000 galaxies results        
        normed_Vsun.append(15.20)   #16 sources 1000 gaalxies 
        normed_errVsun.append(2.54) #16 sources 1000 gaalxies
    
        #normed_Wsun[0]= 8.7 #1000 galaxies results
        #normed_errWsun[0]=  0.5  #1000 galaxies results        
        normed_Wsun.append(8.13)  #16 sources 1000 gaalxies 
        normed_errWsun.append(0.77) #16 sources 1000 gaalxies
        
        #normed_dTd[0]= -0.2 #1000 galaxies results
        #normed_errdTd[0]=  1.1  #1000 galaxies results    
        normed_dTd.append(-0.03)   #16 sources 1000 gaalxies 
        normed_errdTd.append(3.24) #16 sources 1000 gaalxies                 
                                                                                                
    #Plot
    #N_bri_sources=[103, 150, 200, 300, 400]
    matplotlib.rcParams.update({'font.size': 16})
    fig1 = plt.figure(figsize=(15,10))
    ax = fig1.add_subplot(111)        
    #ax.axvline(0.755,color='k',linestyle='dashed',linewidth=3)
    X=np.linspace(0,1300,100)
                        
    #R0
    ax1 = fig1.add_subplot(241)
    #ax1.set_yscale('log')
    ax1.errorbar(N_bri_sources,normed_R0,yerr=normed_errR0,color=colork2[1],linewidth=2,fmt="^",label='$R_0$ & $\\Delta R_0$',capsize=4, elinewidth=4,markersize=8)
    ax1.errorbar(N_BeSSeL,modval['R0'],yerr=final_unc['R0'],fmt="*",color='k',label='Model A5',markersize=12,capsize=4, elinewidth=4)
    ax1.errorbar(N_BeSSeL_past,8.40,yerr=0.36,fmt="*",color='cadetblue',label='Fit 3',markersize=12,capsize=4, elinewidth=4)
    ax1.set_xticks(x_tickts)     
    Y2=0*X+((mean_complete['R0'][0])+(errors_complete['R0'][0]))
    Y1=0*X+((mean_complete['R0'][0])-(errors_complete['R0'][0]))
    ax1.grid(True)
    plt.setp(ax1.get_xticklabels(), visible=False) 
    if(b_w_n=='b'):
        if(f=='none'):
            dum=ax1.get_ylim()
            dum2=(((dum[1]-dum[0])/5.0)*4.5)+dum[0]
            #ax1.text( xpos_label-20, 9, '$R_0$ $\\rm{(kpc)}$', verticalalignment='center', horizontalalignment='center', color=colork2[1], fontsize=27)                    
            ax1.text( xpos_label-20, dum2, '$R_0$ $\\rm{(kpc)}$', verticalalignment='center', horizontalalignment='center', color=colork2[1], fontsize=27)                                
        if(f=='f01'):
            ax1.set_ylim([8.0,8.79]) #f01
            dum=ax1.get_ylim()
            dum2=(((dum[1]-dum[0])/5.0)*4.5)+dum[0]
            y_tickts=np.arange(8.0,8.8,0.1) #f01
            #ax1.text( xpos_label+150, 8.68, '$R_0$ $\\rm{(kpc)}$', verticalalignment='center', horizontalalignment='center', color=colork2[1], fontsize=27)         #f01        
            ax1.text( xpos_label+150, dum2, '$R_0$ $\\rm{(kpc)}$', verticalalignment='center', horizontalalignment='center', color=colork2[1], fontsize=27)         #f01        
            ax1.set_yticks(y_tickts)  
    if(b_w_n=='w'):
        if(f=='none'):  
            ax1.set_ylim([8.01,9.3]) 
            dum=ax1.get_ylim()
            dum2=(((dum[1]-dum[0])/5.0)*4.5)+dum[0]       
            #ax1.text( xpos_label-20, 9.1, '$R_0$ $\\rm{(kpc)}$', verticalalignment='center', horizontalalignment='center', color=colork2[1], fontsize=27)                 
            ax1.text( xpos_label-20, dum2, '$R_0$ $\\rm{(kpc)}$', verticalalignment='center', horizontalalignment='center', color=colork2[1], fontsize=27)                 
        if(f=='f01'):
            ax1.set_ylim([8.0,8.79]) #f01
            dum=ax1.get_ylim()
            dum2=(((dum[1]-dum[0])/5.0)*4.5)+dum[0]        
            y_tickts=np.arange(8.0,8.8,0.1) #f01
            ax1.text( xpos_label+150, dum2, '$R_0$ $\\rm{(kpc)}$', verticalalignment='center', horizontalalignment='center', color=colork2[1], fontsize=27)         #f01        
            #ax1.text( xpos_label+150, 8.68, '$R_0$ $\\rm{(kpc)}$', verticalalignment='center', horizontalalignment='center', color=colork2[1], fontsize=27)         #f01        
            ax1.set_yticks(y_tickts)  

        ax1.fill_between(X, Y1,Y2,color='grey',alpha=0.3)
    ax1.axhline(y=modval['R0'],color='k',linestyle='dashed',linewidth=3)
    
    #Th0
    ax2 = fig1.add_subplot(242)
    ax2.errorbar(N_bri_sources,normed_Th0,yerr=normed_errTh0,color=colork2[2],linewidth=2,fmt="^",label='$R_0$ & $\\Delta R_0$',capsize=4, elinewidth=4,markersize=8)
    ax2.errorbar(N_BeSSeL,modval['Th0'],yerr=final_unc['Th0'],fmt="*",color='k',label='Current BeSSeL',markersize=12,capsize=4, elinewidth=4)
    ax2.errorbar(N_BeSSeL_past,254,yerr=16,fmt="*",color='cadetblue',label='Fit 3',markersize=12,capsize=4, elinewidth=4)
    ax2.set_xticks(x_tickts)   
    Y2=0*X+((mean_complete['Th0'][0])+(errors_complete['Th0'][0]))
    Y1=0*X+((mean_complete['Th0'][0])-(errors_complete['Th0'][0]))    
    ax2.grid(True)
    plt.setp(ax2.get_xticklabels(), visible=False) 
    #plt.setp(ax2.get_yticklabels(), visible=False) 
    if(b_w_n=='b'):
        if(f=='none'):        
            ax2.set_ylim([230,296])
            dum=ax2.get_ylim()
            dum2=(((dum[1]-dum[0])/5.0)*4.5)+dum[0]
            y_tickts=np.arange(230,300,10)
            ax2.text( xpos_label-20, dum2, '$\\Theta_0$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[2], fontsize=27)                         
            #ax2.text( xpos_label-20, 286, '$\\Theta_0$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[2], fontsize=27)                         
            ax2.set_yticks(y_tickts)    
        if(f=='f01'):   
            ax2.set_ylim([230,270.3]) # f01            
            dum=ax2.get_ylim()
            dum2=(((dum[1]-dum[0])/5.0)*4.5)+dum[0]
            y_tickts=np.arange(230,270,5) #f01
            #ax2.text( xpos_label+140, 265, '$\\Theta_0$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[2], fontsize=27)                 #f01
            ax2.text( xpos_label+140, dum2, '$\\Theta_0$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[2], fontsize=27)                 #f01
            ax2.set_yticks(y_tickts)    
    if(b_w_n=='w'):
        if(f=='none'):        
            ax2.set_ylim([225,296])
            dum=ax2.get_ylim()
            dum2=(((dum[1]-dum[0])/5.0)*4.5)+dum[0]            
            y_tickts=np.arange(230,300,10)
            #ax2.text( xpos_label-20, 285, '$\\Theta_0$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[2], fontsize=27)                                 
            ax2.text( xpos_label-20, dum2, '$\\Theta_0$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[2], fontsize=27)                                 
            ax2.set_yticks(y_tickts) 
        if(f=='f01'):                    
            ax2.set_ylim([230,270.3]) # f01
            dum=ax2.get_ylim()
            dum2=(((dum[1]-dum[0])/5.0)*4.5)+dum[0]            
            y_tickts=np.arange(230,270,5) #f01
            ax2.text( xpos_label+140, dum2, '$\\Theta_0$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[2], fontsize=27)                 #f01
            #ax2.text( xpos_label+140, 265, '$\\Theta_0$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[2], fontsize=27)                 #f01
            ax2.set_yticks(y_tickts)    
        ax2.fill_between(X, Y1,Y2,color='grey',alpha=0.3)
    ax2.axhline(y=modval['Th0'],color='k',linestyle='dashed',linewidth=3)

    #dTd    
    ax3 = fig1.add_subplot(243)
    normed_dTd = [x for x in means['dTd']]
    normed_errdTd = [x for x in errors['dTd']]
    ax3.errorbar(N_bri_sources,normed_dTd,yerr=normed_errdTd,color=colork2[8],linewidth=2,fmt="^",capsize=4, elinewidth=4,markersize=8)
    ax3.errorbar(N_BeSSeL,modval['dTd'],yerr=final_unc['dTd'],fmt="*",color='k',markersize=12,capsize=4, elinewidth=4)
    ax3.errorbar(N_BeSSeL_past,2.3,yerr=2.1,fmt="*",color='cadetblue',label='Fit 3',markersize=12,capsize=4, elinewidth=4)
    #ax3.set_ylim(ylim_8)
    #ax3.set_yticks(y_tickts_8)    
    #ax3.yaxis.tick_right()
    #ax3.yaxis.set_label_position('right')    
    Y2=0*X+((mean_complete['dTd'][0])+(errors_complete['dTd'][0]))
    Y1=0*X+((mean_complete['dTd'][0])-(errors_complete['dTd'][0]))
    ax3.grid(True)
    #ax3.tick_params(axis='y',labelsize=14)
    plt.setp(ax3.get_xticklabels(), visible=False)     
    if(b_w_n=='b'):
        if(f=='none'):  
            dum=ax3.get_ylim()
            dum2=(((dum[1]-dum[0])/5.0)*4.5)+dum[0]                  
            y_tickts=np.arange(0,30,5)
            ax3.set_yticks(y_tickts)    
            #ax3.text( xpos_label+5, 25, '$\\frac{d\\Theta}{dR}$ $\\rm{(km/s \\cdot kpc)}$', verticalalignment='center', horizontalalignment='center', color=colork2[8], fontsize=27)                     
            ax3.text( xpos_label+5, dum2, '$\\frac{d\\Theta}{dR}$ $\\rm{(km/s \\cdot kpc)}$', verticalalignment='center', horizontalalignment='center', color=colork2[8], fontsize=27)                     
        if(f=='f01'):     
            ax3.set_ylim([-3.5,5.0]) #f01
            dum=ax3.get_ylim()
            dum2=(((dum[1]-dum[0])/5.0)*4.5)+dum[0]
            y_tickts=np.arange(-3,5.0,1) #f01
            ax3.set_yticks(y_tickts)    
            #ax3.text( xpos_label+50, 4, '$\\frac{d\\Theta}{dR}$ $\\rm{(km/s \\cdot kpc)}$', verticalalignment='center', horizontalalignment='center', color=colork2[8], fontsize=27)                     #f01
            ax3.text( xpos_label+50, dum2, '$\\frac{d\\Theta}{dR}$ $\\rm{(km/s \\cdot kpc)}$', verticalalignment='center', horizontalalignment='center', color=colork2[8], fontsize=27)                     #f01
    if(b_w_n=='w'):
        if(f=='none'):        
            ax3.set_ylim([-4,32])
            dum=ax3.get_ylim()
            dum2=(((dum[1]-dum[0])/5.0)*4.5)+dum[0]            
            #ax3.text( xpos_label+20, 27.5, '$\\frac{d\\Theta}{dR}$ $\\rm{(km/s \\cdot kpc)}$', verticalalignment='center', horizontalalignment='center', color=colork2[8], fontsize=27)                     
            ax3.text( xpos_label+20, dum2, '$\\frac{d\\Theta}{dR}$ $\\rm{(km/s \\cdot kpc)}$', verticalalignment='center', horizontalalignment='center', color=colork2[8], fontsize=27)                     
        if(f=='f01'):                    
            ax3.set_ylim([-3.5,5.0]) #f01
            dum=ax3.get_ylim()
            dum2=(((dum[1]-dum[0])/5.0)*4.5)+dum[0]            
            y_tickts=np.arange(-3,5.0,1) #f01
            ax3.set_yticks(y_tickts)    
            ax3.text( xpos_label+50, dum2, '$\\frac{d\\Theta}{dR}$ $\\rm{(km/s \\cdot kpc)}$', verticalalignment='center', horizontalalignment='center', color=colork2[8], fontsize=27)                     #f01
            #ax3.text( xpos_label+50, 4, '$\\frac{d\\Theta}{dR}$ $\\rm{(km/s \\cdot kpc)}$', verticalalignment='center', horizontalalignment='center', color=colork2[8], fontsize=27)                     #f01
        ax3.fill_between(X, Y1,Y2,color='grey',alpha=0.3)
    ax3.axhline(y=modval['dTd'],color='k',linestyle='dashed',linewidth=3)   
    
    #U_sun
    ax4 = fig1.add_subplot(244)#,sharey=ax1)
    ax4.errorbar(N_bri_sources,normed_Usun,yerr=normed_errUsun,color=colork2[5],linewidth=2,fmt="^",label='$R_0$ & $\\Delta R_0$',capsize=4, elinewidth=4,markersize=8)
    ax4.errorbar(N_BeSSeL,modval['U0'],yerr=final_unc['U0'],fmt="*",color='k',label='Current BeSSeL',markersize=12,capsize=4, elinewidth=4)
    #ax4.set_yticks(y_tickts)
    #ax4.set_ylim(ylim)
    Y2=0*X+((mean_complete['U0'][0])+(errors_complete['U0'][0]))
    Y1=0*X+((mean_complete['U0'][0])-(errors_complete['U0'][0]))
    ax4.grid(True)
    #y_tickts=np.arange(230,290,10)
    #ax4.set_yticks(y_tickts)        
    plt.setp(ax4.get_xticklabels(), visible=False) 
    if(b_w_n=='b'):
        if(f=='none'):
            dum=ax4.get_ylim()
            dum2=(((dum[1]-dum[0])/5.0)*4.5)+dum[0]                    
            y_tickts=np.arange(8,15,1)
            ax4.set_yticks(y_tickts)  
            #ax4.text( xpos_label-35, 14, '$U_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[5], fontsize=27)                         
            ax4.text( xpos_label-35, dum2, '$U_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[5], fontsize=27)                         
        if(f=='f01'):                    
            ax4.set_ylim([8.7,12.7]) #f01
            dum=ax4.get_ylim()
            dum2=(((dum[1]-dum[0])/5.0)*4.5)+dum[0]            
            y_tickts=np.arange(9,13,0.5) #f01
            ax4.set_yticks(y_tickts)    
            #ax4.text( xpos_label+125, 12.2, '$U_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[5], fontsize=27)                 #f01
            ax4.text( xpos_label+125, dum2, '$U_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[5], fontsize=27)                 #f01
    if(b_w_n=='w'):
        if(f=='none'):        
            ax4.set_ylim([8.5,14.2]) #f01
            dum=ax4.get_ylim()
            dum2=(((dum[1]-dum[0])/5.0)*4.5)+dum[0]            
            #ax4.text( xpos_label-15, 13.5, '$U_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[5], fontsize=27)                         
            ax4.text( xpos_label-15, dum2, '$U_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[5], fontsize=27)                         
        if(f=='f01'):                    
            ax4.set_ylim([8.7,12.7]) #f01
            dum=ax4.get_ylim()
            dum2=(((dum[1]-dum[0])/5.0)*4.5)+dum[0]            
            y_tickts=np.arange(9,13,0.5) #f01
            ax4.set_yticks(y_tickts)    
            #ax4.text( xpos_label+125, 12.2, '$U_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[5], fontsize=27)                 #f01
            ax4.text( xpos_label+125, dum2, '$U_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[5], fontsize=27)                 #f01
        ax4.fill_between(X, Y1,Y2,color='grey',alpha=0.3)
    ax4.axhline(y=modval['U0'],color='k',linestyle='dashed',linewidth=3)     

    #V_sun
    ax5 = fig1.add_subplot(245,sharex=ax1)
    normed_Vsun = [x for x in means['V0']]
    normed_errVsun = [x  for x in errors['V0']]                
    #ax5.set_yscale('log',nonposy='clip')
    ax5.errorbar(N_bri_sources,normed_Vsun,yerr=normed_errVsun,color=colork2[6],linewidth=2,fmt="^",capsize=4, elinewidth=4,markersize=8)
    ax5.errorbar(N_BeSSeL,modval['V0'],yerr=final_unc['V0'],fmt="*",color='k',markersize=12,capsize=4, elinewidth=4)
    ax5.set_xlim(xlim)
    Y2=0*X+((mean_complete['V0'][0])+(errors_complete['V0'][0]))
    Y1=0*X+((mean_complete['V0'][0])-(errors_complete['V0'][0]))
    ax5.grid(True)
    y_tickts=np.arange(8,24,2)
    ax5.set_yticks(y_tickts)  
    ax5.axhline(y=modval['V0'],color='k',linestyle='dashed',linewidth=3,label='Initial Values')                               
    #plt.setp(ax5.get_yticklabels(), visible=False) 
    if(b_w_n=='b'):
        if(f=='none'):  
            dum=ax5.get_ylim()
            dum2=(((dum[1]-dum[0])/5.0)*4.5)+dum[0]      
            #ax5.text( xpos_label+40,10, '$V_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[6], fontsize=27)                      
            ax5.text( xpos_label+40,dum2, '$V_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[6], fontsize=27)                      
        if(f=='f01'):                    
            ax5.set_ylim([8.1,23]) 
            dum=ax5.get_ylim()
            dum2=(((dum[1]-dum[0])/5.0)*4.5)+dum[0]            
            #ax5.text( xpos_label+125,21.2, '$V_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[6], fontsize=27)                 #f01
            ax5.text( xpos_label+125,dum2, '$V_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[6], fontsize=27)                 #f01
        #p1 = Rectangle((0, 0), 1, 1, fc="grey",alpha=0.3)
        #ax5.legend([p1],['Complete Simulation'],bbox_to_anchor=(1.3, -0.17),numpoints=1,ncol=2)                    
        #ax5.legend(bbox_to_anchor=(0.2, -0.17),numpoints=1,ncol=2) #f01        
    if(b_w_n=='w'):
        if(f=='none'):        
            ax5.set_ylim([7.0,23])
            dum=ax5.get_ylim()
            dum2=(((dum[1]-dum[0])/5.0)*4.5)+dum[0]
            #ax5.text( xpos_label+20,9, '$V_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[6], fontsize=27)                                    
            ax5.text( xpos_label+20,dum2, '$V_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[6], fontsize=27)                                        
        if(f=='f01'):                    
            ax5.set_ylim([8.1,23]) 
            dum=ax5.get_ylim()
            dum2=(((dum[1]-dum[0])/5.0)*4.5)+dum[0]            
            #ax5.text( xpos_label+125,21.2, '$V_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[6], fontsize=27)                 #f01
            ax5.text( xpos_label+125,dum2, '$V_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[6], fontsize=27)                 #f01
        ax5.fill_between(X, Y1,Y2,color='grey',alpha=0.3)
        #p1 = Rectangle((0, 0), 1, 1, fc="grey",alpha=0.3)
        #ax5.legend([p1],['Complete Simulation'],bbox_to_anchor=(0.1, -0.17),numpoints=1,ncol=2)        
        #p1 = Rectangle((0, 0), 1, 1, fc="grey",alpha=0.3)
        #ax5.legend([p1,],['Complete Simulation'],bbox_to_anchor=(1.3, -0.17),numpoints=1,ncol=2)        
        #ax5.legend(bbox_to_anchor=(2.4, -0.17),numpoints=1,ncol=2)                


    #W_sun
    ax6 = fig1.add_subplot(246,sharex=ax2)    
    normed_Wsun = [x for x in means['W0']]
    normed_errWsun = [x  for x in errors['W0']]          
    ax6.errorbar(N_bri_sources,normed_Wsun,yerr=normed_errWsun,color=colork2[7],linewidth=2,fmt="^",capsize=4, elinewidth=4,markersize=8)
    ax6.errorbar(N_BeSSeL, modval['W0'],yerr=final_unc['W0'],fmt="*",color='k',markersize=12,capsize=4, elinewidth=4)
    #ax6.set_ylim()   
    ax6.set_xlim(xlim)
    ax6.set_ylim([7.9,10.0])
    Y2=0*X+((mean_complete['W0'][0])+(errors_complete['W0'][0]))
    Y1=0*X+((mean_complete['W0'][0])-(errors_complete['W0'][0]))
    ax6.grid(True)
    y_tickts=np.arange(7.4,10.0,0.4)
    ax6.set_yticks(y_tickts)        
    ax6.axhline(y=modval['W0'],color='k',linestyle='dashed',linewidth=3,label='Initial Values')    
    #plt.setp(ax6.get_yticklabels(), visible=False) 
    if(b_w_n=='b'):
        if(f=='none'):        
            ax6.set_ylim([7.3,10])
            dum=ax6.get_ylim()
            dum2=(((dum[1]-dum[0])/5.0)*4.5)+dum[0]
            y_tickts=np.arange(7.4,10.0,0.4)
            ax6.set_yticks(y_tickts)        
            #ax6.text( xpos_label+50,7.7, '$W_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[7], fontsize=27)               
            ax6.text( xpos_label+50,dum2, '$W_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[7], fontsize=27)               
        if(f=='f01'):                    
            ax6.set_ylim([7.3,10])
            dum=ax6.get_ylim()
            dum2=(((dum[1]-dum[0])/5.0)*4.5)+dum[0]                        
            #ax6.text( xpos_label+120,9.65, '$W_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[7], fontsize=27)               #f01
            ax6.text( xpos_label+120,dum2, '$W_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[7], fontsize=27)               #f01
        #ax6.legend(bbox_to_anchor=(0.35, -0.17),numpoints=1,ncol=2)            
        #ax6.legend(bbox_to_anchor=(0.0, -0.17),numpoints=1,ncol=2)
    if(b_w_n=='w'):
        if(f=='none'):
            ax6.set_ylim([7.3,10])                    
            dum=ax6.get_ylim()
            dum2=(((dum[1]-dum[0])/5.0)*4.5)+dum[0]            
            #ax6.text( xpos_label+30,7.6, '$W_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[7], fontsize=27)               
            ax6.text( xpos_label+30,dum2, '$W_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[7], fontsize=27)               
        if(f=='f01'):                    
            ax6.set_ylim([7.3,10])
            dum=ax6.get_ylim()
            dum2=(((dum[1]-dum[0])/5.0)*4.5)+dum[0]            
            #ax6.text( xpos_label+120,9.72, '$W_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[7], fontsize=27)               #f01
            ax6.text( xpos_label+120,dum2, '$W_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[7], fontsize=27)               #f01
        #ax6.legend(bbox_to_anchor=(0.0, -0.17),numpoints=1,ncol=2)
        ax6.fill_between(X, Y1,Y2,color='grey',alpha=0.3)
            
    #U_s
    ax7 = fig1.add_subplot(247,sharex=ax3)
    ax7.errorbar(N_bri_sources,normed_US,yerr=normed_errUS,color=colork2[3],linewidth=2,fmt="^",capsize=4, elinewidth=4,markersize=8)
    p2=ax7.errorbar(N_BeSSeL, modval['US'],yerr=final_unc['US'],fmt="*",color='k',label='Model A5',markersize=12,capsize=4, elinewidth=4)
    p3=ax7.errorbar(N_BeSSeL_past,2.3,yerr=2.1,fmt="*",color='cadetblue',label='Fit 3',markersize=12,capsize=4, elinewidth=4)
    #ax7.set_ylim(ylim)   
    ax7.set_xticks(x_tickts) 
    ax7.set_xlim(xlim)          
    Y2=0*X+((mean_complete['US'][0])+(errors_complete['US'][0]))
    Y1=0*X+((mean_complete['US'][0])-(errors_complete['US'][0]))
    ax7.grid(True)
    p4=ax7.axhline(y=modval['US'],color='k',linestyle='dashed',linewidth=3,label='Initial Values')            
    #plt.setp(ax7.get_yticklabels(), visible=False) 
    if(b_w_n=='b'):
        if(f=='none'):        
            ax7.set_ylim([-8,7])
            dum=ax7.get_ylim()
            dum2=(((dum[1]-dum[0])/5.0)*4.5)+dum[0]            
            y_tickts=np.arange(-6,8,2)
            ax7.text( xpos_label-20, dum2, '$\\bar{U}_s$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[3], fontsize=27)               
            #ax7.text( xpos_label-20, -6, '$\\bar{U}_s$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[3], fontsize=27)               
            ax7.set_yticks(y_tickts)
        if(f=='f01'):                    
            ax7.set_ylim([-1,7]) #f01
            dum=ax7.get_ylim()
            dum2=(((dum[1]-dum[0])/5.0)*4.5)+dum[0]              
            y_tickts=np.arange(0,7,1) #f01
            #ax7.text( xpos_label+140, 6, '$\\bar{U}_s$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[3], fontsize=27)               #f01
            ax7.text( xpos_label+140, dum2, '$\\bar{U}_s$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[3], fontsize=27)               #f01
            ax7.set_yticks(y_tickts)        
        p1 = Rectangle((0, 0), 1, 1, fc="grey",alpha=0.3)
        ax7.legend([p1,p2,p3,p4],['Complete Simulation','Model A5','Fit 3','Initial Values'],numpoints=1,ncol=4,bbox_to_anchor=(1.6, -0.17))                                
    if(b_w_n=='w'):
        if(f=='none'):        
            ax7.set_ylim([-1,7]) #f01 
            dum=ax7.get_ylim()
            dum2=(((dum[1]-dum[0])/5.0)*4.5)+dum[0]            
            y_tickts=np.arange(0,7,1) #f01
            #ax7.text( xpos_label+90, 0, '$\\bar{U}_s$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[3], fontsize=27)               #f01
            ax7.text( xpos_label+90, dum2, '$\\bar{U}_s$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[3], fontsize=27)               #f01
            ax7.set_yticks(y_tickts)         
        if(f=='f01'):                    
            ax7.set_ylim([-1,7]) #f01
            dum=ax7.get_ylim()
            dum2=(((dum[1]-dum[0])/5.0)*4.5)+dum[0]            
            y_tickts=np.arange(0,7,1) #f01
            #ax7.text( xpos_label+140, 6, '$\\bar{U}_s$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[3], fontsize=27)               #f01
            ax7.text( xpos_label+140, dum2, '$\\bar{U}_s$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[3], fontsize=27)               #f01
            ax7.set_yticks(y_tickts)    
        ax7.fill_between(X, Y1,Y2,color='grey',alpha=0.3)
            
    #Vs
    ax8 = fig1.add_subplot(248,sharex=ax4)#,sharey=ax5)
    normed_VS = [x for x in means['VS']]
    normed_errVS = [x for x in errors['VS']]
    ax8.errorbar(N_bri_sources,normed_VS,yerr=normed_errVS,color=colork2[4],linewidth=2,fmt="^",capsize=4, elinewidth=4,markersize=8)
    ax8.errorbar(N_BeSSeL, modval['VS']*-1,yerr=final_unc['VS'],fmt="*",color='k',markersize=12,capsize=4, elinewidth=4)
    ax8.errorbar(N_BeSSeL_past,-6.5,yerr=1.8,fmt="*",color='cadetblue',label='Fit 3',markersize=12,capsize=4, elinewidth=4)
    if(b_w_n=='b'):
        ax8.arrow(45, -5.3, 0, -1.5, head_width=15, head_length=1, fc='k', ec='k')
        ax8.text(120,-6.5, '(-14.7)', verticalalignment='center', horizontalalignment='center', color='k', fontsize=14)                       #f01
    if(b_w_n=='w'):
        ax8.arrow(100, -6.0, 0, -1.5, head_width=15, head_length=1, fc='k', ec='k')
        ax8.text(250,-7.2, '(-14.7)', verticalalignment='center', horizontalalignment='center', color='k', fontsize=14)                       #f01
    #ax8.set_ylim(ylim_4)
    #ax8.set_yticks(y_tickts_4)
    ax8.set_xlim(xlim)
    ax8.set_xticks(x_tickts)       
    #ax8.yaxis.tick_right()
    #ax8.yaxis.set_label_position('right') 
    Y2=0*X+(-1*(mean_complete['VS'][0])+(errors_complete['VS'][0]))
    Y1=0*X+(-1*(mean_complete['VS'][0])-(errors_complete['VS'][0]))
    ax8.grid(True)
    #ax8.set_yscale("log", nonposy='clip')
    #ax8.tick_params(axis='y',labelsize=14)    
    if(b_w_n=='b'):
        if(f=='none'):        
            ax8.set_ylim([-60,13])  
            dum=ax8.get_ylim()
            dum2=(((dum[1]-dum[0])/5.0)*4.5)+dum[0]                      
            y_tickts=np.arange(-50,30,10)
            ax8.set_yticks(y_tickts)  
            #ax8.text(xpos_label-20,-50, '$\\bar{V}_s$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[4], fontsize=27)                       
            ax8.text(xpos_label-20,dum2, '$\\bar{V}_s$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[4], fontsize=27)                       
        if(f=='f01'):                    
            ax8.set_ylim([-8.9,11.5]) #f01
            dum=ax8.get_ylim()
            dum2=(((dum[1]-dum[0])/5.0)*4.5)+dum[0]            
            y_tickts=np.arange(-6,12,3) #f01
            ax8.set_yticks(y_tickts)    
            #ax8.text(xpos_label+150,9, '$\\bar{V}_s$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[4], fontsize=27)                       #f01
            ax8.text(xpos_label+150,dum2, '$\\bar{V}_s$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[4], fontsize=27)                       #f01
        #ax8.legend(bbox_to_anchor=(0.2, -0.17),numpoints=1,ncol=2)                    
    if(b_w_n=='w'):
        if(f=='none'):        
            ax8.set_ylim([-67,15])
            dum=ax8.get_ylim()
            dum2=(((dum[1]-dum[0])/5.0)*4.5)+dum[0]            
            y_tickts=np.arange(-60,20,10)
            ax8.set_yticks(y_tickts)            
            #ax8.text(xpos_label-60,-55, '$\\bar{V}_s$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[4], fontsize=27)                       
            ax8.text(xpos_label-60,dum2, '$\\bar{V}_s$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[4], fontsize=27)                       
        if(f=='f01'):         
            ax8.set_ylim([-11.8,11.5]) #f01
            dum=ax8.get_ylim()
            dum2=(((dum[1]-dum[0])/5.0)*4.5)+dum[0]            
            y_tickts=np.arange(-9,12,3) #f01  
            ax8.set_yticks(y_tickts)            
            #ax8.text(xpos_label+150,9, '$\\bar{V}_s$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[4], fontsize=27)                       #f01
            ax8.text(xpos_label+150,dum2, '$\\bar{V}_s$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[4], fontsize=27)                       #f01
        #ax8.legend(bbox_to_anchor=(0.7, -0.17),numpoints=1,ncol=2)        
        ax8.fill_between(X, Y1,Y2,color='grey',alpha=0.3)
    ax8.axhline(y=modval['VS'],color='k',linestyle='dashed',linewidth=3)    

    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.setp(ax.get_yticklabels(), visible=False)  
    plt.setp(ax.get_xticklabels(), visible=False)          
    ax.set_xlabel('Number of Sources',labelpad=25)
    if(b_w_n=='b'):
        ax.set_ylabel('SET A',labelpad=35,fontweight='bold',fontsize=33)    
    if(b_w_n=='w'):
        ax.set_ylabel('SET B',labelpad=35,fontweight='bold',fontsize=33)    
    #ax.set_ylabel('Normed Parameter',labelpad=40)
    ax.tick_params(
    axis='both',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off',right='off', left='off', labelleft='off')
    #fig1.patch.set_visible(False)
    fig1.tight_layout() 
    
    #plt.savefig(path_outputs+'/whole_sigma15%.pdf',bbox_inches='tight')      
    #plt.savefig(path_outputs+'/bessel_nogauss.pdf',bbox_inches='tight')          
    #plt.savefig(path_outputs+'/bessel_gauss_distance.pdf',bbox_inches='tight')   
    #plt.savefig(path_outputs+'/bessel_gdist_betterunc.pdf',bbox_inches='tight')   
    #plt.savefig(path_outputs+'/bessel_gdist_nicess10.pdf',bbox_inches='tight')  
    #plt.savefig(path_outputs+'/bessel_gpara_-100%new.pdf',bbox_inches='tight')                     
    #plt.savefig(path_outputs+'/bessel_gpara_20%new.pdf',bbox_inches='tight')                     
    plt.savefig(path_outputs+'/evo_'+name_ouputfolder+'.pdf',bbox_inches='tight')
    #Printing Results
    modval2=modval.copy()
    parameters=modval2.keys()
    for i in range(len(parameters)):
        modval2[parameters[i]]=[]
        modval2[parameters[i]].append(modval[parameters[i]])        
        #modval2[parameters[i]].append(means[i])
        #modval2[parameters[i]].append(sigmas[i])
    para=['$R_0 \, \\rm{(kpc)}$','$U_{\sun} \, \\rm{(km \, s^{-1})}$','$U_s  \, \\rm{(km \, s^{-1})} $','d$\\Theta$/d$R  \, \\rm{(km \, s^{-1} \, kpc^{-1})}$','$V_{\sun}  \, \\rm{(km \, s^{-1})} $','$V_s  \, \\rm{(km \, s^{-1})} $','$W_{\sun}  \, \\rm{(km \, s^{-1})} $','$\\Theta_0  \, \\rm{(km \, s^{-1})} $']        
    if (latex_tables=='yes'):
        print('---------------------------------------------------------------')
        print('------------Final Values---------------------------------------------')
        print('---------------------------------------------------------------')
        print('Parameter & Initial & Similar &'+str (N_bri_sources[1]) +'&'+ str (N_bri_sources[2]) +'&'+str (N_bri_sources[3])+'&'+str (N_bri_sources[4])+'&'+str (N_bri_sources[5])+' & Complete &&')
        print(' & Value & BeSSeL & Sample & Sample & Sample & Sample & Sample & Sample &&')    
        print('\\hline')
        for j in range(len(parameters)):
            if (parameters[j]=='R0'):
                round_val=2
            else:
                round_val=1        
            #print ('%s' %(round_val))
            print('%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s' % (para[j]+'&',
            str(modval2[parameters[j]][0])+'  &$',
            str(round(means[parameters[j]][0],round_val))+' \pm',str(round(errors[parameters[j]][0],round_val))+'  $&$',
            str(round(means[parameters[j]][1],round_val))+' \pm',str(round(errors[parameters[j]][1],round_val))+'  $&$',
            str(round(means[parameters[j]][2],round_val))+' \pm',str(round(errors[parameters[j]][2],round_val))+'  $&$',
            str(round(means[parameters[j]][3],round_val))+' \pm',str(round(errors[parameters[j]][3],round_val))+'  $&$',
            str(round(means[parameters[j]][4],round_val))+' \pm',str(round(errors[parameters[j]][4],round_val))+'  $&$',
            str(round(means[parameters[j]][5],round_val))+' \pm',str(round(errors[parameters[j]][5],round_val))+'  $&$',
            str(round(mean_complete[parameters[j]][0],round_val))+' \pm',str(round(errors[parameters[j]][0],round_val))+'$  &&'))
        print('\\hline')


  

    

##########################################################################################################################
###################################     Pearson Coefficents Plots   ######################################################
##########################################################################################################################
if (plots_confirmation_nsources=='yes'):                
    pearso_coeffs=['R_0-T_0','R_0-d_T','R_0-U_0','R_0-V_0','R_0-W_0','R_0-V_s','R_0-U_s',
                                'T_0-d_T','T_0-U_0','T_0-V_0','T_0-W_0','T_0-V_s','T_0-U_s',
                                        'd_T-U_0','d_T-V_0','d_T-W_0','d_T-V_s','d_T-U_s',
                                                    'U_0-V_0','U_0-W_0','U_0-V_s','U_0-U_s',
                                                            'V_0-W_0','V_0-V_s','V_0-U_s',
                                                                        'W_0-V_s','W_0-U_s',
                                                                                'V_s-U_s']
    pearso_coeffs_labels=['$R_0$-$\\Theta_0$','$R_0$-$d\\Theta/dR$',       '$R_0$-$U_{\\odot}$',        '$R_0$-$V_{\\odot}$',        '$R_0$-$W_{\\odot}$',        '$R_0$-$\\bar{V}_{s}$',        '$R_0$-$\\bar{U}_{s}$',
                                                '$\\Theta_0$-$d\\Theta/dR$', '$\\Theta_0$-$U_{\\odot}$',  '$\\Theta_0$-$V_{\\odot}$',  '$\\Theta_0$-$W_{\\odot}$',  '$\\Theta_0$-$\\bar{V}_{s}$',  '$\\Theta_0$-$\\bar{U}_{s}$',
                                                                            '$d\\Theta/dR$-$U_{\\odot}$','$d\\Theta/dR$-$V_{\\odot}$','$d\\Theta/dR$-$W_{\\odot}$','$d\\Theta/dR$-$\\bar{V}_{s}$','$d\\Theta/dR$-$\\bar{U}_{s}$',                 
                                                                                                        '$U_{\\odot}$-$V_{\\odot}$', '$U_{\\odot}$-$W_{\\odot}$', '$U_{\\odot}$-$\\bar{V}_{s}$', '$U_{\\odot}$-$\\bar{U}_{s}$', 
                                                                                                                                        '$V_{\\odot}$-$W_{\\odot}$', '$V_{\\odot}$-$\\bar{V}_{s}$', '$V_{\\odot}$-$\\bar{U}_{s}$',
                                                                                                                                                                    '$W_{\\odot}$-$\\bar{V}_{s}$', '$W_{\\odot}$-$\\bar{U}_{s}$',                                                                                                                                                                                                '$\\bar{V}_{s}$-$\\bar{U}_{s}$']
    def find_corr(k):
        if(k[0:3]=='R_0'):
            colu=0
        if(k[0:3]=='T_0'):
            colu=1  
        if(k[0:3]=='d_T'):
            colu=2                                  
        if(k[0:3]=='U_0'):
            colu=3
        if(k[0:3]=='V_0'):
            colu=4
        if(k[0:3]=='W_0'):
            colu=5
        if(k[0:3]=='V_s'):
            colu=6            
        if(k[0:3]=='U_s'):
            colu=7
        if(k[4:]=='R_0'):
            lin=0
        if(k[4:]=='T_0'):
            lin=1            
        if(k[4:]=='d_T'):
            lin=2      
        if(k[4:]=='U_0'):
            lin=3
        if(k[4:]=='V_0'):
            lin=4
        if(k[4:]=='W_0'):
            lin=5
        if(k[4:]=='V_s'):
            lin=6            
        if(k[4:]=='U_s'):
            lin=7        
        z=(colu,lin)
        return(z)
    #------------------------------------------------------------------------------------------------
    #Pearson
    Order_Pearson=OrderedDict()
    for val in pearso_coeffs:
        Order_Pearson[val]=[] 
    for k in pearso_coeffs:
        z= find_corr(k)
        for i in dummy_7:
            Order_Pearson[k].append(round(final_coeefs[i][z[0]][z[1]],2))
    
    #Complete Sample
    Order_Pearson_com=OrderedDict()
    for val in pearso_coeffs:
        Order_Pearson_com[val]=[]
    for k in pearso_coeffs: 
        z= find_corr(k)
        Order_Pearson_com[k].append(final_coeefs['p_com'][z[0]][z[1]])
    #------------------------------------------------------------------------------
    #p-line
    Order_pline=OrderedDict()
    for val in pearso_coeffs:
        Order_pline[val]=[] 
    for k in pearso_coeffs:
        z= find_corr(k)
        for i in dummy_7:
            Order_pline[k].append(final_coeefs2[i][z[0]][z[1]])  
    
    Order_pline_com=OrderedDict()
    for val in pearso_coeffs:
        Order_pline_com[val]=[]
    for k in pearso_coeffs: 
        z= find_corr(k)
        Order_pline_com[k].append(final_coeefs2['p_com'][z[0]][z[1]])            
    
    #Order_Pearson['R_0-T_0'][0]=0.48
    #Order_Pearson['T_0-V_0'][0]=-0.73
    #Order_Pearson['T_0-V_s'][0]=0.62 
    #Order_Pearson['U_0-U_s'][0]=0.51
    #Order_Pearson['V_0-V_s'][0]=-0.68
    colork3=['g','c','r','b','m']
                    
    #Plot                                
    matplotlib.rcParams.update({'font.size': 26})        
    fig1 = plt.figure(figsize=(12,12))
    ax1 = fig1.add_subplot(111)
    Big_Pearsons=['R_0-T_0','T_0-V_0','T_0-V_s','U_0-U_s','V_0-V_s'] #Big ones Pearsons
    #Plotting the evolutions
    w=0
    #subtracting last number
    if(pconf=='yes'):
        N_bri_sources.remove(16)
    for k in Big_Pearsons:
        if k[-3:]=='V_s':
            tmp=np.arange(1,len(Order_Pearson[k]))
            for i in tmp:
                Order_Pearson[k][i]=((Order_Pearson[k][i])*1)
        ax1.scatter(N_bri_sources,Order_Pearson[k],color=colork3[w],s=60) #        
        label_num=pearso_coeffs.index(k)
        if(b_w_n=='b'):
            setname='SET A'
            ax1.plot(N_bri_sources,Order_Pearson[k],color=colork3[w], linewidth=3,label=pearso_coeffs_labels[label_num])# bessel
        if(b_w_n=='w'):
            setname='SET B'
            ax1.plot(N_bri_sources,Order_Pearson[k],color=colork3[w], linewidth=3)# whole
        w+=1
    #Plotting Reid Values
    w=0
    a=[]
    for k in Big_Pearsons:
        i,j=find_corr(k)
        ax1.plot(103, Reid_pearson[i][j],'*',color=colork3[w],markersize=13)    
        w+=1
    
    #ax1.plot(103, 10.990,'*',color='k',markersize=13,label='BeSSeL Sample')    
    ax1.grid(True)
    ax1.set_xlabel('Number of Sources')
    ax1.set_ylabel('Pearson Coefficient')
    if(b_w_n=='b'):    
        ax1.annotate(setname,(300,0),verticalalignment='center', horizontalalignment='center',fontweight='bold',fontsize=33)
    if(b_w_n=='w'):    
        ax1.annotate(setname,(550,0),verticalalignment='center', horizontalalignment='center',fontweight='bold',fontsize=33)
    ax1.set_ylim([-1.05,1.05])#1.49])   
    if(b_w_n=='w'):
        ax1.set_xlim([80,1010]) 
    if(b_w_n=='b'):
        ax1.set_xlim([80,510])
    #Plotting complete sample lines
    w=0    
    for i in Big_Pearsons:    
        if(b_w_n=='w'):
            if (i[-3:]=='V_s'):
                Order_Pearson_com[i][0]=((Order_Pearson_com[i][0])*-1)
            ax1.axhline(y=Order_Pearson_com[i],linewidth=3,color=colork3[w],ls='dashed') #whole
        w+=1
    #Label for complete samples
    if(b_w_n=='w'):
        ax1.axhline(y=-50,color='k',linewidth=3,ls='dashed',label='Complete Simulation')
        #ax1.plot(103, 5.0,'*',color='k',markersize=13,label='Current BeSSeL')
        ax1.legend(scatterpoints=1,bbox_to_anchor=(0.8,1.13),ncol=1,numpoints=1,prop={'size':27}) # whole
    if(b_w_n=='b'):
        ax1.plot(103, 5.0,'*',color='k',markersize=13,label='Model A5')
        #ax1.axhline(y=-50,color='k',linewidth=3,ls='dashed',label='Complete Simulation')
        ax1.legend(scatterpoints=1,bbox_to_anchor=(0.99,-0.088),ncol=3,numpoints=1,prop={'size':27})
        '''
        handles,labels = ax1.get_legend_handles_labels()
        handles = [handles[0], handles[1], handles[6], handles[3], handles[4], handles[2], handles[5]]
        labels = [labels[0], labels[1], labels[6],labels[3], labels[4], labels[2],labels[5]]
        ax1.legend(handles,labels,scatterpoints=1,bbox_to_anchor=(0.95, 0.13),ncol=3,numpoints=1,prop={'size':22})
        '''
        #ax1.legend() #bessel
    #Saving
    plt.savefig(path_outputs+'/Pearson_evo_'+name_ouputfolder+'.pdf',bbox_inches='tight')    
    #plt.savefig(path_outputs+'/Pearson_evol.pdf',bbox_inches='tight')  
    
##########################################################################################################################
###################################     PDFs for 1000 sources     ##########################################
##########################################################################################################################    
if (plots_confirmation=='yes'):    
    #Reid results:
    Th0V0_Reid=255.2
    Th0V0_sigma_Reid=5.1
    Th0V0R0_Reid=30.57
    Th0V0R0_sigma_Reid=0.43
    V0VS_Reid=17.1
    V0VS_sigma_Reid=1.0

    fig1 = plt.figure(figsize=(15,4))
    matplotlib.rcParams.update({'font.size': 16})
    zz=0 #which sample?
          
    ax1 = fig1.add_subplot(131)
    normed_Th0 = [x for x in dic[dic_table_names[zz][2:]]['Th0']]    
    normed_V0 = [x for x in dic[dic_table_names[zz][2:]]['V0']]        
    add1=[normed_Th0[x]+normed_V0[x] for x in range(len(normed_Th0))]        
    a,bins_data,c=ax1.hist(add1,bins=10,color=colork2[1],normed=True) #   
    (mu_Th0V0, sigma_Th0V0) = norm.fit(add1)
    area = sum(np.diff(bins_data)*a)
    xmin,xmax=ax1.get_xlim()
    x = np.linspace(xmin/2, xmax*2, 1000)
    y=mlab.normpdf(x, mu_Th0V0, sigma_Th0V0)*area        
    y1=mlab.normpdf(x, Th0V0_Reid, Th0V0_sigma_Reid)*area        
    ax1.plot(x,y1,color='grey',linewidth=2.5)    
    ax1.plot(x,y, linestyle='--',color='k',linewidth=2.5)
    ax1.fill_between(x, 0, y1, facecolor='grey', alpha=0.4)
    ax1.grid(True)
    ax1.set_xlabel('$\\Theta_0+V_{\\odot}$ (km/s)')
    ax1.set_ylabel('Normed Count Number')
    y_tickts=np.arange(0,0.20,0.04)
    ax1.set_yticks(y_tickts) 
    ax1.set_xticks([246,250,254,258,262])     
    ax1.set_xlim([xmin,264])  

    ax2 = fig1.add_subplot(132)
    normed_R0 = [x for x in dic[dic_table_names[zz][2:]]['R0']]    
    add2=[((normed_Th0[x]+normed_V0[x])/normed_R0[x]) for x in range(len(normed_Th0))]        
    a,bins_data,c=ax2.hist(add2,bins=10,color=colork2[2],normed=True) #       
    (mu_Th0V0R0, sigma_Th0V0R0) = norm.fit(add2)
    area = sum(np.diff(bins_data)*a)
    xmin,xmax=ax2.get_xlim()
    x = np.linspace(xmin/2, xmax*2, 1000)
    y=mlab.normpdf(x, mu_Th0V0R0, sigma_Th0V0R0)*area        
    y1=mlab.normpdf(x, Th0V0R0_Reid, Th0V0R0_sigma_Reid)*area        
    ax2.plot(x,y1,color='grey',linewidth=2.5,label='Results Model A5')    
    ax2.plot(x,y, linestyle='--',color='k',linewidth=2.5,label='Gaussian Fitting')
    ax2.fill_between(x, 0, y1, facecolor='grey', alpha=0.4)    
    #plt.setp(ax2.get_yticklabels(), visible=False) 
    ax2.grid(True)
    ax2.set_xlabel('$(\\Theta_0+V_{\\odot})/R_0$ (km/s$\\cdot$kpc)')
    x_tickts=np.arange(30.0,32,0.3)
    y_tickts=np.arange(0,2.0,0.4)
    ax2.set_xticks(x_tickts)     
    ax2.set_yticks(y_tickts) 
    ax2.set_xlim([29.9,31.3]) 
    ax2.legend(bbox_to_anchor=(1.3, 1.23),numpoints=1,ncol=2,prop={'size':15})
                                                                                                        
    ax3 = fig1.add_subplot(133)
    normed_VS = [x for x in dic[dic_table_names[zz][2:]]['VS']]    
    add3=[normed_V0[x]-normed_VS[x] for x in range(len(normed_Th0))]        
    a,bins_data,c=ax3.hist(add3,bins=8,color=colork2[3],normed=True) #        
    (mu_V0VS, sigma_V0VS) = norm.fit(add3)
    area = sum(np.diff(bins_data)*a)
    xmin,xmax=ax3.get_xlim()
    x = np.linspace(xmin/2, xmax*2, 1000)
    y=mlab.normpdf(x, mu_V0VS, sigma_V0VS)*area        
    y1=mlab.normpdf(x, V0VS_Reid, V0VS_sigma_Reid)*area        
    ax3.plot(x,y1,color='grey',linewidth=2.5)    
    ax3.plot(x,y, linestyle='--',color='k',linewidth=2.5)
    ax3.fill_between(x, 0, y1, facecolor='grey', alpha=0.4)    
    #plt.setp(ax3.get_yticklabels(), visible=False) 
    ax3.grid(True)
    ax3.set_xlabel('$V_{\\odot}+\\mid \\bar{V_s} \\mid$ (km/s)')   
    y_tickts=np.arange(0,0.5,0.1)
    ax3.set_yticks(y_tickts) 
    ax3.set_xticks([8,12,16,20,24]) 
    ax3.set_xlim([7.9,25.4])
    ax3.set_ylim([0,0.42])         
    #plt.savefig(path_outputs+'/2PDFs.pdf',bbox_inches='tight')   
    
    plt.savefig(path_outputs+'/PDFs_'+name_ouputfolder+'2.pdf',bbox_inches='tight')    
        
    if (latex_tables=='yes'):
        print('Th0+V0: %s +/- %s' %(np.round(mu_Th0V0,1),np.round(sigma_Th0V0,1)))
        print('(Th0+V0)/R0: %s +/- %s' %(np.round(mu_Th0V0R0,2),np.round(sigma_Th0V0R0,2)))
        print('V0+VS: %s +/- %s' %(np.round(mu_V0VS,1),np.round(sigma_V0VS,1)))
        print('--------------------------------------------------')
        print('Reid Results')    
        print('Th0+V0: %s +/- %s' %(Th0V0_Reid,Th0V0_sigma_Reid))
        print('(Th0+V0)/R0: %s +/- %s' %(Th0V0R0_Reid,Th0V0R0_sigma_Reid))
        print('V0+VS: %s +/- %s' %(V0VS_Reid,V0VS_sigma_Reid))    

##########################################################################################################################
###################################     Errors for paramters     ##########################################
##########################################################################################################################    
if (plots_confirmation_nsources=='yes'):    
    if(b_w_n=='b'):
        x_tickts=[0,100,200,300,400,500]
        #y_tickts=[0.4,0.7,1.0,1.3,1.6]
        xpos_label=250
        #ylim=[0.5,1.5]
        xlim=[0,550]
        #ylim_4=[-3.5,5.5]        
        #ylim_8=[-7,8] 
        #y_tickts_4=[-3.0,-1.0,1.0,3.0,5.0]      
        #y_tickts_8=[7,4,1,-2,-5]
        N_BeSSeL=110
        N_BeSSeL_past=30
    elif(b_w_n=='w'):
        x_tickts=[0,250,500,750,1000]
        #y_tickts=[-1.8,-1.1,-0.4,0.3,1.0,1.7]
        xpos_label=550
        #ylim=[0.25,1.75]
        xlim=[0,1050]
        #ylim_4=[-43,7]        
        #ylim_8=[-150,11]
        #y_tickts_4=[1.0,-10,-20,-30,-40]
        #y_tickts_8=[-120,-80,-40,1.0]
        N_BeSSeL=120
        N_BeSSeL_past=30
    else:
        x_tickts=[]
        #y_tickts=[]
        xpos_label=440
        ylim=[]
        xlim=[]
        ylim_4=[]        
        ylim_8=[]
        y_tickts_4=[]
        y_tickts_8=[]
    #Plot
    #N_bri_sources=[103, 150, 200, 300, 400, 800]
    matplotlib.rcParams.update({'font.size': 16})
    fig1 = plt.figure(figsize=(15,10))
    ax = fig1.add_subplot(111)        
    #ax.axvline(0.755,color='k',linestyle='dashed',linewidth=3)
    X=np.linspace(0,1300,100)
    if(pconf=='yes'):
        N_bri_sources.append(16)
                                
    #R0
    ax1 = fig1.add_subplot(241)
    #ax1.set_yscale('log')
    ax1.plot(N_bri_sources,normed_errR0,'o-',color=colork2[1],linewidth=2,markersize=10)
    if (errors_conf=='yes' and b_w_n=='w'):
        normed_errR0_north=[0.072218509400291592, 0.05732922959913779, 0.052195876753302173, 0.047979124402721098, 0.044183074782344248, 0.044067533203469511]#, 0.24]
        ax1.plot(N_bri_sources,normed_errR0_north,'*-',color=colork2[1],linewidth=2,markersize=12)
    ax1.set_xticks(x_tickts)
    ax1.grid(True)
    plt.setp(ax1.get_xticklabels(), visible=False) 
    dum=ax1.get_ylim()
    ax1.set_ylim(0.0,dum[1])    
    dum2=(dum[1]/5.0)
    ax1.text( xpos_label, dum2, '$R_0$ $\\rm{(kpc)}$', verticalalignment='center', horizontalalignment='center', color=colork2[1], fontsize=27)                                    
    if(b_w_n=='w'):
        ax1.axhline(errors_complete['R0'],color=colork2[1],lw=2,ls='dashed')

    #Th0
    ax2 = fig1.add_subplot(242)
    ax2.plot(N_bri_sources,normed_errTh0,'o-',color=colork2[2],linewidth=2,markersize=10)
    if (errors_conf=='yes' and b_w_n=='w'):
        normed_errTh0_north=[3.2287550215369403, 3.0113503521611453, 2.958463158288974, 2.9453714996285103, 2.6429657056506919, 2.6377043249375531]#, 9.11]
        ax2.plot(N_bri_sources,normed_errTh0_north,'*-',color=colork2[2],linewidth=2,markersize=12)
    ax2.set_xticks(x_tickts)   
    ax2.grid(True)
    plt.setp(ax2.get_xticklabels(), visible=False) 
    dum=ax2.get_ylim()
    ax2.set_ylim(0.0,dum[1])    
    dum2=(dum[1]/5.0)
    ax2.text( xpos_label, dum2, '$\\Theta_0$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[2], fontsize=27)                         
    if(b_w_n=='w'):
        ax2.axhline(errors_complete['Th0'],color=colork2[2],lw=2,ls='dashed')

    #dTd    
    ax3 = fig1.add_subplot(243)
    normed_dTd = [x for x in means['dTd']]
    normed_errdTd = [x for x in errors['dTd']]
    ax3.plot(N_bri_sources,normed_errdTd,'o-',color=colork2[8],linewidth=2,markersize=10)    
    if (errors_conf=='yes' and b_w_n=='w'):
        normed_errdTd_north=[1.0148304395809185, 0.86480930126301281, 0.76527552430426893, 0.65843141112440406, 0.6115427010932406, 0.60967848345673492]#, 3.24]
        ax3.plot(N_bri_sources,normed_errdTd_north,'*-',color=colork2[8],linewidth=2,markersize=12)    
    ax3.grid(True)
    plt.setp(ax3.get_xticklabels(), visible=False)   
    dum=ax3.get_ylim()
    ax3.set_ylim(0.0,dum[1])      
    dum2=(dum[1]/5.0)
    if(b_w_n=='b'):
        ax3.text( xpos_label+50, dum2, '$\\frac{d\\Theta}{dR}$ $\\rm{(km/s \\cdot kpc)}$', verticalalignment='center', horizontalalignment='center', color=colork2[8], fontsize=27)                             
    if(b_w_n=='w'):
        ax3.text( xpos_label-50, dum2, '$\\frac{d\\Theta}{dR}$ $\\rm{(km/s \\cdot kpc)}$', verticalalignment='center', horizontalalignment='center', color=colork2[8], fontsize=27)                     
        ax3.axhline(errors_complete['dTd'],color=colork2[8],lw=2,ls='dashed')
        
    #U_sun
    ax4 = fig1.add_subplot(244)#,sharey=ax1)
    ax4.plot(N_bri_sources,normed_errUsun,'o-',color=colork2[5],linewidth=2,markersize=10)    
    if (errors_conf=='yes' and b_w_n=='w'):
        normed_errUsun_north=[1.0002256950808652, 0.87479978336451025, 0.81300695479186702, 0.72258451051087447, 0.76294779855241224, 0.76432012272345673]#, 1.58]
        ax4.plot(N_bri_sources,normed_errUsun_north,'*-',color=colork2[5],linewidth=2,markersize=12)    
    ax4.grid(True)  
    plt.setp(ax4.get_xticklabels(), visible=False) 
    dum=ax4.get_ylim()
    ax4.set_ylim(0.0,dum[1])    
    dum2=(dum[1]/5.0)
    ax4.text( xpos_label, dum2, '$U_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[5], fontsize=27)                         
    if(b_w_n=='w'):
        ax4.axhline(errors_complete['U0'],color=colork2[5],lw=2,ls='dashed')
        
    #V_sun
    ax5 = fig1.add_subplot(245,sharex=ax1)
    normed_Vsun = [x for x in means['V0']]
    normed_errVsun = [x  for x in errors['V0']]                
    #ax5.set_yscale('log',nonposy='clip')
    ax5.plot(N_bri_sources,normed_errVsun,'o-',color=colork2[6],linewidth=2,markersize=10)    
    if (errors_conf=='yes' and b_w_n=='w'):
        normed_errVsun_north=[2.1561942587809657, 2.1009328418629489, 2.3608767546632148, 2.1742338090736211, 1.9559574015459351, 1.9546563560834944]#, 2.54]
        ax5.plot(N_bri_sources,normed_errVsun_north,'*-',color=colork2[6],linewidth=2,markersize=12)    
    ax5.set_xlim(xlim)
    ax5.grid(True)
    dum=ax5.get_ylim()
    ax5.set_ylim(0.0,dum[1])    
    dum2=(dum[1]/5.0)
    ax5.text( xpos_label,dum2, '$V_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[6], fontsize=27)                      
    if(b_w_n=='w'):
        ax5.axhline(errors_complete['V0'],color=colork2[6],lw=2,ls='dashed')
                
    #W_sun
    ax6 = fig1.add_subplot(246,sharex=ax2)    
    normed_Wsun = [x for x in means['W0']]
    normed_errWsun = [x  for x in errors['W0']]
    ax6.plot(N_bri_sources,normed_errWsun,'o-',color=colork2[7],linewidth=2,markersize=10,label='SET B')              
    if (errors_conf=='yes' and b_w_n=='w'):
        normed_errWsun_north=[0.41142068931933895, 0.41026141999505611, 0.36242484291021498, 0.3001786130806432, 0.27501590369890672, 0.27237873751323821]#, 0.77]
        ax6.plot(N_bri_sources,normed_errWsun_north,'*-',color=colork2[7],linewidth=2,markersize=12,label='SET A')              
    ax6.set_xlim(xlim)
    ax6.grid(True)      
    dum=ax6.get_ylim()
    ax6.set_ylim(0.0,dum[1])
    dum2=(dum[1]/5.0)
    ax6.text( xpos_label,dum2, '$W_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[7], fontsize=27)               
    if(b_w_n=='w'):
        ax6.axhline(errors_complete['W0'],color=colork2[7],lw=2,ls='dashed',label='Complete Simulation')
        ax6.legend(scatterpoints=1,bbox_to_anchor=(3.7,-0.13),ncol=3,numpoints=1,prop={'size':23})
        
    #U_s
    ax7 = fig1.add_subplot(247,sharex=ax3)
    ax7.plot(N_bri_sources,normed_errUS,'o-',color=colork2[3],linewidth=2,markersize=10)        
    if (errors_conf=='yes' and b_w_n=='w'):
        normed_errUS_north=[1.5520671855303172, 1.2665994203490258, 1.071265916098133, 0.97563063307099451, 0.90799586405958121, 0.90375556682433345]#, 3.66]        
        ax7.plot(N_bri_sources,normed_errUS_north,'*-',color=colork2[3],linewidth=2,markersize=12)        
    ax7.set_xticks(x_tickts) 
    ax7.set_xlim(xlim)          
    ax7.grid(True)
    dum=ax7.get_ylim()
    ax7.set_ylim(0.0,dum[1])    
    dum2=(dum[1]/5.0)
    ax7.text( xpos_label,dum2, '$\\bar{U}_s$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[3], fontsize=27)               
    if(b_w_n=='w'):
        ax7.axhline(errors_complete['US'],color=colork2[3],lw=2,ls='dashed')
        
    #Vs
    ax8 = fig1.add_subplot(248,sharex=ax4)#,sharey=ax5)
    normed_VS = [-x for x in means['VS']]
    normed_errVS = [x for x in errors['VS']]
    ax8.plot(N_bri_sources,normed_errVS,'o-',color=colork2[4],linewidth=2,markersize=10)        
    if (errors_conf=='yes' and b_w_n=='w'):
        normed_errVS_north=[3.5908735490963757, 3.1656443455399379, 3.1349308631910948, 2.8169148744930745, 2.3403510256320703, 2.3678431327801244]#, 9.78]
        ax8.plot(N_bri_sources,normed_errVS_north,'*-',color=colork2[4],linewidth=2,markersize=12)        
    ax8.set_xlim(xlim)
    ax8.set_xticks(x_tickts)       
    ax8.grid(True)
    dum=ax8.get_ylim()
    ax8.set_ylim(0.0,dum[1])
    dum2=(dum[1]/5.0)
    ax8.text(xpos_label,dum2, '$\\bar{V}_s$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[4], fontsize=27)                       
    if(b_w_n=='w'):
        ax8.axhline(errors_complete['VS'],color=colork2[4],lw=2,ls='dashed')
        
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.setp(ax.get_yticklabels(), visible=False)  
    plt.setp(ax.get_xticklabels(), visible=False)          
    ax.set_xlabel('Number of Sources',labelpad=25)
    #ax.set_ylabel('Normed Parameter',labelpad=40)
    ax.tick_params(
    axis='both',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off',right='off', left='off', labelleft='off')
    #fig1.patch.set_visible(False)
    fig1.tight_layout()     
    plt.savefig(path_outputs+'/err_'+name_ouputfolder+'.pdf',bbox_inches='tight')
os.chdir(path_origin)
print('Done with simulation, fitting and analysis')
print('Check the folder %s for the results' % ('output_'+name_ouputfolder))
print('Thanks for using me!, please remember to reference: M. Reid (CfA), H. van Langevelde (JIVE) and L.H. Quiroga-Nunez (Leiden)')

