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
plots_confirmation = raw_input("Do you want plots of results distributions and pearson coefficents? (yes or no):")
plots_function_confirmation= raw_input("Do you want plots of pearson coefficents comparison? (yes or no):")
clue= raw_input("Why there are different samples? (nsources,unc_pi_muxy,random,declination,fudge,all,none):")
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
    else:
        prtfiles.append(dummy)
os.chdir('../../galaxy_digest/')
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

#message_note="""     
#    <p> <b> Important Note </b>: These specifications DO NOT include the current stage 
#        of BeSSeL porject: 100 sources in -30 < dec <70, which were included in ALL simulations.
#    </p>                  
#    """

message_note="""     
    <p> <b> Important Note </b>.
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
Reid_pearson=[[ 1.000, 0.465, 0.103, 0.452, 0.023,-0.003, 0.002, 0.517],
              [ 0.465, 1.000, 0.136, 0.243,-0.796,-0.009, 0.809, 0.171],
              [ 0.103, 0.136, 1.000,-0.124,-0.009, 0.025, 0.018,-0.094],
              [ 0.452, 0.243,-0.124, 1.000,-0.014,-0.017,-0.025, 0.839],
              [ 0.023,-0.796,-0.009,-0.014, 1.000, 0.011,-0.990,-0.006],
              [-0.003,-0.009, 0.025,-0.017, 0.011, 1.000,-0.010,-0.002],
              [ 0.002, 0.809, 0.018,-0.025,-0.990,-0.010, 1.000,-0.028],
              [ 0.517, 0.171,-0.094, 0.839,-0.006,-0.002,-0.028, 1.000]]
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
fhu.write(message_plots)
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
    num_sources_per_galaxy['p_bri_'+str(i)]=N_bri_sources[i]+100
    N_bri_sources[i]=N_bri_sources[i]+100
dummy_7=num_sources_per_galaxy.keys()

if (plots_function_confirmation=='yes'):
    #------------------------------------------------------------------------------------------------
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
                                                                                                                                                                  '$W_{\\odot}$-$\\bar{V}_{s}$', '$W_{\\odot}$-$\\bar{U}_{s}$',
                                                                                                                                                                                                 '$\\bar{V}_{s}$-$\\bar{U}_{s}$']
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
    #------------------------------------------------------------------------------------------------
    colork2=['g','k','r','b','y','m','c','orange','pink','sienna','peru','dimgrey',
               'peachpuff','darkgoldenrod','firebrick','silver','rosybrown','crimson',
               'cadetblue','tomato','chocolate']
    '''
    k=np.arange(len(pearso_coeffs))
    k1,k2,k3=np.split(k, 3)
    matplotlib.rcParams.update({'font.size': 22})        
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    for k in k1:
        ax1.scatter(N_bri_sources,Order_Pearson[pearso_coeffs[k]],color=colork2[k],s=10) #        
        ax1.plot(N_bri_sources,Order_Pearson[pearso_coeffs[k]],color=colork2[k], label=pearso_coeffs_labels[k]) #        
    ax1.grid(True)
    ax1.legend(scatterpoints=1,loc=8,ncol=3)
    ax1.set_xlabel('Number of Sources')
    ax1.set_ylabel('Pearson Coefficient')
    ax1.set_ylim([-1.0,1.0])    
    
    matplotlib.rcParams.update({'font.size': 22})        
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    for k in k2:
        ax1.scatter(N_bri_sources,Order_Pearson[pearso_coeffs[k]],color=colork2[k],s=10) #        
        ax1.plot(N_bri_sources,Order_Pearson[pearso_coeffs[k]],color=colork2[k], label=pearso_coeffs_labels[k]) #        
    ax1.grid(True)
    ax1.legend(scatterpoints=1,loc=8,ncol=3)
    ax1.set_xlabel('Number of Sources')
    ax1.set_ylabel('Pearson Coefficient')
    ax1.set_ylim([-1.0,1.0])    
    
    matplotlib.rcParams.update({'font.size': 22})        
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    for k in k3:
        ax1.scatter(N_bri_sources,Order_Pearson[pearso_coeffs[k]],color=colork2[k],s=10) #        
        ax1.plot(N_bri_sources,Order_Pearson[pearso_coeffs[k]],color=colork2[k], label=pearso_coeffs_labels[k]) #        
    ax1.grid(True)
    ax1.legend(scatterpoints=1,loc=8,ncol=3)
    ax1.set_xlabel('Number of Sources')
    ax1.set_ylabel('Pearson Coefficient')
    ax1.set_ylim([-1.0,1.0])    
    
    matplotlib.rcParams.update({'font.size': 22})        
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    for k in k3:
        ax1.scatter(N_bri_sources,Order_Pearson[pearso_coeffs[k]],color=colork2[k],s=10) #        
        ax1.plot(N_bri_sources,Order_Pearson[pearso_coeffs[k]],color=colork2[k], label=pearso_coeffs_labels[k]) #        
    ax1.grid(True)
    ax1.legend(scatterpoints=1,loc=8,ncol=3)
    ax1.set_xlabel('Number of Sources')
    ax1.set_ylabel('Pearson Coefficient')
    ax1.set_ylim([-1.0,1.0])        
    '''
    #------------------------------------------------------------------------------------------------
    if(clue=='fudge'):
        colork3=['g','c','r','b','m']
        matplotlib.rcParams.update({'font.size': 22})        
        fig1 = plt.figure(figsize=(20,8))
        ax1 = fig1.add_subplot(111)
        k4=[0,7,9,14,16]
        w=0
        for k in k4:
            ax1.scatter(fudge_fact,Order_Pearson[pearso_coeffs[k]],color=colork3[w],s=50) #        
            ax1.plot(fudge_fact,Order_Pearson[pearso_coeffs[k]],color=colork3[w], linewidth=2,label=pearso_coeffs_labels[k]) #        
            w+=1
        ax1.set_xlim([0.003,1.3]) 
        ax1.grid(True)
        ax1.invert_xaxis() 
        ax1.set_xscale('log')
        ax1.set_xlabel('Fudge Factor')
        ax1.set_ylabel('Pearson Coefficient')
        ax1.set_ylim([-1.5,1.01])   
        ax1.plot(0.1, 0.465,'*',color=colork3[0],markersize=20)
        ax1.plot(0.1,-0.796,'*',color=colork3[1],markersize=20)
        ax1.plot(0.1, 0.809,'*',color=colork3[2],markersize=20)
        ax1.plot(0.1, 0.839,'*',color=colork3[3],markersize=20)
        ax1.plot(0.1,-0.990,'*',color=colork3[4],markersize=20)
        ax1.plot(0.1, 10.990,'*',color='k',markersize=20,label='Reid, 2014')    
        w=0
        for i in k4:    
            ax1.axhline(y=Order_Pearson[pearso_coeffs[i]][-1],linewidth=2,color=colork3[w],ls='dashed')
            w+=1
        ax1.axhline(y=-50,color='k',linewidth=2,ls='dashed',label='$f$=0')   
        ax1.legend(scatterpoints=1,loc=4,ncol=2,numpoints=1)
        plt.savefig(path_outputs+'/Pearson_fudge.png',bbox_inches='tight')     
        #------------------------------------------------------------------------------------------------
        colork3=['g','c','r','b','m']
        matplotlib.rcParams.update({'font.size': 22})        
        fig1 = plt.figure(figsize=(20,8))
        ax1 = fig1.add_subplot(111)
        k4=[0,7,9,14,16]
        w=0
        for k in k4:
            ax1.scatter(fudge_fact,Order_pline[pearso_coeffs[k]],color=colork3[w],s=50) #        
            widt=[6,2,3.5,2,2]
            ax1.plot(fudge_fact,Order_pline[pearso_coeffs[k]],color=colork3[w], linewidth=widt[w],label=pearso_coeffs_labels[k]) #        
            w+=1
        ax1.grid(True)
        ax1.axhline(y=0.05,color='k',linewidth=2,ls='dashed',label='5% Prob.')
        ax1.set_xscale('log')
        ax1.set_xlim([0.003,1.3]) 
        ax1.set_ylim([-0.1,1.3])
        ax1.invert_xaxis() 
        ax1.set_xlabel('Fudge Factor')
        ax1.set_ylabel('p-value')
        ax1.legend(scatterpoints=1,loc=9,ncol=3,numpoints=1)     
        plt.savefig(path_outputs+'/pline_fudge.png',bbox_inches='tight')     
        #------------------------------------------------------------------------------------------------
    if(clue=='nsources'):
        colork3=['g','c','r','b','m']
        matplotlib.rcParams.update({'font.size': 22})        
        fig1 = plt.figure(figsize=(20,8))
        ax1 = fig1.add_subplot(111)
        k4=[0,7,9,14,16]
        w=0
        for k in k4:
            ax1.scatter(N_bri_sources,Order_Pearson[pearso_coeffs[k]],color=colork3[w],s=50) #        
            ax1.plot(N_bri_sources,Order_Pearson[pearso_coeffs[k]],color=colork3[w], linewidth=2,label=pearso_coeffs_labels[k]) #        
            w+=1
        ax1.plot(103, 0.465,'*',color=colork3[0],markersize=10)
        ax1.plot(103,-0.796,'*',color=colork3[1],markersize=10)
        ax1.plot(103, 0.809,'*',color=colork3[2],markersize=10)
        ax1.plot(103, 0.839,'*',color=colork3[3],markersize=10)
        ax1.plot(103,-0.990,'*',color=colork3[4],markersize=10)
        ax1.plot(103, 10.990,'*',color='k',label='Reid et al. 2014')    
        ax1.grid(True)
        ax1.set_xlabel('Number of Sources')
        ax1.set_ylabel('Pearson Coefficient')
        ax1.set_ylim([-1.05,1.05])   
        ax1.set_xlim([80,510]) 
        w=0    
        for i in k4:    
            ax1.axhline(y=Order_Pearson_com[pearso_coeffs[i]],linewidth=2,color=colork3[w],ls='dashed')
            w+=1
        ax1.axhline(y=-50,color='k',linewidth=2,ls='dashed',label='Complete Sample')
        ax1.legend(scatterpoints=1,loc=8,ncol=3,numpoints=1,prop={'size':22})
        #plt.savefig(path_outputs+'/Pearson_nsources.png',bbox_inches='tight')     
        #------------------------------------------------------------------------------------------------
        colork3=['g','c','r','b','m']
        matplotlib.rcParams.update({'font.size': 22})        
        fig1 = plt.figure(figsize=(20,8))
        ax1 = fig1.add_subplot(111)
        k4=[0,7,9,14,16]
        w=0
        for k in k4:
            ax1.scatter(N_bri_sources,Order_pline[pearso_coeffs[k]],color=colork3[w],s=50) #        
            ax1.plot(N_bri_sources,Order_pline[pearso_coeffs[k]],color=colork3[w], linewidth=2,label=pearso_coeffs_labels[k]) #        
            w+=1
        ax1.grid(True)
        ax1.axhline(y=0.05,color='k',linewidth=2,ls='dashed',label='5% Probability')
        ax1.set_xlabel('Number of Sources')
        ax1.set_ylabel('p line')
        ax1.set_ylim([-0.01,0.6])   
        ax1.legend(scatterpoints=1,loc=0,ncol=3,numpoints=1)     
        #plt.savefig(path_outputs+'/pline_nsources.png',bbox_inches='tight')       
    
    
    keys_cov=cov_cov.keys()
    cov_mean=OrderedDict()
    if(initial_values['sam_bri_dec']=='1'):
        for k in range(n_bri):
            cov_mean['bri_'+str(k)]=OrderedDict()    
    if(initial_values['samp_rand']=='1'):    
        for k in range(n_ran):
            cov_mean['ran_'+str(k)]=OrderedDict()   
    cov_mean['com']=OrderedDict()           
    ssamples=cov_mean.keys()
    for sample in ssamples:
        for ipar in fitpars:
            for jpar in fitpars:
                cov_mean[sample][ipar,jpar]=[]
    #------------------------------------------------------------------------------------------------
    cov_std=OrderedDict()
    if(initial_values['sam_bri_dec']=='1'):
        for k in range(n_bri):
            cov_std['bri_'+str(k)]=OrderedDict()    
    if(initial_values['samp_rand']=='1'):    
        for k in range(n_ran):
            cov_std['ran_'+str(k)]=OrderedDict()   
    cov_std['com']=OrderedDict()           
    ssamples=cov_std.keys()
    for sample in ssamples:
        for ipar in fitpars:
            for jpar in fitpars:
                cov_std[sample][ipar,jpar]=[]
                
    for ipar in fitpars:
        for jpar in fitpars:                    
            for ll in keys_cov:
                cov_mean[ll][ipar,jpar]=mean(cov_cov[ll][ipar,jpar])
                cov_std[ll][ipar,jpar]=std(cov_cov[ll][ipar,jpar])
    #------------------------------------------------------------------------------------------------
    del keys_cov[-1]
    cov_mean_all=OrderedDict()
    for ipar in fitpars:
        for jpar in fitpars:
            cov_mean_all[ipar+'-'+jpar]=[]
    
    cov_std_all=OrderedDict()
    for ipar in fitpars:
        for jpar in fitpars:
            cov_std_all[ipar+'-'+jpar]=[]
                
    
    for ll in keys_cov:
        for ipar in fitpars:
            for jpar in fitpars:
                cov_mean_all[ipar+'-'+jpar].append(cov_mean[ll][ipar,jpar])
                cov_std_all[ipar+'-'+jpar].append(cov_std[ll][ipar,jpar])    
#------------------------------------------------------------------------------------------------
    #For the complete sample        
    cov_mean_com=OrderedDict()
    for ipar in fitpars:
        for jpar in fitpars:
            cov_mean_com[ipar+'-'+jpar]=[]

    cov_std_com=OrderedDict()
    for ipar in fitpars:
        for jpar in fitpars:
            cov_std_com[ipar+'-'+jpar]=[]
                                    
    for ipar in fitpars:
        for jpar in fitpars:
            cov_mean_com[ipar+'-'+jpar].append(cov_mean['com'][ipar,jpar])
            cov_std_com[ipar+'-'+jpar].append(cov_std['com'][ipar,jpar])        
            
    #------------------------------------------------------------------------------------------------            
    if(clue=='nsources'):
        colork3=['g','c','r','b','m']
        matplotlib.rcParams.update({'font.size': 22})        
        fig1 = plt.figure(figsize=(20,8))
        ax1 = fig1.add_subplot(111)
        k5=['R0-Th0','Th0-V0','Th0-VS','V0-VS','U0-US']
        w=0
        k4=[0,7,9,14,16]
        for k in range(len(k4)):
            ax1.errorbar(N_bri_sources,cov_mean_all[k5[k]],yerr=cov_std_all[k5[k]],color=colork3[w])
            #scatter(N_bri_sources,cov_mean_all[k5[k]],color=colork3[w],s=50) #        
            ax1.plot(N_bri_sources,cov_mean_all[k5[k]],color=colork3[w], linewidth=2,label=pearso_coeffs_labels[k4[k]]) #        
            w+=1
        ax1.plot(103, 0.465,'*',color=colork3[0],markersize=10)
        ax1.plot(103,-0.796,'*',color=colork3[1],markersize=10)
        ax1.plot(103, 0.809,'*',color=colork3[2],markersize=10)
        ax1.plot(103,-0.839,'*',color=colork3[3],markersize=10)
        ax1.plot(103,-0.990,'*',color=colork3[4],markersize=10)
        ax1.plot(103, 10.990,'*',color='k',label='Reid et al. 2014')    
        ax1.grid(True)
        '''
        low_cov_com=[]
        high_cov_com=[]
        for i in k5:    
            low_cov_com.append(cov_mean_com[i][0]-cov_std_com[i][0])
            high_cov_com.append(cov_mean_com[i][0]+cov_std_com[i][0])    
        aa=numpy.linspace(0,N_bri_sources[-1],len(k5))
        '''
        w=0
        for i in k5:    
            ax1.axhline(y=cov_mean_com[i],linewidth=2,color=colork3[w],ls='dashed')
            #ax1.fill_between(aa, low_cov_com, high_cov_com)          
            w+=1
        ax1.axhline(y=-50,color='k',linewidth=2,ls='dashed',label='Complete Sample')
        ax1.legend(scatterpoints=1,loc=8,ncol=3,numpoints=1)
        ax1.set_xlabel('Number of Sources')
        ax1.set_ylabel('Pearson Coefficient')
        ax1.set_ylim([-2.3,1.0])    
        plt.savefig(path_outputs+'/meanpears_reid_nsources.png',bbox_inches='tight')     
    
    if(float(initial_values['Unc_confirmation_spec'])==2):
        #------------------------------------------------------------------------------------------------        
        colork3=['g','c','r','b','m']
        matplotlib.rcParams.update({'font.size': 22})        
        fig1 = plt.figure(figsize=(20,8))
        ax1 = fig1.add_subplot(111)
        k4=[0,7,9,14,16]
        w=0
        for k in k4:
            ax1.scatter(err_muxy[1:],Order_Pearson[pearso_coeffs[k]][1:],color=colork3[w],s=50) #        
            ax1.plot(err_muxy[1:],Order_Pearson[pearso_coeffs[k]][1:],color=colork3[w], linewidth=2,label=pearso_coeffs_labels[k]) #        
            w+=1
        ax1.grid(True)
        ax1.set_xlabel('Uncertinty in Parallax and proper motion (mas/[yr])')
        ax1.set_ylabel('Pearson Coefficient')
        ax1.set_ylim([-1.5,1.0])   
        ax1.set_xlim([0.03,11.0])
        ax1.invert_xaxis() 
        ax1.set_xscale('log')
        x=[0.07,3.3]  
        ax1.fill_between(x, -10, 10, facecolor='yellow', alpha=0.2)
        ax1.plot([], [], color='yellow',  alpha=0.3, linewidth=20,label='Range $\\Delta \\mu_{x,y}$ BeSSeL')
        ax1.axvline(x=0.6,color='k',linewidth=2,ls='dashed',label=' $\\overline{\\Delta \\mu_{x,y}}$ BeSSeL')
        #ax1.axvline(x=0.03,color='k',linewidth=2,ls='dashed',label='$\\bar{\\pi}$ in BeSSeL')
        ax1.legend(scatterpoints=1,loc=8,ncol=3,numpoints=1)        
        plt.savefig(path_outputs+'/Pearson_erros_muxy.png',bbox_inches='tight')     
        #------------------------------------------------------------------------------------------------
        colork3=['g','r','c','b','m']
        matplotlib.rcParams.update({'font.size': 22})        
        fig1 = plt.figure(figsize=(20,8))
        ax1 = fig1.add_subplot(111)
        k4=[0,7,9,14,16]
        w=0
        for k in k4:
            ax1.scatter(err_muxy[1:],Order_pline[pearso_coeffs[k]][1:],color=colork3[w],s=50) #        
            ax1.plot(err_muxy[1:],Order_pline[pearso_coeffs[k]][1:],color=colork3[w], linewidth=2,label=pearso_coeffs_labels[k]) #        
            w+=1
        ax1.grid(True)
        ax1.axhline(y=0.05,color='k',linewidth=2,ls='dashed',label='5% Prob.')
        ax1.set_xlabel('Uncertinty in Parallax and proper motion (mas/[yr])')
        ax1.set_ylabel('p value')
        ax1.set_ylim([-0.01,0.8])   
        ax1.set_xlim([0.03,11])   
        ax1.invert_xaxis() 
        ax1.set_xscale('log')
        ax1.legend(scatterpoints=1,loc=0,ncol=2,numpoints=1)     
        plt.savefig(path_outputs+'/pline_erros_muxy.png',bbox_inches='tight')             
        
        #------------------------------------------------------------------------------------------------------------------------
        means={}
        for keypar in fitpars:
            means[keypar]=[]
        for i in range(len(dic_table_names)-1):#putting out teh last = complete sample
            for keypar in fitpars:
                if keypar != 'void':
                    means[keypar].append(np.mean(dic[dic_table_names[i][2:]][keypar]))
        errors={}
        for keypar in fitpars:
            errors[keypar]=[]
        for i in range(len(dic_table_names)-1):#putting out teh last = complete sample
            for keypar in fitpars:
                if keypar != 'void':
                    errors[keypar].append(np.std(dic[dic_table_names[i][2:]][keypar]))        
        
        jj=3 # 0=R0,1=Th0, 3=sun, 6=uvs
        colork3=['g','m','c','b','r','k','y']
        matplotlib.rcParams.update({'font.size': 22})        
        fig1 = plt.figure(figsize=(20,8))
        ax1 = fig1.add_subplot(111)
        #for keypar in fitpars:
        if(fitpars[jj]=='R0'):
            ax1.errorbar(N_bri_sources,means[fitpars[jj]],yerr=errors[fitpars[jj]],color=colork3[0],linewidth=2,fmt="^",label='$R_0$ & $\\Delta R_0$',capsize=4, elinewidth=4,markersize=8)
        if(fitpars[jj]=='Th0'):
            ax1.errorbar(N_bri_sources,means[fitpars[jj]],yerr=errors[fitpars[jj]],color=colork3[0],linewidth=2,fmt="^",label='$\\Theta_0$ & $\\Delta \\Theta_0$',capsize=4, elinewidth=4,markersize=8)
        if(fitpars[jj]=='U0'):
            zz=[3,4,5]
            label_sun=['$U_{\\odot}$','$V_{\\odot}$','$W_{\\odot}$']
            w=0
            for kk in zz:                    
                ax1.errorbar(N_bri_sources,means[fitpars[kk]],yerr=errors[fitpars[kk]],color=colork3[w],linewidth=2,fmt="^")#,label=label_sun[w])
                w+=1
            ax1.errorbar(-5,-5,yerr=2,color='k',linewidth=2,fmt="^",label='($U_{\\odot}$,$V_{\\odot}$,$W_{\\odot}$)')
        if(fitpars[jj]=='VS'):
            zz=[6,7]
            label_uvs=['$\\bar{U}_s$ & $\\Delta U_s$','$\\bar{V}_s$ & $\\Delta V_s$']
            w=0
            for kk in zz:                    
                ax1.errorbar(N_bri_sources,means[fitpars[kk]],yerr=errors[fitpars[kk]],color=colork3[w],linewidth=2,fmt="^")#,label=label_uvs[w])
                w+=1                
            ax1.errorbar(-5,-5,yerr=2,color='k',linewidth=2,fmt="^",label='($\\bar{U}_s$, $\\bar{V}_s$)')
        ax1.grid(True)
        ax1.set_xlabel('Number of Sources')
        if(fitpars[jj]=='R0'):
            ax1.errorbar(106,8.34,yerr=0.16,fmt="*",color='b',label='Current BeSSeL',markersize=12,capsize=4, elinewidth=4)
            ax1.set_ylabel('$R_0$ [$kpc$]')
        if(fitpars[jj]=='Th0'):
            ax1.errorbar(106,240,yerr=8,fmt="*",color='b',label='Current BeSSeL',markersize=12,capsize=4, elinewidth=4)
            ax1.set_ylabel('$\\Theta_0$ [$km/s$]')            
        if(fitpars[jj]=='U0'):
            ax1.errorbar(106,11.1,yerr=3.9,fmt="*",color=colork3[0],markersize=12,capsize=4, elinewidth=4)
            ax1.errorbar(109,14.6,yerr=5,fmt="*",color=colork3[1],markersize=12,capsize=4, elinewidth=4)
            ax1.errorbar(100,7.2,yerr=1,fmt="*",color=colork3[2],markersize=12,capsize=4, elinewidth=4)
            ax1.errorbar(0,240,yerr=5,fmt="*",color='k',label='Current BeSSeL errors',markersize=12,capsize=4, elinewidth=4)
            ax1.set_ylabel('$U_{\\odot}$,$V_{\\odot}$,$W_{\\odot}$ [$km/s$]')            
        if(fitpars[jj]=='VS'):
            ax1.errorbar(106,-1.5,yerr=6.8,fmt="*",color=colork3[1],markersize=12,capsize=4, elinewidth=4)
            ax1.errorbar(106,2.8,yerr=2,fmt="*",color=colork3[2],markersize=12,capsize=4, elinewidth=4)
            ax1.errorbar(0,240,yerr=5,fmt="*",color='k',label='Current BeSSeL errors',markersize=12)
            ax1.set_ylabel('$\\bar{U}_{s}$,$\\bar{V}_{s}$ [$km/s$]')            
        if(fitpars[jj]=='R0'):
            ax1.set_ylim([8.0,8.55])
        if(fitpars[jj]=='Th0'):
            ax1.set_ylim([225.0,250.]) 
        if(fitpars[jj]=='U0'):
            ax1.set_ylim([6.0,20.]) 
        if(fitpars[jj]=='VS'):
            ax1.set_ylim([-15.0,40.0]) 
        
        ax1.set_xlim([90,520])
        #Perfect data= complete 1500 sources
        '''
        xx=[-1,2000]  
        if(fitpars[jj]=='U0' or fitpars[jj]=='VS'):
            if(fitpars[jj]=='VS'):
                w=0
                for kk in zz:
                    ax1.fill_between(xx,np.mean(dic[dic_table_names[-1][2:]][fitpars[kk]])+np.std(dic[dic_table_names[-1][2:]][fitpars[kk]]),np.mean(dic[dic_table_names[-1][2:]][fitpars[kk]])-np.std(dic[dic_table_names[-1][2:]][fitpars[kk]]), facecolor=colork3[w], alpha=0.3)
                    w+=1                
            else:
                w=0
                for kk in zz:
                    ax1.fill_between(xx,np.mean(dic[dic_table_names[-1][2:]][fitpars[kk]])+np.std(dic[dic_table_names[-1][2:]][fitpars[kk]]),np.mean(dic[dic_table_names[-1][2:]][fitpars[kk]])-np.std(dic[dic_table_names[-1][2:]][fitpars[kk]]), facecolor=colork3[w], alpha=0.3)
                    w+=1
        else:
            ax1.fill_between(xx,np.mean(dic[dic_table_names[-1][2:]][fitpars[jj]])+np.std(dic[dic_table_names[-1][2:]][fitpars[jj]]),np.mean(dic[dic_table_names[-1][2:]][fitpars[jj]])-np.std(dic[dic_table_names[-1][2:]][fitpars[jj]]), facecolor='yellow', alpha=0.2)
        '''
        if(fitpars[jj]=='R0'):
            #ax1.plot([], [], color='yellow',  alpha=0.3, linewidth=15,label='$R_0$ complete sample')
            ax1.axhline(y=modval[fitpars[jj]],color='k',linewidth=2,ls='dashed',label='$R_0$ model')            
        if(fitpars[jj]=='Th0'):
            #ax1.plot([], [], color='yellow',  alpha=0.3, linewidth=15,label='$\\Theta_0$ complete sample')
            ax1.axhline(y=modval[fitpars[jj]],color='k',linewidth=2,ls='dashed',label='$\\Theta_0$ model')
        if(fitpars[jj]=='U0'):
            w=0
            #ax1.plot([], [], color='k', alpha=0.3, linewidth=15,label='($U_{\\odot}$,$V_{\\odot}$,$W_{\\odot}$) complete sample')
            ax1.axhline(y=-240,color='k',linewidth=2,ls='dashed',label='($U_{\\odot}$,$V_{\\odot}$,$W_{\\odot}$) model')     
            ax1.text(0.5, 0.4, '$U_{\\odot}$',verticalalignment='bottom', horizontalalignment='left',transform=ax1.transAxes,color=colork3[0],fontsize=29)
            ax1.text(0.75, 0.65, '$V_{\\odot}$',verticalalignment='bottom', horizontalalignment='left',transform=ax1.transAxes,color=colork3[1],fontsize=29)            
            ax1.text(0.3, 0.09, '$W_{\\odot}$',verticalalignment='bottom', horizontalalignment='left',transform=ax1.transAxes,color=colork3[2],fontsize=29)
            for kk in zz:   
                #ax1.plot([], [], color=colork3[w],  alpha=0.3, linewidth=15,label=label_sun[w]+' complete sample')
                ax1.axhline(y=modval[fitpars[kk]],color=colork3[w],linewidth=2,ls='dashed')#,label=label_sun[w]+' model')
                w+=1
        if(fitpars[jj]=='VS'):
            w=0
            ax1.axhline(y=modval[fitpars[kk]],color=colork3[w],linewidth=2,ls='dashed',label='$\\bar{U}_s$ & $\\bar{V}_s$ model')
            #ax1.plot([], [], color='k', alpha=0.3, linewidth=15,label='($\\bar{U}_s$, $\\bar{V}_s$) complete sample')
            pp=[0.3,0.3]
            for kk in zz:   
                ax1.plot([], [], color=colork3[w],  alpha=0.3, linewidth=15)#,label=label_uvs[w]+'complete sample')
                w+=1
        #ax1.plot([], [], color='blue',  alpha=0.2, linewidth=15,label='($\\Delta \\pi$,$\\Delta \\mu$) BeSSeL')
        #ax1.axvline(x=103,color='r',linewidth=2,ls='dashed',label='Current BeSSeL')
        ax1.legend(numpoints=1,loc=0,ncol=2,prop={'size':28}) 
        plt.savefig(path_outputs+'/'+fitpars[jj]+'_nsources.png',bbox_inches='tight')  
        #------------------------------------------------------------------------------------------------------------------------
    
    if(float(initial_values['Fudge_Factor'])==1):
        Values_Fudge=initial_values['Values_Fudge']
        fudgeval=Values_Fudge.replace (",", " ")
        fudge_values = [float(x) for x in fudgeval.split()]                          
        means={}
        for keypar in fitpars:
            means[keypar]=[]
        for i in range(len(dic_table_names)-1):
            for keypar in fitpars:
                if keypar != 'void':
                    means[keypar].append(np.mean(dic[dic_table_names[i][2:]][keypar]))
        errors={}
        for keypar in fitpars:
            errors[keypar]=[]
        for i in range(len(dic_table_names)-1):
            for keypar in fitpars:
                if keypar != 'void':
                    errors[keypar].append(np.std(dic[dic_table_names[i][2:]][keypar]))

        var_csv = raw_input("Where are the .csv files or fitting files (path)?: ")
        if(var_csv[-1]!='/'):
            var_csv=var_csv+'/'
        list_csvs=os.listdir(var_csv)
        means_unc={}
        means_unc['unc_para']=[]
        means_unc['unc_mux']=[]
        means_unc['unc_muy']=[]
        csv_sort={}
        os.chdir(var_csv)
        for l in range(sam_bri):
            csv_sort['csv_'+str(l)]=[]            
        for filename in list_csvs:        
            for kk in range(sam_bri):
                if(filename[-6:]=='_'+str(kk)+'.csv'):
                    csv_sort['csv_'+str(kk)].append(filename)
        for kk in range(sam_bri):
            files_csv=csv_sort['csv_'+str(kk)]
            mean_unc_galaxy_para=[]
            mean_unc_galaxy_mux=[]            
            mean_unc_galaxy_muy=[]
            for filename in files_csv:        
                with open(filename) as csvfile:
                    #Opening the file
                    reader = csv.DictReader(csvfile)
                    row_count = sum(1 for row in reader)
                with open(filename) as csvfile:
                    #Opening the file
                    reader = csv.DictReader(csvfile)        
                    #Declaring variables
                    unc_mud=numpy.zeros(row_count)      
                    unc_para=numpy.zeros(row_count)        
                    unc_mux=numpy.zeros(row_count)   
                    yy=0             
                    for row in reader:
                        unc_mud[yy]=row['err_mud(mas/yr)']
                        unc_mux[yy]=row['err_mu_x(mas/yr)']            
                        unc_para[yy]=row['err_parallax(mas)']          
                        yy+=1
                mean_unc_galaxy_para.append(mean(unc_para))
                mean_unc_galaxy_mux.append(mean(unc_para))
                mean_unc_galaxy_muy.append(mean(unc_para))                                
            means_unc['unc_para'].append(mean(mean_unc_galaxy_para))
            means_unc['unc_mux'].append(mean(mean_unc_galaxy_mux))
            means_unc['unc_muy'].append(mean(mean_unc_galaxy_muy))                  
        resolution=[a*b for a,b in zip(means_unc['unc_para'],fudge_values)]    
  
        #Plot
        jj=1 # 0=R0,1=Th0, 3=sun, 6=uvs
        sources=100
        colork3=['g','r','m','b','m','g','c']
        matplotlib.rcParams.update({'font.size': 22})        
        fig1 = plt.figure(figsize=(20,8))
        ax1 = fig1.add_subplot(111)
        #for keypar in fitpars:
        if(sources==100):
            if(fitpars[jj]=='R0'):
                ax1.errorbar(resolution[:-1],means[fitpars[jj]][:-1],yerr=errors[fitpars[jj]][:-1],color=colork3[0],linewidth=2,fmt="^",label='$R_0$ & $\\Delta R_0$ 100 sources',capsize=4, elinewidth=4,markersize=8)
            if(fitpars[jj]=='Th0'):
                ax1.errorbar(resolution[:-1],means[fitpars[jj]][:-1],yerr=errors[fitpars[jj]][:-1],color=colork3[0],linewidth=2,fmt="^",label='$\\Theta_0$ & $\\Delta \\Theta_0$ 100 sources',capsize=4, elinewidth=4,markersize=8)
            if(fitpars[jj]=='U0'):
                zz=[3,4,5]
                label_sun=['$U_{\\odot}$','$V_{\\odot}$','$W_{\\odot}$']
                w=0
                for kk in zz:                    
                    ax1.errorbar(resolution[:-1],means[fitpars[kk]][:-1],yerr=errors[fitpars[kk]][:-1],color=colork3[w],linewidth=2,fmt="^")#,label=label_sun[w])
                    w+=1
            if(fitpars[jj]=='VS'):
                zz=[6,7]
                label_uvs=['$\\bar{U}_s$ & $\\Delta U_s$','$\\bar{V}_s$ & $\\Delta V_s$']
                w=0
                for kk in zz:                    
                    ax1.errorbar(resolution[:-1],means[fitpars[kk]][:-1],yerr=errors[fitpars[kk]][:-1],color=colork3[w],linewidth=2,fmt="^",label=label_uvs[w])
                    w+=1                
        else:
            if(fitpars[jj]=='R0'):
                ax1.errorbar(resolution[1:],means[fitpars[jj]][1:],yerr=errors[fitpars[jj]][1:],color=colork3[0],linewidth=2,fmt="^",label='$R_0$ & $\\Delta R_0$ for 200 sources')
            if(fitpars[jj]=='Th0'):
                ax1.errorbar(resolution[1:],means[fitpars[jj]][1:],yerr=errors[fitpars[jj]][1:],color=colork3[0],linewidth=2,fmt="^",label='$\\Theta_0$ & $\\Delta \\Theta_0$ for 200 sources')
            if(fitpars[jj]=='U0'):
                zz=[2,3,4]
                label_sun=['$U_{\\odot}$','$V_{\\odot}$','$W_{\\odot}$']
                w=0
                for kk in zz:                    
                    ax1.errorbar(resolution[:-1],means[fitpars[kk]][:-1],yerr=errors[fitpars[kk]][:-1],color=colork3[w],linewidth=2,fmt="^",label=(label_sun[w]))
                    w+=1
            if(fitpars[jj]=='VS'):
                zz=[6,7]
                label_uvs=['$\\bar{U}_s$ & $\\Delta U_s$','$\\bar{V}_s$ & $\\Delta V_s$']
                w=0
                for kk in zz:                    
                    ax1.errorbar(resolution[:-1],means[fitpars[kk]][:-1],yerr=errors[fitpars[kk]][:-1],color=colork3[w],linewidth=2,fmt="^",label=label_uvs[w])
                    w+=1  
        #For overplot
        #kkk=np.array(aaa)
        #kk=np.array(resolution[:-1])
        #lll=kk-abs(kkk-kk)
        res=[]
        for s in aaa: res.append(s+s*0.05)
        res_=[]
        for s in aaa: res_.append(s-s*0.03)
        ax1.errorbar(res,bbb,ccc,color=colork3[3],linewidth=2,fmt="*",label='$R_0$ & $\\Delta \\Theta_0$ 200 sources',capsize=4, elinewidth=4,markersize=8) 
        ax1.errorbar(res_,bb,cc,color=colork3[2],linewidth=2,fmt="s",label='$R_0$ & $\\Delta \\Theta_0$ 400 sources',capsize=4, elinewidth=4,markersize=8)
        #Labels
        ax1.grid(True)
        ax1.set_xlabel('Resolution [$mas$]')
        if(fitpars[jj]=='R0'):
            ax1.set_ylabel('$R_0$ [$kpc$]')
        if(fitpars[jj]=='Th0'):
            ax1.set_ylabel('$\\Theta_0$ [$kpc$]')            
        if(fitpars[jj]=='U0'):
            ax1.set_ylabel('$U_{\\odot}$,$V_{\\odot}$,$W_{\\odot}$ [$km/s$]')            
        if(fitpars[jj]=='VS'):
            ax1.set_ylabel('$\\bar{U}_{s}$,$\\bar{V}_{s}$ [$km/s$]')            
        if(fitpars[jj]=='R0'):
            ax1.set_ylim([8.0,11.5])
        if(fitpars[jj]=='Th0'):
            ax1.set_ylim([230.0,380.]) 
        if(fitpars[jj]=='U0'):
            ax1.set_ylim([6.0,25.]) 
        if(fitpars[jj]=='VS'):
            ax1.set_ylim([-15.0,40.0]) 
        ax1.set_xlim([0.00005,0.2])
        ax1.invert_xaxis()
        ax1.set_xscale('log')
        xx=[0.00004,1.16]  
        if(fitpars[jj]=='U0' or fitpars[jj]=='VS'):
            if(fitpars[jj]=='VS'):
                w=0
                for kk in zz:
                    ax1.fill_between(xx,means[fitpars[kk]][-1]+errors[fitpars[kk]][-1],means[fitpars[kk]][-1]-errors[fitpars[kk]][-1], facecolor=colork3[w], alpha=0.3)
                    w+=1                
            else:
                w=0
                for kk in zz:
                    ax1.fill_between(xx,means[fitpars[kk]][-1]+errors[fitpars[kk]][-1],means[fitpars[kk]][-1]-errors[fitpars[kk]][-1], facecolor=colork3[w], alpha=0.3)
                    w+=1
        else:
            ax1.fill_between(xx,means[fitpars[jj]][-1]+errors[fitpars[jj]][-1],means[fitpars[jj]][-1]-errors[fitpars[jj]][-1], facecolor='yellow', alpha=0.2)
        x=[0.005,0.153]  
        ax1.fill_between(x,-500,500, facecolor='blue', alpha=0.1)
        if(fitpars[jj]=='R0'):
            ax1.plot([], [], color='yellow',  alpha=0.3, linewidth=15,label='$R_0$ perfect data')
            ax1.axhline(y=modval[fitpars[jj]],color='k',linewidth=2,ls='dashed',label='$R_0$ model')            
        if(fitpars[jj]=='Th0'):
            ax1.plot([], [], color='yellow',  alpha=0.3, linewidth=15,label='$\\Theta_0$ perfect data')
            ax1.axhline(y=modval[fitpars[jj]],color='k',linewidth=2,ls='dashed',label='$\\Theta_0$ model')
        if(fitpars[jj]=='U0'):
            w=0
            for kk in zz:   
                ax1.plot([], [], color=colork3[w],  alpha=0.3, linewidth=15,label=label_sun[w]+' perfect data')
                ax1.axhline(y=modval[fitpars[kk]],color=colork3[w],linewidth=2,ls='dashed',label=label_sun[w]+' model')
                w+=1
        if(fitpars[jj]=='VS'):
            w=0
            ax1.axhline(y=modval[fitpars[kk]],color=colork3[w],linewidth=2,ls='dashed',label='$\\bar{U}_s$ & $\\bar{V}_s$ model')
            pp=[0.3,0.3]
            for kk in zz:   
                ax1.plot([], [], color=colork3[w],  alpha=0.3, linewidth=15,label=label_uvs[w]+'perfect data')
                w+=1
        ax1.plot([], [], color='blue',  alpha=0.2, linewidth=15,label='($\\Delta \\pi$,$\\Delta \\mu$) BeSSeL')
        #ax1.axvline(x=0.03,color='r',linewidth=2,ls='dashed',label='$\\bar{\\pi}$ in BeSSeL')
        ax1.legend(numpoints=1,loc=0,ncol=1) 
        ax1.get_xaxis().tick_bottom()
        ax2 = ax1.twiny()
        ax2.axes.get_xaxis().set_ticks([])
        ax2.invert_xaxis()
        ax2.set_xscale('log')
        ax2.set_xlim(ax1.get_xlim())
        if(sources==100):        
            ax2.set_xticks(resolution[:-1])
        else:
            ax2.set_xticks(resolution[1:])
        if(sources==100):        
            ax2.set_xticklabels(fudge_values[:-1])
        else:
            ax2.set_xticklabels(fudge_values[1:])
        ax2.set_xlabel("Fudge Factor")
        #fig1.clf()
        #plt.savefig(path_outputs+'/'+fitpars[jj]+'_fudge.png',bbox_inches='tight')  
        '''
        aa=resolution[1:]
        bb=means['Th0'][1:]
        cc=errors['Th0'][1:]
        ax1.errorbar(aa,bb,cc,color=colork3[0],linewidth=2,fmt="*",label='$R_0$ result and its uncertainty for 400 sources') 
        aaa=resolution[1:]
        bbb=means['Th0'][1:]
        ccc=errors['Th0'][1:]        
        ax1.errorbar(aaa,bbb,ccc,color=colork3[0],linewidth=2,fmt="*",label='$R_0$ result and its uncertainty for 200 sources') 
        plt.savefig(path_outputs+'/Pearson_erros_muxy.png',bbox_inches='tight')  
        '''
        #------------------------------------------------------------------------------------------------            
    if(clue=='nsources'):
        colork3=['g','c','r','b','m']
        matplotlib.rcParams.update({'font.size': 22})        
        fig1 = plt.figure(figsize=(20,8))
        ax1 = fig1.add_subplot(111)
        k5=['R0-Th0','Th0-V0','Th0-VS','V0-VS','U0-US']
        w=0
        k4=[0,7,9,14,16]
        for k in range(len(k4)):
            ax1.errorbar(err_muxy[1:],cov_mean_all[k5[k]][1:],yerr=cov_std_all[k5[k]][1:],color=colork3[w])
            #scatter(errmuxy,cov_mean_all[k5[k]],color=colork3[w],s=50) #        
            ax1.plot(err_muxy[1:],cov_mean_all[k5[k]][1:],color=colork3[w], linewidth=2,label=pearso_coeffs_labels[k4[k]]) #        
            w+=1
        ax1.legend(scatterpoints=1,loc=8,ncol=3,numpoints=1)
        ax1.set_xlabel('Number of Sources')
        ax1.set_ylabel('Pearson Coefficient')
        ax1.set_ylim([-2.3,1.0])    
        plt.savefig(path_outputs+'/meanpears_reid_erros_muxy.png',bbox_inches='tight')   
else:
    os.chdir(path_origin)
    print('No plots!')
    print('Done with simulation, fitting and analysis')
    print('Thanks for using me!, please remember to reference us')  
        #------------------------------------------------------------------------------------------------
#Reading the files                         
if (plots_confirmation=='yes'):
    from scipy.stats import norm as norm
    colork2=['g','k','r','b','y','m','c','orange','pink','sienna','peru','dimgrey',
               'peachpuff','darkgoldenrod','firebrick','silver','rosybrown','crimson',
               'cadetblue','tomato','chocolate']


    # Seven Plots
    fig1 = plt.figure(figsize=(15,10))
    matplotlib.rcParams.update({'font.size': 16})
    zz=0 #which sample?
    ax = fig1.add_subplot(111)
            
    ax1 = fig1.add_subplot(241)
    mean_R0 = modval['R0'] # initial values
    (mu_R0, sigma_R0) = norm.fit(dic[dic_table_names[zz][2:]]['R0'])
    normed_R0 = [x / mean_R0 for x in dic[dic_table_names[zz][2:]]['R0']]    
    (mu_normed_R0, sigma_normed_R0) = norm.fit(normed_R0)
    a,bins_data,c=ax1.hist(normed_R0,bins=8,color=colork2[0], alpha=1.0) #     
    area = sum(numpy.diff(bins_data)*a)
    #xmin,xmax=ax5.get_xlim()
    xmin,xmax=ax1.get_xlim()
    x = np.linspace(xmin, xmax, 100)
    y=mlab.normpdf(x, mu_normed_R0, sigma_normed_R0)*area
    ax1.plot(x,y, linestyle='--',color='k',linewidth=2.0)
    ax1.grid(True)
    plt.setp(ax1.get_xticklabels(), visible=False) 
    ax1.set_xticks([0.96,0.98,1.0,1.02,1.04])
    ax1.text( 1.02, 28, '$R_0$', verticalalignment='bottom', horizontalalignment='left', color=colork2[0], fontsize=27)

    ax2 = fig1.add_subplot(242,sharey=ax1)
    mean_US = modval['US'] # initial values    
    (mu_US, sigmaUS) = norm.fit(dic[dic_table_names[zz][2:]]['US'])
    normed_US = [x / mean_US for x in dic[dic_table_names[zz][2:]]['US']]
    (mu_normed_US, sigma_normed_US) = norm.fit(normed_US)
    a,bins_data,c=ax2.hist(normed_US,bins=8,color=colork2[2], alpha=1.0) # 
    area = sum(numpy.diff(bins_data)*a)
    #xmin,xmax=ax7.get_xlim()
    xmin,xmax=ax2.get_xlim()
    x = np.linspace(xmin, xmax, 100)
    y=mlab.normpdf(x, mu_normed_US, sigma_normed_US)*area
    ax2.plot(x,y, linestyle='--',color='k',linewidth=2.0)    
    ax2.grid(True)
    plt.setp(ax2.get_yticklabels(), visible=False)    
    plt.setp(ax2.get_xticklabels(), visible=False)    
    ax2.set_xticks([-0.5,0.0,0.5,1.0,1.5,2.0,2.5])
    ax2.text( 1.9,  28, '$U_s$', verticalalignment='bottom', horizontalalignment='left', color=colork2[2], fontsize=27)               
      
    mean_U0 = modval['U0'] # initial values    
    ax3 = fig1.add_subplot(243,sharey=ax1)
    (mu_U0, sigmaU0) = norm.fit(dic[dic_table_names[zz][2:]]['U0'])
    normed_U0 = [x / mean_U0 for x in dic[dic_table_names[zz][2:]]['U0']]
    (mu_normed_U0, sigma_normed_U0) = norm.fit(normed_U0)
    a,bins_data,c=ax3.hist(normed_U0,bins=8,color=colork2[3], alpha=1.0) # 
    area = sum(numpy.diff(bins_data)*a)
    xmin,xmax=ax3.get_xlim()
    #xmin,xmax=ax7.get_xlim()    
    x = np.linspace(xmin, xmax, 100)
    y=mlab.normpdf(x, mu_normed_U0, sigma_normed_U0)*area
    ax3.plot(x,y, linestyle='--',color='k',linewidth=2.0)    
    ax3.grid(True)
    plt.setp(ax3.get_yticklabels(), visible=False)
    plt.setp(ax3.get_xticklabels(), visible=False)
    ax3.set_xticks([0.55,0.7,0.85,1.0,1.15,1.3,1.45])    
    ax3.text( 1.13,  28, '$U_{\\odot}$', verticalalignment='bottom', horizontalalignment='left', color=colork2[3], fontsize=27)         
    
    ax4 = fig1.add_subplot(244,sharey=ax1)
    mean_dTd = modval['dTd'] # initial values
    (mu_mean_dTd, sigma_mean_dTd) = norm.fit(dic[dic_table_names[zz][2:]]['dTd'])
    normed_dTd = [x / mean_dTd for x in dic[dic_table_names[zz][2:]]['dTd']]    
    (mu_normed_dTd, sigma_normed_dTd) = norm.fit(normed_dTd)
    a,bins_data,c=ax4.hist(normed_dTd,bins=8,color=colork2[9], alpha=1.0) #     
    area = sum(numpy.diff(bins_data)*a)
    #xmin,xmax=ax8.get_xlim()
    xmin,xmax=ax4.get_xlim()    
    x = np.linspace(xmin, xmax, 100)
    y=mlab.normpdf(x, mu_normed_dTd, sigma_normed_dTd)*area
    ax4.plot(x,y, linestyle='--',color='k',linewidth=2.0)
    ax4.grid(True)
    plt.setp(ax4.get_xticklabels(), visible=False)
    plt.setp(ax4.get_yticklabels(), visible=False)     
    ax4.set_xticks([-11,-7,-3,1,5,9,13])
    ax4.text( 8, 26.7, '$\\frac{d \\theta}{dR}$', verticalalignment='bottom', horizontalalignment='left', color=colork2[9], fontsize=27)
          
    ax5 = fig1.add_subplot(245,sharex=ax1)
    mean_Th0 = modval['Th0'] # initial values    
    (mu_Th0, sigmaTh0) = norm.fit(dic[dic_table_names[zz][2:]]['Th0'])
    normed_Th0 = [x / mean_Th0 for x in dic[dic_table_names[zz][2:]]['Th0']]
    (mu_normed_Th0, sigma_normed_Th0) = norm.fit(normed_Th0)
    a,bins_data,c=ax5.hist(normed_Th0,bins=8,color=colork2[5], alpha=1.0) # 
    area = sum(numpy.diff(bins_data)*a)
    xmin,xmax=ax5.get_xlim()
    x = np.linspace(xmin, xmax, 100)
    y=mlab.normpdf(x, mu_normed_Th0, sigma_normed_Th0)*area
    ax5.plot(x,y, linestyle='--',color='k',linewidth=2.0)    
    ax5.grid(True)
    ax5.set_yticks([0,5,10,15,20,25,30])  
    ax5.text( 1.02,  23.5, '$\\Theta_0$', verticalalignment='bottom', horizontalalignment='left', color=colork2[5], fontsize=27)         
  
    ax6 = fig1.add_subplot(246,sharey=ax5,sharex=ax2)
    mean_V0 = modval['V0'] # initial values    
    (mu_V0, sigmaV0) = norm.fit(dic[dic_table_names[zz][2:]]['V0'])
    normed_V0 = [x / mean_V0 for x in dic[dic_table_names[zz][2:]]['V0']]
    (mu_normed_V0, sigma_normed_V0) = norm.fit(normed_V0)
    a,bins_data,c=ax6.hist(normed_V0,bins=8,color=colork2[6], alpha=1.0) #     
    area = sum(numpy.diff(bins_data)*a)
    #xmin,xmax=ax6.get_xlim()
    xmin,xmax=ax6.get_xlim()    
    x = np.linspace(xmin, xmax, 100)
    y=mlab.normpdf(x, mu_normed_V0, sigma_normed_V0)*area
    ax6.plot(x,y, linestyle='--',color='k',linewidth=2.0)    
    ax6.grid(True)
    plt.setp(ax6.get_yticklabels(), visible=False) 
    ax6.text( 1.9,  23, '$V_{\\odot}$', verticalalignment='bottom', horizontalalignment='left', color=colork2[6], fontsize=27)        
 
    mean_W0 = modval['W0'] # initial values    
    ax7 = fig1.add_subplot(247,sharey=ax5,sharex=ax3)
    (mu_W0, sigmaW0) = norm.fit(dic[dic_table_names[zz][2:]]['W0'])
    normed_W0 = [x / mean_W0 for x in dic[dic_table_names[zz][2:]]['W0']]
    (mu_normed_W0, sigma_normed_W0) = norm.fit(normed_W0)
    a,bins_data,c=ax7.hist(normed_W0,bins=8,color=colork2[7], alpha=1.0) # 
    area = sum(numpy.diff(bins_data)*a)
    xmin,xmax=ax7.get_xlim()
    x = np.linspace(xmin, xmax, 100)
    y=mlab.normpdf(x, mu_normed_W0, sigma_normed_W0)*area
    ax7.plot(x,y, linestyle='--',color='k',linewidth=2.0)    
    ax7.grid(True)
    plt.setp(ax7.get_yticklabels(), visible=False) 
    ax7.text( 1.13,  23, '$W_{\\odot}$', verticalalignment='bottom', horizontalalignment='left', color=colork2[7], fontsize=27)         

    mean_VS = modval['VS'] # initial values    
    ax8 = fig1.add_subplot(248,sharey=ax5,sharex=ax4)
    (mu_VS, sigmaVS) = norm.fit(dic[dic_table_names[zz][2:]]['VS'])
    normed_VS = [x / mean_VS for x in dic[dic_table_names[zz][2:]]['VS']]
    (mu_normed_VS, sigma_normed_VS) = norm.fit(normed_VS)
    a,bins_data,c=ax8.hist(normed_VS,bins=8,color=colork2[14], alpha=1.0) # 
    area = sum(numpy.diff(bins_data)*a)
    xmin,xmax=ax8.get_xlim()
    x = np.linspace(xmin, xmax, 100)
    y=mlab.normpdf(x, mu_normed_VS, sigma_normed_VS)*area
    ax8.plot(x,y, linestyle='--',color='k',linewidth=2.0)    
    ax8.grid(True)
    plt.setp(ax8.get_yticklabels(), visible=False)  
    ax8.text( 8,  23, '$V_s$', verticalalignment='bottom', horizontalalignment='left', color=colork2[14], fontsize=27)            
                
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.setp(ax.get_yticklabels(), visible=False)  
    plt.setp(ax.get_xticklabels(), visible=False)          
    ax.set_ylabel('Counts',labelpad=34)
    ax.set_xlabel('Normed Parameter',labelpad=40)
    ax.tick_params(
    axis='both',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off',right='off', left='off', labelleft='off')
    #fig1.patch.set_visible(False)
    fig1.tight_layout()
    
    
    
    plt.savefig('../../plots/seven.pdf',bbox_inches='tight')
    
    modval2=modval.copy()
    parameters=modval2.keys()
    means=[mu_R0,mu_U0,mu_US,mu_mean_dTd,mu_V0,mu_VS,mu_W0,mu_Th0]
    sigmas=[sigma_R0,sigmaU0,sigmaUS,sigma_mean_dTd,sigmaV0,sigmaVS,sigmaW0,sigmaTh0]
    for i in range(len(parameters)):
        modval2[parameters[i]]=[]
        modval2[parameters[i]].append(modval[parameters[i]])        
        modval2[parameters[i]].append(means[i])
        modval2[parameters[i]].append(sigmas[i])

    print('---------------------------------------------------------------')
    print('------------Final Values---------------------------------------------')
    print('---------------------------------------------------------------')
    print('Parameter & Value & Mean & Sigma&')
    for p in parameters:
        print('%s %s %s %s' % (p+'       &', str(modval2[p][0])+'  &',str(round(modval2[p][1],2))+' &',str(round(modval2[p][2],2))+' &'))    

    #Plots R_0
    mean_r0, sigma_r0 = modval['R0'], 0.16 # mean and standard deviation
    x = np.linspace(modval['R0']-0.5,modval['R0']+0.5,1000) #Linear Space

    fig1 = plt.figure(figsize=(15,15))
    ax1 = fig1.add_subplot(111)
    matplotlib.rcParams.update({'font.size': 23})
    nbins=10
    #plt.plot(x, (norm.pdf(x, mean_r0, sigma_r0)),'b--',label='Reid et al. 2014')#*(len(sample_100['R0'])/10.0)))
    #only for 100 sources fudge = 0.1
    #bins_data= [8.163 ,  8.2019,  8.2408,  8.2797,  8.3186,  8.3575,  8.3964,
    #8.4353,  8.4742,  8.5131,  8.552]
    zz=0 #which sample?
    (mu, sigma) = norm.fit(dic[dic_table_names[zz][2:]]['R0'])
    #a,bins_data,c=ax1.hist(dic[dic_table_names[zz][2:]]['R0'],bins_data,color=colork2[0], alpha=1.0)#,normed=True) # 
    normed_R0 = [x / mean_r0 for x in dic[dic_table_names[zz][2:]]['R0']]
    a,bins_data,c=ax1.hist(normed_R0,bins=8,color=colork2[0], alpha=1.0,normed=True) # 
    #ax1.axvline(x=mean_r0,color='k',ls='dashed',label='Initial $R_0$=8.34 kpc')
    #ax1.axvline(x=np.mean(dic[dic_table_names[zz][2:]]['R0']),color='g',ls='dashed',label='Mean value $R_0$=8.35 $\\pm$ 0.08 kpc')
    ax1.grid(True)
    ax1.legend(loc=9)
    #text( 7.85,  1.5, '$\\bar{R}=%s$ $kpc$' % mean_r0, verticalalignment='bottom', horizontalalignment='left', color='b', fontsize=23)
    #text( 7.85,  1.1, '$\\sigma=%s$ $kpc$'  % sigma_r0, verticalalignment='bottom', horizontalalignment='left', color='b', fontsize=23)
    #text( 7.85,  2.5, '$\\bar{R}=%s$ $kpc$' % round(mu,3), verticalalignment='bottom', horizontalalignment='left', color='g', fontsize=23)
    #text( 7.85,  2.1, '$\\sigma=%s$ $kpc$' % round(sigma,3), verticalalignment='bottom', horizontalalignment='left', color='g', fontsize=23)
    xlabel('Distance Center Galaxy $R_0$ (kpc)')
    ylabel('N1ormed Counts') 
    #plt.savefig(path_outputs+'/R_0_100.pdf',bbox_inches='tight')
    
    #Histogram for theta0
    mean_theta0, sigma_theta0 = modval['Th0'], 8 # mean and standard deviation
    x = np.linspace(modval['Th0']-15.5,modval['Th0']+15.5,1000) #Linear Space
    #Histogram for Theta0
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    matplotlib.rcParams.update({'font.size': 17})
    nbins=10
    zz=3 #which sample?
    plt.plot(x, (norm.pdf(x, mean_theta0, sigma_theta0)),'b--',label='Reid et al. 2014')#*(len(sample_100['R0'])/10.0)))
    a,bins_data,c=ax1.hist(dic[dic_table_names[zz][2:]]['Th0'],nbins,color='r', alpha=0.5, normed=True) # ) #       
    (mu, sigma) = norm.fit(dic[dic_table_names[zz][2:]]['Th0'])
    y = mlab.normpdf(bins_data, mu, sigma)
    #ax1.plot(bins_data, y, 'r--', linewidth=2,label='Gaussian Fitting')     
    ax1.axvline(x=mean_theta0,color='k',ls='dashed',label='Initial $\\Theta_0=240$ $km/s$')
    ax1.grid(True)
    ax1.legend(loc=2)
    text( 243,  0.09, '$\\bar{\\theta}_0= \ %s$ $km/s$' % mean_theta0,     verticalalignment='bottom', horizontalalignment='left',     color='b', fontsize=23)
    text( 243,  0.08, '$\\sigma= \ %s$ $km/s$' %sigma_theta0,     verticalalignment='bottom', horizontalalignment='left',     color='b', fontsize=23) 
    text( 243,  0.12, '$\\bar{\\theta}_0= \ %s$ $km/s$' % round(mu,3),     verticalalignment='bottom', horizontalalignment='left',     color='r', fontsize=23)
    text( 243,  0.11, '$\\sigma= \ %s$ $km/s$ ' % round(sigma,3),  verticalalignment='bottom', horizontalalignment='left',     color='r', fontsize=23)
    xlabel('Circular Rotation Speed of the Sun $\\Theta_0$ ($km/s$)')
    ylabel('Normed Counts')
    xlim([225,255])
    #title("Results for $\\Theta_0$")#Parallaxes")
    plt.savefig(path_outputs+'/Theta_0_100.pdf',bbox_inches='tight')
    
    
    from scipy.stats import norm as norm
    #UVW sun
    zz=3
    mean_U0, sigma_U0 = modval['U0'], aprior_err['U0'] # mean and standard deviation
    x_U0 = np.linspace(modval['U0']-5.5,modval['U0']+5.5,1000) #Linear Space
    mean_V0, sigma_V0 = modval['V0'], aprior_err['V0'] # mean and standard deviation
    x_V0 = np.linspace(modval['V0']-5.5,modval['V0']+5.5,1000) #Linear Space
    mean_W0, sigma_W0 = modval['W0'], aprior_err['W0'] # mean and standard deviation
    x_W0 = np.linspace(modval['W0']-5.5,modval['W0']+5.5,1000) #Linear Space
    #plot UVW_sun
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    matplotlib.rcParams.update({'font.size': 17})
    nbins=10
    m1,=plt.plot(x_U0, norm.pdf(x_U0, mean_U0, sigma_U0),'b',label='Reid et al.')#*(len(sample_100['R0'])/10.0)))
    m2,=plt.plot(x_V0, norm.pdf(x_V0, mean_V0, sigma_V0),'g')#,label='Reid 2014')#*(len(sample_100['R0'])/10.0)))
    m3,=plt.plot(x_W0, norm.pdf(x_W0, mean_W0, sigma_W0),'r')#,label='Reid 2014')#*(len(sample_100['R0'])/10.0)))
    a,bins_data_U0,c=ax1.hist(dic[dic_table_names[zz][2:]]['U0'],nbins,color='b', alpha=0.5,normed=True)#, label="100 sources",normed=True) #       
    a,bins_data_V0,c=ax1.hist(dic[dic_table_names[zz][2:]]['V0'],nbins,color='g', alpha=0.5,normed=True)#, label="100 sources",normed=True) #       
    a,bins_data_W0,c=ax1.hist(dic[dic_table_names[zz][2:]]['W0'],nbins,color='r', alpha=0.5,normed=True)#, label="100 sources",normed=True) #       
    (mu_U0, sigm_U0) = norm.fit(dic[dic_table_names[zz][2:]]['U0'])
    (mu_V0, sigm_V0) = norm.fit(dic[dic_table_names[zz][2:]]['V0'])
    (mu_W0, sigm_W0) = norm.fit(dic[dic_table_names[zz][2:]]['W0'])
    y_U0 = mlab.normpdf(bins_data_U0, mu_U0, sigm_U0)
    y_V0 = mlab.normpdf(bins_data_V0, mu_V0, sigm_V0)
    y_W0 = mlab.normpdf(bins_data_W0, mu_W0, sigm_W0)   
    ax1.axvline(x=mean_U0,color='b',label='Initial $U_{\\odot}= \ %s$ $km/s$' % mean_U0)
    ax1.axvline(x=mean_V0,color='g',label='Initial $V_{\\odot}= \ %s$ $km/s$' % mean_V0)
    ax1.axvline(x=mean_W0,color='r',label='Initial $W_{\\odot}= \ %s$ $km/s$' % mean_W0)
    ax1.legend(loc=(0.31,0.6))
    ax1.grid(True)
    legend = ax1.get_legend()
    legend.legendHandles[0].set_color('k')
    #U0
    text( 15.2,  1.08, '$\\bar{U}_{\\odot}= \ %s$ $km/s$' % mean_U0,     verticalalignment='bottom', horizontalalignment='left',     color='k', fontsize=23)
    text( 15.2,  1.00, '$\\sigma= \ %s$ $km/s$' %sigma_U0,     verticalalignment='bottom', horizontalalignment='left',     color='k', fontsize=23) 
    text( 15.2,  0.88, '$\\bar{U}_{\\odot}= \ %s$ $km/s$' % round(mu_U0,3),     verticalalignment='bottom', horizontalalignment='left',     color='b', fontsize=23)
    text( 15.2,  0.80, '$\\sigma= \ %s$ $km/s$ ' % round(sigm_U0,3),  verticalalignment='bottom', horizontalalignment='left',     color='b', fontsize=23)
    #V0
    text( 15.2,  0.48, '$\\bar{V}_{\\odot}= \ %s$ $km/s$' % mean_V0,     verticalalignment='bottom', horizontalalignment='left',     color='k', fontsize=23)
    text( 15.2,  0.40, '$\\sigma= \ %s$ $km/s$' %sigma_V0,     verticalalignment='bottom', horizontalalignment='left',     color='k', fontsize=23) 
    text( 15.2,  0.28, '$\\bar{V}_{\\odot}= \ %s$ $km/s$' % round(mu_V0,3),     verticalalignment='bottom', horizontalalignment='left',     color='g', fontsize=23)
    text( 15.2,  0.20, '$\\sigma= \ %s$ $km/s$ ' % round(sigm_V0,3),  verticalalignment='bottom', horizontalalignment='left',     color='g', fontsize=23)
    #W0
    text( 3.7,  1.08, '$\\bar{W}_{\\odot}= \ %s$ $km/s$' % mean_W0,     verticalalignment='bottom', horizontalalignment='left',     color='k', fontsize=23)
    text( 3.7,  1.00, '$\\sigma= \ %s$ $km/s$' %sigma_W0,     verticalalignment='bottom', horizontalalignment='left',     color='k', fontsize=23) 
    text( 3.7,  0.88, '$\\bar{W}_{\\odot}= \ %s$ $km/s$' % round(mu_W0,3),     verticalalignment='bottom', horizontalalignment='left',     color='r', fontsize=23)
    text( 3.7,  0.80, '$\\sigma= \ %s$ $km/s$ ' % round(sigm_W0,3),  verticalalignment='bottom', horizontalalignment='left',     color='r', fontsize=23)
    #title("Results for $U_{\\odot}$ ")
    xlabel('Solar Velocity Components ($U_{\\odot},V_{\\odot},W_{\\odot}$) [$km/s$]')
    ylabel('Normed Counts')
    xlim([3.5,18.6])
    plt.savefig(path_outputs+'/UVW_100.pdf',bbox_inches='tight')
    

    if(clue=='nsources'):
        bins_data= [8.163 ,  8.2019,  8.2408,  8.2797,  8.3186,  8.3575,  8.3964,
        8.4353,  8.4742,  8.5131,  8.552]
        zz=np.arange(0,len(dic_table_names)-1,2)
        for i in zz:
            a,bins_data,c=ax1.hist(dic[dic_table_names[i][2:]]['R0'],bins_data,color=colork2[i], alpha=0.2, label=str(N_bri_sources[i]) +" masers",normed=True) # 
        #np.histogram(sample_100['R0'], bins=10, normed=True)
        #b,bins_datab,cb=ax1.hist(sample_200['R0'],bins_data,color='b', alpha=0.5, label="200 sources") #       
        #c,bins_datac,cc=ax1.hist(sample_complete['R0'],bins_data,color='r', alpha=0.5, label="Complete") #   
        #(mu, sigma) = norm.fit(sample_100['R0'])
        #y = mlab.normpdf(bins_data, mu, sigma)
        #ax1.plot(bins_data, y, 'r--', linewidth=2,label='Gaussian Fitting')
        ax1.axvline(x=mean_r0,color='k',ls='dashed',label='Initial $R_0=8.34$ $kpc$')
        ax1.grid(True)
        ax1.legend(loc=0)
        #text( 7.82,  7.0, '$\\bar{R}=8.34$ $kpc$, $\\sigma=0.16$ $kpc$',     verticalalignment='bottom', horizontalalignment='left',     color='b', fontsize=21)
        text( 7.85,  7.0, '$\\bar{R}=%s$ $kpc$' % mean_r0,     verticalalignment='bottom', horizontalalignment='left',     color='b', fontsize=23)
        text( 7.85,  6.5, '$\\sigma=%s$ $kpc$' %sigma_r0,     verticalalignment='bottom', horizontalalignment='left',     color='b', fontsize=23) 
        #text( 7.82,  6.0, '$\\bar{R}=%s$ $kpc$, $\\sigma=%s$ $kpc$' % (round(mu,3), round(sigma,3)),     verticalalignment='bottom', horizontalalignment='left',     color='r', fontsize=21)
        #text( 7.85,  5.5, '$\\bar{R}=%s$ $kpc$' % round(mu,3),     verticalalignment='bottom', horizontalalignment='left',     color='r', fontsize=23)
        #text( 7.85,  5.0, '$\\sigma=%s$ $kpc$ ' % round(sigma,3),  verticalalignment='bottom', horizontalalignment='left',     color='r', fontsize=23)
        xlabel('Distance Center Galaxy $R_0$ ($kpc$)')
        ylabel('Normed Counts')
        #title("Results for $R_0$ ")# Parallaxes")
        plt.savefig(path_outputs+'/R_0_100.pdf',bbox_inches='tight')
    
    #Histogram for theta0
    mean_theta0, sigma_theta0 = modval['Th0'], aprior_err['Th0'] # mean and standard deviation
    x = np.linspace(modval['Th0']-15.5,modval['Th0']+15.5,1000) #Linear Space
    #Histogram for Theta0
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    matplotlib.rcParams.update({'font.size': 17})
    nbins=10
    plt.plot(x, (norm.pdf(x, mean_theta0, sigma_theta0)),'b--',label='Reid 2014')#*(len(sample_100['R0'])/10.0)))
    a,bins_data,c=ax1.hist(sample_100['Th0'],nbins,color='g', alpha=0.5, label="100 sources",normed=True) # ) #       
    #b,bins_datab,cb=ax1.hist(sample_200['Th0'],bins_data,color='b', alpha=0.5, label="200 sources") #       
    #c,bins_datac,cc=ax1.hist(sample_complete['Th0'],bins_data,color='r', alpha=0.5, label="Complete") #  
    (mu, sigma) = norm.fit(sample_100['Th0'])
    y = mlab.normpdf(bins_data, mu, sigma)
    ax1.plot(bins_data, y, 'r--', linewidth=2,label='Gaussian Fitting')     
    ax1.axvline(x=mean_theta0,color='k',ls='dashed',label='Initial $\\Theta_0=240$ $km/s$')
    ax1.grid(True)
    ax1.legend(loc=0)
    text( 226,  0.22, '$\\bar{\\theta}_0= \ %s$ $km/s$' % mean_theta0,     verticalalignment='bottom', horizontalalignment='left',     color='b', fontsize=23)
    text( 226,  0.20, '$\\sigma= \ %s$ $km/s$' %sigma_theta0,     verticalalignment='bottom', horizontalalignment='left',     color='b', fontsize=23) 
    text( 226,  0.17, '$\\bar{\\theta}_0= \ %s$ $km/s$' % round(mu,3),     verticalalignment='bottom', horizontalalignment='left',     color='r', fontsize=23)
    text( 226,  0.15, '$\\sigma= \ %s$ $km/s$ ' % round(sigma,3),  verticalalignment='bottom', horizontalalignment='left',     color='r', fontsize=23)
    xlabel('Circular Rotation Speed of the Sun $\\Theta_0$ ($km/s$)')
    ylabel('Normed Counts')
    xlim([225,255])
    #title("Results for $\\Theta_0$")#Parallaxes")
    plt.savefig(path_outputs+'/Theta_0_100.pdf',bbox_inches='tight')
    
    from scipy.stats import norm as norm
    #UVW sun
    mean_U0, sigma_U0 = modval['U0'], aprior_err['U0'] # mean and standard deviation
    x_U0 = np.linspace(modval['U0']-5.5,modval['U0']+5.5,1000) #Linear Space
    mean_V0, sigma_V0 = modval['V0'], aprior_err['V0'] # mean and standard deviation
    x_V0 = np.linspace(modval['V0']-5.5,modval['V0']+5.5,1000) #Linear Space
    mean_W0, sigma_W0 = modval['W0'], aprior_err['W0'] # mean and standard deviation
    x_W0 = np.linspace(modval['W0']-5.5,modval['W0']+5.5,1000) #Linear Space
    #plot UVW_sun
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    matplotlib.rcParams.update({'font.size': 17})
    nbins=10
    m1,=plt.plot(x_U0, norm.pdf(x_U0, mean_U0, sigma_U0),'b',label='Reid 2014')#*(len(sample_100['R0'])/10.0)))
    m2,=plt.plot(x_V0, norm.pdf(x_V0, mean_V0, sigma_V0),'g')#,label='Reid 2014')#*(len(sample_100['R0'])/10.0)))
    m3,=plt.plot(x_W0, norm.pdf(x_W0, mean_W0, sigma_W0),'r')#,label='Reid 2014')#*(len(sample_100['R0'])/10.0)))
    a,bins_data_U0,c=ax1.hist(sample_100['U0'],nbins,color='b', alpha=0.5,normed=True)#, label="100 sources",normed=True) #       
    a,bins_data_V0,c=ax1.hist(sample_100['V0'],nbins,color='g', alpha=0.5,normed=True)#, label="100 sources",normed=True) #       
    a,bins_data_W0,c=ax1.hist(sample_100['W0'],nbins,color='r', alpha=0.5,normed=True)#, label="100 sources",normed=True) #       
    (mu_U0, sigm_U0) = norm.fit(sample_100['U0'])
    (mu_V0, sigm_V0) = norm.fit(sample_100['V0'])
    (mu_W0, sigm_W0) = norm.fit(sample_100['W0'])
    y_U0 = mlab.normpdf(bins_data_U0, mu_U0, sigm_U0)
    y_V0 = mlab.normpdf(bins_data_V0, mu_V0, sigm_V0)
    y_W0 = mlab.normpdf(bins_data_W0, mu_W0, sigm_W0)
    n1,=ax1.plot(bins_data_U0, y_U0, 'k--', linewidth=2,label='Gaussian Fitting')     
    n2=ax1.plot(bins_data_V0, y_V0, 'k--', linewidth=2)#,label='Gaussian Fitting')     
    n3=ax1.plot(bins_data_W0, y_W0, 'k--', linewidth=2)#,label='Gaussian Fitting')     
    #ax1.axvline(x=mean_U0,color='b',label='Initial $U_{\\odot}= \ %s$ $km/s$' % mean_U0)
    #ax1.axvline(x=mean_V0,color='g',label='Initial $V_{\\odot}= \ %s$ $km/s$' % mean_V0)
    #ax1.axvline(x=mean_W0,color='r',label='Initial $W_{\\odot}= \ %s$ $km/s$' % mean_W0)
    ax1.legend(loc=0)
    ax1.grid(True)
    legend = ax1.get_legend()
    legend.legendHandles[0].set_color('k')
    #U0
    text( 8.3,  1.22, '$\\bar{U}_{\\odot}= \ %s$ $km/s$' % mean_U0,     verticalalignment='bottom', horizontalalignment='left',     color='b', fontsize=23)
    text( 8.3,  1.14, '$\\sigma= \ %s$ $km/s$' %sigma_U0,     verticalalignment='bottom', horizontalalignment='left',     color='b', fontsize=23) 
    text( 8.3,  1.02, '$\\bar{U}_{\\odot}= \ %s$ $km/s$' % round(mu_U0,3),     verticalalignment='bottom', horizontalalignment='left',     color='k', fontsize=23)
    text( 8.3,  0.94, '$\\sigma= \ %s$ $km/s$ ' % round(sigm_U0,3),  verticalalignment='bottom', horizontalalignment='left',     color='k', fontsize=23)
    #V0
    text( 16,  0.48, '$\\bar{V}_{\\odot}= \ %s$ $km/s$' % mean_V0,     verticalalignment='bottom', horizontalalignment='left',     color='g', fontsize=23)
    text( 16,  0.40, '$\\sigma= \ %s$ $km/s$' %sigma_V0,     verticalalignment='bottom', horizontalalignment='left',     color='g', fontsize=23) 
    text( 16,  0.28, '$\\bar{V}_{\\odot}= \ %s$ $km/s$' % round(mu_V0,3),     verticalalignment='bottom', horizontalalignment='left',     color='k', fontsize=23)
    text( 16,  0.20, '$\\sigma= \ %s$ $km/s$ ' % round(sigm_V0,3),  verticalalignment='bottom', horizontalalignment='left',     color='k', fontsize=23)
    #W0
    text( 0.5,  0.48, '$\\bar{W}_{\\odot}= \ %s$ $km/s$' % mean_W0,     verticalalignment='bottom', horizontalalignment='left',     color='r', fontsize=23)
    text( 0.5,  0.40, '$\\sigma= \ %s$ $km/s$' %sigma_W0,     verticalalignment='bottom', horizontalalignment='left',     color='r', fontsize=23) 
    text( 0.5,  0.28, '$\\bar{W}_{\\odot}= \ %s$ $km/s$' % round(mu_W0,3),     verticalalignment='bottom', horizontalalignment='left',     color='k', fontsize=23)
    text( 0.5,  0.20, '$\\sigma= \ %s$ $km/s$ ' % round(sigm_W0,3),  verticalalignment='bottom', horizontalalignment='left',     color='k', fontsize=23)
    #title("Results for $U_{\\odot}$ ")
    xlabel('Solar Velocity Components ($U_{\\odot},V_{\\odot},W_{\\odot}$) [$km/s$]')
    ylabel('Normed Counts')
    xlim([0.0,22])
    plt.savefig(path_outputs+'/UVW_100.pdf',bbox_inches='tight')
    
    '''
    
    #Histogram for U_sun
    mean_U0, sigma_U0 = modval['U0'], aprior_err['U0'] # mean and standard deviation
    x_U0 = np.linspace(modval['U0']-5.5,modval['U0']+5.5,1000) #Linear Space
    #plot U_sun
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    matplotlib.rcParams.update({'font.size': 17})
    nbins=10
    plt.plot(x_U0, (norm.pdf(x_U0, mean_U0, sigma_U0)),'b--',label='Reid 2014')#*(len(sample_100['R0'])/10.0)))
    a,bins_data_U0,c=ax1.hist(sample_100['U0'],nbins,color='g', alpha=0.5, label="100 sources",normed=True) #       
    #b,bins_datab,cb=ax1.hist(sample_200['U0'],bins_data,color='b', alpha=0.5, label="200 sources") #       
    #c,bins_datac,cc=ax1.hist(sample_complete['U0'],bins_data,color='r', alpha=0.5, label="Complete") # 
    (mu_U0, sigm_U0) = norm.fit(sample_100['U0'])
    y = mlab.normpdf(bins_data_U0, mu_U0, sigm_U0)
    ax1.plot(bins_data_U0, y, 'r--', linewidth=2,label='Gaussian Fitting')     
    ax1.axvline(x=mean_U0,color='k',ls='dashed',label='$U_{\\odot}=11.1$ $km/s$')
    ax1.legend(loc=0)
    ax1.grid(True)
    text( 6.2,  0.72, '$\\bar{U}_{\\odot}= \ %s$ $km/s$' % mean_U0,     verticalalignment='bottom', horizontalalignment='left',     color='b', fontsize=23)
    text( 6.2,  0.68, '$\\sigma= \ %s$ $km/s$' %sigma_U0,     verticalalignment='bottom', horizontalalignment='left',     color='b', fontsize=23) 
    text( 6.2,  0.58, '$\\bar{U}_{\\odot}= \ %s$ $km/s$' % round(mu_U0,3),     verticalalignment='bottom', horizontalalignment='left',     color='r', fontsize=23)
    text( 6.2,  0.54, '$\\sigma= \ %s$ $km/s$ ' % round(sigm_U0,3),  verticalalignment='bottom', horizontalalignment='left',     color='r', fontsize=23)
    #title("Results for $U_{\\odot}$ ")
    xlabel('$U_{\\odot}$ ($km/s$)')
    ylabel('Normed Counts')
    xlim([6,16])
    plt.savefig(path_plots+'Theta_0_100.pdf',bbox_inches='tight')
    
    #plot V_sun
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    matplotlib.rcParams.update({'font.size': 17})
    nbins=6
    a,bins_data,c=ax1.hist(sample_100['V0'],nbins,color='g', alpha=0.5, label="100 sources") #       
    b,bins_datab,cb=ax1.hist(sample_200['V0'],bins_data,color='b', alpha=0.5, label="200 sources") #       
    c,bins_datac,cc=ax1.hist(sample_complete['V0'],bins_data,color='r', alpha=0.5, label="Complete") #       
    ax1.axvline(x=14.6,color='k',ls='dashed',label='$V_{\\odot}=14.6$ $km/s$')
    ax1.legend(loc=0)
    ax1.grid(True)
    title("Results for $V_{\\odot}$ ")
    xlabel('$V_{\\odot}$ ($km/s$)')
    ylabel('Counts ($N$)')
    
    #plot W_sun
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    matplotlib.rcParams.update({'font.size': 17})
    nbins=6
    a,bins_data,c=ax1.hist(sample_100['W0'],nbins,color='g', alpha=0.5, label="100 sources") #       
    b,bins_datab,cb=ax1.hist(sample_200['W0'],bins_data,color='b', alpha=0.5, label="200 sources") #       
    c,bins_datac,cc=ax1.hist(sample_complete['W0'],bins_data,color='r', alpha=0.5, label="Complete") #       
    ax1.axvline(x=7.2,color='k',ls='dashed',label='$W_{\\odot}=7.2$ $km/s$')
    ax1.legend(loc=0)
    ax1.grid(True)
    title("Results for $W_{\\odot}$ ")
    xlabel('$W_{\\odot}$ ($km/s$)')
    ylabel('Counts ($N$)')
    '''
    
    #plot U_s and V_s
    mean_US, sigma_US = modval['US'], aprior_err['US'] # mean and standard deviation
    mean_VS, sigma_VS = modval['VS'], aprior_err['VS'] # mean and standard deviation
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    matplotlib.rcParams.update({'font.size': 17})
    nbins=10
    a,bins_data_us,c=ax1.hist(sample_100['US'],nbins,color='g', alpha=0.5,normed=True) #       
    a,bins_data_vs,c=ax1.hist(sample_100['VS'],nbins,color='b', alpha=0.5,normed=True) #       
    (mu_US, sigm_US) = norm.fit(sample_100['US'])
    y_us = mlab.normpdf(bins_data_us, mu_US, sigm_US)
    ax1.plot(bins_data_us, y_us, 'g--', linewidth=2,label='Gaussian Fitting')     
    (mu_VS, sigm_VS) = norm.fit(sample_100['VS'])
    y_vs = mlab.normpdf(bins_data_vs, mu_VS, sigm_VS)
    ax1.plot(bins_data_vs, y_vs, 'b--', linewidth=2)#,label='Gaussian Fitting')     
    ax1.legend(loc=0)
    ax1.grid(True)
    #Vs
    text( -4.7,  0.28, 'Initial Value',     verticalalignment='bottom', horizontalalignment='left',     color='k', fontsize=23)
    text( -4.7,  0.25, '$\\bar{U}_{s}= \ %s$ $km/s$' % mean_US,     verticalalignment='bottom', horizontalalignment='left',     color='k', fontsize=23)
    text( -4.7,  0.17, '$\\bar{U}_{s}= \ %s$ $km/s$' % round(mu_US,3),     verticalalignment='bottom', horizontalalignment='left',     color='g', fontsize=23)
    text( -4.7,  0.15, '$\\sigma= \ %s$ $km/s$ ' % round(sigm_US,3),  verticalalignment='bottom', horizontalalignment='left',     color='g', fontsize=23)
    #Us
    text( 2,  0.28, 'Initial Value',     verticalalignment='bottom', horizontalalignment='left',     color='k', fontsize=23)
    text( 2,  0.25, '$\\bar{V}_{s}= \ %s$ $km/s$' % mean_VS,     verticalalignment='bottom', horizontalalignment='left',     color='k', fontsize=23)
    text( 2,  0.17, '$\\bar{V}_{s}= \ %s$ $km/s$' % round(mu_VS,3),     verticalalignment='bottom', horizontalalignment='left',     color='b', fontsize=23)
    text( 2,  0.15, '$\\sigma= \ %s$ $km/s$ ' % round(sigm_VS,3),  verticalalignment='bottom', horizontalalignment='left',     color='b', fontsize=23)
    #
    legend = ax1.get_legend()
    legend.legendHandles[0].set_color('k')
    #title("Results of Velocity Source for 1000 Siumulations")
    xlabel('Mean Velocity Components Source ($U_s$,$V_s$) ($km/s$)')
    ylabel('Normed Counts')
    plt.savefig(path_outputs+'/UsVs_100.pdf',bbox_inches='tight')
    
    '''
    #plot U_s
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    matplotlib.rcParams.update({'font.size': 17})
    nbins=10
    a,bins_data,c=ax1.hist(sample_100['US'],nbins,color='g', alpha=0.5, label="100 sources") #       
    b,bins_datab,cb=ax1.hist(sample_200['US'],bins_data,color='b', alpha=0.5, label="200 sources") #       
    c,bins_datac,cc=ax1.hist(sample_complete['US'],bins_data,color='r', alpha=0.5, label="Complete") #       
    ax1.axvline(x=0.0,color='k',ls='dashed',label='$U_s=0.0$ $km/s$')
    ax1.legend(loc=0)
    ax1.grid(True)
    #title("Results of Velocity Source for 1000 Siumulations")
    xlabel('Mean Velocity Components Source $U_s$ ($km/s$)')
    ylabel('Counts ($N$)')
    
    #plot V_s
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    matplotlib.rcParams.update({'font.size': 17})
    nbins=6
    a,bins_data,c=ax1.hist(sample_100['VS'],nbins,color='g', alpha=0.5, label="100 sources") #       
    b,bins_datab,cb=ax1.hist(sample_200['VS'],bins_data,color='b', alpha=0.5, label="200 sources") #       
    c,bins_datac,cc=ax1.hist(sample_complete['VS'],bins_data,color='r', alpha=0.5, label="Complete") #       
    ax1.axvline(x=0.0,color='k',ls='dashed',label='$V_s=0.0$ $km/s$')
    ax1.legend(loc=0)
    ax1.grid(True)
    #title("Results of Velocity Source for 1000 Siumulations")
    xlabel('Mean Velocity Components Source $V_s$ ($km/s$)')
    ylabel('Counts ($N$)')
    '''
    sample_a=['R0', 'Th0', 'V0', 'VS']
    samples_a=[]
    for pair in itertools.combinations(sample_a,2):
        samples_a.append(pair)
    sample_b=['Th0', 'V0', 'VS']
    samples_b=[]
    for pair in itertools.combinations(sample_b,2):
        samples_b.append(pair)
    sample_c=['U0', 'US', 'VS']
    samples_c=[]
    for pair in itertools.combinations(sample_c,2):
        samples_c.append(pair)    

    for i in range(len(samples_a)):
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        matplotlib.rcParams.update({'font.size': 17})
        ax1.scatter(sample_100[samples_a[i][0]],sample_100[samples_a[i][1]],color='g', alpha=0.5, label='Sample 100') #        
        ax1.scatter(sample_200[samples_a[i][0]],sample_200[samples_a[i][1]],color='b', alpha=0.5, label='Sample 200') #        
        ax1.scatter(sample_complete[samples_a[i][0]],sample_complete[samples_a[i][1]],color='r', alpha=0.5, label='Complete Sample') #        
        ax1.grid(True)
        ax1.legend(scatterpoints=1,loc=4)
        #title("Correlation between")
        #ax1.set_xlabel('$\\theta_0 (km \ s^{-1})$')
        #ax1.set_ylabel('$R_{\\odot} (kpc)$')
        #ax1.set_ylim([8.0,8.6])
        #ax1.set_xlim([230.,250.])
        plt.savefig(path_outputs+'/Pearson_'+samples_a[i][0]+'_'+samples_a[i][1]+'.pdf',bbox_inches='tight')    
        
    for i in range(len(samples_b)):
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        matplotlib.rcParams.update({'font.size': 17})
        ax1.scatter(sample_100[samples_b[i][0]],sample_100[samples_b[i][1]],color='g', alpha=0.5, label='Sample 100') #        
        ax1.scatter(sample_200[samples_b[i][0]],sample_200[samples_b[i][1]],color='b', alpha=0.5, label='Sample 200') #        
        ax1.scatter(sample_complete[samples_b[i][0]],sample_complete[samples_b[i][1]],color='r', alpha=0.5, label='Complete Sample') #        
        ax1.grid(True)
        ax1.legend(scatterpoints=1,loc=4)
        plt.savefig(path_outputs+'/Pearson_'+samples_b[i][0]+'_'+samples_b[i][1]+'.pdf',bbox_inches='tight')    
        
    for i in range(len(samples_c)):
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        matplotlib.rcParams.update({'font.size': 17})
        ax1.scatter(sample_100[samples_c[i][0]],sample_100[samples_c[i][1]],color='g', alpha=0.5, label='Sample 100') #        
        ax1.scatter(sample_200[samples_c[i][0]],sample_200[samples_c[i][1]],color='b', alpha=0.5, label='Sample 200') #        
        ax1.scatter(sample_complete[samples_c[i][0]],sample_complete[samples_c[i][1]],color='r', alpha=0.5, label='Complete Sample') #        
        ax1.grid(True)
        ax1.legend(scatterpoints=1,loc=4)        
        plt.savefig(path_outputs+'/Pearson_'+samples_c[i][0]+'_'+samples_c[i][1]+'.pdf',bbox_inches='tight')    
    
    '''
    #Correlation theta_0 vs R0
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    matplotlib.rcParams.update({'font.size': 17})
    ax1.scatter(sample_100['Th0'],sample_100['R0'],color='g', alpha=0.5, label='Sample 100') #        
    ax1.scatter(sample_200['Th0'],sample_200['R0'],color='b', alpha=0.5, label='Sample 200') #        
    ax1.scatter(sample_complete['Th0'],sample_complete['R0'],color='r', alpha=0.5, label='Complete Sample') #        
    ax1.grid(True)
    ax1.legend(scatterpoints=1,loc=4)
    #title("Correlation between")
    ax1.set_xlabel('$\\theta_0 (km \ s^{-1})$')
    ax1.set_ylabel('$R_{\\odot} (kpc)$')
    ax1.set_ylim([8.0,8.6])
    ax1.set_xlim([230.,250.])
    
    #Correlation theta_0 vs v0
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    matplotlib.rcParams.update({'font.size': 23})
    ax1.scatter(sample_100['Th0'],sample_100['V0'],color='g', alpha=0.5, label='Sample 100') #        
    ax1.scatter(sample_200['Th0'],sample_200['V0'],color='b', alpha=0.5, label='Sample 200') #        
    ax1.scatter(sample_complete['Th0'],sample_complete['V0'],color='r', alpha=0.5, label='Complete Sample') #        
    ax1.grid(True)
    ax1.legend(scatterpoints=1)
    #title("Correlation between")
    ax1.set_xlabel('$\\theta_0 (km \ s^{-1})$')
    ax1.set_ylabel('$V_{\\odot} (km \ s^{-1})$')
    ax1.set_ylim([10.,20.])
    ax1.set_xlim([230.,255.])
    
    
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    matplotlib.rcParams.update({'font.size': 23})
    ax1.scatter(sample_100['Th0'],sample_100['VS'],color='g', alpha=0.5, label='Sample 100') #        
    #ax1.scatter(sample_200['Th0'],sample_200['V0'],color='b', alpha=0.5, label='Sample 200') #        
    ax1.scatter(sample_complete['Th0'],sample_complete['VS'],color='r', alpha=0.5, label='Complete Sample') #        
    ax1.grid(True)
    ax1.legend(scatterpoints=1)
    #title("Correlation between")
    ax1.set_xlabel('$\\theta_0 (km \ s^{-1})$')
    ax1.set_ylabel('$V_{s} (km \ s^{-1})$')
    
    
    
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    matplotlib.rcParams.update({'font.size': 23})
    ax1.scatter(sample_100['U0'],sample_100['VS'],color='g', alpha=0.5, label='Sample 100') #        
    ax1.scatter(sample_200['U0'],sample_200['V0'],color='b', alpha=0.5, label='Sample 200') #        
    ax1.scatter(sample_complete['U0'],sample_complete['VS'],color='r', alpha=0.5, label='Complete Sample') #        
    ax1.grid(True)
    ax1.legend(scatterpoints=1)
    #title("Correlation between")
    ax1.set_xlabel('$U_s (km \ s^{-1})$')
    ax1.set_ylabel('$V_{s} (km \ s^{-1})$')
    '''
    os.chdir(path_origin)
    print('Plots saved!')
    print('Done with simulation, fitting and analysis')
    print('Thanks for using me!, please remember to reference us')
    plt.close('all')
else:
    os.chdir(path_origin)
    print('No plots!')
    print('Done with simulation, fitting and analysis')
    print('Thanks for using me!, please remember to reference us')

