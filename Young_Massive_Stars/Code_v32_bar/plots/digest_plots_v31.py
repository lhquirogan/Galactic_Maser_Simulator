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
#Colors
colork2=['k','g','r','b','sienna','m','c','orange','chocolate','firebrick','rosybrown','crimson',
               'cadetblue','tomato']
#------------------------------------------------------------------------------------
#Question Plots
plots_confirmation = raw_input("Do you want plots of parameters distributions? (yes or no):")
plots_confirmation_nsources = raw_input("Do you want plots of different number of sources? (yes or no):")
plots_confirmation_merged = raw_input("Do you want merged plots sources and parameters? (yes or no):")
plots_confirmation_pearson = raw_input("Do you want Pearson evolution plots? (yes or no):")


##########################################################################################################################
###################################     Galac. Paramters Histograms BeSSeL      ##########################################
##########################################################################################################################
if (plots_confirmation=='yes'):
    # Eight Plots
    fig1 = plt.figure(figsize=(15,10))
    matplotlib.rcParams.update({'font.size': 16})
    zz=0 #which sample?
    label_height=450
    ax = fig1.add_subplot(111)
          
    ax1 = fig1.add_subplot(241)
    mean_R0,sigmas_R0_Reid = modval['R0'],final_unc['R0']  # initial values
    (mu_R0, sigmaR0) = norm.fit(dic[dic_table_names[zz][2:]]['R0'])
    normed_R0 = [x for x in dic[dic_table_names[zz][2:]]['R0']]    
    (mu_normed_R0, sigma_normed_R0) = norm.fit(normed_R0)
    a,bins_data,c=ax1.hist(normed_R0,bins=7,color=colork2[1], alpha=0.7) #     
    area = sum(numpy.diff(bins_data)*a)
    xmin,xmax=ax1.get_xlim()
    #xmin,xmax=ax5.get_xlim()
    x = np.linspace(xmin, xmax, 100)
    y=mlab.normpdf(x, mu_normed_R0, sigma_normed_R0)*area
    y1=mlab.normpdf(x, mean_R0, sigmas_R0_Reid)*area
    ax1.plot(x,y1, color='grey',linewidth=2.5)    
    ax1.plot(x,y, linestyle='--',color='k',linewidth=2.5)
    ax1.fill_between(x, 0, y1, facecolor='grey', alpha=0.4)
    ax1.set_yticks([0,100,200,300,400])
    #ax1.set_xticks([0.96,0.98,1.0,1.02,1.04])
    #ax1.set_xticks([0.7,0.85,1.0,1.15,1.3])
    ax1.set_xlim([8.05,8.64]) 
    ax1.set_xticks([8.1,8.2,8.3,8.4,8.5,8.6])
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
    a,bins_data,c=ax2.hist(normed_Th0,bins=10,color=colork2[2], alpha=0.7) # 
    area = sum(numpy.diff(bins_data)*a) 
    xmin,xmax=ax2.get_xlim()
    #xmin,xmax=ax6.get_xlim()   
    x = np.linspace(xmin, xmax, 100)
    y=mlab.normpdf(x, mu_normed_Th0, sigma_normed_Th0)*area
    y1=mlab.normpdf(x, mean_Th0, sigmas_Th0_Reid)*area
    ax2.plot(x,y1, color='grey',linewidth=2.5) 
    ax2.plot(x,y, linestyle='--',color='k',linewidth=2.0)     
    ax2.fill_between(x, 0, y1, facecolor='grey', alpha=0.4)    
    #ax2.set_xticks([0.6,0.8,1.0,1.2,1.4])
    ax2.set_xlim([226,255]) 
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
    a,bins_data,c=ax3.hist(normed_dTd,bins=8,color=colork2[8], alpha=0.7) #     
    area = sum(numpy.diff(bins_data)*a)
    xmin,xmax=ax3.get_xlim()    
    x = np.linspace(xmin, xmax, 100)
    y1=mlab.normpdf(x, mean_dTd, sigmas_dTd_Reid)*area
    ax3.plot(x,y1, color='grey',linewidth=2.5)     
    ax3.fill_between(x, 0, y1, facecolor='grey', alpha=0.4)
    y=mlab.normpdf(x, mu_normed_dTd, sigma_normed_dTd)*area
    ax3.plot(x,y, linestyle='--',color='k',linewidth=2.0)
    #ax3.set_xlim([-12,15]) 
    ax3.set_xlim([-5,4])     
    ax3.set_xticks([-4,-2,0,2,4])  
    ax3.grid(True)
    plt.setp(ax3.get_yticklabels(), visible=False)     
    ax3.text(mu_mean_dTd, label_height+8, '$\\frac{d \\Theta}{dR}$ $\\rm{(km/s \\cdot kpc)}$', verticalalignment='center', horizontalalignment='center', color=colork2[8], fontsize=25)
    #ax3.axvline(x=mean_dTd,linewidth=3,color='k',ls='dotted')
          
    ax4 = fig1.add_subplot(244,sharey=ax1)
    mean_U0,sigmas_U0_Reid = modval['U0'],final_unc['U0'] # initial values    
    (mu_U0, sigmaU0) = norm.fit(dic[dic_table_names[zz][2:]]['U0'])
    normed_U0 = [x for x in dic[dic_table_names[zz][2:]]['U0']]
    (mu_normed_U0, sigma_normed_U0) = norm.fit(normed_U0)
    a,bins_data,c=ax4.hist(normed_U0,bins=6,color=colork2[5], alpha=0.7) # 
    area = sum(numpy.diff(bins_data)*a)
    xmin,xmax=ax4.get_xlim()
    x = np.linspace(xmin, xmax, 100)
    y=mlab.normpdf(x, mu_normed_U0, sigma_normed_U0)*area
    y1=mlab.normpdf(x, mean_U0, sigmas_U0_Reid)*area
    ax4.plot(x,y1, color='grey',linewidth=2.5)  
    ax4.plot(x,y, linestyle='--',color='k',linewidth=2.0)    
    ax4.fill_between(x, 0, y1, facecolor='grey', alpha=0.4)    
    ax4.set_yticks([0,100,200,300,400])
    ax4.set_ylim([0,499])        
    ax4.set_xlim([6.5,14.01])     
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
    a,bins_data,c=ax5.hist(normed_V0,bins=12,color=colork2[6], alpha=0.7) #     
    area = sum(numpy.diff(bins_data)*a)
    xmin,xmax=ax5.get_xlim()    
    x = np.linspace(xmin, xmax, 100)
    y=mlab.normpdf(x, mu_normed_V0, sigma_normed_V0)*area
    y1=mlab.normpdf(x, mean_V0, sigmas_V0_Reid)*area
    ax5.plot(x,y1, color='grey',linewidth=2.5)     
    ax5.plot(x,y, linestyle='--',color='k',linewidth=2.0)    
    ax5.fill_between(x, 0, y1, facecolor='grey', alpha=0.4)    
    ax5.set_xlim([5,26.5935]) 
    ax5.set_xticks([6,9,12,15,18,21,24])  
    ax5.grid(True)
    ax5.text( mu_V0,  label_height, '$V_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[6], fontsize=25)        
    #ax5.axvline(x=mean_V0,linewidth=3,color='k',ls='dotted')

        
    ax6 = fig1.add_subplot(246,sharey=ax5)    
    mean_W0,sigmas_W0_Reid = modval['W0'],final_unc['W0'] # initial values
    (mu_W0, sigmaW0) = norm.fit(dic[dic_table_names[zz][2:]]['W0'])
    normed_W0 = [x for x in dic[dic_table_names[zz][2:]]['W0']]
    (mu_normed_W0, sigma_normed_W0) = norm.fit(normed_W0)
    a,bins_data,c=ax6.hist(normed_W0,bins=7,color=colork2[7], alpha=0.7) # 
    area = sum(numpy.diff(bins_data)*a)
    xmin,xmax=ax6.get_xlim()
    x = np.linspace(xmin, xmax, 100)
    y=mlab.normpdf(x, mu_normed_W0, sigma_normed_W0)*area
    y1=mlab.normpdf(x, mean_W0, sigmas_W0_Reid)*area
    ax6.plot(x,y1, color='grey',linewidth=2.5,label='Results Model A5')     
    ax6.plot(x,y, linestyle='--',color='k',linewidth=2.0)    
    ax6.fill_between(x, 0, y1, facecolor='grey', alpha=0.4)    
    ax6.set_xlim([6.5,10.5])     
    ax6.set_xticks([7,8,9,10])
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
    a,bins_data,c=ax7.hist(normed_US,bins=6,color=colork2[3], alpha=0.7) # 
    area = sum(numpy.diff(bins_data)*a)
    xmin,xmax=ax7.get_xlim()
    #xmin,xmax=ax7.get_xlim()
    x = np.linspace(xmin, xmax, 100)
    y=mlab.normpdf(x, mu_normed_US, sigma_normed_US)*area
    y1=mlab.normpdf(x, mean_US, sigmas_US_Reid)*area
    ax7.plot(x,y1, color='grey',linewidth=2.5) 
    ax7.plot(x,y, linestyle='--',color='k',linewidth=2.0,label='Gaussian Fitting')     
    ax7.fill_between(x, 0, y1, facecolor='grey', alpha=0.4)    
    #ax7.set_xticks([-0.5,0.0,0.5,1.0,1.5,2.0,2.5])
    #ax7.set_xticks([-1.0,0.0,1.0,2.0,3.0])
    ax7.set_xlim([-1.9,7.9])     
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
    a,bins_data,c=ax8.hist(normed_VS,bins=9,color=colork2[4], alpha=0.7) # 
    area = sum(numpy.diff(bins_data)*a)
    xmin,xmax=ax8.get_xlim()
    #xmin,xmax=ax8.get_xlim()
    x = np.linspace(xmin, xmax, 100)
    y=mlab.normpdf(x, mu_normed_VS, sigma_normed_VS)*area
    y1=mlab.normpdf(x, mean_VS, sigmas_VS_Reid)*area
    ax8.plot(x,y1, color='grey',linewidth=2.5)     
    ax8.plot(x,y, linestyle='--',color='k',linewidth=2.0)    
    ax8.fill_between(x, 0, y1, facecolor='grey', alpha=0.4)    
    #ax8.set_xticks([-14,-9,-4,1,6,11,16])
    ax8.set_ylim([0,499])    
    ax8.set_xlim([-15.54688889,15.9])     
    ax8.set_xticks([-12,-8,-4,0,4,8,12])
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
        
    #plt.savefig(path_outputs+'/eight.pdf',bbox_inches='tight')
    
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
    b_w_n = raw_input("BeSSeL or whole galaxy? (b or w or none):")
    f = raw_input("Quality of data (none, f01, f00):")
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
                        
    #Run only once if you want to include 16 and 103 sources better results
    N_bri_sources.append(16)
    normed_R0[0]= 8.34 #1000 galaxies results
    normed_errR0[0]=  0.08  #1000 galaxies results
    normed_R0.append(8.34)   #16 sources 1000 gaalxies 
    normed_errR0.append(0.24) #16 sources 1000 gaalxies 

    normed_Th0[0]= 239.8 #1000 galaxies results
    normed_errTh0[0]=  4.2  #1000 galaxies results
    normed_Th0.append(240.32)   #16 sources 1000 gaalxies 
    normed_errTh0.append(9.11) #16 sources 1000 gaalxies     
    
    normed_US[0]= 2.9 #1000 galaxies results
    normed_errUS[0]=  1.6  #1000 galaxies results    
    normed_US.append(2.99)  #16 sources 1000 gaalxies 
    normed_errUS.append(3.66) #16 sources 1000 gaalxies     

    normed_VS[0]= 1.0 #1000 galaxies results
    normed_errVS[0]=  4.9  #1000 galaxies results    
    normed_VS.append(1.21)   #16 sources 1000 gaalxies 
    normed_errVS.append(9.78) #16 sources 1000 gaalxies         

    normed_Usun[0]= 10.5 #1000 galaxies results
    normed_errUsun[0]=  1.2  #1000 galaxies results    
    normed_Usun.append(10.79)   #16 sources 1000 gaalxies 
    normed_errUsun.append(1.58) #16 sources 1000 gaalxies             
       
    normed_Vsun[0]= 15.5 #1000 galaxies results
    normed_errVsun[0]=  3.1  #1000 galaxies results        
    normed_Vsun.append(15.20)   #16 sources 1000 gaalxies 
    normed_errVsun.append(2.54) #16 sources 1000 gaalxies

    normed_Wsun[0]= 8.7 #1000 galaxies results
    normed_errWsun[0]=  0.5  #1000 galaxies results        
    normed_Wsun.append(8.13)  #16 sources 1000 gaalxies 
    normed_errWsun.append(0.77) #16 sources 1000 gaalxies
    
    normed_dTd[0]= -0.2 #1000 galaxies results
    normed_errdTd[0]=  1.1  #1000 galaxies results    
    normed_dTd.append(-0.03)   #16 sources 1000 gaalxies 
    normed_errdTd.append(3.24) #16 sources 1000 gaalxies                 
                                                                                                
    #Plot
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
            ax1.text( xpos_label-20, 9, '$R_0$ $\\rm{(kpc)}$', verticalalignment='center', horizontalalignment='center', color=colork2[1], fontsize=27)                    
        if(f=='f01'):
            y_tickts=np.arange(8.0,8.8,0.1) #f01
            ax1.text( xpos_label+150, 8.68, '$R_0$ $\\rm{(kpc)}$', verticalalignment='center', horizontalalignment='center', color=colork2[1], fontsize=27)         #f01        
            ax1.set_yticks(y_tickts)  
            ax1.set_ylim([8.0,8.79]) #f01
    if(b_w_n=='w'):
        if(f=='none'):        
            ax1.text( xpos_label-20, 9.1, '$R_0$ $\\rm{(kpc)}$', verticalalignment='center', horizontalalignment='center', color=colork2[1], fontsize=27)                
            ax1.set_ylim([8.01,9.3])  
        if(f=='f01'):            
            y_tickts=np.arange(8.0,8.8,0.1) #f01
            ax1.text( xpos_label+150, 8.68, '$R_0$ $\\rm{(kpc)}$', verticalalignment='center', horizontalalignment='center', color=colork2[1], fontsize=27)         #f01        
            ax1.set_yticks(y_tickts)  
            ax1.set_ylim([8.0,8.79]) #f01
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
            y_tickts=np.arange(230,300,10)
            ax2.text( xpos_label-20, 286, '$\\Theta_0$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[2], fontsize=27)                         
            ax2.set_yticks(y_tickts)    
            ax2.set_ylim([230,296])
        if(f=='f01'):            
            y_tickts=np.arange(230,270,5) #f01
            ax2.text( xpos_label+140, 265, '$\\Theta_0$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[2], fontsize=27)                 #f01
            ax2.set_yticks(y_tickts)    
            ax2.set_ylim([230,270.3]) # f01
    if(b_w_n=='w'):
        if(f=='none'):        
            y_tickts=np.arange(230,300,10)
            ax2.text( xpos_label-20, 285, '$\\Theta_0$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[2], fontsize=27)                                 
            ax2.set_yticks(y_tickts) 
            ax2.set_ylim([225,296])
        if(f=='f01'):                    
            y_tickts=np.arange(230,270,5) #f01
            ax2.text( xpos_label+140, 265, '$\\Theta_0$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[2], fontsize=27)                 #f01
            ax2.set_yticks(y_tickts)    
            ax2.set_ylim([230,270.3]) # f01
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
            y_tickts=np.arange(0,30,5)
            ax3.set_yticks(y_tickts)    
            ax3.text( xpos_label+5, 25, '$\\frac{d\\Theta}{dR}$ $\\rm{(km/s \\cdot kpc)}$', verticalalignment='center', horizontalalignment='center', color=colork2[8], fontsize=27)                     
        if(f=='f01'):                    
            y_tickts=np.arange(-3,5.0,1) #f01
            ax3.set_yticks(y_tickts)    
            ax3.text( xpos_label+50, 4, '$\\frac{d\\Theta}{dR}$ $\\rm{(km/s \\cdot kpc)}$', verticalalignment='center', horizontalalignment='center', color=colork2[8], fontsize=27)                     #f01
            ax3.set_ylim([-3.5,5.0]) #f01
    if(b_w_n=='w'):
        if(f=='none'):        
            ax3.text( xpos_label+20, 27.5, '$\\frac{d\\Theta}{dR}$ $\\rm{(km/s \\cdot kpc)}$', verticalalignment='center', horizontalalignment='center', color=colork2[8], fontsize=27)                     
            ax3.set_ylim([-4,32])
        if(f=='f01'):                    
            y_tickts=np.arange(-3,5.0,1) #f01
            ax3.set_yticks(y_tickts)    
            ax3.text( xpos_label+50, 4, '$\\frac{d\\Theta}{dR}$ $\\rm{(km/s \\cdot kpc)}$', verticalalignment='center', horizontalalignment='center', color=colork2[8], fontsize=27)                     #f01
            ax3.set_ylim([-3.5,5.0]) #f01
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
            y_tickts=np.arange(8,15,1)
            ax4.set_yticks(y_tickts)  
            ax4.text( xpos_label-35, 14, '$U_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[5], fontsize=27)                         
        if(f=='f01'):                    
            y_tickts=np.arange(9,13,0.5) #f01
            ax4.set_yticks(y_tickts)    
            ax4.text( xpos_label+125, 12.2, '$U_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[5], fontsize=27)                 #f01
            ax4.set_ylim([8.7,12.7]) #f01
    if(b_w_n=='w'):
        if(f=='none'):        
            ax4.text( xpos_label-15, 13.5, '$U_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[5], fontsize=27)                         
            ax4.set_ylim([8.5,14.2]) #f01
        if(f=='f01'):                    
            y_tickts=np.arange(9,13,0.5) #f01
            ax4.set_yticks(y_tickts)    
            ax4.text( xpos_label+125, 12.2, '$U_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[5], fontsize=27)                 #f01
            ax4.set_ylim([8.7,12.7]) #f01
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
            ax5.text( xpos_label+40,10, '$V_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[6], fontsize=27)                      
        if(f=='f01'):                    
            ax5.text( xpos_label+125,21.2, '$V_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[6], fontsize=27)                 #f01
            ax5.set_ylim([8.1,23]) 
        #p1 = Rectangle((0, 0), 1, 1, fc="grey",alpha=0.3)
        #ax5.legend([p1],['Complete Simulation'],bbox_to_anchor=(1.3, -0.17),numpoints=1,ncol=2)                    
        #ax5.legend(bbox_to_anchor=(0.2, -0.17),numpoints=1,ncol=2) #f01        
    if(b_w_n=='w'):
        if(f=='none'):        
            ax5.text( xpos_label+20,9, '$V_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[6], fontsize=27)                                    
            ax5.set_ylim([7.0,23])    
        if(f=='f01'):                    
            ax5.text( xpos_label+125,21.2, '$V_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[6], fontsize=27)                 #f01
            ax5.set_ylim([8.1,23]) 
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
            y_tickts=np.arange(7.4,10.0,0.4)
            ax6.set_yticks(y_tickts)        
            ax6.text( xpos_label+50,7.7, '$W_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[7], fontsize=27)               
            ax6.set_ylim([7.3,10])
        if(f=='f01'):                    
            ax6.text( xpos_label+120,9.65, '$W_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[7], fontsize=27)               #f01
            ax6.set_ylim([7.3,10])            
        #ax6.legend(bbox_to_anchor=(0.35, -0.17),numpoints=1,ncol=2)            
        #ax6.legend(bbox_to_anchor=(0.0, -0.17),numpoints=1,ncol=2)
    if(b_w_n=='w'):
        if(f=='none'):
            ax6.set_ylim([7.3,10])                    
            ax6.text( xpos_label+30,7.6, '$W_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[7], fontsize=27)               
        if(f=='f01'):                    
            ax6.text( xpos_label+120,9.72, '$W_{\\odot}$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[7], fontsize=27)               #f01
            ax6.set_ylim([7.3,10])
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
            y_tickts=np.arange(-6,8,2)
            ax7.text( xpos_label-20, -6, '$\\bar{U}_s$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[3], fontsize=27)               
            ax7.set_yticks(y_tickts)
            ax7.set_ylim([-8,7])
        if(f=='f01'):                    
            y_tickts=np.arange(0,7,1) #f01
            ax7.text( xpos_label+140, 6, '$\\bar{U}_s$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[3], fontsize=27)               #f01
            ax7.set_yticks(y_tickts)  
            ax7.set_ylim([-1,7]) #f01        
        p1 = Rectangle((0, 0), 1, 1, fc="grey",alpha=0.3)
        ax7.legend([p1,p2,p3,p4],['Complete Simulation','Model A5','Fit 3','Initial Values'],numpoints=1,ncol=4,bbox_to_anchor=(1.6, -0.17))                                
    if(b_w_n=='w'):
        if(f=='none'):        
            y_tickts=np.arange(0,7,1) #f01
            ax7.text( xpos_label+90, 0, '$\\bar{U}_s$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[3], fontsize=27)               #f01
            ax7.set_yticks(y_tickts)  
            ax7.set_ylim([-1,7]) #f01        
        if(f=='f01'):                    
            y_tickts=np.arange(0,7,1) #f01
            ax7.text( xpos_label+140, 6, '$\\bar{U}_s$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[3], fontsize=27)               #f01
            ax7.set_yticks(y_tickts)  
            ax7.set_ylim([-1,7]) #f01  
        ax7.fill_between(X, Y1,Y2,color='grey',alpha=0.3)
            
    #Vs
    ax8 = fig1.add_subplot(248,sharex=ax4)#,sharey=ax5)
    normed_VS = [x for x in means['VS']]
    normed_errVS = [x for x in errors['VS']]
    ax8.errorbar(N_bri_sources,normed_VS,yerr=normed_errVS,color=colork2[4],linewidth=2,fmt="^",capsize=4, elinewidth=4,markersize=8)
    ax8.errorbar(N_BeSSeL, modval['VS'],yerr=final_unc['VS'],fmt="*",color='k',markersize=12,capsize=4, elinewidth=4)
    ax8.errorbar(N_BeSSeL_past,-14.7,yerr=1.8,fmt="*",color='cadetblue',label='Fit 3',markersize=12,capsize=4, elinewidth=4)
    #ax8.set_ylim(ylim_4)
    #ax8.set_yticks(y_tickts_4)
    ax8.set_xlim(xlim)
    ax8.set_xticks(x_tickts)       
    #ax8.yaxis.tick_right()
    #ax8.yaxis.set_label_position('right') 
    Y2=0*X+((mean_complete['VS'][0])+(errors_complete['VS'][0]))
    Y1=0*X+((mean_complete['VS'][0])-(errors_complete['VS'][0]))
    ax8.grid(True)
    #ax8.set_yscale("log", nonposy='clip')
    #ax8.tick_params(axis='y',labelsize=14)    
    if(b_w_n=='b'):
        if(f=='none'):        
            y_tickts=np.arange(-50,30,10)
            ax8.set_yticks(y_tickts)  
            ax8.text(xpos_label-20,-50, '$\\bar{V}_s$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[4], fontsize=27)                       
            ax8.set_ylim([-60,13])            
        if(f=='f01'):                    
            y_tickts=np.arange(-9,12,3) #f01
            ax8.set_yticks(y_tickts)    
            ax8.text(xpos_label+150,9, '$\\bar{V}_s$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[4], fontsize=27)                       #f01
            ax8.set_ylim([-8.9,11.5]) #f01
        #ax8.legend(bbox_to_anchor=(0.2, -0.17),numpoints=1,ncol=2)                    
    if(b_w_n=='w'):
        if(f=='none'):        
            y_tickts=np.arange(-60,20,10)
            ax8.set_yticks(y_tickts)            
            ax8.text(xpos_label-60,-55, '$\\bar{V}_s$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[4], fontsize=27)                       
            ax8.set_ylim([-67,15])
        if(f=='f01'):         
            y_tickts=np.arange(-12,12,3) #f01  
            ax8.set_yticks(y_tickts)            
            ax8.text(xpos_label+150,9, '$\\bar{V}_s$ $\\rm{(km/s)}$', verticalalignment='center', horizontalalignment='center', color=colork2[4], fontsize=27)                       #f01
            ax8.set_ylim([-11.8,11.5]) #f01
        #ax8.legend(bbox_to_anchor=(0.7, -0.17),numpoints=1,ncol=2)        
        ax8.fill_between(X, Y1,Y2,color='grey',alpha=0.3)
    ax8.axhline(y=modval['VS'],color='k',linestyle='dashed',linewidth=3)    

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
    
    #plt.savefig(path_outputs+'/whole_sigma15%.pdf',bbox_inches='tight')      
    #plt.savefig(path_outputs+'/bessel_nogauss.pdf',bbox_inches='tight')          
    #plt.savefig(path_outputs+'/bessel_gauss_distance.pdf',bbox_inches='tight')   
    #plt.savefig(path_outputs+'/bessel_gdist_betterunc.pdf',bbox_inches='tight')   
    #plt.savefig(path_outputs+'/bessel_gdist_nicess10.pdf',bbox_inches='tight')  
    #plt.savefig(path_outputs+'/bessel_gpara_-100%new.pdf',bbox_inches='tight')                     
    #plt.savefig(path_outputs+'/bessel_gpara_20%new.pdf',bbox_inches='tight')                     
    #Printing Results
    modval2=modval.copy()
    parameters=modval2.keys()
    for i in range(len(parameters)):
        modval2[parameters[i]]=[]
        modval2[parameters[i]].append(modval[parameters[i]])        
        #modval2[parameters[i]].append(means[i])
        #modval2[parameters[i]].append(sigmas[i])
    para=['$R_0 \, \\rm{(kpc)}$','$U_{\sun} \, \\rm{(km \, s^{-1})}$','$U_s  \, \\rm{(km \, s^{-1})} $','d$\\Theta$/d$R  \, \\rm{(km \, s^{-1} \, kpc^{-1})}$','$V_{\sun}  \, \\rm{(km \, s^{-1})} $','$V_s  \, \\rm{(km \, s^{-1})} $','$W_{\sun}  \, \\rm{(km \, s^{-1})} $','$\\Theta_0  \, \\rm{(km \, s^{-1})} $']        
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
###################################     Merged Plot with subplot      ####################################################
##########################################################################################################################
if (plots_confirmation_merged=='yes'):     
    from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
    from mpl_toolkits.axes_grid1.inset_locator import mark_inset 
    matplotlib.rcParams.update({'font.size': 16})
    fig1 = plt.figure(figsize=(20,20))
    ax = fig1.add_subplot(111)   
    #fig1,ax1 = plt.subplots(figsize=(10,10))
    X=numpy.linspace(0,1300,100)
    
    ax1 = fig1.add_subplot(221)    
    jj=0 #R0
    normed_R0 = [x for x in means[fitpars[jj]]]
    normed_errR0 = [x for x in errors[fitpars[jj]]]
    ax1.errorbar(N_bri_sources,normed_R0,yerr=normed_errR0,color=colork2[1],linewidth=2,fmt="^",label='$R_0$ & $\\Delta R_0$',capsize=4, elinewidth=4,markersize=8)
    ax1.errorbar(107,8.34,yerr=0.16,fmt="*",color='k',label='Current BeSSeL',markersize=12,capsize=4, elinewidth=4)
    ax1.set_xlim([90,520])
    ax1.set_ylim([8.1,9.2])
    Y2=0*X+((mean_complete[fitpars[jj]][0])+(errors_complete[fitpars[jj]][0]))
    Y1=0*X+((mean_complete[fitpars[jj]][0])-(errors_complete[fitpars[jj]][0]))
    ax1.fill_between(X, Y1,Y2,color='grey',alpha=0.3)
    ax1.grid(True)
    ax1.set_ylabel('Galactocentric \n Distance ($\\rm{kpc}$)')
    plt.setp(ax1.get_xticklabels(), visible=False) 
    ax1.text( 450, 8.2, '$R_0$', verticalalignment='bottom', horizontalalignment='left', color=colork2[1], fontsize=27)        
    
    jj=1 #Th0
    ax2 = fig1.add_subplot(222)
    normed_Th0 = [x for x in means[fitpars[jj]]]
    normed_errTh0 = [x  for x in errors[fitpars[jj]]]
    ax2.errorbar(N_bri_sources,normed_Th0,yerr=normed_errTh0,color=colork2[2],linewidth=2,fmt="^",label='$R_0$ & $\\Delta R_0$',capsize=4, elinewidth=4,markersize=8)
    ax2.errorbar(107,240,yerr=8,fmt="*",color='k',markersize=12,capsize=4, elinewidth=4)
    ax2.set_xlim([90,520])
    #ax1.set_ylim([8.1,9.2])  
    Y2=0*X+((mean_complete[fitpars[jj]][0])+(errors_complete[fitpars[jj]][0]))
    Y1=0*X+((mean_complete[fitpars[jj]][0])-(errors_complete[fitpars[jj]][0]))
    ax2.fill_between(X, Y1,Y2,color='grey',alpha=0.3)
    ax2.grid(True)
    plt.setp(ax2.get_xticklabels(), visible=False) 
    ax2.set_ylabel('Solar rotation \n speed ($\\rm{km \, s^{-1}}$)')
    ax2.text( 450, 232, '$\\Theta_0$', verticalalignment='bottom', horizontalalignment='left', color=colork2[2], fontsize=27)        
    
    jj=7 #U_s
    ax3 = fig1.add_subplot(223,sharex=ax1)
    normed_US = [x for x in means[fitpars[jj]]]
    normed_errUS = [x for x in errors[fitpars[jj]]]
    ax3.errorbar(N_bri_sources,normed_US,yerr=normed_errUS,color=colork2[3],linewidth=2,fmt="^",label='$R_0$ & $\\Delta R_0$',capsize=4, elinewidth=4,markersize=8)
    ax3.errorbar(107, 2.9,yerr=2.1,fmt="*",color='k',markersize=12,capsize=4, elinewidth=4)
    ax3.set_xlim([90,520]) 
    #ax1.set_ylim([8.1,9.2])
    Y2=0*X+((mean_complete[fitpars[jj]][0]/modval[fitpars[jj]])+(errors_complete[fitpars[jj]][0]/modval[fitpars[jj]]))
    Y1=0*X+((mean_complete[fitpars[jj]][0]/modval[fitpars[jj]])-(errors_complete[fitpars[jj]][0]/modval[fitpars[jj]]))
    ax3.fill_between(X, Y1,Y2,color='grey',alpha=0.3)
    ax3.grid(True)
    #plt.setp(ax3.get_xticklabels(), visible=False) 
    ax3.set_ylabel('Radial average \n peculiar motion ($\\rm{km \, s^{-1}}$)')
    ax3.set_xlabel('Number of Sources')
    ax3.text( 450, 3, '$U_s$', verticalalignment='bottom', horizontalalignment='left', color=colork2[3], fontsize=27)        
    
    jj=6 #V_s
    ax4 = fig1.add_subplot(224,sharex=ax2)
    normed_VS = [x  for x in means[fitpars[jj]]]
    normed_errVS = [x  for x in errors[fitpars[jj]]]
    ax4.errorbar(N_bri_sources,normed_VS,yerr=normed_errVS,color=colork2[4],linewidth=2,fmt="^",label='$R_0$ & $\\Delta R_0$',capsize=4, elinewidth=4,markersize=8)
    ax4.errorbar(107, 1.6, yerr=6.8,fmt="*",color='k',markersize=12,capsize=4, elinewidth=4)
    Y2=0*X+((mean_complete[fitpars[jj]][0]/modval[fitpars[jj]])+(errors_complete[fitpars[jj]][0]/modval[fitpars[jj]]))
    Y1=0*X+((mean_complete[fitpars[jj]][0]/modval[fitpars[jj]])-(errors_complete[fitpars[jj]][0]/modval[fitpars[jj]]))
    ax4.fill_between(X, Y1,Y2,color='grey',alpha=0.3)
    ax4.set_xlim([90,520]) 
    ax4.grid(True)
    #ax4.tick_params(axis='y',labelsize=14)    
    ax4.set_ylabel('Tangential average \n peculiar motion ($\\rm{km \, s^{-1}}$)')
    ax4.set_xlabel('Number of Sources')
    ax4.text( 450, -3, '$V_s$', verticalalignment='bottom', horizontalalignment='left', color=colork2[4], fontsize=27)        
    
    
    plt.setp(ax.get_yticklabels(), visible=False)  
    plt.setp(ax.get_xticklabels(), visible=False)          
    
    
    zz=0 #which sample?
    mean_R0 = modval['R0'] # initial values
    (mu_R0, sigma_R0) = norm.fit(dic[dic_table_names[zz][2:]]['R0'])
    normed_R0 = [x for x in dic[dic_table_names[zz][2:]]['R0']]    
    # connecting lines between the bbox and the inset axes area
    #ax1.arrow(N_bri_sources[0],normed_R0[0], 55, 0.23, head_length=10, head_width=0.03, fc='g', ec='g')
    #ax1.arrow(107,8.34, 75, 0.23, head_length=10, head_width=0.03, fc='k', ec='k')
    # this is an inset axes over the main axes
    l1 = plt.axes([.17, .7, .17, .17])
    a,bins_data,c=plt.hist(normed_R0,bins=8,color=colork2[1], alpha=1.0) #     
    area = sum(numpy.diff(bins_data)*a)
    x = np.linspace(bins_data[0], bins_data[-1], 100)
    y=mlab.normpdf(x, mu_R0, sigma_R0)*area
    plt.plot(x,y, linestyle='--',color='k',linewidth=2.0)
    plt.grid(True)
    #plt.xticks([8.1,8.2,8.3,8.4,8.5])
    plt.yticks([0,10,20,30])
    plt.xlim([8.08,8.54])
    
    mean_Th0 = modval['Th0'] # initial values
    (mu_Th0, sigma_Th0) = norm.fit(dic[dic_table_names[zz][2:]]['Th0'])
    normed_Th0 = [x for x in dic[dic_table_names[zz][2:]]['Th0']]    
    # connecting lines between the bbox and the inset axes area
    #ax1.arrow(N_bri_sources[0],normed_R0[0], 55, 0.23, head_length=10, head_width=0.03, fc='g', ec='g')
    #ax1.arrow(107,8.34, 75, 0.23, head_length=10, head_width=0.03, fc='k', ec='k')
    # this is an inset axes over the main axes
    l2 = plt.axes([.6, .7, .17, .17])
    a,bins_data,c=plt.hist(normed_Th0,bins=8,color=colork2[2], alpha=1.0) #     
    area = sum(numpy.diff(bins_data)*a)
    x = np.linspace(bins_data[0], bins_data[-1], 100)
    y=mlab.normpdf(x, mu_Th0, sigma_Th0)*area
    plt.plot(x,y, linestyle='--',color='k',linewidth=2.0)
    plt.grid(True)
    plt.xticks([230,235,240,245,250])
    plt.yticks([0,10,20,30])
    plt.xlim([225,253])
    
    mean_US = modval['US'] # initial values
    (mu_US, sigma_US) = norm.fit(dic[dic_table_names[zz][2:]]['US'])
    normed_US = [x for x in dic[dic_table_names[zz][2:]]['US']]    
    # connecting lines between the bbox and the inset axes area
    #ax1.arrow(N_bri_sources[0],normed_R0[0], 55, 0.23, head_length=10, head_width=0.03, fc='g', ec='g')
    #ax1.arrow(107,8.34, 75, 0.23, head_length=10, head_width=0.03, fc='k', ec='k')
    # this is an inset axes over the main axes
    l3 = plt.axes([.17, .14, .17, .17])
    a,bins_data,c=plt.hist(normed_US,bins=8,color=colork2[3], alpha=1.0) #     
    area = sum(numpy.diff(bins_data)*a)
    x = np.linspace(bins_data[0], bins_data[-1], 100)
    y=mlab.normpdf(x, mu_US, sigma_US)*area
    plt.plot(x,y, linestyle='--',color='k',linewidth=2.0)
    plt.grid(True)
    #plt.xticks([230,235,240,245,250])
    plt.yticks([0,10,20,30])
    plt.xlim([-3,9])
    
    mean_VS = modval['VS'] # initial values
    (mu_VS, sigma_VS) = norm.fit(dic[dic_table_names[zz][2:]]['VS'])
    normed_VS = [x for x in dic[dic_table_names[zz][2:]]['VS']]    
    # connecting lines between the bbox and the inset axes area
    #ax1.arrow(N_bri_sources[0],normed_R0[0], 55, 0.23, head_length=10, head_width=0.03, fc='g', ec='g')
    #ax1.arrow(107,8.34, 75, 0.23, head_length=10, head_width=0.03, fc='k', ec='k')
    # this is an inset axes over the main axes
    l4 = plt.axes([.6, .14, .17, .17])
    a,bins_data,c=plt.hist(normed_VS,bins=8,color=colork2[4], alpha=1.0) #     
    area = sum(numpy.diff(bins_data)*a)
    x = np.linspace(bins_data[0], bins_data[-1], 100)
    y=mlab.normpdf(x, mu_VS, sigma_VS)*area
    plt.plot(x,y, linestyle='--',color='k',linewidth=2.0)
    plt.grid(True)
    #plt.xticks([230,235,240,245,250])
    plt.yticks([0,10,20,30])
    #plt.setp(plt.get_xticklabels(), visible=False) 
    #plt.xlim([-3,9])
    plt.savefig(path_outputs+'/draft_plot.pdf',bbox_inches='tight')      
    

##########################################################################################################################
###################################     Pearson Coefficents Plots   ######################################################
##########################################################################################################################
if (plots_confirmation_pearson=='yes'):                
    b_w_n = raw_input("BeSSeL or whole galaxy? (b or w or none):")
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
    
    Order_Pearson['R_0-T_0'][0]=0.48
    Order_Pearson['T_0-V_0'][0]=-0.73
    Order_Pearson['T_0-V_s'][0]=0.62 
    Order_Pearson['U_0-U_s'][0]=0.51
    Order_Pearson['V_0-V_s'][0]=-0.68
    colork3=['g','c','r','b','m']
                    
    #Plot                                
    matplotlib.rcParams.update({'font.size': 26})        
    fig1 = plt.figure(figsize=(12,12))
    ax1 = fig1.add_subplot(111)
    Big_Pearsons=['R_0-T_0','T_0-V_0','T_0-V_s','U_0-U_s','V_0-V_s'] #Big ones Pearsons
    #Plotting the evolutions
    w=0
    for k in Big_Pearsons:
        ax1.scatter(N_bri_sources,Order_Pearson[k],color=colork3[w],s=60) #        
        label_num=pearso_coeffs.index(k)
        if(b_w_n=='b'):
            ax1.plot(N_bri_sources,Order_Pearson[k],color=colork3[w], linewidth=3,label=pearso_coeffs_labels[label_num])# bessel
        if(b_w_n=='w'):
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
    ax1.set_ylim([-1.05,1.05])#1.49])   
    if(b_w_n=='w'):
        ax1.set_xlim([80,1010]) 
    if(b_w_n=='b'):
        ax1.set_xlim([80,510])
    #Plotting complete sample lines
    w=0    
    for i in Big_Pearsons:    
        if(b_w_n=='w'):
            ax1.axhline(y=Order_Pearson_com[i],linewidth=3,color=colork3[w],ls='dashed') #whole
        w+=1
    #Label for complete samples
    if(b_w_n=='w'):
        ax1.axhline(y=-50,color='k',linewidth=3,ls='dashed',label='Complete Simulation')
        #ax1.plot(103, 5.0,'*',color='k',markersize=13,label='Current BeSSeL')
        ax1.legend(scatterpoints=1,bbox_to_anchor=(0.8,1.17),ncol=1,numpoints=1,prop={'size':27}) # whole
    if(b_w_n=='b'):
        ax1.plot(103, 5.0,'*',color='k',markersize=13,label='Model A5')
        #ax1.axhline(y=-50,color='k',linewidth=3,ls='dashed',label='Complete Simulation')
        ax1.legend(scatterpoints=1,bbox_to_anchor=(0.99,-0.13),ncol=3,numpoints=1,prop={'size':27})
        '''
        handles,labels = ax1.get_legend_handles_labels()
        handles = [handles[0], handles[1], handles[6], handles[3], handles[4], handles[2], handles[5]]
        labels = [labels[0], labels[1], labels[6],labels[3], labels[4], labels[2],labels[5]]
        ax1.legend(handles,labels,scatterpoints=1,bbox_to_anchor=(0.95, 0.13),ncol=3,numpoints=1,prop={'size':22})
        '''
        #ax1.legend() #bessel
    #Saving
    #plt.savefig(path_outputs+'/Pearson_evol.pdf',bbox_inches='tight')  
    
##########################################################################################################################
###################################     PDFs for 1000 sources     ##########################################
##########################################################################################################################    
if (plots_confirmation_PDFs=='yes'):    
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
    a,bins_data,c=ax1.hist(add1,bins=12,color=colork2[1],normed=True) #   
    (mu_Th0V0, sigma_Th0V0) = norm.fit(add1)
    area = sum(numpy.diff(bins_data)*a)
    xmin,xmax=ax1.get_xlim()
    x = np.linspace(xmin, xmax, 100)
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
    ax1.set_xlim([245,268])  

    ax2 = fig1.add_subplot(132)
    normed_R0 = [x for x in dic[dic_table_names[zz][2:]]['R0']]    
    add2=[((normed_Th0[x]+normed_V0[x])/normed_R0[x]) for x in range(len(normed_Th0))]        
    a,bins_data,c=ax2.hist(add2,bins=12,color=colork2[2],normed=True) #       
    (mu_Th0V0R0, sigma_Th0V0R0) = norm.fit(add2)
    area = sum(numpy.diff(bins_data)*a)
    xmin,xmax=ax2.get_xlim()
    x = np.linspace(xmin, xmax, 100)
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
    ax2.set_xlim([29.8,31.4]) 
    ax2.legend(bbox_to_anchor=(1.3, 1.23),numpoints=1,ncol=2,prop={'size':15})
                                                                                                        
    ax3 = fig1.add_subplot(133)
    normed_VS = [x for x in dic[dic_table_names[zz][2:]]['VS']]    
    add3=[normed_V0[x]+normed_VS[x] for x in range(len(normed_Th0))]        
    a,bins_data,c=ax3.hist(add3,bins=12,color=colork2[3],normed=True) #        
    (mu_V0VS, sigma_V0VS) = norm.fit(add3)
    area = sum(numpy.diff(bins_data)*a)
    xmin,xmax=ax3.get_xlim()
    x = np.linspace(xmin, xmax, 100)
    y=mlab.normpdf(x, mu_V0VS, sigma_V0VS)*area        
    y1=mlab.normpdf(x, V0VS_Reid, V0VS_sigma_Reid)*area        
    ax3.plot(x,y1,color='grey',linewidth=2.5)    
    ax3.plot(x,y, linestyle='--',color='k',linewidth=2.5)
    ax3.fill_between(x, 0, y1, facecolor='grey', alpha=0.4)    
    #plt.setp(ax3.get_yticklabels(), visible=False) 
    ax3.grid(True)
    ax3.set_xlabel('$V_{\\odot}+\\bar{V_s}$ (km/s)')   
    y_tickts=np.arange(0,0.5,0.1)
    ax3.set_yticks(y_tickts) 
    x_tickts=np.arange(9,29,4)
    ax3.set_xticks(x_tickts) 
    ax3.set_xlim([8,25.1])
    ax3.set_ylim([0,0.42])         
    #plt.savefig(path_outputs+'/PDFs.pdf',bbox_inches='tight')   
    
    print('Th0+V0: %s +/- %s' %(np.round(mu_Th0V0,1),np.round(sigma_Th0V0,1)))
    print('(Th0+V0)/R0: %s +/- %s' %(np.round(mu_Th0V0R0,2),np.round(sigma_Th0V0R0,2)))
    print('V0+VS: %s +/- %s' %(np.round(mu_V0VS,1),np.round(sigma_V0VS,1)))
    print('--------------------------------------------------')
    print('Reid Results')    
    print('Th0+V0: %s +/- %s' %(Th0V0_Reid,Th0V0_sigma_Reid))
    print('(Th0+V0)/R0: %s +/- %s' %(Th0V0R0_Reid,Th0V0R0_sigma_Reid))
    print('V0+VS: %s +/- %s' %(V0VS_Reid,V0VS_sigma_Reid))    