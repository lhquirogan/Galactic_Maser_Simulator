# -*- coding: utf-8 -*-
################################################################################
#Importing Libraries
################################################################################
global debug, interactive
debug = 0
interactive = 1
from pylab import *
from math import *
import matplotlib
import sys
from random import *
import csv
import ephem
import numpy
from aux.vectorastrometry31 import astrometryToPhaseSpace
from aux.coordinates31 import CoordinateTransformation
'''
RA_= 4.10260447470627
DEC_=-0.9717806719519921
p_=0.11372493164013868
mua_cosd_=-5.9094577681165354
mud_=-4.59847062429938
v_lsr_=-82.8531650512696
'''

def obtain_vr_va_vz(RA_,DEC_,p_,mua_cosd_,mud_,v_lsr_):
    r0=8.34
    th0=240
    u_sun_=10.7
    v_sun_=15.6
    w_sun_=8.9
    trans= CoordinateTransformation(1) # equa to galac
    l_,b_=trans.transformSkyCoordinates(RA_,DEC_)
    mu_l_, mu_b_= trans.transformProperMotions(RA_, DEC_, mua_cosd_, mud_)
    #x,y,z,vx,vy,vz=astrometryToPhaseSpace(l, b, p, mu_l, mu_b, v_lsr)
    #x,y,z,vx,vy,vz=astrometryToPhaseSpace(RA, DEC, p, mua_cosd, mud, v_lsr)
    x_baric,y_baric,z_baric,vx_baric,vy_baric,vz_baric=astrometryToPhaseSpace(l_, b_, p_, mu_l_, mu_b_, v_lsr_)
    
    
    x_=y_baric/1000
    y_=-1*((x_baric/1000)-r0)
    z_=z_baric/1000
    
    vx_new_=vx_baric+u_sun_
    vy_new_=vy_baric+v_sun_+th0
    vz_new_=vz_baric+w_sun_
    
    vx_=vy_new_
    vy_=-vx_new_
    vz_=vz_new_
    
    th_ = math.atan2(y_,x_) 
    va_= vx_*sin(th_) - vy_*cos(th_)
    vr_=-vx_*cos(th_) - vy_*sin(th_)
    return(vr_)
    


import numpy as np
import astropy
from astropy.io import ascii
from astropy.table import Table, Column
import scipy as scipy
from scipy.stats import kde
from scipy.stats.stats import pearsonr
%cd "/Users/luishenryquiroganunez/Documents/Leiden_University/3_Semester/MSc_thesis/Codes/Reid_code/fit_galaxy_bayesian/"
#Reading all the files names
t = ascii.read('/Users/luishenryquiroganunez/Documents/Leiden_University/3_Semester/MSc_thesis/Codes/Reid_code/fit_galaxy_bayesian/source_files.txt',format='basic')

ra_=['' for x in range(len(t))]
dc_=['' for x in range(len(t))]
strs=['' for x in range(len(t))]
data=['' for x in range(len(t))]
#Physical Paramaters
para=['' for x in range(len(t))]
para_unc=['' for x in range(len(t))]
mux=['' for x in range(len(t))]
mux_unc=['' for x in range(len(t))]
muy=['' for x in range(len(t))]
muy_unc=['' for x in range(len(t))]
vlsr=['' for x in range(len(t))]
vlsr_unc=['' for x in range(len(t))]

for i in range(len(t)):
        strs[i]=str(t[i][0])
        #strs[i]=str(t[i])
        with open(strs[i]) as f:
            content = f.readlines()
        matching = [s for s in content if "RA  (hhmmss.s)" in s]
        puto=matching[-1].index('.')
        ra_h=float(matching[-1][0:puto-4])
        ra_m=float(matching[-1][puto-4:puto-2])/60
        ra_s=float(matching[-1][puto-2:puto+6])/3600
        ra_[i]=(ra_h+ra_m+ra_s)*15
        kk=content.index(matching[-1])
        #Error in parallax
        pui=content[kk+1]
        puto2=pui.index('.')        
        dc_h=float(pui[0:puto2-4])
        dc_m=float(pui[puto2-4:puto2-2])/60
        dc_s=float(pui[puto2-2:puto2+6])/3600
        dc_[i]=(dc_h+dc_m+dc_s)       
        data=content[kk+2]
        s=list(data)
        a=s[8:16]
        k=''.join(a)
        para_unc[i]=float(k)
        b=s[0:7]
        l=''.join(b)
        para[i]=float(l)
        #Error in prop x
        data=content[kk+3]
        s=list(data)
        a=s[8:16]
        k=''.join(a)
        mux_unc[i]=float(k)
        b=s[0:7]
        l=''.join(b)
        mux[i]=float(l)
        #Error in prop y        
        data=content[kk+4]
        s=list(data)
        a=s[8:16]
        k=''.join(a)
        muy_unc[i]=float(k)
        b=s[0:7]
        l=''.join(b)
        muy[i]=float(l)
        #Error in v_lsr        
        data=content[kk+5]
        s=list(data)
        a=s[8:16]
        k=''.join(a)
        vlsr_unc[i]=float(k)
        b=s[0:7]
        l=''.join(b)
        vlsr[i]=float(l)
        

vr_=['' for x in range(len(t))]
for i in range(len(t)):
    RA_=ra_[i]*pi/180
    DEC_=dc_[i]*pi/180
    vr_[i]=obtain_vr_va_vz(RA_,DEC_,para[i],mux[i],muy[i],vlsr[i])        

plt.hist(vr_)
plt.savefig('../../Code_v31/plots/radial_dist.png')