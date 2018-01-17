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
#Importing auxiliary Libraries
from aux.vectorastrometry31 import phaseSpaceToAstrometry
from aux.coordinates31 import CoordinateTransformation
from aux.gamePlotUtils31 import *
from random import *
from astropy.io import ascii
import os
################################################################################
#Creating Functions
################################################################################
def mygauss( dist, sig ):
    """Evaluates normalized gaussian at dist over sig with mean=0"""
    p = math.exp((-0.5*math.pow(dist,2))/(math.pow(sig,2)))
    return p
################################################################################
def dist_seg( xp, yp, x1, y1, x2, y2):
    """Determines the distance of a point (xp, yp) from a line segment (x1,y1,x2,y2) to the closest point)"""
    xst, yst, xfi, yfi = x1, y1, x2, y2
    #Cannot have start to the right of final
    if (x1>x2):
        xst, xfi, yst, yfi = x2, x1, y2, y1
        if debug: print 'swapped!, change the order of x1 and x2'
    dist = 1e99
    if ( (xfi-xst) != 0 ):
        #not vertical
        a1 = (yfi - yst)/(xfi -xst)
        b1 = yst - a1*xst  
        if (a1 != 0.):
            #not horizontal
            a2 = -1./a1
            b2 = yp - a2*xp
            xs = (b2-b1)/(a1-a2)
            ys = xs*a1+b1
            dist = math.sqrt(math.pow((xs-xp),2) + math.pow((ys-yp),2))
            #but could be outside line segment
            if ( xs > xfi ):
                dist = math.sqrt(math.pow((xfi-xp),2) + math.pow((yfi-yp),2))
            if ( xs < xst ):
                dist = math.sqrt(math.pow((xst-xp),2) + math.pow((yst-yp),2))
        else:
            #horizontal, so that's easy
            dist = abs(yp-yst)
            #outside line segment
            if (xp > xfi):
                dist = math.sqrt(math.pow((xfi-xp),2) + math.pow((yfi-yp),2))
            if (xp < xst):
                dist = math.sqrt(math.pow((xst-xp),2) + math.pow((yst-yp),2))
    else:
        #vertical
        dist = abs(xp-xst)
        #but yfi could be smaller than yst
        if (yfi < yst):
            # inverting yfi and yst, then always:  yfi > yst
            tmp = yfi
            yfi = yst
            yst = tmp
        #outside line segment    
        if (yp > yfi):
            dist = math.sqrt(math.pow((xfi-xp),2) + math.pow((yfi-yp),2))
        if (yp < yst):
            dist = math.sqrt(math.pow((xst-xp),2) + math.pow((yst-yp),2))
    if debug: print 'evaluated: ', xp, yp, xst, yst, xfi, yfi, dist ##????
    return dist
################################################################################
def translbeq(sample,par ):
    """Return sample with Galactic coordinates added"""
    d_sun=par.r0 
    vatr0=par.v0t
    U_sun=par.usun
    V_sun=par.vsun
    W_sun=par.wsun
    out = []
    for i in range( len(sample) ):
        y = sample[i]['y']
        yn = y - d_sun
        x = sample[i]['x']
        z = sample[i]['z']
        l = math.atan2(x,(-1*yn)) 
        d = math.sqrt(math.pow(x,2) + math.pow(yn,2) + math.pow(z,2))
        dg= math.sqrt(math.pow(x,2) + math.pow(y,2) + math.pow(z,2))
        b = math.asin( z / d )
        tmp = ephem.Galactic( l,b, epoch='2000') 
        equ = ephem.Equatorial(tmp)
        ra = float(equ.ra)
        dc = float(equ.dec)
        vr = sample[i]['vr']#U
        vz = sample[i]['vz']#W
        va=sample[i]['va']#V
        th = math.atan2(y,x) 
        vx= va*sin(th) - vr*cos(th)
        vy=-va*cos(th) - vr*sin(th)
        # subtracting the velocity respect to the sun
        vx_sun = vx - vatr0 
        #Velocities vx and vy seen from the LSR, vl (towards the sun) and vt (tangential)
        vl= vx_sun*sin(l)-vy*cos(l)#directly to the LSR
        vt= vx_sun*cos(l)+vy*sin(l)
        f = sample[i]['f']
        f_W=f*3.846*1e26
        d_m=d*1000*3.0857*1e16
        fa =f_W/(4*pi*(d_m**2))
        f_jy_MMB= (fa/1953.125)*1e26#Observed  Flux density in Jy with 0.09 km/s as channel width MMB.
        f_jy_Arecibo= (fa/1525.879)*1e26#Observed  Flux density in Jy with 0.14 km/s as channel width Arecibo.
        f_jy_VLBA= (fa/40000)*1e26#Observed  Flux density in Jy with 1 km/s as channel width VLBA.
        #Dummy Luminosity
        d_dummy=10#kpc
        d_m_dummy=d_dummy*1000*3.0857*1e16
        fa_dummy =f_W/(4*pi*(d_m_dummy**2))
        f_jy_dummy= (fa_dummy/4448.757013)*1e26#Observed  Flux density in Jy with 0.2km/s as channel width
        #Spiral Arm
        ic = sample[i]['ic'] #Identification Color
        #New velocity transfomation. New system sun position to the left
        position_sun=[-d_sun*1000,0,0] #pc
        x_new=-y*1000 #pc
        y_new=x*1000
        z_new=z*1000
        x_barycentric= x_new - position_sun[0]
        y_barycentric= y_new - position_sun[1]
        z_barycentric= z_new - position_sun[2]
        vel_sun=[U_sun,V_sun+vatr0,W_sun]
        vx_new=-vy
        vy_new=vx
        vz_new=vz
        #Heliocentric                
        vx_barycentric = vx_new - vel_sun[0] # km/s
        vy_barycentric = vy_new - vel_sun[1] 
        vz_barycentric = vz_new - vel_sun[2] 
        #Important= vl = vel_{LSR} and vrad_sun = vel_{Helio}
        #velocity transfomation from barycentric to l and b        
        l_b, b_b, parallax, mul_cosb, mub, vrad_sun= phaseSpaceToAstrometry(x_barycentric, y_barycentric, z_barycentric, \
            vx_barycentric,vy_barycentric,vz_barycentric)
        #velocity transfomation from l and b to alpha and delta
        trans= CoordinateTransformation(0) #gal to eq.
        mua, mud= trans.transformProperMotions(l_b, b_b, mul_cosb, mub)
        mu_x=mua
        #Uncertainties
        if(par.Unc_confirmation==1.0):
            if(par.Unc_confirmation_spec==0.0):
                '''
                Calculation:
                EVN calculator, set U band, data rate=32 Mbits/s---> 2 polarization, 1 suband, 4 MHz Total Bandwidth, int_time= 120 minutes
                EVN Calculator= 0.341 mJy/beam noise in the map!
                But we are interested in spectral line of 1 km/s = 50 kHz--> noise per channel map = 0.341*sqrt(80)
                ---> 0.003 JY
                snr=(f_jy_VLBA/0.003)*math.sqrt(2))#Sp [Jy]/sigma Noise Image =  Sp [Jy]/Jy/hr * sqrt(2hr) 
                err_parallax=MATH.SQRT(math.sqrt(4/pi))*(1.013/(math.sqrt(8*(numpy.log(2)))))*(1/snr)/math.sqrt(4) #1.013=res. vlba for 12 GHz masers, eq 1 Reid et al 1988, 330:809ApJ
                '''
                snr=(f_jy_VLBA)/(0.003*math.sqrt(2))#Sp [Jy]/sigma Noise Image =  Sp [Jy]/Jy/hr * sqrt(2hr)     
                err_parallax=math.sqrt(4/pi)*(1.013/(math.sqrt(8*(numpy.log(2)))))*(1/snr)/math.sqrt(4) #1.013=res. vlba for 12 GHz masers, eq 1 Reid et al 1988, 330:809ApJ
                err_mu_x=err_parallax 
                err_mud=err_parallax 
                err_vlsr=par.sigvr
                if(par.Push_betterunc==1):
                    perce_para=err_parallax*100/parallax
                    if(perce_para>20):
                        err_parallax=parallax*0.2
                        err_mu_x=err_parallax 
                        err_mud=err_parallax                                     
            elif(par.Unc_confirmation_spec==1.0):                
                err_parallax=parallax*par.Unc_confirmation_percent/100
                err_mu_x=err_parallax # mas/1yr
                err_mud=err_parallax # mas/1yr
                err_vlsr=par.sigvr# 5km/s                                               
        elif(par.Unc_confirmation==0.0):
            err_parallax=err_mu_x=err_mud=err_vlsr=0.0                    
        #Adding errors and calling them: _observables
        #NOT INCLUDING FUDGE FACTOR check export files
        if(par.Perfect_data==1.0):
            parallax_obs=gauss(parallax,0.0)            
            mux_obs=gauss(mu_x,0.0)
            mud_obs=gauss(mud,0.0) 
            vl_obs=gauss(vl,0.0)
        elif(par.Perfect_data==0.0):           
            if(par.Perfect_parallax==1.0):
                parallax_obs=gauss(parallax,0.0)
            elif(par.Perfect_parallax==0.0):
                '''
                if(par.Mark_trick==1.0):
                    parallax_obs1=gauss(parallax,err_parallax)            
                    d_news=bias_factor(parallax_obs1,err_parallax)
                    #d_news=(1/parallax_obs1)/abs(bias_factors)
                    parallax_obs=1/d_news
                elif(par.Mark_trick==0.0):
                '''
            parallax_obs=gauss(parallax,err_parallax)            
            mux_obs=gauss(mu_x,err_mu_x)
            mud_obs=gauss(mud,err_mud) 
            vl_obs=gauss(vl,err_vlsr)                                                                                
        if(par.Distance_spread==1.0):            
            new_distance=gauss(d,err_parallax/(parallax**2))                                    
            parallax_obs=1/abs(new_distance)
        if(par.Only_nice_sources==1.0):   
            perce_para_nices=err_parallax*100/parallax_obs
            if(perce_para_nices<=20):                                                                
                if(par.Noerrorbutobs_para==1.0):                          
                    #err_parallax=err_mu_x=err_mud=err_vlsr=0.0                                                         
                    err_parallax=0.0   
                point = { 'x(kpc)':x, 'y(kpc)':y, 'z(kpc)':z, 'va(km/s)':va, 'vr(km/s)':vr,\
                        'vz(km/s)':vz, 'f(Lsun)':f, 'ic(Arm)':ic,'l(rad)':l, 'b(rad)':b,\
                        'dsun(kpc)':d,'dg(kpc)':dg, 'vl(km/s)':vl, 'vt(km/s)':vt, 'fa(W/m2)':fa,
                        'ra(rad)':ra,'dc(rad)':dc,'f_MMB(Jy)':f_jy_MMB,'f_Arecibo(Jy)':f_jy_Arecibo,\
                        'f_dummy(Jy)':f_jy_dummy,'parallax(mas)':parallax,'vrad_sun(km/s)':vrad_sun,\
                        'mua(mas/yr)':mua,'mud(mas/yr)': mud,'err_parallax(mas)':err_parallax,\
                        'mux(mas/yr)':mu_x,'err_mu_x(mas/yr)':err_mu_x, 'err_mud(mas/yr)':err_mud,\
                        'err_vls(km/s)':err_vlsr,'f_jy_VLBA(Jy)':f_jy_VLBA,'parallax_obs(mas)':parallax_obs,\
                        'mux_obs(mas/yr)':mux_obs,'mud_obs(mas/yr)':mud_obs,'vl_obs(km/s)':vl_obs}
                    #'parallax_obs1':parallax_obs1,'d_news':d_news}
                out.append(point)
        elif(par.Only_nice_sources==0.0):                                                                                  
                if(par.Noerrorbutobs_para==1.0):                          
                    err_parallax=err_mu_x=err_mud=err_vlsr=0.0                                                         
                    #err_parallax=0.0   
                point = { 'x(kpc)':x, 'y(kpc)':y, 'z(kpc)':z, 'va(km/s)':va, 'vr(km/s)':vr,\
                        'vz(km/s)':vz, 'f(Lsun)':f, 'ic(Arm)':ic,'l(rad)':l, 'b(rad)':b,\
                        'dsun(kpc)':d,'dg(kpc)':dg, 'vl(km/s)':vl, 'vt(km/s)':vt, 'fa(W/m2)':fa,
                        'ra(rad)':ra,'dc(rad)':dc,'f_MMB(Jy)':f_jy_MMB,'f_Arecibo(Jy)':f_jy_Arecibo,\
                        'f_dummy(Jy)':f_jy_dummy,'parallax(mas)':parallax,'vrad_sun(km/s)':vrad_sun,\
                        'mua(mas/yr)':mua,'mud(mas/yr)': mud,'err_parallax(mas)':err_parallax,\
                        'mux(mas/yr)':mu_x,'err_mu_x(mas/yr)':err_mu_x, 'err_mud(mas/yr)':err_mud,\
                        'err_vls(km/s)':err_vlsr,'f_jy_VLBA(Jy)':f_jy_VLBA,'parallax_obs(mas)':parallax_obs,\
                        'mux_obs(mas/yr)':mux_obs,'mud_obs(mas/yr)':mud_obs,'vl_obs(km/s)':vl_obs}
                #'parallax_obs1':parallax_obs1,'d_news':d_news}            
                out.append(point)
    return out
################################################################################
def Comparison_MMB(sample,par):
    '''Returns the number of sources and the slope of
        the observed flux density function using the limitations
        of the MMB survey (Green et al 2008, MNRAS 392, 783).
    '''
    nbar = 8
    MMB_lim_flux=0.17
    Total_min_flux_MMB=MMB_lim_flux*par.MMB_sigma
    t_Jy = ascii.read('../surveys/MMB_Dist_Flux_edited.csv')    
    t_Jy.remove_row(0)
    kk=[]
    for i in range(len(t_Jy)):
        if(t_Jy[i][5]=='-'):
            kk.append(i)
    t_Jy.remove_rows(kk)
    z=numpy.zeros(len(t_Jy))
    for i in range(len(t_Jy)):            
         z[i]=t_Jy[i][5]
    t=where(z>Total_min_flux_MMB)
    Jy_MMB=z[t[0]]
    min_data_MMB_obs=numpy.log10(min(Jy_MMB))
    max_data_MMB_obs=numpy.log10(max(Jy_MMB))
    max_long=60
    min_long=-174.0
    max_lati=2.0
    min_lati=-2.0
    lp = [ (sample[i]['l(rad)']*180/math.pi) for i in range(len(sample)) ]
    bp = [ (sample[i]['b(rad)']*180/math.pi) for i in range(len(sample)) ] 
    f_Jyp = [ sample[i]['f_MMB(Jy)'] for i in range(len(sample))]
    k=0
    for i in range(len(f_Jyp)):
        if(lp[i]>min_long and lp[i]<max_long and bp[i]>=min_lati and bp[i]<=max_lati and f_Jyp[i]>(Total_min_flux_MMB)):
            #limits of the Arecibo Survey see: Green et al 2008, MNRAS 392, 783.
            k+=1
    f_Jy_MMB=np.zeros(k)
    w=0
    for i in range( len(f_Jyp) ):
        if(lp[i]>min_long and lp[i]<max_long and bp[i]>=min_lati and bp[i]<=max_lati and f_Jyp[i]>(Total_min_flux_MMB)):
            f_Jy_MMB[w]=f_Jyp[i]
            w+=1
    print 'Sources that satisfy MMB limits= %s (Simulation) vs %s (Observed)' % (len(f_Jy_MMB),len(Jy_MMB))
    #Histogram Plot
    bins_data=10**(numpy.linspace(min_data_MMB_obs,max_data_MMB_obs,nbar))
    y_simu,x_simu_aa=numpy.histogram(f_Jy_MMB,bins=bins_data)
    x_simu=numpy.zeros(len(x_simu_aa)-1)
    for i in range(len(x_simu_aa)-1):
        x_simu[i]=(x_simu_aa[i+1]+x_simu_aa[i])/2.0
    remove_zeros_simu=[]
    for i in range(len(y_simu)):
        if (y_simu[i] == 0):
            remove_zeros_simu.append(i)     
    y_simu_new=numpy.delete(y_simu,remove_zeros_simu)
    x_simu_new=numpy.delete(x_simu,remove_zeros_simu)
    log_x_simu=numpy.log10(x_simu_new)
    log_y_simu=numpy.log10(y_simu_new)
    coeffs_simu = numpy.polyfit(log_x_simu,log_y_simu,deg=1)
    #Observations
    y_obs,x_obs_aa=numpy.histogram(Jy_MMB,bins=bins_data)
    remove_zeros_data=[]
    for i in range(len(y_obs)):
        if (y_obs[i] == 0):
            remove_zeros_data.append(i) 
    y_obs_new=numpy.delete(y_obs,remove_zeros_data)
    x_obs_new=numpy.delete(x_obs_aa,remove_zeros_data)
    ysimu_new=numpy.delete(y_simu,remove_zeros_data)        
    #Chi square
    chi_2=[]
    for i in range(len(y_obs_new)):
        chi_2.append(((float(ysimu_new[i]-y_obs_new[i]))**2)/float(y_obs_new[i]))    
    #Kullback–Leibler divergence
    Total_simu=sum(ysimu_new)
    prob_ysimu=ysimu_new/(Total_simu*1.0)
    Total_obs=sum(y_obs_new)
    prob_yobs=y_obs_new/(Total_obs*1.0)        
    KLD=numpy.dot(prob_ysimu,numpy.log10(prob_yobs/prob_ysimu))
    result=[len(f_Jy_MMB),coeffs_simu[0],sum(chi_2),abs(KLD)]
    return(result)
################################################################################
def Comparison_Arecibo(sample,par):
    '''Returns the number of sources and the slope of
        the observed flux density function using the limitations
        of the Arecibo survey (Panadian et al, 2007, ApJ 656, 255).
    '''
    nbar = 8
    Arecibo_lim_flux=0.27
    Total_min_flux_Arecibo=Arecibo_lim_flux*par.Arecibo_sigma
    t_Jy = ascii.read('../surveys/Arecibo_2.txt',format='basic')
    z=numpy.zeros(len(t_Jy))
    for i in range(len(t_Jy)):
         z[i]=t_Jy[i][0]
    t=where(z>Total_min_flux_Arecibo)
    Jy_pan=z[t[0]]
    min_data=numpy.log10(min(Jy_pan))
    max_data=numpy.log10(max(Jy_pan))   
    lp = [ (sample[i]['l(rad)']*180/math.pi) for i in range(len(sample)) ]
    bp = [ (sample[i]['b(rad)']*180/math.pi) for i in range(len(sample)) ] 
    f_Jyp = [ sample[i]['f_Arecibo(Jy)'] for i in range(len(sample)) ] # in Jy
    k=0
    for i in range( len(f_Jyp) ):
        if(lp[i]>=35.2 and lp[i]<=53.7 and bp[i]>=-0.41 and bp[i]<=0.41 and f_Jyp[i]>Total_min_flux_Arecibo):
                #limits of the Arecibo Survey see: Panadian et al, 2007, ApJ 656, 255.
            k+=1 
    f_Jy_Arecibo=np.zeros(k)
    w=0
    for i in range( len(f_Jyp) ):
        if(lp[i]>=35.2 and lp[i]<=53.7 and bp[i]>=-0.41 and bp[i]<=0.41 and f_Jyp[i]>Total_min_flux_Arecibo):
            f_Jy_Arecibo[w]=f_Jyp[i]
            w+=1
    print 'Sources that satisfy Arecibo limits= %s (Simulation) vs %s (Observed)' % (len(f_Jy_Arecibo),len(Jy_pan))
    #Histogram
    bins_data=10**(numpy.linspace(min_data,max_data,nbar))
    y_simu,x_simu_aa=numpy.histogram(f_Jy_Arecibo,bins=bins_data) 
    x_simu=numpy.zeros(len(x_simu_aa)-1)
    for i in range(len(x_simu_aa)-1):
        x_simu[i]=(x_simu_aa[i+1]+x_simu_aa[i])/2.0
    remove_zeros_simu=[]
    for i in range(len(y_simu)):
        if (y_simu[i] == 0):
            remove_zeros_simu.append(i)     
    y_simu_new=numpy.delete(y_simu,remove_zeros_simu)
    x_simu_new=numpy.delete(x_simu,remove_zeros_simu)
    log_x_simu=numpy.log10(x_simu_new)
    log_y_simu=numpy.log10(y_simu_new)
    coeffs_simu = numpy.polyfit(log_x_simu,log_y_simu,deg=1)    
    #Observations
    y_obs,x_obs_aa=numpy.histogram(Jy_pan,bins=bins_data)
    remove_zeros_data=[]
    for i in range(len(y_obs)):
        if (y_obs[i] == 0):
            remove_zeros_data.append(i) 
    y_obs_new=numpy.delete(y_obs,remove_zeros_data)
    x_obs_new=numpy.delete(x_obs_aa,remove_zeros_data)
    ysimu_new=numpy.delete(y_simu,remove_zeros_data)        
    #Chi square
    chi_2=[]
    for i in range(len(y_obs_new)):
        chi_2.append(((float(ysimu_new[i]-y_obs_new[i]))**2)/float(y_obs_new[i]))
    #Kullback–Leibler divergence
    Total_simu=sum(ysimu_new)
    prob_ysimu=ysimu_new/(Total_simu*1.0)
    Total_obs=sum(y_obs_new)
    prob_yobs=y_obs_new/(Total_obs*1.0)        
    KLD=numpy.dot(prob_ysimu,numpy.log10(prob_yobs/prob_ysimu))
    result=[len(f_Jy_Arecibo),coeffs_simu[0],sum(chi_2),abs(KLD)]
    return(result)  
################################################################################
def ellipse_extra(ra,rb,ang,x0,y0,Nb=50):
    '''
    Returns x and y of an ellipse
    ra - major axis length
    rb - minor axis length
    ang - angle
    x0,y0 - position of centre of ellipse
    Nb - No. of points that make an ellipse
    '''
    xpos,ypos=x0,y0
    radm,radn=ra,rb
    an=ang*pi/180.0    
    co,si=numpy.cos(an),numpy.sin(an)
    the=numpy.linspace(0,2*pi,Nb)
    X=radm*numpy.cos(the)*co-si*radn*numpy.sin(the)+xpos
    Y=radm*numpy.cos(the)*si+co*radn*numpy.sin(the)+ypos
    return X,Y    
################################################################################    
def create_array_control(control,parameter):
    lists = control[parameter]
    listsr= lists.replace (",", " ")
    array = [int(x) for x in listsr.split()]               
    return array
################################################################################
def do_plots(sample,par):
    #only_arms_plot(par,save='yes')
    spatial_distribution_plot(sample,par,save='yes')
    velocity_distribution_plot(sample,par,save='yes')
    Luminosity_distribution_plot(sample,save='yes')
    Combined_plots(sample,par,save='yes')
    Survey_Comp(sample,par,save='yes')
    Parrallax_Plots(sample,par,save='yes')
    return     
################################################################################
def create_dir_plots(path,fname_directory,key,i):
    os.chdir(path+'/output/'+fname_directory+'/output_plots')        
    os.makedirs(key+str(i))
    os.chdir(key+str(i))
    return
################################################################################
def gauss_prob(c_d, c_s, d):
    twopi=2*np.pi
    prefactor = np.sqrt(1.0/twopi) / abs(c_s)
    arg = (d - c_d)**2 / (2.0*c_s**2)
    prob = 0.#d0
    if ( arg < 9.99 ):
        prob = prefactor * math.exp(-arg)
    return prob
################################################################################
def bias_factor(parallax,sigma_par):
# Generate a Gaussian parallax pdf
    parallax = 1./15.
    sigma_par= parallax*0.2
    #print('PDFs for parallax = %s +/- %s' %(parallax,sigma_par))
    npts = 200   
    par_min=parallax*0.2
    par_max=parallax*2.5
    #par_min = parallax-(4*sigma_par)
    #par_max = parallax+(4*sigma_par)
    par=np.zeros(npts+1)
    dist=np.zeros(npts+1)
################################################################################
    p_bin_size = (par_max - par_min) / (npts-1)
    p_sum = 0.
    p_int = 0.
    p0    = 0.
    p1    = 0.
#    print ('  Parallax   ProbDen   ProbInt')
    
    for n in range(1,npts+1):
        p = par_min + (n-1)*p_bin_size
        prob=gauss_prob(parallax,sigma_par,p)
        par[n] = prob*p_bin_size
        p_int  = p_int + par[n]
#        print(round(p,5),round(par[n],5),round(p_int,5))
        p_sum = p_sum + p*par[n]
    #print(' Mean parallax = %s' % p_sum)
################################################################################    #     Convert to a distance pdf
    d_sum = 0.
    d_int = 0.
    d_min = 1.0/par_max
    #d_min = (1/parallax)-((1/sigma_par)*4)
    d_max = 1.0/par_min
    #d_max = (1/parallax)+((1/sigma_par)*4)
    d_bin_size = (d_max - d_min) / (npts-1)
#    print ('  Distance   ProbDen   ProbInt')
    for n in range(1,npts+1):   
        par_d = d_min + (n-1)*d_bin_size 
        p = 1/par_d                      
        prob=gauss_prob( parallax, sigma_par, p)
        dist[n] = prob * p**2 * d_bin_size 
        d_int   = d_int + dist[n]
#        print(round(par_d,5),round(dist[n],5),round(d_int,5))
        d_sum = d_sum + par_d*dist[n]
    #print ('Mean distance = %s' % d_sum)
    return(d_sum)    
    ################################################################################