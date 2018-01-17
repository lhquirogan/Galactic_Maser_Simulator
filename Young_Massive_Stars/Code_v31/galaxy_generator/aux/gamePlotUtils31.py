import csv
import numpy
from pylab import *
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import astropy
from astropy.io import ascii
import matplotlib.lines as mlines
import heapq
import numpy as np
from matplotlib.ticker import NullFormatter
import math
import os
from aux.gameModelPars31 import *
import scipy.stats as stats
#---------------------------------------------------------------------------------------------------------
def roundup(x):
    return int(math.ceil(x / 100.0)) * 100
#---------------------------------------------------------------------------------------------------------
def readcsv(fname):
    sample = []

    fcsv = open(fname)
    for row in csv.DictReader( fcsv):
        for keys in row:
            row[keys] = float(row[keys])
        sample.append(row)
    return sample
#---------------------------------------------------------------------------------------------------------
def only_arms_plot(par,save='no'):
    plt.clf()
    if(save=='yes'):
        #plt.ioff()
        pass
    else:
        print('The plots were NOT saved!')
        
    figure(figsize=(10,10))
    #matplotlib.rcParams.update({'font.size': 27}) 
    segs = par.gs.segs_spiral(par.rtrunc, par.nspirplot) 
    for j in range( len(segs) ):
        name_arm=['Norma Arm','Carina-Sagitarius Arm','Perseus Arm','Crux-Scutum Arm','Local Arm']
        line_arm_color=['r-','c-','g-','m-','b-']
        #for k in range(len(segs[0]['x'])):
        #    segs[j]['y'][k]=segs[j]['y'][k]*-1
        plot( segs[j]['x'], segs[j]['y'], line_arm_color[j],label=name_arm[j], linewidth=3)
        #plot(segs[j]['y'], segs[j]['x'],  line_arm_color[j],label=name_arm[j])
    norma = mlines.Line2D([], [], color='r')
    carina = mlines.Line2D([], [], color='c')
    perseus= mlines.Line2D([], [], color='g')
    crux = mlines.Line2D([], [], color='m')
    local = mlines.Line2D([], [], color='b')    
    #legend1=plt.legend([perseus],['Perseus'],prop={'size':22},loc=4,ncol=1)           
    #legend1=plt.legend([norma,carina,perseus,crux,local],['Norma','Carina-Sagitarius',\
     #           'Perseus','Crux-Scutum','Local'],prop={'size':22},loc=4,ncol=1)
    #plt.gca().add_artist(legend1)
    text( 12.0,  12.0, 'Norma',     verticalalignment='bottom', horizontalalignment='right',     color='red', fontsize=25)
    text( -8.0, -14.0, 'Carina \n Sagitarius',     verticalalignment='bottom', horizontalalignment='right',     color='c', fontsize=25)
    text( -2.0,  11.0, 'Perseus',     verticalalignment='bottom', horizontalalignment='right',     color='green', fontsize=25)
    text( 13.0, -13.0, 'Crux-Scutum',     verticalalignment='bottom', horizontalalignment='right',     color='m', fontsize=25)
    text( -8.0, 5.5, 'Local',     verticalalignment='bottom', horizontalalignment='right',     color='b', fontsize=25)
    xlabel('X $[kpc]$')
    ylabel('Y $[kpc]$')
    ring = par.cr.central_ring(par.nspirplot) 
    plot(ring[0]['x'], ring[0]['y'], 'k--', linewidth = 3)
    text(0.0,0.0, 'GC', fontsize=15)
    plot(0.0,0.0,'o',color='g',markersize=10)
    axis([-par.rtrunc-5.0,par.rtrunc+5.0,-par.rtrunc-5.0,par.rtrunc+5.0])
    axis('equal')
    plot( 0, par.r0, '*', color ='y',markersize=13)
    #plot( -par.r0, 0, '*', color ='y',markersize=10)
    plt.grid(True)
    name_plot_1="arms.pdf"
    if(save=='yes'):
        plt.savefig(name_plot_1,bbox_inches='tight')
        #print(name_plot_1)
    plt.close()
    plt.ion()
    return
    
#---------------------------------------------------------------------------------------------------------    
def only_arms_plot3d(par,save='no'):
    plt.clf()
    if(save=='yes'):
        #plt.ioff()
        pass
    else:
        print('The plots were NOT saved!')
        

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111,projection='3d')        
    #matplotlib.rcParams.update({'font.size': 27}) 
    segs = par.gs.segs_spiral(par.rtrunc, par.nspirplot) 
    for j in range( len(segs) ):
        name_arm=['Norma Arm','Carina-Sagitarius Arm','Perseus Arm','Crux-Scutum Arm','Local Arm']
        line_arm_color=['r-','c-','g-','m-','b-']
        #for k in range(len(segs[0]['x'])):
        #    segs[j]['y'][k]=segs[j]['y'][k]*-1
        ax.plot( segs[j]['x'], segs[j]['y'], 0,line_arm_color[j],label=name_arm[j], linewidth=3)
        #plot(segs[j]['y'], segs[j]['x'],  line_arm_color[j],label=name_arm[j])
    norma = mlines.Line2D([], [], color='r')
    carina = mlines.Line2D([], [], color='c')
    perseus= mlines.Line2D([], [], color='g')
    crux = mlines.Line2D([], [], color='m')
    local = mlines.Line2D([], [], color='b')    
    #legend1=plt.legend([perseus],['Perseus'],prop={'size':22},loc=4,ncol=1)           
    #legend1=plt.legend([norma,carina,perseus,crux,local],['Norma','Carina-Sagitarius',\
     #           'Perseus','Crux-Scutum','Local'],prop={'size':22},loc=4,ncol=1)
    #plt.gca().add_artist(legend1)
    #ax.text( 12.0,  12.0, 0,'Norma',     verticalalignment='bottom', horizontalalignment='right',     color='red', fontsize=25)
    #ax.text( -8.0, -14.0,0, 'Carina \n Sagitarius',     verticalalignment='bottom', horizontalalignment='right',     color='c', fontsize=25)
    #ax.text( -2.0,  11.0,0, 'Perseus',     verticalalignment='bottom', horizontalalignment='right',     color='green', fontsize=25)
    #ax.text( 13.0, -13.0,0, 'Crux-Scutum',     verticalalignment='bottom', horizontalalignment='right',     color='m', fontsize=25)
    #ax.text( -8.0, 5.5,0, 'Local',     verticalalignment='bottom', horizontalalignment='right',     color='b', fontsize=25)
    #ax.set_xlabel('X $[kpc]$')
    #ax.set_ylabel('Y $[kpc]$')
    ring = par.cr.central_ring(par.nspirplot) 
    ax.plot(ring[0]['x'], ring[0]['y'], 'w--', linewidth = 3)
    #ax.text(0.0,0.0,0, 'GC', fontsize=15)
    ax.scatter(0.0,0.0,0,marker='+',color='k',s=10)
    axis([-par.rtrunc-5.0,par.rtrunc+5.0,-par.rtrunc-5.0,par.rtrunc+5.0])
    ax.axis('equal')
    ax.axis('off')
    
    x, y = ogrid[0:img.shape[0], 0:img.shape[1]]
    ax = gca(projection='3d')
    ax.plot_surface(x, y, 10, rstride=5, cstride=5, facecolors=img)    
    #ax.plot( 0, par.r0,0, '*', color ='y',markersize=13)
    #plot( -par.r0, 0, '*', color ='y',markersize=10)
    #ax.grid(True)
    name_plot_1="arms.pdf"
    if(save=='yes'):
        plt.savefig(name_plot_1,bbox_inches='tight')
        #print(name_plot_1)
    plt.close()
    plt.ion()
    return    
#---------------------------------------------------------------------------------------------------------
def spatial_distribution_plot(sample,par,save='no'):
    if(save=='yes'):
        #plt.ioff()
        pass
    else:
        print('The plots were NOT saved!')
    
    xp = [ sample[i]['x(kpc)'] for i in range(len(sample)) ]
    yp = [ sample[i]['y(kpc)'] for i in range(len(sample)) ]
    zp = [ sample[i]['z(kpc)'] for i in range(len(sample)) ]
    ic = [ int(sample[i]['ic(Arm)']) for i in range(len(sample)) ]
    rp = []
    for i in range(len(xp)):
        rp.append( sqrt(xp[i]**2+yp[i]**2) )
    #Only sources (exponential decay)
    plt.clf()
    matplotlib.rcParams.update({'font.size': 17})
    figure()
    aa=1-(1/(exp(par.rmintrunc/par.halfmr)))
    bb=1-(1/(exp(par.rtrunc/par.halfmr)))
    x=numpy.zeros(len(sample))
    y=numpy.zeros(len(sample))
    for i in range(len(sample)):
        r_o=-par.halfmr*numpy.log(1-numpy.random.uniform(aa,bb))
        tethas=numpy.random.uniform(0,2*pi)
        x[i]=r_o*cos(tethas) 
        y[i]=r_o*sin(tethas) 
        plot(x[i], y[i], '.', color='b', alpha = 0.4, markersize = 10)
    xlabel('X $[kpc]$')
    ylabel('Y $[kpc]$')
    axis([-par.rtrunc-5,par.rtrunc+5,-par.rtrunc-5,par.rtrunc+5])
    axis('equal')
    grid(True)
    name_plot_1="spatial_distr_decayxy.pdf"
    if(save=='yes'):
        plt.savefig(name_plot_1,bbox_inches='tight')
        #print(name_plot_1)
    plt.close()

    #Histograms N vs radial (exp decay)
    matplotlib.rcParams.update({'font.size': 27})        
    fig3 = plt.figure(figsize=(12,12))  
    ax1 = fig3.add_subplot(2,1,1)
    #title("Radial maser distribution")
    nbar_d=np.arange(0.0,16.0,1.0)
    number,bins,rectangles=ax1.hist(rp,nbar_d,color ='r',width=1.0,alpha=0.5)
    #ax1.set_xlabel('r $[kpc]$')
    ax1.set_ylabel('Number \n $[counts]$')
    ax1.set_yticks([0,100,200,300])
    plt.setp( ax1.get_xticklabels(), visible=False)
    #ax1.axes.get_xaxis().set_ticks([]) 
    ax1.grid(True)
    ax1.set_xlim([2,16])

    ax2 = fig3.add_subplot(2,1,2,sharex=ax1)
    #title("Radial density distribution of masers")
    surf_dens=numpy.zeros(len(number))
    for i in range(len(number)):
        surf_dens[i]=number[i]/(pi*(numpy.power(bins[i+1],2)-numpy.power(bins[i],2)))
    aa=np.delete(bins,-1) # erase the last value of bins, to match the size of surf_dens
    ax2.bar(aa,surf_dens,color='b',width=1.0,alpha=0.5)
    ax2.set_yticks([0,2,4,6,8,10])
    ax2.set_xlabel('R $[kpc]$')
    ax2.set_ylabel('Surface density \n $[counts \ kpc^{-2}]$') 
    ax2.grid(True)
    ax2.set_xlim([2,16])
    ax2.set_ylim([0,12.5])
    fig3.subplots_adjust(hspace=0.02)
    name_plot_2="spatial_distr_decay.pdf"
    if(save=='yes'):
        plt.savefig(name_plot_2,bbox_inches='tight')
        #print(name_plot_2)
    plt.close()
    
    #Only arms + sources
    matplotlib.rcParams.update({'font.size': 27}) 
    figure(figsize=(12,12))
    segs = par.gs.segs_spiral(par.rtrunc, par.nspirplot) # why different than par.nspirseg? OJO
    for j in range( len(segs) ):
        name_arm=['Norma Arm','Carina-Sagitarius Arm','Perseus Arm','Crux-Scutum Arm','Local Arm']
        line_arm_color=['r-','c-','g-','m-','b-']
        #for k in range(len(segs[0]['x'])):
        #    segs[j]['y'][k]=segs[j]['y'][k]*-1
        plot( segs[j]['x'], segs[j]['y'], line_arm_color[j],label=name_arm[j], linewidth=3)
        #plot(segs[j]['y'], segs[j]['x'],  line_arm_color[j],label=name_arm[j])
    norma = mlines.Line2D([], [], color='r')
    carina = mlines.Line2D([], [], color='c')
    perseus= mlines.Line2D([], [], color='g')
    crux = mlines.Line2D([], [], color='m')
    local = mlines.Line2D([], [], color='b')    
    text( 12.0,  12.0, 'Norma',     verticalalignment='bottom', horizontalalignment='right',     color='red', fontsize=25)
    text( -8.0, -14.0, 'Carina \n Sagitarius',     verticalalignment='bottom', horizontalalignment='right',     color='c', fontsize=25)
    text( -2.0,  11.0, 'Perseus',     verticalalignment='bottom', horizontalalignment='right',     color='green', fontsize=25)
    text( 13.0, -13.0, 'Crux-Scutum',     verticalalignment='bottom', horizontalalignment='right',     color='m', fontsize=25)
    text( -8.0, 5.5, 'Local',     verticalalignment='bottom', horizontalalignment='right',     color='b', fontsize=25)
    for i in range(len(xp) ):
        plot( xp[i], yp[i], '.', alpha = 0.4 , color = par.colors[ic[i]], markersize =8.0)
    xlabel('X $[kpc]$')
    ylabel('Y $[kpc]$')
    axis([-par.rtrunc-5.0,par.rtrunc+5.0,-par.rtrunc-5.0,par.rtrunc+5.0])
    axis('equal')
    plot( 0, par.r0, '*', color ='y',markersize=13)
    ring = par.cr.central_ring(par.nspirplot) 
    plot(ring[0]['x'], ring[0]['y'], 'k--', linewidth = 3)
    text(0.0,0.0, 'GC', fontsize=15)
    plot(0.0,0.0,'o',color='g',markersize=10)
    #plot( -par.r0, 0, '*', color ='y',markersize=10)
    plt.grid(True)
    name_plot_3="spatial_distr_xyarms.pdf"
    if(save=='yes'):
        plt.savefig(name_plot_3,bbox_inches='tight')
        #print(name_plot_3)
    plt.close()
    #z distribution vs r 
    x=rp
    y=zp
    nullfmt = NullFormatter()         # no labels
    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left+width+0.02
    rect_scatter = [left, bottom, width, height]
    #rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]
    # start with a rectangular Figure
    matplotlib.rcParams.update({'font.size': 40}) 
    figure(figsize=(20,8))
    axScatter = plt.axes(rect_scatter)
    #axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)
    # no labels
    #axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)
    # the scatter plot:
    axScatter.scatter(x, y)
    # now determine nice limits by hand:
    nbins=20
    xmax=np.max(np.fabs(x))+5.0
    xmin=0.0
    ymax=np.max(np.fabs(y))*2.0
    ymin=ymax*-1
    axScatter.set_xlim( (xmin, xmax) )
    axScatter.set_ylim( (ymin, ymax) )
    axScatter.set_xlabel('R $[kpc]$')
    axScatter.set_ylabel('Z $[kpc]$')
    axScatter.grid(True)
    a2,b2,c2=axHisty.hist(y, nbins,color='g', orientation='horizontal')
    #axHistx.set_xlim( axScatter.get_xlim() )
    #axHistx.set_yticks(np.arange(min(a1), max(a1), 200))
    #axHistx.set_ylabel('N $[counts]$')
    axHisty.set_ylim( axScatter.get_ylim() )
    axHisty.set_xticks(np.linspace(0, roundup(max(a2)), 3))
    axHisty.set_xlabel('N $[counts]$')
    name_plot_4="spatial_distr_z.pdf"
    if(save=='yes'):
        plt.savefig(name_plot_4,bbox_inches='tight')
        #print(name_plot_4)
    plt.close()
    plt.ion()
    return
    
def velocity_distribution_plot(sample,par,save='no'):
    
    if(save=='yes'):
        #plt.ioff()
        pass
    else:
        print('The plots were NOT saved!')
        
    vpr = [ sample[i]['vr(km/s)'] for i in range(len(sample)) ]
    vpz = [ sample[i]['vz(km/s)'] for i in range(len(sample)) ]
    
    #Plot vz vs vr with histograms
    plt.clf()
    matplotlib.rcParams.update({'font.size': 23}) 
    x=vpr
    y=vpz
    nullfmt = NullFormatter()         # no labels
    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left+width+0.02
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]
    # start with a rectangular Figure
    plt.figure(figsize=(12,12))
    axScatter = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)
    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)
    # the scatter plot:
    axScatter.scatter(x, y)
    # now determine nice limits by hand:
    binwidth = 1
    xymax = np.max( [np.max(np.fabs(x)), np.max(np.fabs(y))] ) + 1
    xymin = xymax*-1
    axScatter.set_xlim( (xymin-2.0, xymax+2.0) )
    axScatter.set_ylim( (xymin-2.0, xymax+2.0) )
    axScatter.set_xlabel('U $[km \ s^{-1}]$')
    axScatter.set_ylabel('W $[km \ s^{-1}]$')
    bins = np.arange(xymin, xymax, binwidth)
    a1,b1,c1=axHistx.hist(x, bins=bins,color='r')
    a2,b2,c2=axHisty.hist(y, bins=bins,color='g', orientation='horizontal')
    axHistx.set_xlim( axScatter.get_xlim() )
    axHistx.set_yticks(np.arange(min(a2), max(a2), 200))
    axHistx.set_ylabel('N $[counts]$')
    axHisty.set_ylim( axScatter.get_ylim() )
    axHisty.set_xticks(np.arange(min(a2), max(a2), 200))
    axHisty.set_xlabel('N $[counts]$')
    name_plot_1="velocity_vr-vz-histo.pdf"
    if(save=='yes'):
        plt.savefig(name_plot_1,bbox_inches='tight')
        #print(name_plot_1)
    plt.close()
    '''    
    #Plot gaussian 
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    matplotlib.rcParams.update({'font.size': 23}) 
    nbins=20
    mu_v0r, sigma_v0r = par.v0r, par.sigvr # mean and standard deviation
    x = numpy.linspace(-20.0,20.0,1000)
    y=mlab.normpdf(x,mu_v0r,sigma_v0r)*9200
    a,bins_data,c=ax1.hist(vpr,nbins,color='g', alpha=0.7) 
    ax1.plot(x,y,linewidth=3)       
    #ax1.fill_between(x, y.min(), y, facecolor='blue', alpha=0.2)
    ax1.grid(True)
    #ax1.legend(loc=0)
    #title("Radial Velocity Distribution")
    xlabel('Radial Velocity $U$ $[km\ s^{-1}]$')
    ylabel('N $[Counts]$')
    name_plot_2="velocity_distr_vr-gaussian.pdf"
    if(save=='yes'):
        plt.savefig(name_plot_2,bbox_inches='tight')
        #print(name_plot_2)
    plt.close()
    '''
    plt.ion()
    return
    
def Luminosity_distribution_plot(sample,save='no'):
    
    if(save=='yes'):
        #plt.ioff()
        pass
    else:
        print('The plots were NOT saved!')
    plt.clf()   
    #Takin data
    fp = [ sample[i]['f(Lsun)'] for i in range(len(sample)) ]
    #f_Jy_dummyp = [ sample[i]['f_dummy(Jy)'] for i in range(len(sample)) ]
    f_Jyp = [ sample[i]['f_MMB(Jy)'] for i in range(len(sample)) ]   

    #Parameters to plot
    nbar = 20
    min_dat=numpy.log10(min(fp))
    max_dat=numpy.log10(max(fp))    
    bins_dat=10**(numpy.linspace(min_dat,max_dat,nbar))    
    min_data=numpy.log10(min(f_Jyp))
    max_data=numpy.log10(max(f_Jyp))    
    bins_data=10**(numpy.linspace(min_data,max_data,nbar))
    
    #Plot
    matplotlib.rcParams.update({'font.size': 25})
    matplotlib.rcParams['xtick.major.pad']='8'     
    fig3 = plt.figure(figsize=(20,8))
    #Subplot(lineas,columnas,numero)
    ax1 = fig3.add_subplot(121)
    ax1.hist(fp,bins=bins_dat,log=True,color='cadetblue',alpha=0.6) 
    #ax1.set_aspect('equal', adjustable='box')
    ax1.set_xscale("log")
    ax1.set_xlim([8E-9,2E-3])    
    ax1.set_ylim([1,5E2])    
    ax1.set_yticks([1E0,1E1,1E2])    
    ax1.set_xlabel('Luminosity Peak $\\rm{L_p} \, \\rm{(L_{\\odot})}$')
    ax1.set_ylabel('Source Counts')
    ax1.grid(True)
    z=numpy.array([1E-8,1E-7,1E-6,1E-5,1E-4,1E-3])
    ax1.set_xticks(z)    
    ax2 = fig3.add_subplot(122,sharey=ax1)
    plt.setp( ax2.get_yticklabels(), visible=False)
    ax2.hist(f_Jyp,bins=bins_data,color ='firebrick',log=True,alpha=0.6)#,weights=weights)     
    #ax2.hist(f_Jy_dummyp,bins=bins_datas,log=True,color='g') 
    #ax2.set_aspect('equal', adjustable='box')
    ax2.set_xscale("log")
    #ax2.axes.get_yaxis().set_ticks([]) 
    ax2.set_xlim([3E-2,8E4])    
    ax2.set_ylim([1,5E2])    
    ax2.set_yticks([1E0,1E1,1E2])    
    ax2.set_xlabel('Flux density Peak $\\rm{S_p \, (Jy)}$')
    zz=numpy.array([1E-1,1E0,1E1,1E2,1E3,1E4])
    ax2.set_xticks(zz)    
    ax2.grid(True)
    fig3.subplots_adjust(wspace=0.01)
    name_plot_1="luminosity_distr.pdf"
    if(save=='yes'):
        plt.savefig(name_plot_1,bbox_inches='tight')
        #print(name_plot_1)
    plt.close()

    matplotlib.rcParams.update({'font.size': 25})
    matplotlib.rcParams['xtick.major.pad']='8'     
    fig4 = plt.figure(figsize=(12,8))
    #Subplot(lineas,columnas,numero)
    ax3 = fig4.add_subplot(111)
    ax3.hist(fp,bins=bins_dat,log=True,color='cadetblue',alpha=0.6) 
    #ax1.set_aspect('equal', adjustable='box')
    ax3.set_xscale("log")
    ax3.set_xlim([8E-9,2E-3])    
    ax3.set_ylim([1,5E2])    
    ax3.set_yticks([1E0,1E1,1E2])    
    ax3.set_xlabel('Peak luminosity $\\rm{L_p}$ $\\rm{(L_{\\odot})}$')
    ax3.set_ylabel('Source Counts')
    ax3.grid(True)
    z=numpy.array([1E-8,1E-7,1E-6,1E-5,1E-4,1E-3])
    ax3.set_xticks(z)        
    name_plot_2="luminosity_distribution_only.pdf"
    if(save=='yes'):
        plt.savefig(name_plot_2,bbox_inches='tight')
    plt.close()    

    matplotlib.rcParams.update({'font.size': 25})
    matplotlib.rcParams['xtick.major.pad']='8'     
    fig5 = plt.figure(figsize=(12,8))
    #Subplot(lineas,columnas,numero)
    ax4 = fig5.add_subplot(111)
    #plt.setp( ax2.get_yticklabels(), visible=False)
    ax4.hist(f_Jyp,bins=bins_data,color ='firebrick',log=True,alpha=0.6)#,weights=weights)     
    #ax2.hist(f_Jy_dummyp,bins=bins_datas,log=True,color='g') 
    #ax2.set_aspect('equal', adjustable='box')
    ax4.set_xscale("log")
    #ax2.axes.get_yaxis().set_ticks([]) 
    ax4.set_xlim([3E-2,8E4])    
    ax4.set_ylim([1,5E2])    
    ax4.set_yticks([1E0,1E1,1E2])  
    ax4.set_ylabel('Source Counts')      
    ax4.set_xlabel('Peak flux density $\\rm{S_p}$ (Jy)')
    zz=numpy.array([1E-1,1E0,1E1,1E2,1E3,1E4])
    ax4.set_xticks(zz)    
    ax4.grid(True)
    name_plot_3="flux_density_distribution_only.pdf"
    if(save=='yes'):
        plt.savefig(name_plot_3,bbox_inches='tight')
    plt.close()    
    plt.ion()
    return
    
def Combined_plots(sample,par,save='no'):
    if(save=='yes'):
        #plt.ioff()
        pass
    else:
        print('The plots were NOT saved!')
    #plt.clf()  
    xp = [ sample[i]['x(kpc)'] for i in range(len(sample)) ]
    yp = [ sample[i]['y(kpc)'] for i in range(len(sample)) ]
    zp = [ sample[i]['z(kpc)'] for i in range(len(sample)) ]
    vpr = [ sample[i]['vr(km/s)'] for i in range(len(sample)) ]
    vpa = [ sample[i]['va(km/s)'] for i in range(len(sample)) ]
    vpz = [ sample[i]['vz(km/s)'] for i in range(len(sample)) ]
    ic = [ int(sample[i]['ic(Arm)']) for i in range(len(sample)) ]
    dgp = [ sample[i]['dg(kpc)'] for i in range(len(sample)) ]   
    lp = [ (sample[i]['l(rad)']*180/math.pi) for i in range(len(sample)) ]
    bp = [ (sample[i]['b(rad)']*180/math.pi) for i in range(len(sample)) ]
    vl = [ sample[i]['vl(km/s)'] for i in range(len(sample)) ]
    vt = [ sample[i]['vt(km/s)'] for i in range(len(sample)) ]
    fp = [ sample[i]['f(Lsun)'] for i in range(len(sample)) ]
    f_Jyp = [ sample[i]['f_MMB(Jy)'] for i in range(len(sample)) ]   
    ra = [ sample[i]['ra(rad)']*180/math.pi for i in range(len(sample)) ]
    dc = [ sample[i]['dc(rad)']*180/math.pi for i in range(len(sample)) ]
    rp = []
    muap=[ sample[i]['mux_obs(mas/yr)'] for i in range(len(sample)) ]#OJO these are observables not exact
    mudp=[ sample[i]['mud_obs(mas/yr)'] for i in range(len(sample)) ]#OJO these are observables not exact
    for i in range(len(xp)):
        rp.append( sqrt(xp[i]**2+yp[i]**2) ) 
    
    #Plot xy with lumosities 
    matplotlib.rcParams.update({'font.size': 27}) 
    figure(figsize=(12,12))
    segs = par.gs.segs_spiral(par.rtrunc, par.nspirplot) # why different than par.nspirseg? OJO
    for j in range( len(segs) ):
        name_arm=['Norma Arm','Carina-Sagitarius Arm','Perseus Arm','Crux-Scutum Arm','Local Arm']
        line_arm_color=['r-','c-','g-','m-','b-']
        #for k in range(len(segs[0]['x'])):
        #    segs[j]['y'][k]=segs[j]['y'][k]*-1
        plot( segs[j]['x'], segs[j]['y'], line_arm_color[j],label=name_arm[j], linewidth=3)
        #plot(segs[j]['y'], segs[j]['x'],  line_arm_color[j],label=name_arm[j])
    norma = mlines.Line2D([], [], color='r')
    carina = mlines.Line2D([], [], color='c')
    perseus= mlines.Line2D([], [], color='g')
    crux = mlines.Line2D([], [], color='m')
    local = mlines.Line2D([], [], color='b')    
    text( 12.0,  12.0, 'Norma',     verticalalignment='bottom', horizontalalignment='right',     color='red', fontsize=25)
    text( -8.0, -14.0, 'Carina \n Sagitarius',     verticalalignment='bottom', horizontalalignment='right',     color='c', fontsize=25)
    text( -2.0,  11.0, 'Perseus',     verticalalignment='bottom', horizontalalignment='right',     color='green', fontsize=25)
    text( 14.0, -12.0, 'Crux-Scutum',     verticalalignment='bottom', horizontalalignment='right',     color='m', fontsize=25)
    text( -8.0, 5.5, 'Local',     verticalalignment='bottom', horizontalalignment='right',     color='b', fontsize=25)
    text( -0.8, -2.1, 'Central \n Ring', verticalalignment='bottom', horizontalalignment='center', color='k', fontsize=16)    
    for i in range(len(xp) ):
        plot( xp[i], yp[i], '.', alpha = 0.4 , color = par.colors[ic[i]], markersize = (numpy.log10(fp[i])*4)+37)
        #plot( -yp[i], xp[i], '.', alpha = 0.4 , color = par.colors[ic[i]], markersize = (numpy.log10(fp[i])*4)+37)
        z=i
    #small_dot, =plot( -yp[z], xp[z],  '.', alpha = 0.4 , color = par.colors[ic[z]], markersize = 5)
    #big_dot,   =plot( -yp[z], xp[z], '.', alpha = 0.4 , color = par.colors[ic[z]], markersize = 25)
    small_dot, =plot( xp[z], yp[z], '.', alpha = 0.4 , color = par.colors[ic[z]], markersize = 5)
    big_dot,   =plot( xp[z], yp[z], '.', alpha = 0.4 , color = par.colors[ic[z]], markersize = 25)
    legend2=legend([small_dot, big_dot],['$10^{-8} L_{\odot}$','$10^{-3} L_{\odot}$'],loc=4,ncol=2,numpoints=1,prop={'size':23},borderpad=0.2)
    plt.gca().add_artist(legend2) 
    ring = par.cr.central_ring(par.nspirplot) 
    plot(ring[0]['x'], ring[0]['y'], 'k--', linewidth = 3)
    text(0.0,0.0, 'GC', fontsize=15)  
    plot(0.0,0.0,'o',color='g',markersize=10)
    xlabel('X $[kpc]$')
    ylabel('Y $[kpc]$')
    axis([-par.rtrunc-10.0,par.rtrunc+10.0,-par.rtrunc-10.0,par.rtrunc+10.0])
    axis('equal')
    plot( 0, par.r0, '*', color ='y',markersize=13)
    #plot( -par.r0, 0, '*', color ='y',markersize=10)
    plt.grid(True)
    name_plot_0="samplot.pdf"
    if(save=='yes'):
        plt.savefig(name_plot_0,bbox_inches='tight')
        #print(name_plot_0)
    plt.close()
    
    #Rotation Curve Plot
    #Keplerian radius
    aaa=numpy.arange(0.1,20.0,0.1)
    M_galaxy=8e10 #solar masses
    M_ratio=M_galaxy/(8e10)
    r_ratio=(aaa/8.34)
    v_keplerian=150*np.sqrt(M_ratio/r_ratio)
    #solid body
    solid_v=100*aaa
    #Halo
    halo_v=240.0*np.sqrt((1-((1.0/aaa)*arctan(aaa/1.0))))
    aa=numpy.arange(0.0,25.0,0.1)
    cte=(par.v0t)-(par.rot_curve*par.r0)
    bb=(par.rot_curve*aa)+cte
    #bb=numpy.ones(len(aa))*par.v0t
    #Plot   
    x=rp
    y=vpa    
    nullfmt = NullFormatter()         # no labels
    
    matplotlib.rcParams.update({'font.size': 17}) 
    figure(figsize=(7,10))
    # definitions for the axes
    rect_scatter = [0.1, 0.1, 0.89, 0.6]
    rect_histx = [0.1, 0.72, 0.89, 0.19]
    #rect_histy = [left_h, bottom, width, height]
    # start with a rectangular Figure
    axScatter = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    #axHisty = plt.axes(rect_histy)
    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    #axHisty.yaxis.set_major_formatter(nullfmt)
    # the scatter plot:
    #axScatter.scatter(x, y)
    axScatter.plot( rp, vpa, '.', alpha =0.5, color = 'orange', label='Masers', markersize=15)
    axScatter.plot(aa , bb , '-' , color = 'k',label='$d\Theta/dR$',linewidth=3)
    axScatter.plot(par.r0,par.v0t, '*', color ='y',markersize=13,label='Sun')
    axScatter.plot(aaa,v_keplerian,linestyle='dashed',color='k', label='Keplerian velocity',linewidth=3)
    axScatter.plot(aaa,solid_v,linestyle='dotted',color='k', label='Solid-Body rotation',linewidth=3)
    axScatter.plot(aaa,halo_v,linestyle='dashdot',color='k',label='Halo rotation',linewidth=3)
    # now determine nice limits by hand:
    nbins=12
    axScatter.set_xlabel('Galactocentric Distance $R$ (kpc)',fontsize=18)
    axScatter.set_ylabel('Rotation Speed $\\Theta(R)$ (km/s)',fontsize=18)
    axScatter.grid(True)
    axScatter.set_yticks([0,50,100,150,200,250,300])
    axScatter.set_ylim([0,350])
    axScatter.set_xlim([0,17])
    axScatter.set_xticks([0,3,6,9,12,15])
    axScatter.legend(loc=8,prop={'size':17},numpoints=1,ncol=2,bbox_to_anchor = (0.5, 0.0))
    #add_artist(legends)
    a1,b1,c1=axHistx.hist(x, nbins,color='orange',label='Radial Distribution')
    #a2,b2,c2=axHisty.hist(y, nbins,color='g', orientation='horizontal')
    axHistx.set_xlim( axScatter.get_xlim() )
    axHistx.set_xticks(axScatter.get_xticks())
    axHistx.set_yticks([0,100,200,300])
    axHistx.set_ylabel('Source Counts',fontsize=17)
    axHistx.grid(True)
    axHistx.legend(loc=0,prop={'size':18},numpoints=1)#bbox_to_anchor = (0.25, 0.69)    
    #axHisty.set_ylim( axScatter.get_ylim() )
    #axHisty.set_xticks([0,400,800])
    #axHisty.set_xlabel('N $[counts]$')
    name_plot_1="rot_curve.pdf"
    if(save=='yes'):
        plt.savefig(name_plot_1,bbox_inches='tight')
        #print(name_plot_1)
    plt.close()

    from matplotlib import gridspec
    #Rotation Curve Plot PAPER    
    matplotlib.rcParams.update({'font.size': 17}) 
    fig3=figure(figsize=(7,9))
    fig3.subplots_adjust(hspace=0.25)
    gs=gridspec.GridSpec(2,1,height_ratios=[1,2])
    axHistx = fig3.add_subplot(gs[0])    
    axScatter = fig3.add_subplot(gs[1])        
    # definitions for the axes
    #rect_scatter = [0.1, 0.1, 0.89, 0.6]
    #rect_histx = [0.1, 0.72, 0.89, 0.19]
    #rect_histy = [left_h, bottom, width, height]
    # start with a rectangular Figure
    #axScatter.axes(rect_scatter)
    #axHistx.axes(rect_histx)
    #axHisty = plt.axes(rect_histy)
    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    #axHisty.yaxis.set_major_formatter(nullfmt)
    # the scatter plot:
    #axScatter.scatter(x, y)
    axScatter.plot( rp, vpa, '.', alpha =0.5, color = 'orange', label='Masers', markersize=15)
    axScatter.plot(aa , bb , '-' , color = 'k',label='$d\Theta/dR$',linewidth=3)
    axScatter.plot(par.r0,par.v0t, '*', color ='y',markersize=13,label='Sun')
    # now determine nice limits by hand:
    nbins=12
    axScatter.set_xlabel('Galactocentric Distance $R$ (kpc)',fontsize=18)
    axScatter.set_ylabel('Rotation Speed $\\Theta(R)$ (km/s)',fontsize=18)
    axScatter.grid(True)
    axScatter.set_yticks([150,200,250,300])
    axScatter.set_ylim([150,310])
    axScatter.set_xlim([0,17])
    axScatter.set_xticks([0,3,6,9,12,15])
    axScatter.legend(loc=9,prop={'size':17},numpoints=1,ncol=3)
    #add_artist(legends)
    a1,b1,c1=axHistx.hist(x, nbins,color='orange',label='Radial Distribution')
    #a2,b2,c2=axHisty.hist(y, nbins,color='g', orientation='horizontal')
    axHistx.set_xlim( axScatter.get_xlim() )
    axHistx.set_xticks(axScatter.get_xticks())
    axHistx.set_yticks([0,100,200,300])
    axHistx.set_ylim([0,350])    
    axHistx.set_ylabel('Source Counts',fontsize=17)
    axHistx.grid(True)
    axHistx.legend(loc=0,prop={'size':18},numpoints=1)#bbox_to_anchor = (0.25, 0.69)    
    #axHisty.set_ylim( axScatter.get_ylim() )
    #axHisty.set_xticks([0,400,800])
    #axHisty.set_xlabel('N $[counts]$')
    name_plot_14="rot_curve_paper.pdf"
    if(save=='yes'):
        plt.savefig(name_plot_14,bbox_inches='tight')
        #print(name_plot_1)
    plt.close()
    
    
    #Observed Flux Density Function
    matplotlib.rcParams.update({'font.size': 23}) 
    hola=figure()
    ofdf=hola.add_subplot(111)
    nbar = 20
    min_data=numpy.log10(min(f_Jyp))
    max_data=numpy.log10(max(f_Jyp))    
    bins_data=10**(numpy.linspace(min_data,max_data,nbar))
    plt.gca().set_xscale("log")
    #weights=np.ones_like(f_Jyp)/len(f_Jyp)
    ofdf.hist(f_Jyp,bins=bins_data,color ='b',log=True,alpha=0.5)#,weights=weights) 
    ofdf.set_xlabel('$S_p$ (Jy)')
    ofdf.set_ylabel('Source Counts')
    ofdf.set_ylim([0.9,1500])
    ofdf.set_xlim([1E-2,1E6])
    ofdf.set_xticks(numpy.logspace(-1,5,7))
    name_plot_2="observed_flux.pdf"
    if(save=='yes'):
        plt.savefig(name_plot_2,bbox_inches='tight')
        #print(name_plot_2)
    plt.close()
    
    #lat-long- view
    figure()
    matplotlib.rcParams.update({'font.size': 18})
    map = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')
    #title("Profile View")
    tl=numpy.zeros(len(lp))
    tb=numpy.zeros(len(lp))
    for i in range( len(lp) ):
        tl[i], tb[i] = map(lp[i],bp[i])
        map.plot( tl[i], tb[i], 'o', alpha = 0.5 , color = par.colors[ic[i]], markersize = (numpy.log10(fp[i])*2)+21)
        z=i
    map.drawparallels(numpy.arange(-90.,120.,30.),labels=[True])
    map.drawmeridians(np.arange(0.,420.,60.))
    map.drawmapboundary(fill_color='GhostWhite')
    #sun_1, sun_2 = map(360,0.)
    #map.plot( sun_1, sun_2, '*', color ='y',markersize=10)
    small_dot, =map.plot( tl[z], tb[z], 'o', alpha = 0.5 , color = par.colors[ic[z]], markersize = 5)
    big_dot,   =map.plot( tl[z], tb[z], 'o', alpha = 0.5 , color = par.colors[ic[z]], markersize = 15)
    lg=legend([small_dot, big_dot],['$10^{-8} L_{\odot}$','$10^{-3} L_{\odot}$'],loc=8,numpoints=1,prop={'size':16})
    legend1=legend([norma,carina,perseus,crux,local],['Norma Arm','Carina-Sagitarius Arm',\
                'Perseus Arm','Crux-Scutum Arm','Local Arm'],prop={'size':11},loc='upper center',bbox_to_anchor = (0.5, 0.87),ncol=2)
    plt.gca().add_artist(lg)
    plt.gca().add_artist(legend1)
    name_plot_3="l-b.pdf"
    if(save=='yes'):
        plt.savefig(name_plot_3,bbox_inches='tight')
        #print(name_plot_3)
    plt.close()
    
    #V_lsr vs long
    matplotlib.rcParams.update({'font.size': 27})
    vlsr_plot=figure(figsize=(12,8))
    vlsr_ax=vlsr_plot.add_subplot(111)
    min_dat=numpy.log10(min(f_Jyp))
    max_dat=numpy.log10(max(f_Jyp))  
    #title("Galactic Methanol Masers Distribution")
    for i in range( len(lp) ):
        vlsr_ax.plot( lp[i], vl[i], 'o', alpha = 0.5 , color = par.colors[ic[i]], markersize = ((numpy.log10(fp[i])*3.0)+27.0))
        z=i
    #plt.plot( lon, vls, '+', color = 'y',label='Pestalozzi')
    small_dot, =vlsr_ax.plot( lp[z], vl[z], 'o', alpha = 0.4 , color = 'k', markersize = 3)
    big_dot,   =vlsr_ax.plot( lp[z], vl[z], 'o', alpha = 0.4 , color = 'k', markersize = 17.5)
    vlsr_ax.legend([small_dot, big_dot],['$10^{-8} L_{\odot}$','$10^{-3} L_{\odot}$'],loc=4,numpoints=1,prop={'size':32})
    vlsr_ax.axis([180,-180,-200,200])
    vlsr_ax.set_xlabel('Galactic Longitude ')
    vlsr_ax.set_ylabel('$V_{LSR}$ (km/s)')
    vlsr_ax.grid(True)
    name_plot_4="l-v.pdf"
    if(save=='yes'):
        plt.savefig(name_plot_4,bbox_inches='tight')
        #print(name_plot_4)
    plt.close()        

    #Geocentric view in rac and dec
    matplotlib.rcParams.update({'font.size': 18})
    map = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')   
    figure()
    tx=numpy.zeros(len(ra))
    ty=numpy.zeros(len(ra))
    for i in range(len(ra)):
        tx[i], ty[i] = map(ra[i],dc[i])
        map.plot( tx[i], ty[i], 'o', alpha = 0.5 , color = par.colors[ic[i]], markersize = (numpy.log10(f_Jyp[i])*2.0)+5.0)
        z=i
    map.drawparallels(numpy.arange(-90.,120.,30.),labels=[True])
    map.drawmeridians(np.arange(0.,420.,60.))
    map.drawmapboundary(fill_color='GhostWhite')
    small_dot,  =map.plot( tx[z], ty[z], 'o', alpha = 0.5 , color = 'k', markersize = 2.0)
    big_dot,    =map.plot( tx[z], ty[z], 'o', alpha = 0.5 , color = 'k', markersize = 15)
    lg=legend([small_dot, big_dot],['$0.02Jy$','$10^{3} Jy$'],loc=8,numpoints=1,prop={'size':18})
    name_plot_5="rac_dec.pdf"
    if(save=='yes'):
        plt.savefig(name_plot_5,bbox_inches='tight')
        #print(name_plot_5)  
    plt.close()
                        
    #Vt vs long
    matplotlib.rcParams.update({'font.size': 27})
    vlsr_plot2=figure(figsize=(20,8))
    vlsr_ax2=vlsr_plot2.add_subplot(111)
    min_dat2=numpy.log10(min(f_Jyp))
    max_dat2=numpy.log10(max(f_Jyp))  
    #title("Galactic Methanol Masers Distribution")
    for i in range( len(lp) ):
        vlsr_ax2.plot( lp[i], vt[i], 'o', alpha = 0.5 , color = par.colors[ic[i]], markersize = ((numpy.log10(fp[i])*3.0)+27.0))
        z=i
    #plt.plot( lon, vls, '+', color = 'y',label='Pestalozzi')
    small_dot2, =vlsr_ax2.plot( lp[z], vt[z], 'o', alpha = 0.4 , color = 'k', markersize = 3)
    big_dot2,   =vlsr_ax2.plot( lp[z], vt[z], 'o', alpha = 0.4 , color = 'k', markersize = 17.5)
    vlsr_ax2.legend([small_dot2, big_dot2],['$10^{-8} L_{\odot}$','$10^{-3} L_{\odot}$'],loc=4,numpoints=1,prop={'size':32})
    #vlsr_ax2.axis([180,-180,-200,200])
    vlsr_ax2.set_xlabel('Galactic Longitude ($^o$)')
    vlsr_ax2.set_ylabel('$V_{t}$ $(km/s)$')
    vlsr_ax2.grid(True)
    name_plot_6="l-vt.pdf"
    if(save=='yes'):
        plt.savefig(name_plot_6,bbox_inches='tight')
        #print(name_plot_6) 
    plt.close() 
        
    #sqrt(mua^2 + mud^2) vs long                
    matplotlib.rcParams.update({'font.size': 27})
    mu_plot=figure(figsize=(20,10))
    mu_ax=mu_plot.add_subplot(111)  
    #title("Galactic Methanol Masers Distribution")
    for i in range( len(lp) ):
        mu_ax.plot( lp[i],numpy.sqrt((muap[i]**2)+(mudp[i]**2)), 'o', alpha = 0.5 , color = par.colors[ic[i]], markersize = ((numpy.log10(fp[i])*3.0)+27.0))
        z=i
    #plt.plot( lon, vls, '+', color = 'y',label='Pestalozzi')
    small_dot3, =mu_ax.plot( lp[z], numpy.sqrt((muap[i]**2)+(mudp[i]**2)), 'o', alpha = 0.4 , color = 'k', markersize = 3)
    big_dot3,   =mu_ax.plot( lp[z], numpy.sqrt((muap[i]**2)+(mudp[i]**2)), 'o', alpha = 0.4 , color = 'k', markersize = 17.5)
    mu_ax.legend([small_dot3, big_dot3],['$10^{-8} L_{\odot}$','$10^{-3} L_{\odot}$'],loc=1,numpoints=1,prop={'size':32})
    #vlsr_ax2.axis([180,-180,-200,200])
    mu_ax.set_xlabel('Galactic Longitude ($^o$)')
    mu_ax.set_ylabel('$\\mu_{mas/yr}$ $(km/s)$')
    mu_ax.grid(True)
    name_plot_7="l-mu.pdf"
    if(save=='yes'):
        plt.savefig(name_plot_7,bbox_inches='tight')
        #print(name_plot_7)  
    plt.close()
    plt.ion()
    return
    
def files_for_sterre(sample,par,save='no',sources=200):
    import os
    if(save=='yes'):
        os.chdir("/Users/luishenryquiroganunez/Documents/Leiden_University/3_Semester/MSc_thesis/Codes/Code_v16/file_sterre/")
        print('The following file was saved in file_sterre folder:')
    else:
        print('The file was NOT saved!')
    #------------------------------------------------------------------------------------------------------------        
    segs = par.gs.segs_spiral(par.rtrunc, par.nspirplot) # why different than par.nspirseg? OJO
    #------------------------------------------------------------------------------------------------------------
    Norma_x=segs[0]['x']
    Norma_y=segs[0]['y']
    Carina_Sagitarius_x=segs[1]['x']
    Carina_Sagitarius_y=segs[1]['y']
    Perseus_x=segs[2]['x']
    Perseus_y=segs[2]['y']
    Crux_Scutum_x=segs[3]['x']
    Crux_Scutum_y=segs[3]['y']
    Local_x=segs[4]['x']
    Local_y=segs[4]['y']   
    xp = [ sample[i]['x(kpc)'] for i in range(sources) ]
    yp = [ sample[i]['y(kpc)'] for i in range(sources) ]
    zp = [ sample[i]['z(kpc)'] for i in range(sources) ]
    #------------------------------------------------------------------------------------------------------------
    arms = [Norma_x, Norma_y, Carina_Sagitarius_x, Carina_Sagitarius_y, Perseus_x,Perseus_y,Crux_Scutum_x,Crux_Scutum_y,Local_x,Local_y]#Final Table
    masers= [xp,yp,zp]
    if(save=='yes'):
        ascii.write(arms,'arms_sterre.txt',names=['Norma_x', 'Norma_y','Carina-Sagitarius_x','Carina-Sagitarius_y','Perseus_x','Perseus_y','Crux-Scutum_x','Crux-Scutum_y','Local_x','Local_y'])
        ascii.write(masers,'masers_sterre.txt',names=['X[kpc]','Y[kpc]','Z[kpc]'])
        print('File=arms_sterre.txt')  
        print('File=masers_sterre.txt') 
    plt.close() 
    return
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------        
def Survey_Comp(sample,par,save='no'):
    if(save=='yes'):
        #plt.ioff()
        pass
    else:
        print('The plots were NOT saved!')
    #limits of the Arecibo Survey see: Panadian et al, 2007, ApJ 656, 255:
    nbar = 8
    Arecibo_min_long=35.2
    Arecibo_max_long=53.7
    Arecibo_min_lat=-0.41
    Arecibo_max_lat=0.41  
    Arecibo_lim_flux=0.27
    Total_min_flux_Arecibo=Arecibo_lim_flux*par.Arecibo_sigma
    #MMB limits
    MMB_min_long=-180.0
    MMB_max_long=60.0
    MMB_min_lat=-2.0
    MMB_max_lat=2.0  
    MMB_lim_flux=0.17
    Total_min_flux_MMB=MMB_lim_flux*par.MMB_sigma
    
    xp = [ sample[i]['x(kpc)'] for i in range(len(sample)) ]
    yp = [ sample[i]['y(kpc)'] for i in range(len(sample)) ]
    ic = [ int(sample[i]['ic(Arm)']) for i in range(len(sample)) ]
    lp = [ (sample[i]['l(rad)']*180/math.pi) for i in range(len(sample)) ]
    bp = [ (sample[i]['b(rad)']*180/math.pi) for i in range(len(sample)) ] 
    fp = [ sample[i]['f(Lsun)'] for i in range(len(sample)) ] # in solar luminosities
    f_Jyp_Arecibo = [ sample[i]['f_Arecibo(Jy)'] for i in range(len(sample)) ] # in Jy
    f_Jyp_MMB = [ sample[i]['f_MMB(Jy)'] for i in range(len(sample)) ] # in Jy    
    ra = [ sample[i]['ra(rad)']*180/math.pi for i in range(len(sample)) ]
    dc = [ sample[i]['dc(rad)']*180/math.pi for i in range(len(sample)) ]
    #fpa = [ sample[i]['fa(W/m2)'] for i in range(len(sample)) ]#    
    
    #Arecibo Position Plot
    fig=plt.figure(figsize=(20,10))
    matplotlib.rcParams.update({'font.size': 27})
    ax=fig.add_subplot(211)
    #then plot the spiral structure for the points that Arecibo shuold find
    segs = par.gs.segs_spiral(par.rtrunc, par.nspirplot) # why different than par.nspirseg? OJO
    for j in range( len(segs) ):
        name_arm=['Norma Arm','Carina-Sagitarius Arm','Perseus Arm','Crux-Scutum Arm','Local Arm']
        line_arm_color=['r-','c-','g-','m-','b-']
        plot( segs[j]['x'], segs[j]['y'], line_arm_color[j],label=name_arm[j])
        norma = mlines.Line2D([], [], color='r')
    carina = mlines.Line2D([], [], color='c')
    perseus= mlines.Line2D([], [], color='g')
    crux = mlines.Line2D([], [], color='m')
    local = mlines.Line2D([], [], color='b')    
    legend1=plt.legend([norma,carina,perseus,crux,local],['Norma Arm','Carina-Sagitarius Arm',\
                'Perseus Arm','Crux-Scutum Arm','Local Arm'],prop={'size':16},loc=2)
    for i in range(len(xp) ):
        if(lp[i]>=Arecibo_min_long and lp[i]<=Arecibo_max_long and bp[i]>=Arecibo_min_lat and bp[i]<=Arecibo_max_lat and f_Jyp_Arecibo[i]>Total_min_flux_Arecibo):
            plot( xp[i], yp[i], '.', alpha = 0.4 , color = par.colors[ic[i]], markersize = (numpy.log10(fp[i])*4)+37)
            z=i
    title("Top view of the plane of the Milky Way for Arecibo Survey Limits")
    small_dot, =plot( xp[z], yp[z], '.', alpha = 0.4 , color = par.colors[ic[z]], markersize = 5)
    big_dot,   =plot( xp[z], yp[z], '.', alpha = 0.4 , color = par.colors[ic[z]], markersize = 25)
    legend2=legend([small_dot, big_dot],['$10^{-8} L_{\odot}$','$10^{-3} L_{\odot}$'],loc=0,numpoints=1,prop={'size':22})
    plt.gca().add_artist(legend1)
    plt.gca().add_artist(legend2)
    xlabel('X [$kpc$]')
    ylabel('Y [$kpc$]')
    axis([-par.rtrunc-5,par.rtrunc+5,-par.rtrunc-5,par.rtrunc+5])
    axis('equal')
    plot( 0, par.r0, '*', color ='y',markersize=10)
    ring = par.cr.central_ring(par.nspirplot) 
    plot(ring[0]['x'], ring[0]['y'], 'k--', linewidth = 3)
    text(0.0,0.0, 'GC', fontsize=15)
    plot(0.0,0.0,'o',color='y',markersize=10)
    plt.grid(True)
    
    ay=fig.add_subplot(223)
    map = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')
    title("Right Ascention and Declination view")
    tx=numpy.zeros(len(ra))
    ty=numpy.zeros(len(ra))
    for i in range(len(ra)):
        tx[i], ty[i] = map(ra[i],dc[i])
        if(lp[i]>=Arecibo_min_long and lp[i]<=Arecibo_max_long and bp[i]>=Arecibo_min_lat and bp[i]<=Arecibo_max_lat and f_Jyp_Arecibo[i]>Total_min_flux_Arecibo):
            map.plot( tx[i], ty[i], 'o', alpha = 0.5 , color = par.colors[ic[i]], markersize = (numpy.log10(f_Jyp_Arecibo[i])*2)+9)
            z=i
    ra_rc, dc_rc = map(276,-3)
    ra_lc, dc_lc = map(300,23)
    diff_ra=ra_lc-ra_rc
    diff_dc=dc_lc-dc_rc
    map.drawparallels(numpy.arange(-90.,120.,30.),labels=[True])
    map.drawmeridians(np.arange(0.,420.,60.))
    map.drawmapboundary(fill_color='GhostWhite')
    #small_dot,  =map.plot( tx[z], ty[z], 'o', alpha = 0.5 , color = 'k', markersize = 5.6)
    #big_dot,    =map.plot( tx[z], ty[z], 'o', alpha = 0.5 , color = 'k', markersize = 15)
    rect1 = matplotlib.patches.Rectangle((ra_rc,dc_rc), diff_ra, diff_dc, facecolor='grey',alpha=0.5)
    ay.add_patch(rect1)
    #lg=legend([small_dot, big_dot],['$0.02Jy$','$10^{3} Jy$'],loc=8,numpoints=1,prop={'size':20})

    az=fig.add_subplot(224)
    for i in range( len(ra) ):
        if(lp[i]>=Arecibo_min_long and lp[i]<=Arecibo_max_long and bp[i]>=Arecibo_min_lat and bp[i]<=Arecibo_max_lat and f_Jyp_Arecibo[i]>Total_min_flux_Arecibo):
        #limits of the Arecibo Survey see: Panadian et al, 2007, ApJ 656, 255.
            plot( ra[i], dc[i], 'o', alpha = 0.5 , color = par.colors[ic[i]], markersize =(numpy.log10(f_Jyp_Arecibo[i])*2)+9)
            z=i
    small_dot,  =plot( ra[z], dc[z], 'o', alpha = 0.5 , color = 'k', markersize = 5.6)
    big_dot,    =plot( ra[z], dc[z], 'o', alpha = 0.5 , color = 'k', markersize = 15)
    lg=legend([small_dot, big_dot],['$0.02Jy$','$10^{3} Jy$'],loc=4,numpoints=1,prop={'size':20})
    plt.grid(True)
    xlabel('Right Ascension $(^o)$')
    ylabel('Declination $(^o)$')
    transFigure = fig.transFigure.inverted()
    coord1 = transFigure.transform(ay.transData.transform([ra_lc,dc_rc-(diff_dc/4)]))
    coord2 = transFigure.transform(az.transData.transform([282,0]))
    line1 = matplotlib.lines.Line2D((coord1[0],coord2[0]),(coord1[1],coord2[1]),transform=fig.transFigure,color='k')                         
    coord3 = transFigure.transform(ay.transData.transform([ra_lc,dc_lc-(diff_dc/4)]))
    coord4 = transFigure.transform(az.transData.transform([282,19]))
    line2 = matplotlib.lines.Line2D((coord3[0],coord4[0]),(coord3[1],coord4[1]),transform=fig.transFigure,color='k')         
    fig.lines = line1,line2
    if(save=='yes'):
        plt.savefig('Arecibo_positions.pdf',bbox_inches='tight')
    plt.close()
    
    #Arecibo Flux Comparison Plot
    k=0
    for i in range( len(f_Jyp_Arecibo) ):
        if(lp[i]>=Arecibo_min_long and lp[i]<=Arecibo_max_long and bp[i]>=Arecibo_min_lat and bp[i]<=Arecibo_max_lat and f_Jyp_Arecibo[i]>Total_min_flux_Arecibo):
            k+=1 #calculating the number of masers that can be detected by Arecibo in position and sensitivity
    f_Jy_Arecibo=np.zeros(k) # Emission sources detected by Arecibo in Jy
    w=0
    for i in range( len(f_Jyp_Arecibo) ):
        if(lp[i]>=Arecibo_min_long and lp[i]<=Arecibo_max_long and bp[i]>=Arecibo_min_lat and bp[i]<=Arecibo_max_lat and f_Jyp_Arecibo[i]>Total_min_flux_Arecibo):
            f_Jy_Arecibo[w]=f_Jyp_Arecibo[i]
            w+=1
    print 'Sources that satisfy Arecibo limits= %s (Simulation)' % len(f_Jy_Arecibo)    
    #Arecibo Data in Jansky

    actual_directory=os.getcwd()    
    os.chdir(par.root_surveys)    
    t_Jy = ascii.read('Arecibo_2.txt',format='basic')
    os.chdir(actual_directory)
    z=numpy.zeros(len(t_Jy))
    for i in range(len(t_Jy)):
         z[i]=t_Jy[i][0]
    t=where(z>Total_min_flux_Arecibo)
    Jy_pan=z[t[0]]
    print 'Sources that satisfy Arecibo limits= %s (Observations)' % len(Jy_pan)
    #Histogram Plot
    figure(figsize=(20,10))
    matplotlib.rcParams.update({'font.size': 27}) 
    #flux_ax=fig.add_subplot(111)  
    #title("Observed Flux Density Function for Arecibo Survey")
    #min_data=numpy.log10(min(f_Jy_Arecibo))
    #max_data=numpy.log10(max(f_Jy_Arecibo))
    min_data=numpy.log10(min(Jy_pan))
    max_data=numpy.log10(max(Jy_pan))    
    bins_data=10**(numpy.linspace(min_data,max_data,nbar))
    plt.gca().set_xscale("log")
    #weights_sim=np.ones_like(f_Jy_Arecibo)/len(f_Jy_Arecibo)
    #weights_obs=np.ones_like(Jy_pan)/len(Jy_pan)
    #Simulation
    y_simu,x_simu_aa,z_simu=plt.hist(f_Jy_Arecibo,bins_data,alpha = 0.5,color ='g',label='Simulation',log=True)#,weights=weights_sim) # Histograms of sources detected by Arecibo in Jy
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
    poly_simu = numpy.poly1d(coeffs_simu)
    yfit_simu = lambda x_simu_new: 10**(poly_simu(numpy.log10(x_simu_new)))
    simu_line, =plt.loglog(x_simu_new,yfit_simu(x_simu_new),color='g',label='simu')
    #Data, Observation
    y_data,x_aa,z_data=plt.hist(Jy_pan,bins=bins_data,alpha = 0.5,color='b',label='Data',log=True)#,weights=weights_obs)
    '''
    aa=where(y_data==0)
    y_data[aa]=1
    '''
    x_data=numpy.zeros(len(x_aa)-1)
    for i in range(len(x_aa)-1):
        x_data[i]=(x_aa[i+1]+x_aa[i])/2.0
    remove_zeros_data=[]
    for i in range(len(y_data)):
        if (y_data[i] == 0):
            remove_zeros_data.append(i) 
    y_data_new=numpy.delete(y_data,remove_zeros_data)
    x_data_new=numpy.delete(x_data,remove_zeros_data)
    ####
    log_x_data=numpy.log10(x_data_new)
    log_y_data=numpy.log10(y_data_new)
    coeffs_data = numpy.polyfit(log_x_data,log_y_data,deg=1)
    muuu, buuu, r_valueuuu, p_valueuuu, std_erruuu = stats.linregress(log_x_data, log_y_data)
    muuu2, buuu2, r_valueuuu2, p_valueuuu2, std_erruuu2 = stats.linregress(log_x_simu,log_y_simu)
    print 'm= %s' % muuu 
    print 'err= %s' % std_erruuu
    print 'm2= %s' % muuu2 
    print 'err2= %s' % std_erruuu2 
    poly_data = numpy.poly1d(coeffs_data)
    yfit_data = lambda x_data_new: 10**(poly_data(numpy.log10(x_data_new)))
    data_line, =plt.loglog(x_data_new,yfit_data(x_data_new),color='b',label='dat')
    print 'Index factor of the observed flux density function %s (Simulation)' % coeffs_simu[0]
    print 'Index factor of the observed flux density function %s (Observation)' % coeffs_data[0] 
    y_simu,x_simu_aa,z_simu=plt.hist(f_Jy_Arecibo,bins_data,alpha = 0.5,color ='g',label='Simulation',log=True)#,weights=weights_sim) # Histograms of sources detected by Arecibo in Jy
    y_data,x_aa,z_data=plt.hist(Jy_pan,bins=bins_data,alpha = 0.5,color='b',label='Data',log=True)#,weights=weights_obs)
    scatter(x_simu_new,y_simu_new,color='r')
    scatter(x_data_new,y_data_new,color='k')
    xlabel('Flux Density $[Jy]$')
    ylabel('$N$')
    legend([data_line, simu_line],['Observation','Simulation'],loc=0)#,numpoints=1,prop={'size':20})
    if(save=='yes'):
        plt.savefig('Arecibo_flux.pdf',bbox_inches='tight')
    plt.close()

    #Saving variables for the combined plot    
    bins_data_ar=bins_data
    x_simu_new_ar=x_simu_new
    yfit_simu_ar=yfit_simu
    x_data_new_ar=x_data_new
    yfit_data_ar=yfit_data
    y_simu_new_ar=y_simu_new
    y_data_new_ar=y_data_new
    ############################################################################################################################################
    #####################################################          MMB PLOTS       ###################################################################    
    ############################################################################################################################################
    fig=plt.figure(figsize=(20,10))
    matplotlib.rcParams.update({'font.size': 27}) 
    ax=fig.add_subplot(211)
    #then plot the spiral structure for the points that MMB shuold find
    segs = par.gs.segs_spiral(par.rtrunc, par.nspirplot) # why different than par.nspirseg? OJO
    for j in range( len(segs) ):
        name_arm=['Norma Arm','Carina-Sagitarius Arm','Perseus Arm','Crux-Scutum Arm','Local Arm']
        line_arm_color=['r-','c-','g-','m-','b-']
        plot( segs[j]['x'], segs[j]['y'], line_arm_color[j],label=name_arm[j])
    norma = mlines.Line2D([], [], color='r')
    carina = mlines.Line2D([], [], color='b')
    perseus= mlines.Line2D([], [], color='g')
    crux = mlines.Line2D([], [], color='m')
    local = mlines.Line2D([], [], color='c')    
    legend1=plt.legend([norma,carina,perseus,crux,local],['Norma Arm','Carina-Sagitarius Arm',\
                'Perseus Arm','Crux-Scutum Arm','Local Arm'],prop={'size':16},loc=2)
    for i in range(len(xp) ):
        if(lp[i]>MMB_min_long and lp[i]<MMB_max_long and bp[i]>=MMB_min_lat and bp[i]<=MMB_max_lat and f_Jyp_MMB[i]>(Total_min_flux_MMB)):
            #limits of the MMB Survey see: Green et al 2008, MNRAS 392, 783.
            plot( xp[i], yp[i], '.', alpha = 0.2 , color = par.colors[ic[i]], markersize = (numpy.log10(fp[i])*4)+37)
            z=i
    title("Top view of the plane of the Milky Way for MMB Survey Limits")
    small_dot, =plot( xp[z], yp[z], '.', alpha = 0.4 , color = par.colors[ic[z]], markersize = 5)
    big_dot,   =plot( xp[z], yp[z], '.', alpha = 0.4 , color = par.colors[ic[z]], markersize = 25)
    legend2=legend([small_dot, big_dot],['$10^{-8} L_{\odot}$','$10^{-3} L_{\odot}$'],loc=0,numpoints=1,prop={'size':22})
    plt.gca().add_artist(legend1)
    plt.gca().add_artist(legend2)
    xlabel('X [$kpc$]')
    ylabel('Y [$kpc$]')
    axis([-par.rtrunc-5,par.rtrunc+5,-par.rtrunc-5,par.rtrunc+5])
    axis('equal')
    plot( 0, par.r0, '*', color ='y',markersize=10)
    ring = par.cr.central_ring(par.nspirplot) 
    plot(ring[0]['x'], ring[0]['y'], 'k--', linewidth = 3)
    text(0.0,0.0, 'GC', fontsize=15)
    plot(0.0,0.0,'o',color='y',markersize=10)
    plt.grid(True)
        
    ay=fig.add_subplot(223)
    map = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')
    title("Right Ascention and Declination view")
    tx=numpy.zeros(len(ra))
    ty=numpy.zeros(len(ra))
    for i in range(len(ra)):
        tx[i], ty[i] = map(ra[i],dc[i])
        if(lp[i]>MMB_min_long and lp[i]<MMB_max_long and bp[i]>=MMB_min_lat and bp[i]<=MMB_max_lat and f_Jyp_MMB[i]>(Total_min_flux_MMB)):
             #limits of the Arecibo Survey see: Green et al 2008, MNRAS 392, 783.
            map.plot( tx[i], ty[i], 'o', alpha = 0.5 , color = par.colors[ic[i]], markersize = (numpy.log10(f_Jyp_MMB[i])*2)+9)
            #plot( ra[i], dc[i], 'o', alpha = 0.5 , color = par.colors[ic[i]], markersize =(((((numpy.log10(fpa[i])-3)*(1.5-0.1))/(8-3))+3.6)*10))
            z=i    
    map.drawparallels(numpy.arange(-90.,120.,30.),labels=[True])
    map.drawmeridians(np.arange(0.,420.,60.))
    map.drawmapboundary(fill_color='GhostWhite')
    
    az=fig.add_subplot(224)
    for i in range( len(ra) ):
        if(lp[i]>MMB_min_long and lp[i]<MMB_max_long and bp[i]>=MMB_min_lat and bp[i]<=MMB_max_lat and f_Jyp_MMB[i]>(Total_min_flux_MMB)):
             #limits of the Arecibo Survey see: Green et al 2008, MNRAS 392, 783.
            plot( ra[i], dc[i], 'o', alpha = 0.5 , color = par.colors[ic[i]], markersize =(numpy.log10(f_Jyp_MMB[i])*2)+9)
            z=i
    small_dot,  =plot( ra[z], dc[z], 'o', alpha = 0.5 , color = 'k', markersize = 5.6)
    big_dot,    =plot( ra[z], dc[z], 'o', alpha = 0.5 , color = 'k', markersize = 15)
    lg=legend([small_dot, big_dot],['$0.02Jy$','$10^{3} Jy$'],loc=0,numpoints=1,prop={'size':20})
    plt.grid(True)
    xlabel('Right Ascension $(^o)$')
    ylabel('Declination $(^o)$')
    if(save=='yes'):
        plt.savefig('MMB_positions.pdf',bbox_inches='tight')
    plt.close()
    
    #Flux MMB plot
    figure(figsize=(12,12))
    matplotlib.rcParams.update({'font.size': 27}) 
    #title("Observed Flux Density Function for MMB Survey")
    k=0
    for i in range( len(f_Jyp_MMB) ):
        if(lp[i]>MMB_min_long and lp[i]<MMB_max_long and bp[i]>=MMB_min_lat and bp[i]<=MMB_max_lat and f_Jyp_MMB[i]>(Total_min_flux_MMB)):
            #limits of the Arecibo Survey see: Green et al 2008, MNRAS 392, 783.
            k+=1 #calculating the number of masers that can be detected by Arecibo in position adn sensitivity
    f_Jy_MMB=np.zeros(k) # Emission sources detected by Arecibo in Jy
    w=0
    for i in range( len(f_Jyp_MMB) ):
        if(lp[i]>MMB_min_long and lp[i]<MMB_max_long and bp[i]>=MMB_min_lat and bp[i]<=MMB_max_lat and f_Jyp_MMB[i]>(Total_min_flux_MMB)):
            f_Jy_MMB[w]=f_Jyp_MMB[i]
            w+=1
    print 'Sources that satisfy MMB limits= %s (Simulation)' % len(f_Jy_MMB)
    #MMB Data in Jansky
    #Previous catalogue
    #t_Jy = ascii.read('/Users/luishenryquiroganunez/Documents/Leiden_University/3_Semester/MSc_thesis/Surveys_Data/MMB/MMB_Cat_I_to_IV_edited.txt')
    #Current Catalogue
    actual_directory=os.getcwd()    
    os.chdir(par.root_surveys)
    t_Jy = ascii.read('MMB_Dist_Flux_edited.csv')    
    os.chdir(actual_directory)
        

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
    print 'Sources that satisfy MMB limits= %s (Observations)' % len(Jy_MMB)

    #Histogram Plot of Simulation based on Observations limits
    min_data=numpy.log10(min(Jy_MMB))
    max_data=numpy.log10(max(Jy_MMB))    
    bins_data=10**(numpy.linspace(min_data,max_data,nbar))
    plt.gca().set_xscale("log")
    #weights_sim=np.ones_like(f_Jy_MMB)/len(f_Jy_MMB)
    #weights_obs=np.ones_like(Jy_MMB)/len(Jy_MMB)
    y_simu,x_simu_aa,z_simu=plt.hist(f_Jy_MMB,bins_data,alpha = 0.5,color ='g',label='MMB Simulation',log=True)#,weights=weights_sim) # Histograms of sources detected by MMB in Jy
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
    poly_simu = numpy.poly1d(coeffs_simu)
    yfit_simu = lambda x_simu_new: 10**(poly_simu(numpy.log10(x_simu_new)))
    simu_line, =plt.loglog(x_simu_new,yfit_simu(x_simu_new),color='g',label='simu')
    
    #Data, Observation
    y_data,x_aa,z_data=plt.hist(Jy_MMB,bins=bins_data,alpha = 0.5,color='b',label='MMB Data',log=True)#,weights=weights_obs) 
    x_data=numpy.zeros(len(x_aa)-1)
    for i in range(len(x_aa)-1):
        x_data[i]=(x_aa[i+1]+x_aa[i])/2.0
    remove_zeros_data=[]
    for i in range(len(y_data)):
        if (y_data[i] == 0):
            remove_zeros_data.append(i) 
    y_data_new=numpy.delete(y_data,remove_zeros_data)
    x_data_new=numpy.delete(x_data,remove_zeros_data)
    log_x_data=numpy.log10(x_data_new)
    log_y_data=numpy.log10(y_data_new)
    coeffs_data = numpy.polyfit(log_x_data,log_y_data,deg=1)
    muuu, buuu, r_valueuuu, p_valueuuu, std_erruuu = stats.linregress(log_x_data, log_y_data)
    muuu2, buuu2, r_valueuuu2, p_valueuuu2, std_erruuu2 = stats.linregress(log_x_simu,log_y_simu)
    print 'm= %s' % muuu 
    print 'err= %s' % std_erruuu
    print 'm2= %s' % muuu2 
    print 'err2= %s' % std_erruuu2  
      
    poly_data = numpy.poly1d(coeffs_data)
    yfit_data = lambda x_data_new: 10**(poly_data(numpy.log10(x_data_new)))
    data_line, =plt.loglog(x_data_new,yfit_data(x_data_new),color='b',label='dat')
    print 'Index factor of the observed flux density function %s (Simulation)' % coeffs_simu[0]
    print 'Index factor of the observed flux density function %s (Observation)' % coeffs_data[0] 
    y_simu,x_simu_aa,z_simu=plt.hist(f_Jy_MMB,bins_data,alpha = 0.5,color ='g',label='MMB Simulation',log=True)#,weights=weights_sim) # Histograms of sources detected by MMB in Jy
    y_data,x_aa,z_data=plt.hist(Jy_MMB,bins=bins_data,alpha = 0.5,color='b',label='MMB Data',log=True)#,weights=weights_obs) 
    scatter(x_simu_new,y_simu_new,color='r',alpha=1.0)
    scatter(x_data_new,y_data_new,color='k',alpha=1.0)
    xlabel('Flux Density (Jy)')
    ylabel('Source Counts')
    legend([data_line, simu_line],['Observation','Simulation'],loc=0)#,numpoints=1,prop={'size':20})
    if(save=='yes'):
        plt.savefig('MMB_flux.pdf',bbox_inches='tight')
    plt.close()
    ############################################################################################################################################
    #####################################################          COMBINED HISTOHGRAM PLOTS       ###################################################################    
    ############################################################################################################################################    
    
    #Histogram Plot
    matplotlib.rcParams.update({'font.size': 27}) 
    fig5 = plt.figure(figsize=(12,12))
    fig5.subplots_adjust(hspace=0.0)
    ax  = fig5.add_subplot(111)
    ax2 = fig5.add_subplot(211)
    ax1 = fig5.add_subplot(212,sharex=ax2)
    ax1.set_xscale("log")    
    ax.set_xlabel('Peak Flux Density (Jy)')
    ax.set_ylabel('Source Counts')
    ax.xaxis.labelpad = 40
    ax.yaxis.labelpad = 65
    plt.setp( ax.get_xticklabels(), visible=False)
    plt.setp( ax.get_yticklabels(), visible=False)

    #MMB
    #y_simu,x_simu_aa,z_simu=ax2.hist(f_Jy_MMB,bins_data_ar,alpha = 0.5,color ='g',label='MMB Simulation',log=True)#,weights=weights_sim) # Histograms of sources detected by MMB in Jy
    y_simu,x_simu_aa,z_simu=ax2.hist(f_Jy_MMB,bins_data_ar,histtype='step', edgecolor='green',linewidth=5.5,label='MMB Simulation',log=True)#,weights=weights_sim) # Histograms of sources detected by MMB in Jy
    simu_line, =ax2.loglog(x_simu_new-10000,yfit_simu(x_simu_new),color='g',label='simu',linewidth=10)
    y_data,x_aa,z_data=ax2.hist(Jy_MMB,bins=bins_data_ar,alpha = 0.2,color='b',label='MMB Data',log=True)#,weights=weights_obs) 
    data_line, =ax2.loglog(x_data_new-10000,yfit_data(x_data_new),alpha = 0.4,color='b',label='dat',linewidth=10)
    #y_simu,x_simu_aa,z_simu=ax2.hist(f_Jy_MMB,bins_data_ar,alpha = 0.5,color ='g',label='MMB Simulation',log=True)#,weights=weights_sim) # Histograms of sources detected by MMB in Jy
    #y_simu,x_simu_aa,z_simu=ax2.hist(f_Jy_MMB,bins_data_ar,histtype='step', edgecolor='green',linewidth=5.5,label='MMB Simulation',log=True)#,weights=weights_sim) # Histograms of sources detected by MMB in Jy
    y_data,x_aa,z_data=ax2.hist(Jy_MMB,bins=bins_data_ar,alpha = 0.2,color='b',label='MMB Data',log=True)#,weights=weights_obs) 
    #ax2.scatter(x_simu_new,y_simu_new,color='r',alpha=1.0)
    #ax2.scatter(x_data_new,y_data_new,color='k',alpha=1.0)
    ax2.set_xlim([0.3,1000])
    ax2.set_ylim([1,400])
    #ax2.legend([data_line, simu_line],['Observation','Simulation'],loc=1,ncol=1,prop={'size':23})#,numpoints=1,prop={'size':20})
    plt.setp( ax2.get_xticklabels(), visible=False)

    #Arecibo
#    y_simu_ar,x_simu_aa_ar,z_simu_ar=ax1.hist(f_Jy_Arecibo,bins_data_ar,alpha = 0.5,color ='g',label='Simulation',log=True)#,weights=weights_sim) # Histograms of sources detected by Arecibo in Jy
    y_simu_ar,x_simu_aa_ar,z_simu_ar=ax1.hist(f_Jy_Arecibo,bins_data_ar,histtype='step', edgecolor='green',linewidth=5.5,label='Simulation',log=True)#,weights=weights_sim) # Histograms of sources detected by Arecibo in Jy
    #simu_line_ar, =ax1.loglog(x_simu_new_ar,yfit_simu_ar(x_simu_new_ar),color='g',label='simu')
    y_data_ar,x_aa_ar,z_data_ar=ax1.hist(Jy_pan,bins=bins_data_ar,alpha = 0.2,color='b',label='Data',log=True)#,weights=weights_obs)
    #data_line_ar, =ax1.loglog(x_data_new_ar,yfit_data_ar(x_data_new_ar),color='b',label='dat')
    #y_simu_ar,x_simu_aa_ar,z_simu_ar=ax1.hist(f_Jy_Arecibo,bins_data_ar,alpha = 0.5,color ='g',label='Simulation',log=True)#,weights=weights_sim) # Histograms of sources detected by Arecibo in Jy
    y_simu_ar,x_simu_aa_ar,z_simu_ar=ax1.hist(f_Jy_Arecibo,bins_data_ar,histtype='step', edgecolor='green',linewidth=3.5,label='Simulation',log=True)#,weights=weights_sim) # Histograms of sources detected by Arecibo in Jy
    y_data_ar,x_aa_ar,z_data_ar=ax1.hist(Jy_pan,bins=bins_data_ar,alpha = 0.2,color='b',label='Data',log=True)#,weights=weights_obs)
   
    ax1.set_xlim([0.3,1000])
    ax1.set_ylim([0.1,200])
    ax1.legend([data_line, simu_line],['Observation','Simulation'],loc=9,ncol=2,prop={'size':23})#,numpoints=1,prop={'size':20})    
    ax2.text( 300, 120, 'MMB',     verticalalignment='center', horizontalalignment='center',     color='r', fontsize=25)
    ax1.text( 300, 10, 'ARECIBO',     verticalalignment='center', horizontalalignment='center',     color='r', fontsize=25)

    if(save=='yes'):
        plt.savefig('Comp_MMBAr_flux.pdf',bbox_inches='tight')
    plt.close()
    plt.ion()
    return
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
def AVN_plot(sample,par,save='no'):
    if(save=='yes'):
        #plt.ioff()
        pass
    else:
        print('The plots were NOT saved!')
    
    xp = [ sample[i]['x(kpc)'] for i in range(len(sample)) ]
    yp = [ sample[i]['y(kpc)'] for i in range(len(sample)) ]
    zp = [ sample[i]['z(kpc)'] for i in range(len(sample)) ]
    dc = [ sample[i]['dc(rad)']*180/math.pi for i in range(len(sample))]
    d_sun = [ sample[i]['dsun(kpc)'] for i in range(len(sample))]
    err_mu_x = [ sample[i]['err_mu_x(mas/yr)'] for i in range(len(sample))]    
    f_jy_VLBA = [ sample[i]['f_jy_VLBA(Jy)'] for i in range(len(sample))]    

    
    '''
    #Only arms + sources
    matplotlib.rcParams.update({'font.size': 27}) 
    figure(figsize=(12,12))
    segs = par.gs.segs_spiral(par.rtrunc, par.nspirplot) # why different than par.nspirseg? OJO
    for j in range( len(segs) ):
        name_arm=['Norma Arm','Carina-Sagitarius Arm','Perseus Arm','Crux-Scutum Arm','Local Arm']
        line_arm_color=['k-','k-','k-','k-','k-']
        plot( segs[j]['x'], segs[j]['y'], line_arm_color[j], linewidth=3)
    norma = mlines.Line2D([], [], color='r')
    carina = mlines.Line2D([], [], color='c')
    perseus= mlines.Line2D([], [], color='g')
    crux = mlines.Line2D([], [], color='m')
    local = mlines.Line2D([], [], color='b')    
    text( 12.0,  12.0, 'Norma',     verticalalignment='bottom', horizontalalignment='right',     color='k', fontsize=25)
    text( -8.0, -14.0, 'Carina \n Sagitarius',     verticalalignment='bottom', horizontalalignment='right',     color='k', fontsize=25)
    text( -2.0,  11.0, 'Perseus',     verticalalignment='bottom', horizontalalignment='right',     color='k', fontsize=25)
    text( 13.0, -13.0, 'Crux-Scutum',     verticalalignment='bottom', horizontalalignment='right',     color='k', fontsize=25)
    text( -8.0, 5.5, 'Local',     verticalalignment='bottom', horizontalalignment='right',     color='k', fontsize=25)
    for i in range(len(xp) ):
        if(dc[i]>=20):
            plot( xp[i], yp[i], '.', alpha = 0.6 , color ='r', markersize =10.0)
        elif(dc[i]<=-20):
            plot( xp[i], yp[i], '.', alpha = 0.6 , color ='g', markersize =10.0)            
        else:
            plot( xp[i], yp[i], '.', alpha = 0.6 , color ='b', markersize =10.0)                        
    plot(100, 100, '.', color ='g', markersize =20.0,label='$\\delta \\leq -20^{\\degree} $')    
    plot(100, 100, '.', color ='b', markersize =20.0,label='$ -20^{\\degree} \\! < \\delta < 20^{\\degree}$')    
    plot(100, 100, '.', color ='r', markersize =20.0,label='$\\delta \\geq 20^{\\degree} $')    
    legend(loc=9,numpoints=1,ncol=3,columnspacing=0.8,prop={'size':27})
    xlabel('X $[kpc]$')
    ylabel('Y $[kpc]$')
    axis([-par.rtrunc-5.0,par.rtrunc+5.0,-par.rtrunc-5.0,par.rtrunc+5.0])
    axis('equal')
    xlim([-15.2,15.2])
    ylim([-15.5,18.5])
    plot( 0, par.r0, '*', color ='y',markersize=15)
    ring = par.cr.central_ring(par.nspirplot) 
    plot(ring[0]['x'], ring[0]['y'], 'm--', linewidth = 3)
    text(0.0,0.0, 'GC', fontsize=15)
    plot(0.0,0.0,'o',color='g',markersize=10)
    plt.grid(True)
    name_plot_3="../plots/AVN.pdf"
    if(save=='yes'):
        plt.savefig(name_plot_3,bbox_inches='tight')
    plt.close()
    
    '''
    #Histograms N vs radial (exp decay)        
    n_north=[]
    n_south=[]
    n_mid=[]
    nbar_d=np.arange(0.0,25.0,1.0)
    for i in range(len(d_sun)):
        if(dc[i]>=20):
            n_north.append(i)
        elif(dc[i]<=-20):
            n_south.append(i)
        else:
            n_mid.append(i)
    d_north=[sample[k]['dsun(kpc)'] for k in n_north]
    d_south=[sample[k]['dsun(kpc)'] for k in n_south]
    d_mid=[sample[k]['dsun(kpc)'] for k in n_mid]
    
    '''
    
    number_south,bins=numpy.histogram(d_south,nbar_d)
    number_north,bins=numpy.histogram(d_north,nbar_d)
    number_mid,bins=numpy.histogram(d_mid,nbar_d)

    bins_new=numpy.zeros(len(bins)-1)
    for i in range(len(bins)-1):
        bins_new[i]=(bins[i]+bins[i+1])/2.0
    
    matplotlib.rcParams.update({'font.size': 27})       
    fig3 = plt.figure(figsize=(15,8))  
    ax1 = fig3.add_subplot(1,1,1)
    width=1.0
    ax1.bar(bins_new, number_north, width, color='r',bottom=number_south+number_mid,label='$\\delta \\geq 20^{\\degree} $')        
    ax1.bar(bins_new, number_mid, width, color='b',bottom=number_south,label='$ -20^{\\degree} \\! < \\delta < 20^{\\degree}$')
    ax1.bar(bins_new, number_south, width, color='g',label='$\\delta \\leq -20^{\\degree} $')    
    ax1.set_xlabel('Solar Distance $[kpc]$')
    ax1.set_ylabel('Number of Sources')
    ax1.grid(True)
    ax1.legend()
    name_plot_2="../plots/AVN_hist.pdf"
    if(save=='yes'):
        plt.savefig(name_plot_2,bbox_inches='tight')
    plt.close()   
      

    distance_break_ranges=25
    step=25.0/distance_break_ranges
    array=[]
    array.append(0)
   
    for i in range(int(distance_break_ranges)):
        array.append(step*(i+1))
   
    ranges_blocks={}
    ranges_blocks_flux={}
    ranges_blocks_errmux={}
    for i in range(len(array)-1):
        ranges_blocks[i]=[]
        ranges_blocks_flux[i]=[]
        ranges_blocks_errmux[i]=[]
        
    for i in range(len(sample)):
        for j in range(len(array)-1):        
            if(d_sun[i]>=array[j] and d_sun[i]<array[j+1]):
                ranges_blocks[j].append(i)

    for i in range(len(ranges_blocks)):
        a=len(ranges_blocks[i])
        for j in range(a):
            ranges_blocks_flux[i].append(f_jy_VLBA[ranges_blocks[i][j]])
            ranges_blocks_errmux[i].append(err_mu_x[ranges_blocks[i][j]])        
            
    fluxes_median=numpy.zeros(distance_break_ranges)
    errors_median=numpy.zeros(distance_break_ranges)
    fluxes_mean=numpy.zeros(distance_break_ranges)
    errors_mean=numpy.zeros(distance_break_ranges)
        
    for i in range(len(ranges_blocks)):
        k=np.median(ranges_blocks_flux[i])
        kk=np.mean(ranges_blocks_flux[i])
        fluxes_median[i]=k
        fluxes_mean[i]=kk      
        l=np.median(ranges_blocks_errmux[i])
        ll=np.mean(ranges_blocks_errmux[i])
        errors_median[i]=l
        errors_mean[i]=ll
            
            
    where_are_NaNs = np.isnan(fluxes_median)
    fluxes_median[where_are_NaNs] = 0
    errors_median[where_are_NaNs] = 0    
    fluxes_mean[where_are_NaNs] = 0
    errors_mean[where_are_NaNs] = 0  
    
    bins=range(distance_break_ranges+1)
    bins_new=numpy.zeros(len(bins)-1)
    for i in range(len(bins)-1):
        bins_new[i]=(bins[i]+bins[i+1])/2.0  
    
    matplotlib.rcParams.update({'font.size': 27})       
    fig3 = plt.figure(figsize=(15,8))  
    ax1 = fig3.add_subplot(1,1,1)
    width=1.0        
    ax1.bar(bins_new, fluxes_median, width, color='r',alpha=0.5, log=True,label='Median')        
    ax1.bar(bins_new, fluxes_mean, width, color='b',alpha=0.5, log=True,label='Mean')        
    ax1.set_xlabel('Solar Distance $[kpc]$')
    ax1.set_ylabel('Observed Flux $[Jy]$')
    ax1.grid(True)
    ax1.legend()
    name_plot_3="../plots/AVN_hist_fluxes.pdf"
    if(save=='yes'):
        plt.savefig(name_plot_3,bbox_inches='tight')    
    plt.close()   
        
    matplotlib.rcParams.update({'font.size': 27})       
    fig3 = plt.figure(figsize=(15,8))  
    ax1 = fig3.add_subplot(1,1,1)
    width=1.0        
    ax1.bar(bins_new, errors_median, width, color='r',alpha=0.5, log=True,label='Median')        
    ax1.bar(bins_new, errors_mean, width, color='b',alpha=0.5, log=True,label='Mean')        
    ax1.set_xlabel('Solar Distance $[kpc]$')
    ax1.set_ylabel('Error in Proper Motion $[mas/yr]$')
    ax1.grid(True)
    ax1.legend()        
    name_plot_4="../plots/AVN_hist_errmu.pdf"
    if(save=='yes'):
        plt.savefig(name_plot_4,bbox_inches='tight')    
    plt.close()   
    '''
    
    flux_north=[sample[k]['f_jy_VLBA(Jy)'] for k in n_north]
    flux_south=[sample[k]['f_jy_VLBA(Jy)'] for k in n_south]
    flux_mid=[sample[k]['f_jy_VLBA(Jy)'] for k in n_mid]
    
    matplotlib.rcParams.update({'font.size': 27})       
    fig5 = plt.figure(figsize=(15,8))  
    ax5 = fig5.add_subplot(1,1,1)
    ax5.plot( d_north, flux_north, '.' , color='r',label='$\\delta \\geq 20^{\\degree} $')        
    ax5.plot(d_mid, flux_mid,  '.', color='b',label='$ -20^{\\degree} \\! < \\delta < 20^{\\degree}$')
    ax5.plot(d_south, flux_south, '.', color='g',label='$\\delta \\leq -20^{\\degree} $')         
    ax5.set_xlabel('Solar Distance $[kpc]$')
    ax5.set_ylabel('Observed Flux $[Jy]$')
    ax5.set_yscale('log')
    ax5.grid(True)
    lgnd =ax5.legend(numpoints=1)
    lgnd.legendHandles[0]._legmarker.set_markersize(20)
    lgnd.legendHandles[1]._legmarker.set_markersize(20)
    lgnd.legendHandles[2]._legmarker.set_markersize(20)    
    name_plot_5="../plots/AVN_fluxes.pdf"
    if(save=='yes'):
        plt.savefig(name_plot_5,bbox_inches='tight')    
    plt.close() 
    '''
    Vt=10 # km/s dummy velocity
    proper_motion_dummy=numpy.zeros(len(sample))
    for i in range(len(sample)):
        proper_motion_dummy[i]=Vt/ (4.75 * d_sun[i]*1000)
    
    proper_motion_north=[proper_motion_dummy[k]*1000 for k in n_north] #mas / yr
    proper_motion_south=[proper_motion_dummy[k]*1000 for k in n_south] #mas / yr
    proper_motion_mid=[proper_motion_dummy[k]*1000 for k in n_mid] #mas / yr
    
    
    err_mu_north = [ numpy.sqrt((sample[i]['err_mu_x(mas/yr)'])**2+(sample[i]['err_mud(mas/yr)'])**2) for i in n_north]    
    err_mu_south = [ numpy.sqrt((sample[i]['err_mu_x(mas/yr)'])**2+(sample[i]['err_mud(mas/yr)'])**2) for i in n_south]    
    err_mu_mid = [ numpy.sqrt((sample[i]['err_mu_x(mas/yr)'])**2+(sample[i]['err_mud(mas/yr)'])**2) for i in n_mid]    

    mu_north = [ numpy.sqrt((sample[i]['mux(mas/yr)'])**2+(sample[i]['mud(mas/yr)'])**2) for i in n_north]    
    mu_south = [ numpy.sqrt((sample[i]['mux(mas/yr)'])**2+(sample[i]['mud(mas/yr)'])**2) for i in n_south]    
    mu_mid = [ numpy.sqrt((sample[i]['mux(mas/yr)'])**2+(sample[i]['mud(mas/yr)'])**2) for i in n_mid]    
    
            
    snr_north=numpy.zeros(len(d_north))
    snr_south=numpy.zeros(len(d_south))    
    snr_mid=numpy.zeros(len(d_mid))
    
    for i in range(len(d_north)):
        #snr_north[i]=proper_motion_north[i]/err_mu_x_north[i]
        snr_north[i]=mu_north[i]/err_mu_north[i]
        
    for i in range(len(d_south)):
        #snr_south[i]=proper_motion_south[i]/err_mu_x_south[i]
        snr_south[i]=mu_south[i]/err_mu_south[i]
                
    for i in range(len(d_mid)):
        #snr_mid[i]=proper_motion_mid[i]/err_mu_x_mid[i]
        snr_mid[i]=mu_mid[i]/err_mu_mid[i]
                    
    matplotlib.rcParams.update({'font.size': 27})       
    fig3 = plt.figure(figsize=(15,8))  
    ax1 = fig3.add_subplot(1,1,1)
    #ax1.plot( d_sun, proper_motion_dummy/err_mu_x, '.' , color='k')        
    ax1.plot( d_north, snr_north , '.' , color='r',label='$\\delta \\geq 20^{\\degree} $')        
    ax1.plot(d_mid, snr_mid,  '.', color='b',label='$ -20^{\\degree} \\! < \\delta < 20^{\\degree}$')
    ax1.plot(d_south, snr_south, '.', color='g',label='$\\delta \\leq -20^{\\degree} $')         
    ax1.set_xlabel('Solar Distance $[kpc]$')
    ax1.set_ylabel('SNR for proper motion')
    ax1.set_yscale('log')
    ax1.grid(True)
    lgnd =ax1.legend(numpoints=1,loc=0)
    lgnd.legendHandles[0]._legmarker.set_markersize(20)
    lgnd.legendHandles[1]._legmarker.set_markersize(20)
    lgnd.legendHandles[2]._legmarker.set_markersize(20)   
    name_plot_6="../plots/AVN_err.pdf"
    if(save=='yes'):
        plt.savefig(name_plot_6,bbox_inches='tight')    
    plt.close() 
    
    '''
    para_north=[sample[k]['parallax(mas)'] for k in n_north]
    para_south=[sample[k]['parallax(mas)'] for k in n_south]
    para_mid=[sample[k]['parallax(mas)'] for k in n_mid]    

    err_para_north=[sample[k]['err_parallax(mas)'] for k in n_north]
    err_para_south=[sample[k]['err_parallax(mas)'] for k in n_south]    
    err_para_mid=[sample[k]['err_parallax(mas)'] for k in n_mid]    
        
    snr_para_north=[para_north[k]/err_para_north[k] for k in range(len(n_north))]
    snr_para_south=[para_south[k]/err_para_south[k] for k in range(len(n_south))]
    snr_para_mid=[para_mid[k]/err_para_mid[k] for k in range(len(n_mid))]
             
    matplotlib.rcParams.update({'font.size': 27})       
    fig7 = plt.figure(figsize=(15,8))  
    ax7 = fig7.add_subplot(1,1,1)
    #ax1.plot( d_sun, proper_motion_dummy/err_mu_x, '.' , color='k')        
    ax7.plot( d_north, snr_para_north , '.' , color='r',label='$\\delta \\geq 20^{\\degree} $')        
    ax7.plot(d_mid, snr_para_mid,  '.', color='b',label='$ -20^{\\degree} \\! < \\delta < 20^{\\degree}$')
    ax7.plot(d_south, snr_para_south, '.', color='g',label='$\\delta \\leq -20^{\\degree} $')         
    ax7.set_xlabel('Solar Distance $[kpc]$')
    ax7.set_ylabel('SNR for parallax')
    ax7.set_yscale('log')
    ax7.grid(True)
    lgnd =ax7.legend(numpoints=1,loc=0)
    lgnd.legendHandles[0]._legmarker.set_markersize(20)
    lgnd.legendHandles[1]._legmarker.set_markersize(20)
    lgnd.legendHandles[2]._legmarker.set_markersize(20)   
    name_plot_7="../plots/AVN_para.pdf"
    if(save=='yes'):
        plt.savefig(name_plot_7,bbox_inches='tight')    
    plt.close() 
    
    plt.ion()
    return
    
def Parrallax_Plots(sample,par,save='no'):
    '''
    Parallax-Distance histograms (distributions)
    '''    
    if(save=='yes'):
        #plt.ioff()
        pass
    else:
        print('The plots were NOT saved!')
    #plt.clf()  
    dsun = [ sample[i]['dsun(kpc)'] for i in range(len(sample)) ]   
    parralax_obs = [ sample[i]['parallax_obs(mas)'] for i in range(len(sample)) ]   
    parralax = [ sample[i]['parallax(mas)'] for i in range(len(sample)) ]       
    err_parralax = [ sample[i]['err_parallax(mas)'] for i in range(len(sample)) ]       
    sigma_perce_obs = []
    sigma_perce = []
    for i in range(len(parralax)):
        sigma_perce.append(err_parralax[i]*100/parralax[i])    
        sigma_perce_obs.append(err_parralax[i]*100/parralax_obs[i])    
        
    #Histograms 
    matplotlib.rcParams.update({'font.size': 22}) 
    fig5 = plt.figure(figsize=(24,8))
    fig5.subplots_adjust(hspace=0.0)
    ax  = fig5.add_subplot(111)
    ax1 = fig5.add_subplot(131)
    ax2 = fig5.add_subplot(132)
    ax3 = fig5.add_subplot(133)
    #ax1.set_xscale("log")    
    #ax.set_xlabel('Peak Flux Density (Jy)')
    ax.set_ylabel('Source Counts')
    ax.xaxis.labelpad = 40
    ax.yaxis.labelpad = 65
    nbins=100
    plt.setp( ax.get_xticklabels(), visible=False)
    plt.setp( ax.get_yticklabels(), visible=False)
    
    ##Histogram Distance
    a,b,c=ax1.hist(dsun,nbins,color ='g')
    ax1.set_xlim([0,20])
    ax1.set_xlabel('Solar Distance (kpc)')    
    #Parallax Histo
    x,y,z=ax2.hist(parralax,nbins,color ='b')
    ax2.set_xlabel('Parallax (mas)')    
    ax2.set_xlim([-1,2])        
    #Parallax Histo
    x,y,z=ax3.hist(parralax_obs,nbins,color ='r')
    ax3.set_xlabel('Parallax Obs(mas)')    
    ax3.set_xlim([-1,2])    

    if(save=='yes'):
        plt.savefig('dist-para.pdf',bbox_inches='tight')    
    plt.close()
    
    matplotlib.rcParams.update({'font.size': 22}) 
    fig5 = plt.figure(figsize=(24,8))
    fig5.subplots_adjust(hspace=0.0)
    ax  = fig5.add_subplot(111)
    ax1 = fig5.add_subplot(121)
    ax2 = fig5.add_subplot(122)
    #ax1.set_xscale("log")    
    #ax.set_xlabel('Peak Flux Density (Jy)')
    ax.set_ylabel('Source Counts')
    ax.xaxis.labelpad = 40
    ax.yaxis.labelpad = 65
    nbins=100
    plt.setp( ax.get_xticklabels(), visible=False)
    plt.setp( ax.get_yticklabels(), visible=False)
    
    ##Percentage err/parallax
    a,b,c=ax1.hist(sigma_perce,nbins,color ='g')
    #ax1.set_xlim([0,20])
    ax1.set_xlabel('sigma/parallax (%)')    
    ##Percentage err/parallax_obs
    x,y,z=ax2.hist(sigma_perce_obs,nbins,color ='b')
    ax2.set_xlabel('sigma/parallax_obs (%)')    
    #ax2.set_xlim([-1,2])        
         
    if(save=='yes'):
        plt.savefig('sigpara.pdf',bbox_inches='tight')    
    plt.close()
    plt.ion()
    return
    
def low_erros(sample,par,save='no'):
    '''
    Plot only the soruces with err parallax<10%
    '''
    if(save=='yes'):
        #plt.ioff()
        pass
    else:
        print('The plots were NOT saved!')
    #plt.clf()  
    xp = [ sample[i]['x(kpc)'] for i in range(len(sample)) ]
    yp = [ sample[i]['y(kpc)'] for i in range(len(sample)) ]
    ic = [ int(sample[i]['ic(Arm)']) for i in range(len(sample)) ]
    para = [ sample[i]['parallax_obs(mas)'] for i in range(len(sample)) ]
    err_para = [ sample[i]['err_parallax(mas)'] for i in range(len(sample)) ]
    fp = [ sample[i]['f(Lsun)'] for i in range(len(sample)) ]
    rp = []
    for i in range(len(xp)):
        rp.append( sqrt(xp[i]**2+yp[i]**2) )     
    #Plot xy with lumosities 
    matplotlib.rcParams.update({'font.size': 27}) 
    figure(figsize=(10,10))
    segs = par.gs.segs_spiral(par.rtrunc, par.nspirplot) # why different than par.nspirseg? OJO
    for j in range( len(segs) ):
        name_arm=['Norma Arm','Carina-Sagitarius Arm','Perseus Arm','Crux-Scutum Arm','Local Arm']
        line_arm_color=['r-','c-','g-','m-','b-']
        #for k in range(len(segs[0]['x'])):
        #    segs[j]['y'][k]=segs[j]['y'][k]*-1
        plot( segs[j]['x'], segs[j]['y'], line_arm_color[j],label=name_arm[j], linewidth=3)
        #plot(segs[j]['y'], segs[j]['x'],  line_arm_color[j],label=name_arm[j])
    norma = mlines.Line2D([], [], color='r')
    carina = mlines.Line2D([], [], color='c')
    perseus= mlines.Line2D([], [], color='g')
    crux = mlines.Line2D([], [], color='m')
    local = mlines.Line2D([], [], color='b')    
    text( 12.0,  12.0, 'Norma',     verticalalignment='bottom', horizontalalignment='right',     color='red', fontsize=25)
    text( -8.0, -14.0, 'Carina \n Sagitarius',     verticalalignment='bottom', horizontalalignment='right',     color='c', fontsize=25)
    text( -2.0,  11.0, 'Perseus',     verticalalignment='bottom', horizontalalignment='right',     color='green', fontsize=25)
    text( 14.0, -12.0, 'Crux-Scutum',     verticalalignment='bottom', horizontalalignment='right',     color='m', fontsize=25)
    text( -8.0, 5.5, 'Local',     verticalalignment='bottom', horizontalalignment='right',     color='b', fontsize=25)
    text( -0.8, -2.1, 'Central \n Ring', verticalalignment='bottom', horizontalalignment='center', color='k', fontsize=16)    
    kk=0
    for i in range(len(xp) ):
        a=err_para[i]/para[i]
        if(a<=0.1):
            plot( xp[i], yp[i], '.', alpha = 0.4 , color = par.colors[ic[i]], markersize = (numpy.log10(fp[i])*4)+37)
            #plot( -yp[i], xp[i], '.', alpha = 0.4 , color = par.colors[ic[i]], markersize = (numpy.log10(fp[i])*4)+37)
            z=i
            kk+=1
        elif(a>0.1):
            plot( xp[i], yp[i], '.', alpha = 0.8, color = 'k', markersize = (numpy.log10(fp[i])*4)+37)
    print ('Number of sources in the sample:%s' % kk)
    #small_dot, =plot( -yp[z], xp[z],  '.', alpha = 0.4 , color = par.colors[ic[z]], markersize = 5)
    #big_dot,   =plot( -yp[z], xp[z], '.', alpha = 0.4 , color = par.colors[ic[z]], markersize = 25)
    small_dot, =plot( xp[z], yp[z], '.', alpha = 0.4 , color = par.colors[ic[z]], markersize = 5)
    big_dot,   =plot( xp[z], yp[z], '.', alpha = 0.4 , color = par.colors[ic[z]], markersize = 25)
    legend2=legend([small_dot, big_dot],['$10^{-8} L_{\odot}$','$10^{-3} L_{\odot}$'],loc=4,ncol=2,numpoints=1,prop={'size':23},borderpad=0.2)
    plt.gca().add_artist(legend2) 
    ring = par.cr.central_ring(par.nspirplot) 
    plot(ring[0]['x'], ring[0]['y'], 'k--', linewidth = 3)
    text(0.0,0.0, 'GC', fontsize=15)  
    plot(0.0,0.0,'o',color='g',markersize=10)
    xlabel('X $[kpc]$')
    ylabel('Y $[kpc]$')
    axis([-par.rtrunc-10.0,par.rtrunc+10.0,-par.rtrunc-10.0,par.rtrunc+10.0])
    axis('equal')
    plot( 0, par.r0, '*', color ='y',markersize=13)
    #plot( -par.r0, 0, '*', color ='y',markersize=10)
    plt.grid(True)
    name_plot_0="samplot2.pdf"
    if(save=='yes'):
        plt.savefig(name_plot_0,bbox_inches='tight')
        #print(name_plot_0)
    plt.close()    
    plt.ion()
       
def plot_sample_subtraction(root_filename1,root_filename2,par,save='yes'):
    '''
    Plot x and y for extra soruces at different samples
    '''
    from operator import truediv    
    #root_filename1='/Users/luishenryquiroganunez/Documents/Leiden_University/3_Semester/MSc_thesis/Codes/Code_v30/galaxy_generator/output/si30_moresrcs_gpara_100/D1114T203138/gameN1300@D1114T203138_brisample_4.csv'
    #root_filename2='/Users/luishenryquiroganunez/Documents/Leiden_University/3_Semester/MSc_thesis/Codes/Code_v30/galaxy_generator/output/si30_moresrcs_gpara_100/D1114T203138/gameN1300@D1114T203138_brisample_5.csv'    
    with open(root_filename1) as csvfile:
        #Opening the file
        reader = csv.DictReader(csvfile)
    csvfile.close()
    with open(root_filename1) as csvfile:
        #Opening the file
        reader = csv.DictReader(csvfile)        
        #Declaring variables
        x_1=[]
        y_1=[]
        flux_1=[]
        para_1=[]
        err_para_1=[]
        for row in reader:
            x_1.append(row['x(kpc)'])                    
            y_1.append(row['y(kpc)'])      
            flux_1.append(row['f(Lsun)'])    
            para_1.append(row['parallax_obs(mas)'])              
            err_para_1.append(row['err_parallax(mas)'])                                            
    csvfile.close()  
    with open(root_filename2) as csvfile:
        #Opening the file
        reader = csv.DictReader(csvfile)
    csvfile.close()
    with open(root_filename2) as csvfile:
        #Opening the file
        reader = csv.DictReader(csvfile)        
        #Declaring variables
        x_2=[]
        y_2=[]
        flux_2=[]
        para_2=[]
        err_para_2=[]
        for row in reader:
            x_2.append(row['x(kpc)'])                    
            y_2.append(row['y(kpc)'])      
            flux_2.append(row['f(Lsun)'])    
            para_2.append(row['parallax_obs(mas)'])              
            err_para_2.append(row['err_parallax(mas)'])                    
    csvfile.close()  
    if(save=='yes'):
        #plt.ioff()
        pass
    else:
        print('The plots were NOT saved!')
    #plt.clf()  
    list_extra=[]
    xps = list(set(x_2) - set(x_1))
    for i in range(len(xps)):
        list_extra.append(x_2.index(xps[i]))

    yps=[]
    fps=[]
    paras=[]
    err_paras=[]
    for k in list_extra:
        yps.append(y_2[k])
        fps.append(flux_2[k])
        paras.append(para_2[k])
        err_paras.append(err_para_2[k])
        
    '''
    yps = list(set(y_2) - set(y_1))
    fps = list(set(flux_2) - set(flux_1))
    paras = list(set(para_2) - set(para_1))
    err_paras = list(set(err_para_2) - set(err_para_1))
    '''
    xp=[float(xps[i]) for i in range(len(xps))]
    yp=[float(yps[i]) for i in range(len(yps))]
    fp=[float(fps[i]) for i in range(len(fps))]
    para=[float(paras[i]) for i in range(len(paras))]
    err_para=[float(err_paras[i]) for i in range(len(err_paras))]
    sigma_perce_obs=[abs((x*100)/y) for x, y in zip(err_para, para)]
    #Plot xy with lumosities 
    matplotlib.rcParams.update({'font.size': 27}) 
    figure(figsize=(10,10))
    segs = par.gs.segs_spiral(par.rtrunc, par.nspirplot) # why different than par.nspirseg? OJO
    for j in range( len(segs) ):
        name_arm=['Norma Arm','Carina-Sagitarius Arm','Perseus Arm','Crux-Scutum Arm','Local Arm']
        line_arm_color=['r-','c-','g-','m-','b-']
        #for k in range(len(segs[0]['x'])):
        #    segs[j]['y'][k]=segs[j]['y'][k]*-1
        plot( segs[j]['x'], segs[j]['y'], line_arm_color[j],label=name_arm[j], linewidth=3)
        #plot(segs[j]['y'], segs[j]['x'],  line_arm_color[j],label=name_arm[j])
    norma = mlines.Line2D([], [], color='r')
    carina = mlines.Line2D([], [], color='c')
    perseus= mlines.Line2D([], [], color='g')
    crux = mlines.Line2D([], [], color='m')
    local = mlines.Line2D([], [], color='b')    
    text( 12.0,  12.0, 'Norma',     verticalalignment='bottom', horizontalalignment='right',     color='red', fontsize=25)
    text( -8.0, -14.0, 'Carina \n Sagitarius',     verticalalignment='bottom', horizontalalignment='right',     color='c', fontsize=25)
    text( -2.0,  11.0, 'Perseus',     verticalalignment='bottom', horizontalalignment='right',     color='green', fontsize=25)
    text( 14.0, -12.0, 'Crux-Scutum',     verticalalignment='bottom', horizontalalignment='right',     color='m', fontsize=25)
    text( -8.0, 5.5, 'Local',     verticalalignment='bottom', horizontalalignment='right',     color='b', fontsize=25)
    text( -0.8, -2.1, 'Central \n Ring', verticalalignment='bottom', horizontalalignment='center', color='k', fontsize=16)    
    kk=0
    for i in range(len(xp) ):
        plot( xp[i], yp[i], '.', alpha = 0.4 , color ='k', markersize = (numpy.log10(fp[i])*4)+37)
        z=i
        kk+=1
    print ('Number of sources in the sample:%s' % kk)
    small_dot, =plot( xp[z], yp[z], '.', alpha = 0.4 , color = 'k', markersize = 5)
    big_dot,   =plot( xp[z]+100, yp[z]+100, '.', alpha = 0.4 , color = 'k', markersize = 25)
    legend2=legend([small_dot, big_dot],['$10^{-8} L_{\odot}$','$10^{-3} L_{\odot}$'],loc=4,ncol=2,numpoints=1,prop={'size':23},borderpad=0.2)
    plt.gca().add_artist(legend2) 
    ring = par.cr.central_ring(par.nspirplot) 
    plot(ring[0]['x'], ring[0]['y'], 'k--', linewidth = 3)
    text(0.0,0.0, 'GC', fontsize=15)  
    plot(0.0,0.0,'o',color='g',markersize=10)
    xlabel('X $[kpc]$')
    ylabel('Y $[kpc]$')
    axis([-par.rtrunc-10.0,par.rtrunc+10.0,-par.rtrunc-10.0,par.rtrunc+10.0])
    axis('equal')
    xlim([-15,15])
    ylim([-15,15])
    plot( 0, par.r0, '*', color ='y',markersize=13)
    plt.grid(True)
    name_plot_0="samplot_3.pdf"
    if(save=='yes'):
        plt.savefig(name_plot_0,bbox_inches='tight')
        #print(name_plot_0)
    plt.close()    
    
    matplotlib.rcParams.update({'font.size': 22}) 
    fig5 = plt.figure(figsize=(8,8))
    ax = fig5.add_subplot(111)
    #ax1.set_xscale("log")    
    #ax.set_xlabel('Peak Flux Density (Jy)')
    ax.set_ylabel('Source Counts')
    nbins=100
    ##Percentage err/parallax_obs
    x,y,z=ax.hist(sigma_perce_obs,nbins,color ='b')
    ax.set_xlabel('sigma/parallax_obs (%)')    
    ax.set_xticks([100,500,1000,2000])
    #ax.set_xlim([0,114])                 
    if(save=='yes'):
        plt.savefig('sigpara_3.pdf',bbox_inches='tight')    
    plt.close()
    plt.ion()
    return    