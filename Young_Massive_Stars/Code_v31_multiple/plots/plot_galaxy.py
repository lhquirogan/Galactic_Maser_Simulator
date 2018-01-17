#Description: Plot of the galaxy over the model
#------------------------------------------------------------------------------------------------------------
#Importing Libraries
import pyfits
import numpy as np
import astropy
from astropy.io import ascii
from astropy.table import Table, Column
import scipy as scipy
import scipy.ndimage
import matplotlib.pyplot as plt
import matplotlib.patches as patches
#------------------------------------------------------------------------------------------------------------
#Reading galaxy
galaxy_grey=pyfits.getdata('/Users/luishenryquiroganunez/Documents/Research/PhDThesis/first_paper/galaxy_disk/milkywaystructure_8_4_kpc.fits')
#------------------------------------------------------------------------------------------------------------
#Run one galaxy
sample=gsample
#------------------------------------------------------------------------------------------------------------
#Loading simulation data
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
#------------------------------------------------------------------------------------------------------------
#Plot xy with lumosities 
size=23
size_2=11
size_3=20
matplotlib.rcParams.update({'font.size': size}) 
figure(figsize=(12,12))
implot = plt.imshow(galaxy_grey,cmap='Greys_r',origin='lower',extent=[-20,20,-20,20])
segs = par.gs.segs_spiral(par.rtrunc, par.nspirplot) # why different than par.nspirseg? OJO
for j in range( len(segs) ):
    name_arm=['Norma Arm','Carina-Sagitarius Arm','Perseus Arm','Crux-Scutum Arm','Local Arm']
    #line_arm_color=['r-','b-','g-','m-','c-']
    line_arm_color=['r-','c','g-','m-','b-']
    plot( segs[j]['x'], segs[j]['y'], line_arm_color[j],label=name_arm[j], linewidth=3)
norma = mlines.Line2D([], [], color='r')
carina = mlines.Line2D([], [], color='b')
perseus= mlines.Line2D([], [], color='g')
crux = mlines.Line2D([], [], color='k')
local = mlines.Line2D([], [], color='c')    
text( 14.0,  11.0, 'Outer \n Norma',     verticalalignment='bottom', horizontalalignment='right',     color='red', fontsize=size_3)
text( -8.2, -12.0, 'Carina \n Sagitarius',     verticalalignment='bottom', horizontalalignment='right',     color='c', fontsize=size_3)
text(-4.0,  12.0, 'Perseus',     verticalalignment='bottom', horizontalalignment='right',     color='green', fontsize=size_3)
text( 13.0, -12.0, 'Crux-Scutum',     verticalalignment='bottom', horizontalalignment='right',     color='m', fontsize=size_3)
text( -8.0, 4.8, 'Local',     verticalalignment='bottom', horizontalalignment='right',     color='b', fontsize=(size_3-2))
text( -0.8, -2.5, 'Central \n Ring', verticalalignment='bottom', horizontalalignment='center', color='k', fontsize=size_2)    
text( -12.5,  -15.5, '5 kpc',     verticalalignment='bottom', horizontalalignment='center',     color='white', fontsize=(size_3-3))
for i in range(len(xp) ):
    if(f_Jyp[i]>=0.003):
        if(dc[i]<=-20):
            markers='.'
            #plot( xp[i], yp[i], marker=markers, alpha = 0.9 , color = 'orange', markersize = (numpy.log10(fp[i])*4)+37)
            plot( xp[i], yp[i], marker=markers, alpha = 0.9 , color = par.colors[ic[i]], markersize = (numpy.log10(fp[i])*4)+37)
        if(dc[i]>=20):
            markers='.'
            plot( xp[i], yp[i], marker=markers, alpha = 0.9 , color = par.colors[ic[i]], markersize = (numpy.log10(fp[i])*4)+37)
            #plot( xp[i], yp[i], marker=markers, alpha = 0.9 , color = 'lightpink', markersize = (numpy.log10(fp[i])*4)+37)
        if(dc[i]<20 and dc[i]>0):
            num=np.random.normal(10,20)
            #num=np.random.uniform()
            if(num>=0.):
                #rand_color='pink'
                rand_color='orange'                
            elif(num<0):
                rand_color='orange'
            markers='.'
            plot( xp[i], yp[i], marker=markers, alpha = 0.9 , color = par.colors[ic[i]], markersize = (numpy.log10(fp[i])*4)+37)
            #plot( xp[i], yp[i], marker=markers, alpha = 0.9 , color = rand_color, markersize = (numpy.log10(fp[i])*4)+37)
        if(dc[i]>-20 and dc[i]<0):
            num=np.random.normal(-10,20)
            #num=np.random.uniform()
            if(num>0.):
                #rand_color='pink'
                rand_color='orange'
            elif(num<=0.):
                rand_color='orange'
            markers='.'            
            plot( xp[i], yp[i], marker=markers, alpha = 0.9 , color = par.colors[ic[i]], markersize = (numpy.log10(fp[i])*4)+37)
            #plot( xp[i], yp[i], marker=markers, alpha = 0.9 , color = rand_color, markersize = (numpy.log10(fp[i])*4)+37)
            '''
        if(dc[i]>-30 and dc[i]<40):
            markers=(5, 0)      
            plot( xp[i], yp[i], marker=markers, alpha = 0.4 , color = par.colors[ic[i]], markersize = (numpy.log10(fp[i])*4)+37)            
        if(dc[i]>=40):
            markers='.'
            plot( xp[i], yp[i], marker=markers, alpha = 0.4 , color = par.colors[ic[i]], markersize = (numpy.log10(fp[i])*4)+37)            
            '''
            z=i
#small_dot, =plot( xp[z], yp[z], '.', alpha = 0.4 , color = par.colors[ic[z]], markersize = 5)
#big_dot,   =plot( xp[z], yp[z], '.', alpha = 0.4 , color = par.colors[ic[z]], markersize = 25)
small_dot, =plot( xp[z]+100, yp[z]+100, '.' , alpha = 0.8 , color = 'k', markersize = 5)
big_dot,   =plot( xp[z]+100, yp[z]+100, '.' , alpha = 0.8 , color = 'k', markersize = 25)
dot_o,   =plot( xp[z]+100, yp[z]+100, '.' , alpha = 0.9 , color = 'orange', markersize = 18)
dot_p,   =plot( xp[z]+100, yp[z]+100, '.' , alpha = 0.9 , color = 'lightpink', markersize = 18)
#legend2=legend([small_dot, big_dot],['$\\rm{10^{-8} L_{\odot}}$','$\\rm{10^{-3} L_{\odot}}$'],ncol=2,numpoints=1,prop={'size':size_3},borderpad=0.2,bbox_to_anchor=(0.5, -0.12),loc='lower center',labelspacing=2)
legend2=legend([small_dot, big_dot],['$\\rm{10^{-8} L_{\odot}}$','$\\rm{10^{-3} L_{\odot}}$'],ncol=2,numpoints=1,prop={'size':size_3},borderpad=0.2,loc=4,labelspacing=2)
plt.gca().add_artist(legend2) 
#legend3=legend([dot_o,dot_p],['VLBI with $\\rm{SKA1-Mid}$','$\\rm{VLBA}$'],ncol=2,numpoints=1,prop={'size':size},borderpad=0.2,bbox_to_anchor=(0.5, 1.12),loc='upper center')
extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
line4 = Line2D(range(10), range(10), marker='', color="white")
#legend3=legend([extra, dot_o,line4, dot_p], ('                  VLBI with:', '$\\rm{SKA1-Mid}$','','$\\rm{VLBA}$'),ncol=2,numpoints=1,prop={'size':size},borderpad=0.2,bbox_to_anchor=(0.5, 1.18),loc='upper center')
#plt.gca().add_artist(legend3) 
ring = par.cr.central_ring(par.nspirplot) 
plot(ring[0]['x'], ring[0]['y'], 'k--', linewidth = 3)
#text(0.0,0.0, 'GC', fontsize=size_2)  
plot(0.0,0.0,'+',color='k',markersize=10)
#plt.xlabel('X (kpc)',labelpad=-5)
#plt.ylabel('Y (kpc)',labelpad=-27)
#tick_params(axis='both',labelsize=size)
#axis([-par.rtrunc-10.0,par.rtrunc+10.0,-par.rtrunc-10.0,par.rtrunc+10.0])
#axis('equal')
plt.xlim([-16,16])
plt.ylim([-16,16])
plot( 0, par.r0, '*', color ='y',markersize=15)
#plt.xticks([-10,-5,0,5,10])
#plt.yticks([-10,-5,0,5,10])
#plt.grid(True)
currentAxis = plt.gca()
currentAxis.add_patch(patches.Rectangle(
        (-14.9, -14),   # (x,y)
        5,          # width
        0.1,          # height 
        facecolor='white'))
plt.axis('off')
#name_plot_0="../plots/samplot_calendar3.png"
#name_plot_1="../plots/samplot.pdf"
#plt.savefig(name_plot_0,bbox_extra_artists=(legend3,), bbox_inches='tight')
#plt.savefig(name_plot_1,bbox_extra_artists=(lgd,), bbox_inches='tight'))
plt.savefig(name_plot_1, bbox_inches='tight')
#------------------------------------------------------------------------------------------------------------
#Importing Libraries
import pyfits
import numpy as np
import astropy
from astropy.io import ascii
from astropy.table import Table, Column
import scipy as scipy
import scipy.ndimage
import matplotlib.pyplot as plt
#------------------------------------------------------------------------------------------------------------
#Reading l-v
galaxy_lv_cube=pyfits.getdata('/Users/luishenryquiroganunez/Documents/Leiden_University/3_Semester/MSc_thesis/Codes/Code_v27/surveys/COGAL_deep.fits')
where_are_NaNs = np.isnan(galaxy_lv_cube)
galaxy_lv_cube[where_are_NaNs] = 0
galaxy_lv=np.sum(galaxy_lv_cube, axis=0)
galaxy_lv=np.rot90(np.rot90(np.rot90(galaxy_lv)))
#galaxy_lv=np.rot90(galaxy_lv)
#pyfits.writeto('/Users/luishenryquiroganunez/Documents/Leiden_University/3_Semester/MSc_thesis/Codes/Code_v27/surveys/cube_lv.fits',galaxy_lv)
#------------------------------------------------------------------------------------------------------------
#Plot xy with lumosities 
matplotlib.rcParams.update({'font.size': 17})
vlsr_plot=figure(figsize=(12,11))
vlsr_ax=vlsr_plot.add_subplot(111)
min_dat=numpy.log10(min(f_Jyp))
max_dat=numpy.log10(max(f_Jyp))  
par.colors[3]='m'
#title("Galactic Methanol Masers Distribution")
for i in range( len(lp) ):
    vlsr_ax.plot( lp[i], vl[i], 'o', alpha = 0.8 , color = par.colors[ic[i]], markersize = ((numpy.log10(fp[i])*3.0)+27.0))
    z=i
#implot = plt.imshow(galaxy_lv,cmap='gray',origin='lower',extent=[180,-180,-150,150],aspect=0.2)#,vmin=-1, vmax=1)
implot = plt.imshow(galaxy_lv,cmap='gray',origin='lower',extent=[180,-180,-169,169],aspect=0.2,vmin=-5, vmax=12)
small_dot, =vlsr_ax.plot( lp[z]+300, vl[z], 'o', alpha = 0.8 , color = 'k', markersize = 3)
big_dot,   =vlsr_ax.plot( lp[z]+300, vl[z], 'o', alpha = 0.8 , color = 'k', markersize = 17.5)
vlsr_ax.legend([small_dot, big_dot],['$\\rm{10^{-8} L_{\odot}}$','$\\rm{10^{-3} L_{\odot}}$'],loc=2,numpoints=1,prop={'size':19},ncol=2)
vlsr_ax.axis([135,-135,-180,180])
vlsr_ax.set_xlabel('Galactic Longitude (degrees)')
#vlsr_ax.set_xticks([-90,-60,-30,0,30,60,90])
vlsr_ax.set_xticks([-135,-90,-45,0,45,90,135])
vlsr_ax.set_ylabel('$V_{LSR}$ (km/s)',labelpad=-5)
vlsr_ax.grid(True)
tick_params(axis='both',labelsize=15)
vlsr_ax.set_axis_bgcolor('grey')
#name_plot_1="../plots/l-v3.pdf"
#plt.savefig(name_plot_1,bbox_inches='tight')



#------------------------------------------------------------------------------------------------------------
#Plot error uncertianties

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
        kk=content.index(matching[-1])
        #Error in parallax
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

'''
filename='/Users/luishenryquiroganunez/Documents/Leiden_University/3_Semester/MSc_thesis/Codes/Code_v31/galaxy_generator/output/si31_BeSSeLmimic/D0216T235052/gameN1300@D0216T235052_brisample_0.csv'
with open(filename) as csvfile:
    #Opening the file
    reader = csv.DictReader(csvfile)
    row_count = sum(1 for row in reader)
with open(filename) as csvfile:
    #Opening the file
    reader = csv.DictReader(csvfile)        
    #Declaring variables
    err_pi=numpy.zeros(row_count)        
    i=0
    for row in reader:
        err_pi[i]=row['err_parallax(mas)']  
        i+=1
'''        
                             
para=[]
err_parallax=[]
for i in range(len(gsample)):
    para.append(gsample[i]['parallax_obs(mas)'])
    err_parallax.append(gsample[i]['err_parallax(mas)'])

matplotlib.rcParams.update({'font.size': 27})
fig1 = plt.figure()
ax = fig1.add_subplot(111)
nbins=np.arange(0.0,0.25,0.035)
#nbins=10
weights = np.ones_like(para_unc)/len(para_unc)
a,b,c=ax.hist(para_unc,bins=nbins,color='b',alpha=0.4,label='Observational Errors',weights=weights)
weights = np.ones_like(err_parallax)/len(err_parallax)
ax.hist(err_parallax,bins=b,histtype='step', edgecolor='green',linewidth=3.5,label='Simulated Errors',weights=weights)
ax.grid(True)
ax.legend(loc=0,ncol=1)
plt.xlabel('Parallax error $\\Delta \\pi$ (mas)')
plt.ylabel('Normalized Counts')
plt.savefig('/Users/luishenryquiroganunez/Documents/Leiden_University/3_Semester/MSc_thesis/Codes/Code_v31/plots/e_parallax.pdf',bbox_inches='tight')  