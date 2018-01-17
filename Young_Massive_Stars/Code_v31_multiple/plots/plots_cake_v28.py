import numpy as np
import matplotlib.pyplot as plt
import ephem
import random
import math
from math import pi
from operator import truediv
import pylab
from astropy.io import ascii
from matplotlib.patches import Ellipse
#pi=3.141592653589793


def Trans(alpha,delta):
    RA=alpha*pi/180
    DEC=delta*pi/180    
    tmp=ephem.Equatorial(RA,DEC,epoch='2000')
    lb=ephem.Galactic(tmp)
    b=float(lb.lat)*180/pi
    l=float(lb.long)*180/pi
    aa=[l,b]
    return(aa)  
        
def Trans_galaeq(longi,lati):
    l=longi*pi/180
    b=lati*pi/180    
    tmp=ephem.Galactic(l,b,epoch='2000')
    eq=ephem.Equatorial(tmp)
    alpha=float(eq.ra)*180/pi
    delta=float(eq.dec)*180/pi
    aa=[alpha,delta]
    return(aa)
    
def points_generator_gala(num,longi_max,longi_min):
    point=[]
    for i in range(num):
        ran_long=random.uniform(longi_min,longi_max)
        ran_lat=random.uniform(-2.0,2.0)
        z=Trans_galaeq(ran_long,ran_lat)
        ran_dist=random.uniform(0,20)
        k=[ran_long,ran_lat,z[0],z[1],ran_dist]
        point.append(k)
    return(point)    

def points_generator(num,dec_max,dec_min,ra_max,ra_min):
    point=[]
    for i in range(num):
        ran_dec=random.uniform(dec_min,dec_max)
        ran_ra=random.uniform(ra_min,ra_max)
        z=Trans(ran_ra,ran_dec)
        ran_dist=random.uniform(0,20)
        k=[ran_ra,ran_dec,z[0],z[1],ran_dist]
        point.append(k)
    return(point)


#-----------------------------------------------------------------------------------------------------------------------------------------------
#plot 1      
N = 8
thetas = np.arange(0.0, 2*np.pi, 2*np.pi/(N))
radiis = np.ones(N)*(15+par.r0)
widths = np.ones(N)*np.pi/4
colork=['g','k','r','b','y','m','c','orange']
gc=(-1*par.r0)
gc_angle=atan2(gc,0)+pi/2
    
ax = plt.subplot(111, projection='polar')
#ax.scatter(l_1, r_1, color='r')
ax.plot(0,0,'*',color='y',markersize=20)
ax.plot(gc_angle,par.r0,'o',color='g',markersize=10)
segs = par.gs.segs_spiral(par.rtrunc, par.nspirplot)
for j in range( len(segs) ):
    x_array = np.asarray(segs[j]['x'])
    y_array = np.asarray(segs[j]['y'])-par.r0
    r=np.sqrt(x_array**2+y_array**2)
    theta=np.zeros(len(r))
    for i in range(len(r)):
        theta[i]=atan2(y_array[i],x_array[i])+pi/2
    #ax.plot(theta,r, line_arm_color[j],label=name_arm[j], linewidth=3)
    ax.plot(theta,r, 'k', linewidth=3)
ax.set_theta_zero_location("S")
ax.grid(True)
for k in range(N):
    ax.bar(thetas[k], radiis[k], width=widths[k], bottom=0.0,color=colork[k],edgecolor='k',alpha=0.3)
ax.set_rgrids([5,10,15,20],angle=202.5,fontsize=18,label='kpc')
ax.set_xticklabels(['$l=0^{o}$', '$l=45^{o}$', '$l=90^{o}$', '$l=135^{o}$', '$l=180^{o}$', '$l=225^{o}$', '$l=270^{o}$', '$l=315^{o}$'])
ax.set_yticklabels(['5 $kpc$','10 $kpc$','15 $kpc$','20 $kpc$'])
ax.text(gc_angle,par.r0, 'GC', fontsize=15)
for n in range(N):
    dummy_angle=(n*pi/(4))+ 0.39269908169872414
    ax.text(dummy_angle,22.5, str(n+1), fontsize=15,color=colork[n])
plt.show()
#plt.savefig('long_cake.pdf',bbox_inches='tight')
#-----------------------------------------------------------------------------------------------------------------------------------------------

points_1=points_generator_gala(10000,90,-30,0,24)
l_1=[]
b_1=[]
ra_1=[]
dec_1=[]
r_1=[]
for i in range(len(points_1)):
    l_1.append(points_1[i][0]*pi/180)
    b_1.append(points_1[i][1]*pi/180)
    ra_1.append(points_1[i][2]*pi/180)  
    dec_1.append(points_1[i][3]*pi/180)
    r_1.append(points_1[i][4])      
for i in range(len(ra_1)):
    if(ra_1[i]>180):
        ra_1[i]=ra_1[i]-360
points=points_generator(num,dc_max,dc_min,ra_max,ra_min)
points=points_generator(10000,90,-30,315,89)
l=[]
r=[]
b=[]
for i in range(len(points)):
   if(points[i][1]>-2 and points[i][1]< 2):
        b.append(points[i][1])
        l.append(points[i][0])
        r.append(points[i][2])

#------------------------------------------------------------------------------------------------------------------------------    
splits=8
dic={}
for j in range(splits):
    dic['ra_'+str(j)]=[]
    dic['dec_'+str(j)]=[]
    dic['r_'+str(j)]=[]  
    dic['l_'+str(j)]=[]
    dic['b_'+str(j)]=[]        
    mini=j*360/splits
    maxi=(j+1)*360/splits    
    points=points_generator_gala(500,maxi,mini)
    for i in range(len(points)):
        dic['l_'+str(j)].append(points[i][0]*pi/180)
        dic['b_'+str(j)].append(points[i][1]*pi/180)        
        dic['ra_'+str(j)].append(points[i][2]/15)
        dic['dec_'+str(j)].append(points[i][3])
        dic['r_'+str(j)].append(points[i][4])
#------------------------------------------------------------------------------------------------------------------------------
params = {'legend.fontsize': 20}
pylab.rcParams.update(params)
az=plt.subplot(111)
for i in range(len(dic)/5):
    az=plt.subplot(111)    

plt.grid(True)
plt.legend(scatterpoints=1,ncol=8,loc=8)
xlabel('Right Ascension ($hh$)')
ylabel('Declination ($^o$)')
xlim([0,24])
ylim([-80,80])
#plt.savefig('radec_cake.pdf',bbox_inches='tight')        

lp_new=np.zeros(len(lp))
for i in range(len(lp)):
    if (lp[i]<0):
        lp_new[i]=lp[i]+360
    else:
        lp_new[i]=lp[i]

hola=figure()
az=hola.add_subplot(111)
az.scatter(lp_new, bp)
xlabel('longitude')
ylabel('latititde')

xlim([0,24])
ylim([-80,80])

#------------------------------------------------------------------------------------------------------------------------------
colork2=['g','k','r','b','y','m','c','o']
rac=[]
dec=[]
for i in range(len(gsample)):
    dec.append(gsample[i]['dc(rad)']*180/pi)
    rac.append((gsample[i]['ra(rad)']*180/pi)*24/360)
dec_range={}
for p in range(len(colork2)):
    dec_range['range_'+str(p)]=[]
for p in range(len(colork2)):
    min_dummy_dec=min(dic['dec_'+str(p)])
    max_dummy_dec=max(dic['dec_'+str(p)])
    min_dummy_rac=min(dic['ra_'+str(p)])
    max_dummy_rac=max(dic['ra_'+str(p)])
    for j in range(len(dec)):
        if(dec[j]>min_dummy_dec and dec[j]<max_dummy_dec and rac[j]>min_dummy_rac and rac[j]<max_dummy_rac):
            dec_range['range_'+str(p)].append(dec[j])




print('---------------------------------------------------------------')
print('------------RANGES---------------------------------------------')
print('---------------------------------------------------------------')
print('Region|| Color ||    RA   ||     DEC   ')
for p in range(len(colork2)):
    print('%s %s %s %s' % (str(p+1)+'     ||',colork2[p]+'     ||','['+str(int(min(dic['ra_'+str(p)]))).zfill(2)+','+str(int(max(dic['ra_'+str(p)]))).zfill(2)+'] ||','['+str(int(min(dic['dec_'+str(p)]))).zfill(2)+','+str(int(max(dic['dec_'+str(p)]))).zfill(2)+']'+'     ||'))
#------------------------------------------------------------------------------------------------------------------------------
ax = plt.subplot(111, projection='polar')
for m in range(len(colork)):
    ax.scatter(dic['l_'+str(m)],dic['r_'+str(m)], color=colork[m],alpha=0.3)
ax.plot(0,0,'*',color='y',markersize=20)
ax.plot(gc_angle,par.r0,'o',color='g',markersize=10)
segs = par.gs.segs_spiral(par.rtrunc, par.nspirplot)
for j in range( len(segs) ):
    x_array = np.asarray(segs[j]['x'])
    y_array = np.asarray(segs[j]['y'])-par.r0
    r=np.sqrt(x_array**2+y_array**2)
    theta=np.zeros(len(r))
    for i in range(len(r)):
        theta[i]=atan2(y_array[i],x_array[i])+pi/2
    #ax.plot(theta,r, line_arm_color[j],label=name_arm[j], linewidth=3)
    ax.plot(theta,r, 'k', linewidth=3)
ax.set_theta_zero_location("S")
ax.grid(True)
for k in range(N):
    ax.bar(thetas[k], radiis[k], width=widths[k], bottom=0.0,color=colork[k],edgecolor='k',alpha=0.3)
ax.set_rgrids([5,10,15,20],angle=202.5,fontsize=18,label='kpc')
ax.set_xticklabels(['$l=0^{o}$', '$l=45^{o}$', '$l=90^{o}$', '$l=135^{o}$', '$l=180^{o}$', '$l=225^{o}$', '$l=270^{o}$', '$l=315^{o}$'])
ax.set_yticklabels(['5 $kpc$','10 $kpc$','15 $kpc$','20 $kpc$'])
ax.text(gc_angle,par.r0, 'GC', fontsize=15)
plt.show()
#plt.savefig('long_cake+scatter.pdf',bbox_inches='tight')
#------------------------------------------------------------------------------------------------------------------------------        
figure()
matplotlib.rcParams.update({'font.size': 17})
map = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')
title("Right Ascention and Declination view")
k=['r','b','g','k','m','y']
for j in range(len(dic)/5):
    tx=np.zeros(len(dic['ra_'+str(j)]))
    ty=np.zeros(len(dic['ra_'+str(j)]))
    for i in range(len(dic['ra_'+str(j)])):
        tx[i], ty[i] = map(dic['ra_'+str(j)][i]*15,dic['dec_'+str(j)][i])
        map.plot( tx[i], ty[i], 'o', alpha = 0.5,color=colork[j])
#map.drawparallels(numpy.arange(-90.,120.,30.),labels=[True])
#map.drawmeridians(np.arange(0.,420.,60.))
map.drawmapboundary(fill_color='GhostWhite')
map.drawparallels(numpy.arange(-90.,120.,30.),labels=[True])
map.drawmeridians(np.arange(0.,420.,60.))
#plt.savefig('dec_sky.pdf',bbox_inches='tight')

#------------------------------------------------------------------------------------------------------------------------------        
#%cd "/Users/luishenryquiroganunez/Documents/Leiden_University/3_Semester/MSc_thesis/Codes/Reid_code/fit_galaxy_bayesian-2/"
t = ascii.read('source_files.inp')
alphas_reid=[]
deltas_reid=[]
distance_reid=[]
l_reid=[]
b_reid=[]
for j in range(len(t)):
    f = open(t[j][0], 'r')
    m=f.read()
    l=m.split()
    f.close()
    for i in range(len(l)):
        if(l[i]=='RA'):
            a=i    
    for i in range(len(l)):
        if(l[i]=='Dec'):
            b=i
    RA=l[a-2]
    point_ra=RA.find('.')
    deci_ss_ra=RA[point_ra:]
    int_ss_ra=RA[point_ra-2:point_ra]
    ss_ra=float(int_ss_ra+deci_ss_ra)
    k_ra=RA.find(str(int_ss_ra+deci_ss_ra))
    mm_ra=float(RA[k_ra-2:k_ra])
    r_ra=RA.find(RA[k_ra-2:k_ra]+int_ss_ra+deci_ss_ra)
    hh_ra=float(RA[0:r_ra])
    alpha=((hh_ra)+(mm_ra/60.0)+(ss_ra/3600.0))*15
    
    DEC=l[b-2]
    point=DEC.find('.')
    deci_ss=DEC[point:]
    int_ss=DEC[point-2:point]
    ss=float(int_ss+deci_ss)
    k=DEC.find(str(int_ss+deci_ss))
    mm=float(DEC[k-2:k])
    r=DEC.find(DEC[k-2:k]+int_ss+deci_ss)
    dd=float(DEC[0:r])
    if (dd<0):
        delta=(-1)*(abs(dd)+(mm/60.0)+(ss/3600.0))
    else:
        delta=dd+(mm/60.0)+(ss/3600.0)
    distance_reid.append(1/(float(l[b+2])))
    tmp=ephem.Equatorial(alpha*pi/180,delta*pi/180,epoch='2000')
    lb=ephem.Galactic(tmp)
    b_reid.append(float(lb.lat))
    l_reid.append(float(lb.long))    
    alphas_reid.append(alpha*pi/180)
    deltas_reid.append(delta*pi/180)
#--------------------------------------------------------------------------------------------------------------------    
#plot      reid
N = 8
thetas = np.arange(0.0, 2*np.pi, 2*np.pi/(N))
radiis = np.ones(N)*(15+par.r0)
widths = np.ones(N)*np.pi/4
colork=['g','k','r','b','y','m','c','orange']
gc=(-1*par.r0)
gc_angle=atan2(gc,0)+pi/2
    
ax = plt.subplot(111, projection='polar')
ax.scatter(l_reid, distance_reid, color='k')#,alpha=0.8,label='Reid et al. 2014')
ax.plot(0,0,'*',color='y',markersize=20)
ax.plot(gc_angle,par.r0,'o',color='g',markersize=10)
segs = par.gs.segs_spiral(par.rtrunc, par.nspirplot)
for j in range( len(segs) ):
    x_array = np.asarray(segs[j]['x'])
    y_array = np.asarray(segs[j]['y'])-par.r0
    r=np.sqrt(x_array**2+y_array**2)
    theta=np.zeros(len(r))
    for i in range(len(r)):
        theta[i]=atan2(y_array[i],x_array[i])+pi/2
    #ax.plot(theta,r, line_arm_color[j],label=name_arm[j], linewidth=3)
    ax.plot(theta,r, 'k', linewidth=3)
ax.set_theta_zero_location("S")
ax.grid(True)
for k in range(N):
    ax.bar(thetas[k], radiis[k], width=widths[k], bottom=0.0,color=colork[k],edgecolor='k',alpha=0.3)
ax.set_rgrids([5,10,15,20],angle=202.5,fontsize=18,label='kpc')
ax.set_xticklabels(['$l=0^{o}$', '$l=45^{o}$', '$l=90^{o}$', '$l=135^{o}$', '$l=180^{o}$', '$l=225^{o}$', '$l=270^{o}$', '$l=315^{o}$'])
ax.set_yticklabels(['5 $kpc$','10 $kpc$','15 $kpc$','20 $kpc$'])
ax.text(gc_angle,par.r0, 'GC', fontsize=15)
for n in range(N):
    dummy_angle=(n*pi/(4))+ 0.39269908169872414
    ax.text(dummy_angle,22.5, str(n+1), fontsize=15,color=colork[n])
plt.show()
#plt.legend(scatterpoints=1, loc=8)
#plt.savefig('long_cake_reid.pdf',bbox_inches='tight')
#--------------------------------------------------------------------------------------------------------------------
#Reid counting
reids={}
for y in range(N):
   reids['region_'+str(y+1)]=[]        
for i in range(len(l_reid)):
    if(l_reid[i]>=0 and l_reid[i]<=pi/4):
        reids['region_1'].append(l_reid[i])
    if(l_reid[i]>=pi/4 and l_reid[i]<=pi/2):
        reids['region_2'].append(l_reid[i])
    if(l_reid[i]>=pi/2 and l_reid[i]<=3*pi/4):
        reids['region_3'].append(l_reid[i])
    if(l_reid[i]>=3*pi/4 and l_reid[i]<=pi):
        reids['region_4'].append(l_reid[i])        
    if(l_reid[i]>=pi and l_reid[i]<=5*pi/4):
        reids['region_5'].append(l_reid[i])
    if(l_reid[i]>=5*pi/4 and l_reid[i]<=3*pi/2):
        reids['region_6'].append(l_reid[i])
    if(l_reid[i]>=3*pi/2 and l_reid[i]<=7*pi/4):
        reids['region_7'].append(l_reid[i])
    if(l_reid[i]>=7*pi/4 and l_reid[i]<=2*pi):
        reids['region_8'].append(l_reid[i])  
        
        
colork2=['g','k','r','b','y','m','c','o']
print('-----------------------------')
print('-------FILTER   REID-----------')
print('-----------------------------')
print('Region|| Color ||    Number   ')
for p in range(len(colork)):
    print('%s %s %s' % (str(p+1)+'     ||',colork2[p]+'     ||', len(reids['region_'+str(p+1)])))
#------------------------------------------------------------------------------------------------------------------------------        
figure()
matplotlib.rcParams.update({'font.size': 17})
map = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')
title("Right Ascention and Declination view")
k=['r','b','g','brown','m','y']
for j in range(len(dic)/5):
    tx=np.zeros(len(dic['ra_'+str(j)]))
    ty=np.zeros(len(dic['ra_'+str(j)]))
    for i in range(len(dic['ra_'+str(j)])):
        tx[i], ty[i] = map(dic['ra_'+str(j)][i]*15,dic['dec_'+str(j)][i])
        map.plot( tx[i], ty[i], 'o', alpha = 0.5,color=colork[j])
for r in range(len(alphas_reid)):
    tx=np.zeros(len(alphas_reid))
    ty=np.zeros(len(alphas_reid))
    tx[r], ty[r] = map(alphas_reid[r]*180/pi,deltas_reid[r]*180/pi)
    map.plot( tx[r], ty[r], '*',color='c',markersize=10)        
#map.drawparallels(numpy.arange(-90.,120.,30.),labels=[True])
#map.drawmeridians(np.arange(0.,420.,60.))
map.drawmapboundary(fill_color='GhostWhite')
map.drawparallels(numpy.arange(-90.,120.,30.),labels=[True])
map.drawmeridians(np.arange(0.,420.,60.))
#plt.savefig('dec_sky_reid.pdf',bbox_inches='tight')
#------------------------------------------------------------------------------------------------------------------------------        
#Number of sources in Reid Region
deltas_degrees_reid=[]
for i in range(len(deltas_reid)):
    deltas_degrees_reid.append(deltas_reid[i]*180/pi)
total_simu_lives_Reid=[]
for i in range(len(gsample)):
    if((gsample[i]['dc(rad)']*180/pi)<max(deltas_degrees_reid) and (gsample[i]['dc(rad)']*180/pi)>min(deltas_degrees_reid)):
        total_simu_lives_Reid.append(gsample[i]['dc(rad)']*180/pi)
total_simu_lives_Reid=[]
for i in range(len(gsample)):
    if((gsample[i]['dc(rad)']*180/pi)<70 and (gsample[i]['dc(rad)']*180/pi)>-30):
        total_simu_lives_Reid.append(gsample[i]['dc(rad)']*180/pi)    
    
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord


ra_Gaia = data_Gaia[:,1]  #my data in deg, icrs frame
dec_Gaia = data_Gaia[:,2]

ra_Gaia = np.radians(ra_Gaia) -np.pi  #convert in rad to make a plot
dec_Gaia = np.radians(dec_Gaia) 

# definition of galactic plane as point with longitude from 0 to 360 and latitude 0 :
lon = np.linspace(0,360,num=500)
lat = np.zeros(500)
c = SkyCoord(l=lon, b=lat, unit=(u.degree, u.degree), frame='galactic') #convert in icrs frame and in deg to make a plot
c = c.icrs
ra_gal = c.ra.rad 
dec_gal = c.dec.rad 
ra_gal = ra_gal-np.pi

#tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
plt.subplot(111, projection='aitoff', axisbg ='SeaShell') #last parameter (axisbg ='SeaShell') changes color of background
plt.plot(ra_gal, dec_gal, 'ro', markersize=3, alpha=0.7)  #plot galactic plane
plt.plot(ra_Gaia, dec_Gaia, 'bs', markersize=3, alpha=0.7) #plot observations 
plt.grid(True)  #this is what you need to make a grid
#plt.tick_params(ticks=tick_labels, labelsize=10)
#plt.set_xticklabels([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
#plt.xticklabels(tick_labels)
plt.xlabel("RA",fontsize=18)
plt.ylabel("Dec",fontsize=18)
plt.show()


#-------------------

#Values_Fudge=[1.0,0.5,0.25,0.10,0.05,0.025,0.01,0.005,0.0]
para=[]
err_parallax=[]
for i in range(len(gsample)):
    para.append(gsample[i]['parallax_obs(mas)'])
    err_parallax.append(gsample[i]['err_parallax(mas)'])
'''     
map_para={}
Values_Fudge=[1.0,0.1,0.01,0.001,0.0001]
for z in range(len(Values_Fudge)):
    para_obs=[]
    for i in range(len(gsample)):
        unc_dummy=err_parallax[i]*Values_Fudge[z]
        para_obs.append(gauss(para[i],unc_dummy))
    map_para['fudge_'+str(Values_Fudge[z])]=para_obs
'''

#run mine_plot errors to ad the errors in Reid et al. 2014 for para_unc
map_unc={}
Values_Fudge=[5.0,1.0,0.5]
for z in range(len(Values_Fudge)):
    unc_f=[]
    for i in range(len(gsample)):
        unc_dummy=err_parallax[i]*Values_Fudge[z]
        unc_f.append(unc_dummy)
    map_unc['fudge_'+str(Values_Fudge[z])]=unc_f
  
    
      
          
colork3=['g','r','b','c','m','orange','purple','peru','gold']
matplotlib.rcParams.update({'font.size': 22})
ax = plt.subplot(111)
nbins=10
a,b,c=ax.hist(para_unc,alpha=0.3,bins=nbins,color='k',label='$\\Delta\\pi$ Reid data',normed=True)
ax.axvline(np.median(para_unc), color='k', ls='dashed',linewidth=2)
for k in range(len(Values_Fudge)):
    labels='$f$ = %s' %(Values_Fudge[k])
    ax.hist(map_unc['fudge_'+str(Values_Fudge[k])],bins=b,alpha=0.3,color=colork3[k],label=labels,normed=True)
    ax.axvline(np.median(map_unc['fudge_'+str(Values_Fudge[k])]), color=colork3[k], ls='dashed',linewidth=2)
ax.grid(True)
ax.legend(loc=0,ncol=1)
plt.xlabel('$ \\Delta \\pi$ [mas]')
plt.ylabel('Normed Counts')
plt.savefig('fudge_parallax.pdf',bbox_inches='tight')  


err_mux=[]
for i in range(len(gsample)): err_mux.append(gsample[i]['err_mu_x(mas/yr)'])

ax = plt.subplot(111)
a,b,c=ax.hist(err_mux,alpha=0.5,color='b',label='Unc. $\\mu_x$ Simulated = $\\sigma$',normed=True)
ax.hist(mux_unc,alpha=0.5,bins=b,color='g',label='Unc. in $\\mu_x$ Reid data',normed=True)
plt.grid(True)
plt.legend(loc=0)
xlabel('Unc. in $\\mu_x$ [mas/yr]')
ylabel('Normed Counts')
plt.savefig('unc_mux_comparison.pdf',bbox_inches='tight')  


#Jimis data plot
#%cd "/Users/luishenryquiroganunez/Documents/Leiden_University/3_Semester/MSc_thesis/Codes/Reid_code/fit_galaxy_bayesian-2/"
t_Jy = ascii.read('../surveys/MMB_Dist_Flux_edited.csv')    
t_Jy.remove_row(0)
kk=[]
for i in range(len(t_Jy)):
    if(t_Jy[i][5]=='-' or t_Jy[i][0]=='-' or t_Jy[i][1]=='-' or t_Jy[i][2]=='-'):
        kk.append(i)
    
t_Jy.remove_rows(kk)

distance_jimmi=[]
l_jimmi=[]
b_jimmi=[]
for i in range(len(t_Jy)):            
        distance_jimmi.append(float(t_Jy[i][2]))
        l_jimmi.append(float(t_Jy[i][0])*pi/180)
        b_jimmi.append(float(t_Jy[i][1])*pi/180)                

#--------------------------------------------------------------------------------------------------------------------    
#plot jimmi
#plt.figure()
#ax = plt.gca()
N = 8
#In polar coordiantes
r0 = 8.34; th0 = 0#-pi/2.0  # center
a = 4.5; b = 2.2  # semidiameters
th = np.linspace(0, 2*pi, 1000)  # x values of interest
phi=-52.0*pi/180.0
width=10.0
'''
P=numpy.zeros(len(th))
Q=numpy.zeros(len(th))
R=numpy.zeros(len(th))
rr=numpy.zeros(len(th))
for i in range(len(th)):
    P[i]=r0*(((b**2-a**2)*cos(th[i]+th0-2*phi))+((a**2+b**2)*cos(th[i]-th0)))
    R[i]=(b**2-a**2)*cos(2*th[i]-2*phi)+a**2+b**2
for i in range(len(th)):
        Q[i]=numpy.sqrt(2)*a*b*numpy.sqrt(R[i]-(2*r0**2*(sin(th[i]-th0)**2)))
for i in range(len(th)):
        rr[i]=(P[i]+Q[i])/R[i]
'''      
X=a*numpy.cos(th)
Y=b*numpy.sin(th)
x=r0*numpy.cos(th0)+X*numpy.cos(phi)-Y*numpy.sin(phi)
y=r0*numpy.sin(th0)+X*numpy.sin(phi)+Y*numpy.cos(phi)
radii=numpy.sqrt(x**2+y**2)
thetas=numpy.zeros(len(y))
for i in range(len(y)):
    thetas[i]=math.atan2(y[i],x[i])
            
ax = plt.subplot(111, projection='polar')
ax.scatter(thetas,radii,color='r',s=width)
thetas = np.arange(0.0, 2*np.pi, 2*np.pi/(N))
radiis = np.ones(N)*(15+par.r0)
widths = np.ones(N)*np.pi/4
colork=['g','k','r','b','y','m','c','orange']
gc=(-1*par.r0)
gc_angle=atan2(gc,0)+pi/2    
ax = plt.subplot(111, projection='polar')
ax.scatter(l_jimmi, distance_jimmi, color='b',alpha=0.5)#,label='Reid et al. 2014')
ax.plot(0,0,'*',color='y',markersize=20)
ax.plot(gc_angle,par.r0,'o',color='g',markersize=10)
segs = par.gs.segs_spiral(par.rtrunc, par.nspirplot)
for j in range( len(segs) ):
    x_array = np.asarray(segs[j]['x'])
    y_array = np.asarray(segs[j]['y'])-par.r0
    r=np.sqrt(x_array**2+y_array**2)
    theta=np.zeros(len(r))
    for i in range(len(r)):
        theta[i]=atan2(y_array[i],x_array[i])+pi/2
    #ax.plot(theta,r, line_arm_color[j],label=name_arm[j], linewidth=3)
    ax.plot(theta,r, 'k', linewidth=3)
ax.set_theta_zero_location("S")
ax.grid(True)
for k in range(N):
    ax.bar(thetas[k], radiis[k], width=widths[k], bottom=0.0,color=colork[k],edgecolor='k',alpha=0.3)
ax.set_rgrids([5,10,15,20],angle=202.5,fontsize=18,label='kpc')
ax.set_xticklabels(['$l=0^{o}$', '$l=45^{o}$', '$l=90^{o}$', '$l=135^{o}$', '$l=180^{o}$', '$l=225^{o}$', '$l=270^{o}$', '$l=315^{o}$'])
ax.set_yticklabels(['5 $kpc$','10 $kpc$','15 $kpc$','20 $kpc$'])
ax.text(gc_angle,par.r0, 'GC', fontsize=15)
for n in range(N):
    dummy_angle=(n*pi/(4))+ 0.39269908169872414
    ax.text(dummy_angle,22.5, str(n+1), fontsize=15,color=colork[n])
#plt.show()
#plt.legend(scatterpoints=1, loc=8)
#plt.savefig('../plots/long_cake_jimmi.pdf',bbox_inches='tight')

#--------------------------------------------------------------------------------------------------------------------    
#plot jimmi 2
#Jimis data plot
t_Jy = ascii.read('../surveys/MMB_Dist_Flux_edited.csv')    
t_Jy.remove_row(0)
kk=[]
for i in range(len(t_Jy)):
    if(t_Jy[i][5]=='-' or t_Jy[i][0]=='-' or t_Jy[i][1]=='-' or t_Jy[i][2]=='-'):
        kk.append(i)
    
t_Jy.remove_rows(kk)

distance_jimmi=[]
l_jimmi=[]
b_jimmi=[]
for i in range(len(t_Jy)):            
        distance_jimmi.append(float(t_Jy[i][2]))
        l_jimmi.append(float(t_Jy[i][0])*pi/180)
        b_jimmi.append(float(t_Jy[i][1])*pi/180)    

x_jimmi=np.zeros(len(distance_jimmi))
y_jimmi=np.zeros(len(distance_jimmi))
for i in range(len(distance_jimmi)):
    x_jimmi[i]=cos((3*pi/2)+l_jimmi[i])*distance_jimmi[i]
    y_jimmi[i]=sin((3*pi/2)+l_jimmi[i])*distance_jimmi[i]+par.r0
    
figure(figsize=(10,10))
#matplotlib.rcParams.update({'font.size': 27}) 
segs = par.gs.segs_spiral(par.rtrunc, par.nspirplot) # why different than par.nspirseg? OJO
for j in range( len(segs) ):
    line_arm_color=['r-','b-','g-','k-','c-']
    plot( segs[j]['x'], segs[j]['y'], line_arm_color[j], linewidth=3)
norma = mlines.Line2D([], [], color='r')
carina = mlines.Line2D([], [], color='b')
perseus= mlines.Line2D([], [], color='g')
crux = mlines.Line2D([], [], color='k')
local = mlines.Line2D([], [], color='c')    
text( 12.0,  12.0, 'Norma',     verticalalignment='bottom', horizontalalignment='right',     color='red', fontsize=25)
text( -8.0, -14.0, 'Carina \n Sagitarius',     verticalalignment='bottom', horizontalalignment='right',     color='blue', fontsize=25)
text( -2.0,  11.0, 'Perseus',     verticalalignment='bottom', horizontalalignment='right',     color='green', fontsize=25)
text( 13.0, -13.0, 'Crux-Scutum',     verticalalignment='bottom', horizontalalignment='right',     color='black', fontsize=25)
text( -8.0, 5.5, 'Local',     verticalalignment='bottom', horizontalalignment='right',     color='c', fontsize=25)
text(0.0,0.0, 'GC', fontsize=15)
plot(0.0,0.0,'o',color='g',markersize=10)
xlabel('X $[kpc]$')
ylabel('Y $[kpc]$')
axis([-par.rtrunc-5.0,par.rtrunc+5.0,-par.rtrunc-5.0,par.rtrunc+5.0])
axis('equal')
plot( 0, par.r0, '*', color ='y',markersize=13)
ring = par.cr.central_ring(par.nspirplot) 
plot(ring[0]['x'], ring[0]['y'], 'm--', linewidth = 3)
text(0.0,0.0, 'GC', fontsize=15)
#Jimmi's data
scatter(x_jimmi, y_jimmi, color='m',alpha=0.5)
plot(50, 50, '.', color='m',alpha=0.5,markersize =20.0,label='MMB data')
xlim([-20.0,20.0])
ylim([-15.2,15.2])
legend(loc=2,numpoints=1,prop={'size':27})
#Lines
x_line_1=np.arange(0.0,13.1,0.1)
x_line_2=np.arange(-13.0,0.1,0.1)
y_line_1=-1.535*x_line_1+par.r0
y_line_2=1.535*x_line_2+par.r0
plot(x_line_1,y_line_1,'--',color='orange',linewidth = 2)
plot(x_line_2,y_line_2,'--',color='orange',linewidth = 2)
fill_between(x_line_1,min(y_line_1), y_line_1, facecolor='orange',alpha=0.3,interpolate=True)
fill_between(x_line_2,min(y_line_1), y_line_2, facecolor='orange',alpha=0.3,interpolate=True)
plt.grid(True)
#plt.savefig('../plots/out_jimmi.pdf',bbox_inches='tight')

#--------------------------------------------------------------------------------------------------------------------    
#plot jimmi 3
#Jimis data plot
t_Jy = ascii.read('../surveys/MMB_Dist_Flux_edited.csv')    
t_Jy.remove_row(0)
kk=[]
for i in range(len(t_Jy)):
    if(t_Jy[i][5]=='-' or t_Jy[i][0]=='-' or t_Jy[i][2]=='-' ):
    #if(t_Jy[i][5]=='-' or t_Jy[i][0]=='-'):
        kk.append(i)
    
t_Jy.remove_rows(kk)

zz=[]
hh=[]
jj=[]
for i in range(len(t_Jy)):
    if(float(t_Jy[i][0])<32.144 and float(t_Jy[i][0])>0.0):
        zz.append(i)        
    elif(float(t_Jy[i][0])<359.99 and float(t_Jy[i][0])>327.856):
        hh.append(i)        
    else:
        jj.append(i)
        
t_Jy.remove_rows(zz)
t_Jy.remove_rows(hh)
        
distance_jimmi=[]
l_jimmi=[]
b_jimmi=[]
for i in range(len(t_Jy)):            
        distance_jimmi.append(float(t_Jy[i][2]))
        l_jimmi.append(float(t_Jy[i][0])*pi/180)
        b_jimmi.append(float(t_Jy[i][1])*pi/180)    

x_jimmi=np.zeros(len(distance_jimmi))
y_jimmi=np.zeros(len(distance_jimmi))
for i in range(len(distance_jimmi)):
    x_jimmi[i]=cos((3*pi/2)+l_jimmi[i])*distance_jimmi[i]
    y_jimmi[i]=sin((3*pi/2)+l_jimmi[i])*distance_jimmi[i]+par.r0
    
figure(figsize=(10,10))
#matplotlib.rcParams.update({'font.size': 27}) 
segs = par.gs.segs_spiral(par.rtrunc, par.nspirplot) # why different than par.nspirseg? OJO
for j in range( len(segs) ):
    line_arm_color=['r-','b-','g-','k-','c-']
    plot( segs[j]['x'], segs[j]['y'], line_arm_color[j], linewidth=3)
norma = mlines.Line2D([], [], color='r')
carina = mlines.Line2D([], [], color='b')
perseus= mlines.Line2D([], [], color='g')
crux = mlines.Line2D([], [], color='k')
local = mlines.Line2D([], [], color='c')    
text( 12.0,  12.0, 'Norma',     verticalalignment='bottom', horizontalalignment='right',     color='red', fontsize=25)
text( -8.0, -14.0, 'Carina \n Sagitarius',     verticalalignment='bottom', horizontalalignment='right',     color='blue', fontsize=25)
text( -2.0,  11.0, 'Perseus',     verticalalignment='bottom', horizontalalignment='right',     color='green', fontsize=25)
text( 13.0, -13.0, 'Crux-Scutum',     verticalalignment='bottom', horizontalalignment='right',     color='black', fontsize=25)
text( -8.0, 5.5, 'Local',     verticalalignment='bottom', horizontalalignment='right',     color='c', fontsize=25)
text(0.0,0.0, 'GC', fontsize=15)
plot(0.0,0.0,'o',color='g',markersize=10)
xlabel('X $[kpc]$')
ylabel('Y $[kpc]$')
axis([-par.rtrunc-5.0,par.rtrunc+5.0,-par.rtrunc-5.0,par.rtrunc+5.0])
axis('equal')
plot( 0, par.r0, '*', color ='y',markersize=13)
ring = par.cr.central_ring(par.nspirplot) 
plot(ring[0]['x'], ring[0]['y'], 'm--', linewidth = 3)
text(0.0,0.0, 'GC', fontsize=15)
#Jimmi's data
scatter(x_jimmi, y_jimmi, color='m',alpha=0.5)
plot(50, 50, '.', color='m',alpha=0.5,markersize =20.0,label='MMB data')
xlim([-20.0,20.0])
ylim([-15.2,15.2])
legend(loc=2,numpoints=1,prop={'size':27})
#Lines
'''
x_line_1=np.arange(0.0,13.1,0.1)
x_line_2=np.arange(-13.0,0.1,0.1)
y_line_1=-1.535*x_line_1+par.r0
y_line_2=1.535*x_line_2+par.r0
plot(x_line_1,y_line_1,'--',color='orange',linewidth = 2)
plot(x_line_2,y_line_2,'--',color='orange',linewidth = 2)
fill_between(x_line_1,min(y_line_1), y_line_1, facecolor='orange',alpha=0.3,interpolate=True)
fill_between(x_line_2,min(y_line_1), y_line_2, facecolor='orange',alpha=0.3,interpolate=True)
'''
plt.grid(True)
#plt.savefig('../plots/out_jimmi2.pdf',bbox_inches='tight')

#-------------------------------------------------------------------------------------
#plot Jimmi 4
t_Jy = ascii.read('../surveys/MMB_Dist_Flux_edited.csv')    
t_Jy.remove_row(0)
kk=[]
for i in range(len(t_Jy)):
    if(t_Jy[i][5]=='-' or t_Jy[i][0]=='-' or t_Jy[i][1]=='-' ):
    #if(t_Jy[i][5]=='-' or t_Jy[i][0]=='-'):
        kk.append(i)
    
t_Jy.remove_rows(kk)


l_jimmi=[]
b_jimmi=[]
flux_jimmi=[]
dis_confirm=[]
for i in range(len(t_Jy)):            
        l_jimmi.append(float(t_Jy[i][0]))
        b_jimmi.append(float(t_Jy[i][1]))    
        flux_jimmi.append(float(t_Jy[i][5]))
        if(t_Jy[i][2]=='-'):
            dis_confirm.append(0)
        else:
            dis_confirm.append(1)
#lat-long- view
figure(figsize=(25,25))
matplotlib.rcParams.update({'font.size': 18})
map = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')
#map = Basemap(projection='stere',lon_0=0,lat_0=0.0,llcrnrlat=-5,urcrnrlat=5,\
#          llcrnrlon=-150,urcrnrlon=150,resolution='l')
tl=numpy.zeros(len(l_jimmi))
tb=numpy.zeros(len(l_jimmi))
for i in range( len(l_jimmi) ):
    tl[i], tb[i] = map(l_jimmi[i],b_jimmi[i])
    if(dis_confirm[i]==1):
        map.plot( tl[i], tb[i], 'o', alpha = 0.5 , color = 'g')#, markersize = (numpy.log10(flux_jimmi[i])*2)+21)
        z=i
    else:
        map.plot( tl[i], tb[i], 'o', alpha = 0.5 , color = 'r')#, markersize = (numpy.log10(flux_jimmi[i])*2)+21)        
map.drawparallels(numpy.arange(-90.,120.,30.),labels=[True])
map.drawmeridians(np.arange(0.,420.,30.))
map.drawmapboundary(fill_color='GhostWhite')
Dis_dot, =map.plot( tl[z], tb[z], 'o', alpha = 0.5 , color = 'g', markersize = 5)
Nodist_dot,   =map.plot( tl[z], tb[z], 'o', alpha = 0.5 , color = 'r', markersize = 5)
lg=legend([Dis_dot, Nodist_dot],['Distance Info.','No distance Info.'],loc=8,numpoints=1,prop={'size':16})
plt.gca().add_artist(lg)
name_plot_3="l-b_jimmi.pdf"
plt.savefig(name_plot_3,bbox_inches='tight')
#-------------------------------------------------------------------------------------
#Central Ellipse
def ellipse(ra,rb,ang,x0,y0,Nb=50):
    '''ra - major axis length
    rb - minor axis length
    ang - angle
    x0,y0 - position of centre of ellipse
    Nb - No. of points that make an ellipse
    
    based on matlab code ellipse.m written by D.G. Long,
    Brigham Young University, based on the
    CIRCLES.m original
    written by Peter Blattner, Institute of Microtechnology,
    University of
    Neuchatel, Switzerland, blattner@imt.unine.ch
    '''
    xpos,ypos=x0,y0
    radm,radn=ra,rb
    an=ang*pi/180.0
    
    co,si=numpy.cos(an),numpy.sin(an)
    the=numpy.linspace(0,2*pi,Nb)
    X=radm*numpy.cos(the)*co-si*radn*numpy.sin(the)+xpos
    Y=radm*numpy.cos(the)*si+co*radn*numpy.sin(the)+ypos
    return X,Y
    
ra=4.1
rb=2.2
ang=38# (PA=-38)
x0=0.0
y0=0.0

x_ell, y_ell = ellipse(ra,rb,ang,x0,y0)


figure(figsize=(10,10))
#matplotlib.rcParams.update({'font.size': 27}) 
segs = par.gs.segs_spiral(par.rtrunc, par.nspirplot) # why different than par.nspirseg? OJO
for j in range( len(segs) ):
    line_arm_color=['r-','b-','g-','k-','c-']
    plot( segs[j]['x'], segs[j]['y'], line_arm_color[j], linewidth=3)
norma = mlines.Line2D([], [], color='r')
carina = mlines.Line2D([], [], color='b')
perseus= mlines.Line2D([], [], color='g')
crux = mlines.Line2D([], [], color='k')
local = mlines.Line2D([], [], color='c')    
text( 12.0,  12.0, 'Norma',     verticalalignment='bottom', horizontalalignment='right',     color='red', fontsize=25)
text( -8.0, -14.0, 'Carina \n Sagitarius',     verticalalignment='bottom', horizontalalignment='right',     color='blue', fontsize=25)
text( -2.0,  11.0, 'Perseus',     verticalalignment='bottom', horizontalalignment='right',     color='green', fontsize=25)
text( 13.0, -13.0, 'Crux-Scutum',     verticalalignment='bottom', horizontalalignment='right',     color='black', fontsize=25)
text( -8.0, 5.5, 'Local',     verticalalignment='bottom', horizontalalignment='right',     color='c', fontsize=25)

xlabel('X $[kpc]$')
ylabel('Y $[kpc]$')
axis([-par.rtrunc-5.0,par.rtrunc+5.0,-par.rtrunc-5.0,par.rtrunc+5.0])
axis('equal')
plot( 0, par.r0, '*', color ='y',markersize=13)
plot(x_ell, y_ell, 'm--', linewidth = 2)
text(0.0,0.0, 'GC', fontsize=15)
plot(0.0,0.0,'o',color='g',markersize=10)
plt.grid(True)






from mpl_toolkits.basemap import Basemap, cm
# requires netcdf4-python (netcdf4-python.googlecode.com)
from netCDF4 import Dataset as NetCDFFile
import numpy as np
import matplotlib.pyplot as plt

# plot rainfall from NWS using special precipitation
# colormap used by the NWS, and included in basemap.

nc = NetCDFFile('/Users/luishenryquiroganunez/Downloads/nws_precip_20160106_nc/nws_precip_conus_20160106.nc')
# data from http://water.weather.gov/precip/
prcpvar = nc.variables['amountofprecip']
data = 0.01*prcpvar[:]
latcorners = nc.variables['lat'][:]
loncorners = -nc.variables['lon'][:]
lon_0 = -nc.variables['true_lon'].getValue()
lat_0 = nc.variables['true_lat'].getValue()
# create figure and axes instances
fig = plt.figure(figsize=(8,8))
ax = fig.add_axes([0.1,0.1,0.8,0.8])
# create polar stereographic Basemap instance.
m = Basemap(projection='stere',lon_0=lon_0,lat_0=90.,lat_ts=lat_0,\
            llcrnrlat=latcorners[0],urcrnrlat=latcorners[2],\
            llcrnrlon=loncorners[0],urcrnrlon=loncorners[2],\
            rsphere=6371200.,resolution='l',area_thresh=10000)
# draw coastlines, state and country boundaries, edge of map.
m.drawcoastlines()
m.drawstates()
m.drawcountries()
# draw parallels.
parallels = np.arange(0.,90,10.)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
meridians = np.arange(180.,360.,10.)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
ny = data.shape[0]; nx = data.shape[1]
lons, lats = m.makegrid(nx, ny) # get lat/lons of ny by nx evenly space grid.
x, y = m(lons, lats) # compute map proj coordinates.
# draw filled contours.
clevs = [0,1,2.5,5,7.5,10,15,20,30,40,50,70,100,150,200,250,300,400,500,600,750]
cs = m.contourf(x,y,data,clevs,cmap=cm.s3pcpn)
# add colorbar.
cbar = m.colorbar(cs,location='bottom',pad="5%")
cbar.set_label('mm')
# add title
plt.title(prcpvar.long_name+' for period ending '+prcpvar.dateofdata)
plt.show()