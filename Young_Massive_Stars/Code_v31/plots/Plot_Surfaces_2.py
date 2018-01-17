# -*- coding: utf-8 -*-
'''
Arecibo:
Sources that satisfy Arecibo limits= 76 (Observations)
Index factor of the observed flux density function -0.471161157475 (Observation)
################################################################################
MMB:
Sources that satisfy MMB limits= 927 (Observations)
Index factor of the observed flux density function -0.581896115547
 (Observation)
'''
################################################################################
#Path
################################################################################
%cd "/Users/luishenryquiroganunez/Documents/Leiden_University/3_Semester/MSc_thesis/Codes/Code_v31/surface_bestnalpha/Matrix_ind_num_nbar8_sigma4"
################################################################################
#Importing Libraries
################################################################################
import csv
import numpy
from pylab import *
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import astropy
from astropy.io import ascii
import matplotlib.lines as mlines
from math import *
import scipy
from scipy import interpolate
from scipy import ndimage
import os
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from numpy import linspace, meshgrid
from matplotlib.mlab import griddata
################################################################################
#Reading the cube
################################################################################
names_txt_files=[]
num_files=(len([name for name in os.listdir('.') if os.path.isfile(name)]))-1
for i in range(num_files):
    names_txt_files.append('Calc_Num_Ind_'+str(i)+'.txt')
#Creating cubes
data_dummy=numpy.loadtxt(names_txt_files[0])
size=data_dummy.shape
cube_MMB_num=numpy.zeros([int(numpy.sqrt(size[0])),int(numpy.sqrt(size[0])),len(names_txt_files)])
cube_MMB_ind=numpy.zeros([int(numpy.sqrt(size[0])),int(numpy.sqrt(size[0])),len(names_txt_files)])
cube_MMB_chi2=numpy.zeros([int(numpy.sqrt(size[0])),int(numpy.sqrt(size[0])),len(names_txt_files)])
cube_MMB_KLD=numpy.zeros([int(numpy.sqrt(size[0])),int(numpy.sqrt(size[0])),len(names_txt_files)])
cube_Arecibo_num=numpy.zeros([int(numpy.sqrt(size[0])),int(numpy.sqrt(size[0])),len(names_txt_files)])
cube_Arecibo_ind=numpy.zeros([int(numpy.sqrt(size[0])),int(numpy.sqrt(size[0])),len(names_txt_files)])
cube_Arecibo_chi2=numpy.zeros([int(numpy.sqrt(size[0])),int(numpy.sqrt(size[0])),len(names_txt_files)])
cube_Arecibo_KLD=numpy.zeros([int(numpy.sqrt(size[0])),int(numpy.sqrt(size[0])),len(names_txt_files)])
#Constants of the observations
num_MMB_obs=908#837#927
ind_MMB_obs=-0.599716013349 #-0.657976860525#-0.581896115547
num_Arecibo_obs=76
ind_Arecibo_obs=-0.471161157475 
################################################################################
#Getting the data
################################################################################
for i in range(len(names_txt_files)):
    data=numpy.loadtxt(names_txt_files[i])
    #Sorting the data
    num_0=data[:,0]
    ind_0=data[:,1]
    MMB_num=abs(data[:,2]-num_MMB_obs)
    MMB_ind=abs(abs(data[:,3])-abs(ind_MMB_obs))
    MMB_chi2=data[:,4]
    MMB_KLD=data[:,5]
    Arecibo_num=abs(data[:,6]-num_Arecibo_obs)
    Arecibo_ind=abs(abs(data[:,7])-abs(ind_Arecibo_obs))
    Arecibo_chi2=data[:,8]
    Arecibo_KLD=data[:,9]
    data_dims=(int(numpy.sqrt(size[0])),int(numpy.sqrt(size[0])))
    num = np.asarray(num_0).reshape(data_dims)
    ind = np.asarray(ind_0).reshape(data_dims)
    MMB_n = np.asarray(MMB_num).reshape(data_dims)
    MMB_i= np.asarray(MMB_ind).reshape(data_dims)
    MMB_chi = np.asarray(MMB_chi2).reshape(data_dims)
    MMB_KLD2 = np.asarray(MMB_KLD).reshape(data_dims)
    Arecibo_n= np.asarray(Arecibo_num).reshape(data_dims)
    Arecibo_i= np.asarray(Arecibo_ind).reshape(data_dims)
    Arecibo_chi = np.asarray(Arecibo_chi2).reshape(data_dims)
    Arecibo_KLD2 = np.asarray(Arecibo_KLD).reshape(data_dims)    
    cube_MMB_num[:,:,i]=MMB_n
    cube_MMB_ind[:,:,i]=MMB_i
    cube_MMB_chi2[:,:,i]=MMB_chi
    cube_MMB_KLD[:,:,i]=MMB_KLD2
    cube_Arecibo_num[:,:,i]=Arecibo_n
    cube_Arecibo_ind[:,:,i]=Arecibo_i
    cube_Arecibo_chi2[:,:,i]=Arecibo_chi
    cube_Arecibo_KLD[:,:,i]=Arecibo_KLD2
################################################################################
#Colapsing the cubes
final_MMB_num=numpy.mean(cube_MMB_num,axis=2)
final_MMB_ind=numpy.mean(cube_MMB_ind,axis=2)
final_MMB_chi2=numpy.mean(cube_MMB_chi2,axis=2)
final_MMB_KLD=numpy.mean(cube_MMB_KLD,axis=2)
final_Arecibo_num=numpy.mean(cube_Arecibo_num,axis=2)
final_Arecibo_ind=numpy.mean(cube_Arecibo_ind,axis=2)
final_Arecibo_chi2=numpy.mean(cube_Arecibo_chi2,axis=2)
final_Arecibo_KLD=numpy.mean(cube_Arecibo_KLD,axis=2)
################################################################################
#Importing data
################################################################################
Z = ndimage.filters.gaussian_filter(numpy.log(final_MMB_chi2), 2, mode='nearest')
num_masers=numpy.linspace(numpy.min(num[:,0]),numpy.max(num[:,0]),100)
ind_masers=numpy.linspace(numpy.min(ind[0,:]),numpy.max(ind[0,:]),100)
X,Y=numpy.meshgrid(num_masers,ind_masers)
################################################################################
#Getting points as 3d arrays
################################################################################
points=[]
for i in range(100):
    for j in range(100):
        aaa=[ind_masers[i],num_masers[j],Z[i,j]]
        points.append(aaa)
################################################################################
#Defining the central region to map a 2d gaussian
################################################################################
levels_center=[4.9]
fig = plt.figure(figsize=plt.figaspect(0.5))
ax = fig.gca()
cs3=ax.contour(X, Y, Z, zdir='z', offset=4, levels=levels_center,colors='k')
p = cs3.collections[0].get_paths()[0]
v = p.vertices
x_center = v[:,0]
y_center = v[:,1]
N_max,N_min=max(x_center),min(x_center)
al_max,al_min=max(y_center),min(y_center)
################################################################################
#Points inside that region
################################################################################
selected_points=[]
for i in range(len(points)):
    dummy_point=points[i]
    if(dummy_point[1]>N_min and dummy_point[1]<N_max and dummy_point[0]>al_min and dummy_point[0]<al_max):
        selected_points.append(dummy_point)
alphas=[]
Number=[]
for i in range(len(selected_points)):
    alphas.append(selected_points[i][0])
    Number.append(selected_points[i][1])   
################################################################################
#Mean
################################################################################  
al_mean=np.mean(alphas)#-1.43
Nu_mean=np.mean(Number)#1313
################################################################################
#Sigma^2
################################################################################          
alphas=[]
Number=[]
for i in range(len(selected_points)):
    alphas.append(selected_points[i][0]-al_mean)
    Number.append(selected_points[i][1]-Nu_mean)  
sig_2_al=[]
sig_2_Nu=[]
for i in range(len(alphas)):
    sig_2_al.append(alphas[i]*alphas[i])
    sig_2_Nu.append(Number[i]*Number[i])    
np.sqrt(sum(sig_2_al)/len(alphas))#0.18
np.sqrt(sum(sig_2_Nu)/len(alphas))#60              
################################################################################
#Limits for plots
################################################################################  
al_min_=al_mean-np.sqrt(sum(sig_2_al)/len(alphas))
al_max_=al_mean+np.sqrt(sum(sig_2_al)/len(alphas))
Nu_min_=Nu_mean-np.sqrt(sum(sig_2_Nu)/len(alphas))
Nu_max_=Nu_mean+np.sqrt(sum(sig_2_Nu)/len(alphas))
################################################################################################################
#Making the alpha curve
################################################################################################################
alpha_curve=[]
for i in range(len(points)):
    dummy_point=points[i]
    if(dummy_point[1]>Nu_min_ and dummy_point[1]<Nu_max_):
        alpha_curve.append(dummy_point)
x = np.arange(1200)
y= np.split(x,100)
kk=[]
for j in range(len(y)):
    b=[]
    for i in y[j]:
        a=alpha_curve[i][0]
        b.append(alpha_curve[i][2])
    k=(a,np.mean(b))
    kk.append(k)
l=[]
n=[]
for i in range(len(kk)):
    l.append(kk[i][0]) 
    n.append(kk[i][1])        
ll=np.array(l)
nn=np.array(n)
N=nn.reshape([10,10])
num_masers_=numpy.linspace(numpy.min(num[:,0]),numpy.max(num[:,0]),10)
ind_masers_=numpy.linspace(numpy.min(ind[0,:]),numpy.max(ind[0,:]),10)
X_,Y_=numpy.meshgrid(num_masers_,ind_masers_)
################################################################################################################
#Making the N curve
################################################################################################################
N_curve=[]
for i in range(len(points)):
    dummy_point=points[i]
    if(dummy_point[0]>al_min_ and dummy_point[0]<al_max_):
        N_curve.append(dummy_point)
yy=[]
for i in range(100):
    yy.append(np.arange(i,len(N_curve),100))
kkk=[]
for j in range(len(yy)):
    b=[]
    for i in yy[j]:
        a=N_curve[i][1]
        b.append(N_curve[i][2])
    k=(a,np.mean(b))
    kkk.append(k)    
h=[]
jj=[]
for i in range(len(kkk)):
    h.append(kkk[i][0]) 
    jj.append(kkk[i][1])        
hh=np.array(h)
jj=np.array(jj)
alp=jj.reshape([10,10],order='F')
##################################################################################################
#Plot        
##################################################################################################
levels_=numpy.array([4.65,4.7,4.9,5.4,5.7,6.0,7.2,8.4])

matplotlib.rcParams.update({'font.size': 21})
fig = plt.figure(figsize=(8,8))
ax = fig.gca(projection='3d')
surf = ax.plot_wireframe(X, Y, Z, rstride=8, cstride=8, colors='black',antialiased=False)
#,linewidth=0.5
cset = ax.contourf(X, Y, Z, zdir='z', offset=4.0,cmap=cm.rainbow)
#csets = ax.contourf(X_, Y_, N, zdir='x', offset=900,cmap=cm.rainbow)
#csets = ax.contourf(X, Y, Z, zdir='y', offset=-1.1,cmap=cm.rainbow)
csets = ax.contour(X, Y, Z, zdir='z', offset=4, levels=levels_,colors='k',linewidth=2)
csetss = ax.contour(X_, Y_, N, zdir='x', levels=[1313,1314],offset=900,colors='maroon',linewidths=[5,5],linestyles='dashed')
csetsss = ax.contour(X_, Y_, alp, zdir='y',levels=[-1.44,-1.43],offset=-1.1,colors='maroon',linewidths=[5,5],linestyles='dashed')
ax.set_xticks([1000,1200,1400,1600,1800])
ax.set_yticks([-1.2,-1.4,-1.6,-1.8])
ax.set_ylabel('Index Value $(\\alpha)$')
ax.set_xlabel('Number of Sources $(N)$')
ax.set_xlim([900.,1800.])
ax.set_ylim([-2., -1.1])
ax.set_zlim([3.9,  8.6864489])
ax.set_zlabel('log($\\xi^2$)')
zz=[4.0,5.0,6.0,7.0,8.0]
ax.set_zticks(zz)
ax.xaxis.labelpad = 30
ax.yaxis.labelpad = 30
ax.zaxis.labelpad = 12
ax.grid()
cbar=fig.colorbar(cset,shrink=0.7,orientation='horizontal',pad=0.08)
cbar.set_label('log($\\xi^2$)')
tt=[4.2,5.4,6.6,7.8,9.0]
cbar.set_ticks(tt)    

#, aspect=10, pad=0.15