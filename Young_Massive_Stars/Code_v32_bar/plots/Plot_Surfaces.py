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
#------------------------------------------------------------------------
#Path
%cd "/Users/luishenryquiroganunez/Documents/Leiden_University/3_Semester/MSc_thesis/Codes/Code_v31/surface_bestnalpha/Matrix_ind_num_nbar8_sigma4"
#Importing Libraries
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
#------------------------------------------------------------------------
#Reading the cube
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
#------------------------------------------------------------------------
#Getting the data
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
#------------------------------------------------------------------------
#Colapsing the cubes
final_MMB_num=numpy.mean(cube_MMB_num,axis=2)
final_MMB_ind=numpy.mean(cube_MMB_ind,axis=2)
final_MMB_chi2=numpy.mean(cube_MMB_chi2,axis=2)
final_MMB_KLD=numpy.mean(cube_MMB_KLD,axis=2)
final_Arecibo_num=numpy.mean(cube_Arecibo_num,axis=2)
final_Arecibo_ind=numpy.mean(cube_Arecibo_ind,axis=2)
final_Arecibo_chi2=numpy.mean(cube_Arecibo_chi2,axis=2)
final_Arecibo_KLD=numpy.mean(cube_Arecibo_KLD,axis=2)
############################################################################################################
#PLOTSSSS!
############################################################################################################










# One Surface Plot
#h = plt.contourf(num[:,0],ind[0,:],MMB_n)
h=imshow(final_MMB_num, interpolation='none', origin='bottom', extent=[numpy.min(num[:,0]), numpy.max(num[:,0]), numpy.min(ind[0,:]), np.max(ind[0,:])],aspect='auto')
plt.colorbar(h, orientation='vertical')
plt.xlabel("Initial Number of Sources")
plt.ylabel("Initial Index value")
plt.show()
#------------------------------------------------------------------------
# Four Surface Plot
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
a=ax1.imshow(final_MMB_num, interpolation='none', origin='bottom', extent=[numpy.min(num[:,0]), numpy.max(num[:,0]), numpy.min(ind[0,:]), np.max(ind[0,:])],aspect='auto')
#ax1.colorbar(a, orientation='vertical')
ax1.set_title('MMB Survey')
b=ax2.imshow(final_Arecibo_num, interpolation='none', origin='bottom', extent=[numpy.min(num[:,0]), numpy.max(num[:,0]), numpy.min(ind[0,:]), np.max(ind[0,:])],aspect='auto')
ax2.set_title('Arecibo Survey')
c=ax3.imshow(final_MMB_ind, interpolation='none', origin='bottom', extent=[numpy.min(num[:,0]), numpy.max(num[:,0]), numpy.min(ind[0,:]), np.max(ind[0,:])],aspect='auto')
d=ax4.imshow(final_Arecibo_ind, interpolation='none', origin='bottom', extent=[numpy.min(num[:,0]), numpy.max(num[:,0]), numpy.min(ind[0,:]), np.max(ind[0,:])],aspect='auto')
fig.subplots_adjust(right=0.8)
cbar_ax_a = fig.add_axes([0.85, 0.53, 0.05, 0.37])#[left, bottom, width, height] 
fig.colorbar(a, cax=cbar_ax_a)
cbar_ax_b = fig.add_axes([0.85, 0.1, 0.05, 0.37])#[left, bottom, width, height] 
fig.colorbar(b, cax=cbar_ax_b)
plt.show()
#------------------------------------------------------------------------
# Only MMB Surface Plot
#With KLD
matplotlib.rcParams.update({'font.size': 17})
fig, ax1 = plt.subplots(1, 1)
a=ax1.imshow(final_MMB_KLD, interpolation='none', origin='bottom', extent=[numpy.min(num[:,0]), numpy.max(num[:,0]), numpy.min(ind[0,:]), np.max(ind[0,:])],aspect='auto')
#ax1.colorbar(a, orientation='vertical')
ax1.set_title('Initial Conditions $N_i$ and $\\alpha$ compared MMB Survey')
ax1.set_ylabel('Initial Index Value $(\\alpha)$')
ax1.set_xlabel('Initial Number of Sources $(N_i)$')
fig.subplots_adjust(right=0.8)
cbar_ax_a = fig.add_axes([0.85, 0.1, 0.05, 0.8])#[left, bottom, width, height] 
fig.colorbar(a, cax=cbar_ax_a)
plt.ylabel("$KL$ Differentation")
fig.set_label('Initial Index Value')
plt.show()

#With chi 2
matplotlib.rcParams.update({'font.size': 17})
fig, ax1 = plt.subplots(1, 1)
a=ax1.imshow(numpy.log(final_MMB_chi2), interpolation='none', origin='bottom', extent=[numpy.min(num[:,0]), numpy.max(num[:,0]), numpy.min(ind[0,:]), np.max(ind[0,:])],aspect='auto')
#ax1.colorbar(a, orientation='vertical')
ax1.set_title('Initial Conditions $N_i$ and $\\alpha$ compared MMB Survey')
ax1.set_ylabel('Initial Index Value $(\\alpha)$')
ax1.set_xlabel('Initial Number of Sources $(N_i)$')
fig.subplots_adjust(right=0.8)
cbar_ax_a = fig.add_axes([0.85, 0.1, 0.05, 0.8])#[left, bottom, width, height] 
fig.colorbar(a, cax=cbar_ax_a)
plt.ylabel("log($\\chi^2$)")
fig.set_label('Initial Index Value')
plt.show()
#pyplot.savefig("Surface.png",bbox_inches='tight')

#With Num and Ind
matplotlib.rcParams.update({'font.size': 17})
fig, (ax1,ax2) = plt.subplots(2, 1, sharex='col')
a=ax1.imshow(MMB_n, interpolation='none', origin='bottom', extent=[numpy.min(num[:,0]), numpy.max(num[:,0]), numpy.min(ind[0,:]), np.max(ind[0,:])],aspect='auto')
#ax1.colorbar(a, orientation='vertical')
ax1.set_title('Initial Conditions $N_i$ and $\\alpha$ compared MMB Survey')
ax1.set_ylabel('Initial Index Value $(\\alpha)$')
b=ax2.imshow(MMB_i, interpolation='none', origin='bottom', extent=[numpy.min(num[:,0]), numpy.max(num[:,0]), numpy.min(ind[0,:]), np.max(ind[0,:])],aspect='auto')
ax2.set_xlabel('Initial Number of Sources $(N_i)$')
ax2.set_ylabel('Initial Index Value $(\\alpha)$')
fig.subplots_adjust(right=0.8)
cbar_ax_a = fig.add_axes([0.85, 0.53, 0.05, 0.37])#[left, bottom, width, height] 
fig.colorbar(a, cax=cbar_ax_a)
plt.ylabel("|$N-N_{obs}$|")
cbar_ax_b = fig.add_axes([0.85, 0.1, 0.05, 0.37])#[left, bottom, width, height] 
fig.colorbar(b, cax=cbar_ax_b)
plt.ylabel("|$ \\beta-\\beta_{obs}$|")
fig.set_label('Initial Index Value')
plt.show()

#------------------------------------------------------------------------
# Only Arecibo Surface Plot
from matplotlib.colors import LogNorm
matplotlib.rcParams.update({'font.size': 17})
fig, (ax1,ax2) = plt.subplots(2, 1, sharex='col')
a=ax1.imshow(Arecibo_n, interpolation='none', origin='bottom', extent=[numpy.min(num[:,0]), numpy.max(num[:,0]), numpy.min(ind[0,:]), np.max(ind[0,:])],aspect='auto',norm=LogNorm(vmin=1, vmax=100))
#ax1.colorbar(a, orientation='vertical')
ax1.set_title('Initial Conditions $N_i$ and $\\alpha$ compared Arecibo Survey')
ax1.set_ylabel('Initial Index Value $(\\alpha)$')
b=ax2.imshow(Arecibo_i, interpolation='none', origin='bottom', extent=[numpy.min(num[:,0]), numpy.max(num[:,0]), numpy.min(ind[0,:]), np.max(ind[0,:])],aspect='auto',norm=LogNorm(vmin=0.1, vmax=1))
ax2.set_xlabel('Initial Number of Sources $(N_i)$')
ax2.set_ylabel('Initial Index Value $(\\alpha)$')
fig.subplots_adjust(right=0.8)
cbar_ax_a = fig.add_axes([0.85, 0.53, 0.05, 0.37])#[left, bottom, width, height] 
fig.colorbar(a, cax=cbar_ax_a)
plt.ylabel("|$N-N_{obs}$|")
cbar_ax_b = fig.add_axes([0.85, 0.1, 0.05, 0.37])#[left, bottom, width, height] 
fig.colorbar(b, cax=cbar_ax_b)
plt.ylabel("|$ \\beta-\\beta_{obs}$|")
fig.set_label('Initial Index Value')
plt.show()
#pyplot.savefig("Surface_Arecibo.png",bbox_inches='tight')
#------------------------------------------------------------------------
#Calculting the minimum between two plots of the MMB
#Constants of the observations
num_MMB_obs=908 
ind_MMB_obs=-0.599716013349
#Getting the data
data=numpy.loadtxt('Calc_Num_Ind_0.txt')
#Sorting the data
num_0=data[:,0]
ind_0=data[:,1]
MMB_num=abs(data[:,2]-num_MMB_obs)/num_MMB_obs
MMB_ind=abs(data[:,3]-ind_MMB_obs)/abs(ind_MMB_obs)
#data_dims=(5,5)
MMB_n = np.asarray(MMB_num).reshape(data_dims)
MMB_i= np.asarray(MMB_ind).reshape(data_dims)
Total=numpy.log10(((MMB_n+MMB_i)/2.0))
Total=Total*(-1)
Total=Total/(numpy.max(Total))
num = np.asarray(num_0).reshape(data_dims)
ind = np.asarray(ind_0).reshape(data_dims)
img_gaus = ndimage.filters.gaussian_filter(Total, 2, mode='nearest')
num_masers=numpy.linspace(numpy.min(num[:,0]),numpy.max(num[:,0]),100)
ind_masers=numpy.linspace(numpy.min(ind[0,:]),numpy.max(ind[0,:]),100)
X,Y=numpy.meshgrid(num_masers,ind_masers)


#Plot
matplotlib.rcParams.update({'font.size': 17})
#fig, ax = plt.subplots()
imshow(Total, interpolation='none', origin='bottom', extent=[numpy.min(num[:,0]), numpy.max(num[:,0]), numpy.min(ind[0,:]), np.max(ind[0,:])],aspect='auto')
#imshow(img_gaus, interpolation='none', origin='bottom', extent=[numpy.min(num[:,0]), numpy.max(num[:,0]), numpy.min(ind[0,:]), np.max(ind[0,:])],aspect='auto')
#CS = plt.contour(x, y, Total,colors='k')
#heatmap = ax.pcolor(Total)
cbar=plt.colorbar()
levels_=[0.2,0.3,0.4,0.5,0.6]
levels_range=[0.4]
levels_center=[0.4]
cs=plt.contour(X, Y, img_gaus,50, colors='k',levels=levels_)
plt.clabel(cs, inline=1, fontsize=10)
plt.xlabel('Initial Number of Sources $(N_i)$')
plt.ylabel('Initial Index Value $(\\alpha)$')
#cbar.ax.setlabel('Arbitrary Units', rotation=270)
cbar.set_label('Arbitrary Units')
#pyplot.savefig("Surface_both_MMB.png",bbox_inches='tight')


#Calcutaing the center and error
cs2=plt.contour(X, Y, img_gaus,50, colors='k',levels=levels_)
p = cs2.collections[0].get_paths()[0]
v = p.vertices
x = v[:,0]
y = v[:,1]

cs3=plt.contour(X, Y, img_gaus,50, colors='k',levels=levels_center)
p = cs2.collections[0].get_paths()[0]
v = p.vertices
x_center = v[:,0]
y_center = v[:,1]

def fitEllipse(x,y):
    x = x[:,np.newaxis]
    y = y[:,np.newaxis]
    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    S = np.dot(D.T,D)
    C = np.zeros([6,6])
    C[0,2] = C[2,0] = 2; C[1,1] = -1
    E, V =  eig(np.dot(inv(S), C))
    n = np.argmax(np.abs(E))
    a = V[:,n]
    return a
def ellipse_center(a):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    num = b*b-a*c
    x0=(c*d-b*f)/num
    y0=(a*f-b*d)/num
    return np.array([x0,y0])
aa=fitEllipse(x_center,y_center)    
center=ellipse_center(aa)
#array([ 1540.99789764,    -1.56126529])
error_Num=numpy.max([abs(numpy.max(x)-center[0]),abs(numpy.min(x)-center[0])])
error_Ind=numpy.max(abs(numpy.max(y)-center[1]),abs(numpy.min(y)-center[1]))
#(172.83114895668996, 0.22829395805862451)



#------------------------------------------------------------------------
#Calculting the minimum of the MMB (ESAT ES!)
img_gaus = ndimage.filters.gaussian_filter(numpy.log(final_MMB_chi2), 2, mode='nearest')
num_masers=numpy.linspace(numpy.min(num[:,0]),numpy.max(num[:,0]),100)
ind_masers=numpy.linspace(numpy.min(ind[0,:]),numpy.max(ind[0,:]),100)
X,Y=numpy.meshgrid(num_masers,ind_masers)
#Plot
matplotlib.rcParams.update({'font.size': 27})
fig,ax1= plt.subplots(1, 1,figsize=(12,12))
a=ax1.imshow(numpy.log(final_MMB_chi2), interpolation='nearest', origin='bottom', extent=[numpy.min(num[:,0]), numpy.max(num[:,0]), numpy.min(ind[0,:]), np.max(ind[0,:])],aspect='auto')
ax1.set_xticks([1000,1200,1400,1600,1800])
ax1.set_yticks([-1.2,-1.4,-1.6,-1.8,-2.0])
#ax1.set_title('Initial Conditions $N_i$ and $\\alpha$ compared MMB Survey')
ax1.set_ylabel('Initial Index Value $(\\alpha)$')
ax1.set_xlabel('Initial Number of Sources $(N)$')
fig.subplots_adjust(right=0.8)
cbar_ax_a = fig.add_axes([0.85, 0.1, 0.05, 0.8])#[left, bottom, width, height] 
cbar=fig.colorbar(a, cax=cbar_ax_a)
levels_=numpy.array([4.65,4.7,4.9,5.4,5.7,6.0,7.2,8.4])
#levels_=numpy.array([4.65, 4.7, 4.8])
cs=ax1.contour(X, Y, img_gaus,hold='on', levels=levels_, colors='k',origin='image')
#ax1.clabel(cs, inline=5, fontsize=10)
cbar.set_label('log($\\xi^2$)')
cbar.set_ticks([5.0,6.0,7.0,8.0])
#pyplot.savefig("../../plots/Surface_both_MMB.pdf",bbox_inches='tight')

#levels_center=[4.65]
levels_center=[4.9]
cs3=plt.contour(X, Y, img_gaus, colors='k',levels=levels_center)
p = cs3.collections[0].get_paths()[0]
v = p.vertices
x_center = v[:,0]
y_center = v[:,1]

#3D


def fitEllipse(x,y):
    x = x[:,np.newaxis]
    y = y[:,np.newaxis]
    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    S = np.dot(D.T,D)
    C = np.zeros([6,6])
    C[0,2] = C[2,0] = 2; C[1,1] = -1
    E, V =  eig(np.dot(inv(S), C))
    n = np.argmax(np.abs(E))
    a = V[:,n]
    return a
def ellipse_center(a):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    num = b*b-a*c
    x0=(c*d-b*f)/num
    y0=(a*f-b*d)/num
    return np.array([x0,y0])
    
aa=fitEllipse(x_center,y_center)    
center=ellipse_center(aa)
#array([ 1298.10797886,    -1.43666962])

levels_=numpy.array([4.7])
#Calcutaing the center and error
cs2=plt.contour(X, Y, img_gaus, colors='k',levels=levels_)
p = cs2.collections[0].get_paths()[0]
v = p.vertices
x = v[:,0]
y = v[:,1]

error_Num=numpy.max([abs(numpy.max(x)-center[0]),abs(numpy.min(x)-center[0])])
error_Ind=numpy.max(abs(numpy.max(y)-center[1]),abs(numpy.min(y)-center[1]))
#140, 0.33

###############################################################################################
#Calculate  beta and N for one same galaxy
###############################################################################################
%cd "/Users/luishenryquiroganunez/Documents/Leiden_University/3_Semester/MSc_thesis/Codes/Code_v31/surface_bestnalpha/One_same_galaxy"
names_txt_files=os.listdir('.')
data=numpy.loadtxt(names_txt_files[1])
N_MMB=np.mean(data[0])
err_N_MMB=np.std(data[0])
alpha_MMB=np.mean(data[1])
err_alpha_MMB=np.std(data[1])
N_Arecibo=np.mean(data[2])
err_N_Arecibo=np.std(data[2])
alpha_Arecibo=np.mean(data[3])
err_alpha_Arecibo=np.std(data[3])
print ('Number of sources detected (N) in MMB (%s +/- %s) and Arecibo (%s +/- %s)' %(N_MMB, err_N_MMB, N_Arecibo, err_N_Arecibo))
print ('Slope Flux Desnisty Function (alpha) in MMB (%s +/- %s) and Arecibo (%s +/- %s)' %(alpha_MMB, err_alpha_MMB, alpha_Arecibo, err_alpha_Arecibo))
