#-------------------------------------------------------------------------------
#Importing Libraries
import matplotlib
import sys
import numpy
global debug, interactive
debug = 0
interactive = 1
#if (interactive):
#    matplotlib.use('MacOSX')
#else:
#    matplotlib.use('pdf')
from pylab import *
from math import *
from random import *
import csv
from time import clock, time, localtime
import ephem
from aux.gameModelPars31 import *
from aux.gamePlotUtils31 import *
from aux.gameFunctions31 import *
#-------------------------------------------------------------------------------
def full_model(file_name,folder_,conf_size=0.0,conf_alpha=0.0): #(npar,index_):    
    #-------------------------------------------------------------------------------
    #execution 1 then goes to def __init__(self,r0=8.5,nobj=100) in gameModelPars5.py
    #Initializes global paramters...see def __init__(self,r0=8.5,nobj=100) in gameModelPars5.py
    if(conf_size != 0 and conf_alpha != 0):
        par = GaMePars(file_name,folder,conf_size,conf_alpha) 
    else:
        par = GaMePars(file_name,folder_) 
    #-------------------------------------------------------------------------------
    stc = clock()
    stt = time()
    prt1='-------------------------------------------------------------------------------'
    prt2='starting time: '
    prt3= time()-stt 
    print prt1 
    print prt2,
    print clock()-stc,time()-stt 
    #-------------------------------------------------------------------------------
    #Extra Functions
    def segs_spiralxy(xp,yp,segs):
        #execution5.2
        """Finds the closest line segment over all segments of all arms"""
        #xp and yp are a random values between -rtunc(-15kpc) and rtrunc(15kpc)
        # segs are the segments (40 points=nspirseg) of the spiral arms
        dist = 1e9 #beyond Andromeda declaring variable?
        #this is run every time, could be static
        ic = -1
        # j is the number of segments
        for j in range( len(segs) ):
            # i are the points of each segements (nspirseg=40=maxi)
            for i in range( len( segs[j]['x'] )-1):
                d = dist_seg( xp, yp, segs[j]['x'][i], segs[j]['y'][i], \
                                    segs[j]['x'][i+1], segs[j]['y'][i+1] )
                                    #d is the distance from xp,yp, to a tiny segment of an arm
                if debug:print 'segs>>',i,j, d
                if ( d < dist ):
                    dist = d # saving distance as d
                    ic = j # saving the color of the arm as ic
    
        p = mygauss(dist, par.fat) # gaussian distribution
        if debug: print xp, yp, dist, p 
        return p, ic
        
    def draw_sample( nobj, segs):
        #execution 5.1 create disperison points around the arms and for each point give it atributes
        #as intrinsic luminosty (f, power law), tangential velocity (va, gaussian), 
        #radial velocity (vr, gaussian), a-velocity (vz, gaussian), x,y,z(gaussian) position
        """Draws a sample, uses a uniform x,y and a chance to get close to spirals
        all other parameters are drawn from independent Gaussians"""
        sample = []
        ndraw = 0
        cte=(par.v0t)-(par.rot_curve*par.r0)
        while (len(sample) < nobj):
            #if (len(sample) > 0): print sample[len(sample)-1]
            ndraw += 1
            # r_o and thetas: distances and angles to map the xy diagram, rtunc= maximum distance of the spiral arms
            aa=1-(1/(exp(par.rmintrunc/par.halfmr)))
            bb=1-(1/(exp(par.rtrunc/par.halfmr)))
            r_o=-par.halfmr*numpy.log(1-numpy.random.uniform(aa,bb))
            tethas=numpy.random.uniform(0,2*pi)
            x=r_o*cos(tethas) 
            y=r_o*sin(tethas) 
            """This was done if you need really complex functions
            z = random()*2*ztrunc-ztrunc
            va = random()*2*vtrunc-vtrunc
            vr = random()*2*vtrunc-vtrunc
            vz = random()*2*vtrunc-vtrunc
            """
            p = random()
            #p=0.5
            #calculate to which spiral arm is close the point (x,y) and with the distance make a gaussian distribution
            pxy, ic = segs_spiralxy(x, y, segs)
            #could add all the other non simple function
            # and multiply the outcomes
            #do something fast
            if (p < pxy ):
                #now add the functions that can be drawn directly:
                z = gauss( 0, par.hz ) #mean poisition in z, standart deviation of poisition in z
                #r_dis=numpy.sqrt(x**2+y**2+z**2) # in kpc
                #pr=scl_len(r_dis)
                #p_2=random()
                #if (pr > p_2):  
                distan=numpy.sqrt(x**2+y**2+z**2)
                va_rot = (par.rot_curve*distan)+cte        
                #va = gauss( par.v0t, par.sigva ) - gauss( par.v0t_vs, par.sigva_vs)#mean tangential velocity (240), standart deviation of tangential velocity (15) + mean source motion in V
                va = gauss( va_rot,  par.sigva ) - gauss( par.v0t_vs, par.sigva_vs)#  vlsr with rotation curve, standart deviation of tangential velocity (15) + mean source motion in V
                vr = gauss( par.v0r, par.sigvr ) + gauss( par.v0r_us, par.sigvr_us ) #mean radial velocity, standart deviation of radial velocity + mean source motion in U
                vz = gauss( par.v0z, par.sigvz ) #mean z-velocity, standart deviation of z-velocity
                #Luminosity Function in Luminosities
                f=numpy.power((numpy.random.uniform(0,1)*(par.index+1)/par.k)+(numpy.power(par.min_l,par.index+1)),(1/(par.index+1)))# solving the past equation for L
                point = { 'x':x, 'y':y, 'z':z, 'va':va, 'vr':vr, 'vz':vz, 'f':f, 'ic':ic}#,'f_jy':f_jy}
                sample.append(point)
        print ('Sample Created!')       
        print ndraw, ' draws for ',len(sample), ' points'
        return sample
          
    # execution 4 create the 5 spiral arms limitated by a trunctaion radius (15Kpc) and a number of segments of the arms (40)
    #then goes to execution 5 draw sampe, just after these lines
    segs = par.gs.segs_spiral(par.rtrunc,par.nspirseg)
    #segs[each arm=0:5]['x'or'y']= 40 points
    #Check: 
    #plot(segs[0]['x'],segs[0]['y'],label='Norma Arm')
    #plot(segs[1]['x'],segs[1]['y'],label='Carina-Sagitarius Arm')
    #plot(segs[2]['x'],segs[2]['y'],label='Perseus Arm')
    #plot(segs[3]['x'],segs[3]['y'],label='Crux-Scutum Arm')
    #plot(segs[4]['x'],segs[4]['y'],label='Local Arm')
    #legend()
    #-------------------------------------------------------------------------------
    # execution 5 create disperison points around the arms and for each point give it atributes
    #as intrinsic luminosty (f, power law),identification arm and color (ic), tangential velocity (va, gaussian), 
    #radial velocity (vr, gaussian), z-velocity (vz, gaussian), x,y,z(gaussian) position
    #then goes to execution 6 translbeq, just after these lines.
    dsample = draw_sample( par.nobj, segs)
    #dsample[each point=0:npar]-----> each point has 8 atributes
    #Check:
    #aa=np.zeros(npar)
    #for i in range (npar):
    #    aa[i]=dsample[i]['x']
    #bb=np.zeros(npar)    
    #for i in range (npar):
    #    bb[i]=dsample[i]['y']
    #scatter(aa,bb)
    #-------------------------------------------------------------------------------
    #print 'sample: ',clock()-stc,time()-stt 
    #-------------------------------------------------------------------------------
    # execution 6 give more atributes to each point: Return sample with Galactic and equatorial coordinates, observed luminosity
    #added.
    #as galactic latitude (b), disatnce to the sun (d), declination (dc), intrinsic luminosty (f, power law),
    #observed luminosity (fa), identification arm and color (ic), galactic longitude (l), right ascention (ra)
    #tangential velocity (va, gaussian), velocity seen by the sun (radial) (vl), radial velocity GC (vr, gaussian),
    #velocity seen by the sun (tangential) vt, z-velocity (vz, gaussian), x,y,z(gaussian) position
    #gsample = translb( dsample )
    gsample = translbeq( dsample,par )
    #gsample[each point=0:npar]-----> each point has 16 atributes
    #-------------------------------------------------------------------------------
    stamp = 'D%02d%02dT%02d%02d%02d' % localtime()[1:6]
    #stamp=stamp_i+str(int(round(random()*100000,0)))
    #-------------------------------------------------------------------------------
    #Creating output folder
    prev_path=os.getcwd()
    directory='output'
    if not os.path.exists(directory):
        os.makedirs(directory)
    os.chdir(directory)
    #-------------------------------------------------------------------------------
    #Saving parameters in a file
    fname = 'gameN%d@%s.csv' % (len(gsample), stamp)
    fname_directory='si'+file_name[6:-4]+'/'+stamp
    print 'Number or Code of this simulation:'
    print stamp
    '''
    if not os.path.exists(fname_directory):
        os.makedirs(fname_directory)
    '''
    try: 
        os.makedirs(fname_directory)
    except OSError:
        if not os.path.isdir(fname_directory):
            raise
    os.chdir(fname_directory)
    os.chdir(prev_path)
    path_csv='./'+directory+'/'+fname_directory+'/'+fname
    prt7 = 'csv file saved in:'
    prt8 = path_csv
    print ('csv file saved in:') 
    print path_csv
    fcsv = open(path_csv, 'wb')
    wdict = csv.DictWriter(fcsv, gsample[0].keys() )
    wdict.writerow(dict(zip(gsample[0].keys(),gsample[0].keys() )))
    #in python 2.7
    #wdict.writeheader(sample[0].keys() )
    wdict.writerows(gsample)
    fcsv.close()
    comp_MMB=Comparison_MMB(gsample,par)
    comp_Arecibo=Comparison_Arecibo(gsample,par)
    #-------------------------------------------------------------------------------
    if interactive:
        show()
    else:
        savefig('./test.pdf')
    prt9='Simulation done in (sec): '
    prt10=time()-stt 
    print 'Simulation done in (sec): ',clock()-stc,time()-stt 
    #-------------------------------------------------------------------------------
    #saving outputs
    return(gsample,par,fname,fname_directory,comp_MMB,comp_Arecibo,str(prt1)+'\n',str(prt2)+'\n',str(prt3)+'\n',str(prt7)+'\n',str(prt8)+'\n',str(prt9)+'\n',str(prt10)+'\n')