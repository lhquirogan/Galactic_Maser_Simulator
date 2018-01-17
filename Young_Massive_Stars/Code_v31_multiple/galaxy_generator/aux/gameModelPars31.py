import numpy
from math import log, sin, cos, pow
import astropy
from astropy.io import ascii
from astropy.table import Table, Column
from aux.gameModelPars31 import *
import sys
import re
global debug
debug = 0

def removeComments(string):
    '''
    For remove /*COMMENT */ in control txt file
    '''
    string = re.sub(re.compile("/\*.*?\*/",re.DOTALL ) ,"" ,string)
    string = re.sub(re.compile("//.*?\n" ) ,"" ,string) 
    return string
    
def parse_com(name,folder):
    '''
    Reading control file and exporting as dictionary
    '''
    repair = re.compile(r'^(\s*(\S+)\s*=\s*(\S+(?:\s*\,\s*\S+)*)\s*)+$')
    recomment = re.compile(r'^\#')
    control = {}
    name=folder+name
    file = open(name,'r')
    for line in file:
        line = line.rstrip('\n')
        line = removeComments(line)
        yescomment  = recomment.match(line)
        yespair = repair.match(line)
        if yescomment:
            pass
        elif yespair:
            input = yespair.groups()
            control[input[1]]=input[2]
        else:
            print 'Crap, no match:', line
    return control
    
class GaMePars:
    """General parameters for running simulations"""
    class GalSpiral:
        """Define a specific Galactic model"""
        class Spiral:
            """A specific arm description, r=rmin*10^(ext/alpha).
            tetha=alpha*log(r/rmin)+th_min
            """
            def __init__(self,alpha,rmin,thmin,ext):
                self.alpha = alpha # winding angle
                self.rmin = rmin # inner radius of the arm
                self.thmin = thmin # Angle at the inner radius of the arm
                self.ext = ext #extension of the arm
        def __init__(self):
            """Initializes the Wainscoat Galactic spiral model, see Wainscoat et al.,
            ApJ Supplement Series, vol. 83, no. 1, p. 111-146. Pag 115"""
            self.arms = []
            self.arms.append(self.Spiral( 4.25, 3.48, 0.0, 6.0 )) #Norma-Arm
            self.arms.append(self.Spiral( 4.25, 3.48, 3.141, 6.0 )) #Carina-Sagitarius Arm
            self.arms.append(self.Spiral( 4.89, 4.90, 2.525, 6.0 )) #Perseus Arm
            self.arms.append(self.Spiral( 4.89, 4.90, 5.666, 6.0 )) #Crux-Scutum Arm
            self.arms.append(self.Spiral( 4.57, 8.10, 5.847, 0.55 )) #Local Arm 
        def segs_spiral(self,rtrunc,maxi):
            """returns Galactic log spirals as a list of dictionaries,
            specify the number of segs to map it out
            Now about coordinates Wainscoat defines spiral at GC
            with theta 0 at l=0 and progressing counterclockwise.
            This seems consistent with the puting the Sun at pos Y axis
            So X=-Rsin(theta)
            and Y=Rcos(theta)"""
            #return points x and y of each spiral arm
            segs = []
            for j in range(len(self.arms) ):
                rmax = min( rtrunc, self.arms[j].rmin*10**((self.arms[j].ext)/self.arms[j].alpha))
                rstep = (rmax - self.arms[j].rmin)/maxi
                if debug: print j, rmax, rstep
                xp = []
                yp = []
                for i in range(0,int(maxi)):
                    rrun=self.arms[j].rmin+rstep*i
                    theta = self.arms[j].alpha*log(rrun/self.arms[j].rmin)+self.arms[j].thmin
                    xp.append(-1*rrun*sin(theta))
                    yp.append(rrun*cos(theta))
                segs.append( {'x':xp, 'y':yp} )
            return segs

    class CentralRing:
        """Define the central ring in the MW"""
        class EllipseRing:
            """A specific ellipse description
            """
            def __init__(self,major,minor,angle,x_center,y_center):
                self.major = major 
                self.minor = minor 
                self.angle = angle 
                self.x_center = x_center 
                self.y_center = y_center                
        def __init__(self):
            """Initializes the Central Ring model, see Sanna and Green,
            """
            self.ellipse = []
            self.ellipse.append(self.EllipseRing( 4.1, 2.2, 38.0 , 0.0, 0.0))# PA= -38 ==> angle=90-38=52            
        def central_ring(self,maxi):
            """returns central ring as a list of dictionaries,
            specify the number of segs to map it out
            """
            #return points x and y of central ring
            segs_ring = []
            xpos,ypos=self.ellipse[0].x_center,self.ellipse[0].y_center
            radm,radn=self.ellipse[0].major,self.ellipse[0].minor
            an=self.ellipse[0].angle*numpy.pi/180.0
            co,si=numpy.cos(an),numpy.sin(an)
            the=numpy.linspace(0,2*numpy.pi,maxi)
            Xp=radm*numpy.cos(the)*co-si*radn*numpy.sin(the)+xpos
            Yp=radm*numpy.cos(the)*si+co*radn*numpy.sin(the)+ypos
            segs_ring.append( {'x':Xp, 'y':Yp} )
            return segs_ring
                                                                                                                                                                                                              
    def __init__(self,file_name,folder_,conf_size=0.0,conf_alpha=0.0):
        """Initializes global paramters"""
        self.gs = self.GalSpiral()
        self.cr = self.CentralRing()
        self.colors = [ 'r', 'c', 'g', 'm', 'b', 'y']
        self.nspirplot = 80                
        control = parse_com(file_name,folder_)    
        self.r0 = float(control['r0'])        
        self.v0t = float(control['v0t'])
        self.rot_curve=float(control['Rotation_Curve'])
        self.usun= float(control['usun'])
        self.vsun=float(control['vsun'])
        self.wsun=float(control['wsun'])
        self.fat = float(control['fat'])
        self.hz = float(control['hz'])
        self.halfmr=float(control['halfmr']) 
        self.v0r = float(control['v0r'])                                
        self.sigvr = float(control['sigvr'])
        self.sigva = float(control['sigva']) 
        self.v0z = float(control['v0z'])
        self.sigvz = float(control['sigvz'])
        self.v0r_us = float(control['v0r_us'])
        self.sigvr_us = float(control['sigvr_us'])
        self.v0t_vs = float(control['v0t_vs'])                                        
        self.sigva_vs = float(control['sigva_vs'])
        if(conf_alpha==0.0):
            self.index = float(control['index']) 
        else:
            self.index = conf_alpha
        self.min_l = float(control['min_l']) 
        self.max_l = float(control['max_l']) 
        self.k=(self.index+1)/(numpy.power(self.max_l,self.index+1)-numpy.power(self.min_l,self.index+1))        
        self.rtrunc = float(control['rtrunc']) 
        self.rmintrunc = float(control['rmintrunc'])         
        self.nspirseg = float(control['nspirseg']) 
        self.MMB_sigma= float(control['MMB_sigma'])
        self.Arecibo_sigma=float(control['Arecibo_sigma'])        
        if(conf_size==0.0):
            self.nobj = float(control['nobj'])
        else:
            self.nobj=conf_size
        self.Unc_confirmation=float(control['Unc_confirmation'])
        self.Push_betterunc=float(control['Push_betterunc'])
        self.Unc_confirmation_spec=float(control['Unc_confirmation_spec'])
        self.Unc_confirmation_percent=float(control['Unc_confirmation_percent'])
        self.Perfect_data=float(control['Perfect_data'])        
        self.Perfect_parallax=float(control['Perfect_parallax'])          
        self.Distance_spread=float(control['Distance_spread'])  
        self.Noerrorbutobs_para=float(control['Noerrorbutobs_para'])  
        self.Only_nice_sources=float(control['Only_nice_sources']) 
        self.Mark_trick=float(control['Mark_trick']) 
        self.root_surveys=control['Root_Surveys']