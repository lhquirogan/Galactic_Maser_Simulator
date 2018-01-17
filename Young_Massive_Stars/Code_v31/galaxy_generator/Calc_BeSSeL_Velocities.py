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


l=0.67
b=-0.04	
p=0.130
mua_cosd=-1.23
mud=-3.84
v_lsr=61.0

c1=sin()cos()-
c2=
k=astrometryToPhaseSpace(phi, theta, parallax, muphistar, mutheta, vrad)
k=astrometryToPhaseSpace(l, b, p, mua_cosd, mud, v_lsr)
