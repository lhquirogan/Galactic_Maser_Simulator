import sys
var = raw_input("Please enter the parameter's file: ")
sys.argv = ['GaMe_LHQN_v31.py', var]
var2 = raw_input("How many galaxies do you want?: ")
num_gala=int(var2)
for i in range(num_gala):
    print ('#################################################################')
    print('Galaxy simulation number %s of %s  "' % (i+1, num_gala))
    print ('#################################################################')
    execfile('GaMe_LHQN_v31.py')    