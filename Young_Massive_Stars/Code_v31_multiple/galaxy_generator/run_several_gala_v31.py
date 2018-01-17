import sys
from shutil import move
import os
import itertools
import shutil


kvar = raw_input("Do you want to make multiple file parameters (yes or no)?  ")
if (kvar=='yes' or kvar=='y'):
    def replace(file_path,file_path_old_, pattern0, pattern1, pattern2, pattern3, subst0, subst1, subst2,subst3):    
        with open(file_path,'w') as new_file:
            with open(file_path_old_) as old_file:
                for line in old_file:
                    if (line==pattern0):
                        new_file.write(line.replace(pattern0, subst0))
                    elif (line==pattern1):
                        new_file.write(line.replace(pattern1, subst1))
                    elif (line==pattern2):
                        new_file.write(line.replace(pattern2, subst2))
                    elif (line==pattern3):
                        new_file.write(line.replace(pattern3, subst3))                                        
                    else:
                        new_file.write(line)
        return
            
    def seek(filename,parame):
        coin_line=[]
        with open(filename) as f:
            for line in f:
                    if parame in line:
                    #if 'r0=' in line:
                        coin_line.append(line)
        return(coin_line[0])               
    
                                                    
    def add_equal(par):
        for i in range(len(par)):
            par[i]=par[i]+str('=')
        return(par)   
    generator_folder=os.getcwd()+'/'         
    para_folder=os.getcwd()+'/parameters/'
    try: 
        os.makedirs('parameters/')
    except OSError:
        if not os.path.isdir('parameters/'):
            raise                                                                                                                                                        
    varss = raw_input("Please enter the old parameter's file (should be in this directory): ")
    #varss='para_v31_BeSSeLmimic.txt'
    shutil.copy2(generator_folder+varss,para_folder+varss) 
    file_path_old=para_folder+varss
    var_para = raw_input("Parameters to change? (e.g r0, v0t, vsun, v0t_vs): ")
    parameterss=[x.strip() for x in var_para.split(',')]
    parametersss=add_equal(parameterss)
    val=[]
    for i in range(len(parametersss)):
        val_para= raw_input("Parameters values for %s" %parametersss[i])
        values=[x.strip() for x in val_para.split(',')]
        point=parametersss[i],values
        val.append(point)
    '''
    para_change=[]
    for i in range(len(val)):
        for j in range(len(val[i][1])):
            para_change.append(val[i][0]+val[i][1][j])
    '''
    parametersss_values=[]
    for i in range(len(val)):
        parametersss_values.append(val[i][1])
    z=list(itertools.product(parametersss_values[0],parametersss_values[1], parametersss_values[2],parametersss_values[3]))
    
    for m in range(len(z)):
        name_dummy_file='/para_v31_'+str(m)+'.txt'
        file_path=para_folder+name_dummy_file                
        pattern0=seek(file_path_old,parametersss[0])
        subst0=pattern0[0:3]+z[m][0]+pattern0[-63:]
        pattern1=seek(file_path_old,parametersss[1])
        subst1=pattern1[0:4]+z[m][1]+pattern1[-55:]
        pattern2=seek(file_path_old,parametersss[2])
        subst2=pattern2[0:5]+z[m][2]+pattern2[-40:]
        pattern3=seek(file_path_old,parametersss[3])
        subst3=pattern3[0:7]+z[m][3]+pattern3[-98:]
        replace(file_path,file_path_old, pattern0,pattern1,pattern2,pattern3,subst0,subst1,subst2,subst3)

    var2 = raw_input("How many galaxies do you want per set of parameters?: ")

    try:
        os.remove(para_folder+varss)
    except OSError:
        pass
            
    for m in range(len(z)):    
        print ('#################################################################')
        print('Galaxy set number %s of %s  "' % (m+1, len(z)))
        print ('#################################################################')
        var='/para_v31_'+str(m)+'.txt'
        sys.argv = ['GaMe_LHQN_v31.py', var]
        num_gala=int(var2)
        for i in range(num_gala):
            print ('#################################################################')
            print('Galaxy simulation number %s of %s  "' % (i+1, num_gala))
            print ('#################################################################')
            execfile('GaMe_LHQN_v31.py') 

    
else:
    var = raw_input("Please enter the parameter's file: ")
    sys.argv = ['GaMe_LHQN_v31.py', var]
    var2 = raw_input("How many galaxies do you want?: ")
    num_gala=int(var2)
    for i in range(num_gala):
        print ('#################################################################')
        print('Galaxy simulation number %s of %s  "' % (i+1, num_gala))
        print ('#################################################################')
        execfile('GaMe_LHQN_v31.py')    