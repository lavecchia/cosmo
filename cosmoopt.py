# Run optimization of radii and gamma factors to minimize differences between deltaGcalc and deltaGexp
# INPUT: 
#   > experimental deltaG
#   > template geometries in MOPAC input format, previusly optimized with the method to test
#   > initial radii values in a file

from commonfunctions import * # common function
from commonparameters import * # common parameters
import numpy as np
import scipy.optimize
import subprocess
from itertools import count
import ConfigParser # import setting file
from optparse import OptionParser
import sys, os, time
from subprocess import Popen, list2cmdline
import glob
import time #check time of calculation
import datetime
import shutil
import math
import ga_interface
from copy import deepcopy #to deepcopy dictionaries
from pylab import *
#~ import minimize_interface


###########
# PARSE
###########
parser = OptionParser()
parser.add_option("-o", "--outfile", dest="outfilename", default=DEFREPORFILE,
                  help="name of file to write REPORT", metavar="REPORT")
parser.add_option("-c", "--config", dest="configfile", default=CONFIG_FILENAME,
                  help="FILE to setting the program", metavar="CONFIG")
parser.add_option("-i", "--initial", dest="paramfile", default=DEFPARAMFILE,
                  help="FILE to read values of params", metavar="PARAMETERS")
parser.add_option("-r", "--ref", dest="reffilename", default=DEFREFFILE,
                  help="FILE to read values to comparate", metavar="FILE")
parser.add_option("-g", "--gasdir", dest="gasdir", default=DEFGASDIR,
                  help="GASDIR directory where take the gas calculation", metavar="GASDIR")
parser.add_option("-t", "--templatedir", dest="templatedir", default=DEFTEMPLATEDIR,
                  help="TEMPLATEDIR directory where take the template files", metavar="TEMPLATEDIR")
 
(options, args) = parser.parse_args()

outfilename = options.outfilename
paramfile = options.paramfile
reffilename = options.reffilename
gasdir = options.gasdir
templatedir = options.templatedir

cfg = ConfigParser.ConfigParser()  
if not cfg.read([options.configfile]):  
    print "No existe el archivo" 

###########
# VARS
###########
varlist = [] #store name of vars to then print
mcmarklist = [] #store the last mcmark to control the temperature in MC run
#limits
if cfg.has_option("limits", "radiilimit"):  
    radiilimit = float(cfg.get("limits", "radiilimit"))
else:  
    radiilimit = 0.4 #fix maximum displacement from initial values
varlist.append("radiilimit")

if cfg.has_option("limits", "radiirange"):  
    radiirange = float(cfg.get("limits", "radiirange"))
else:  
    radiirange = 0.2 #range of radii
varlist.append("radiirange")


if cfg.has_option("limits", "cosmoradiilimit"):  
    cosmoradiilimit = float(cfg.get("limits", "cosmoradiilimit"))
else:  
    cosmoradiilimit = 0.2 #fix maximum displset_FortiRM12007acement from initial values of cosmo radii. only for internal calculation of MOPAC
varlist.append("cosmoradiilimit")

if cfg.has_option("limits", "cosmoradiirange"):  
    cosmoradiirange = float(cfg.get("limits", "cosmoradiirange"))
else:  
    cosmoradiirange = 0.1 #range of cosmo radii
varlist.append("cosmoradiirange")

if cfg.has_option("limits", "gammalimit"):  
    gammalimit = float(cfg.get("limits", "gammalimit"))
else:  
    gammalimit = 0.099 #fix maximum displacement from initial values
varlist.append("gammalimit")

if cfg.has_option("limits", "gammarange"):  
    gammarange = float(cfg.get("limits", "gammarange"))
else:  
    gammarange = 0.05 #range of gamma
varlist.append("gammarange")

if cfg.has_option("limits", "rsolvlimit"):  
    rsolvlimit = float(cfg.get("limits", "rsolvlimit"))
else:  
    rsolvlimit = 0.2 #fix maximum displacement from initial values
varlist.append("rsolvlimit")

if cfg.has_option("limits", "rsolvrange"):  
    rsolvrange = float(cfg.get("limits", "rsolvrange"))
else:  
    rsolvrange = 0.1 #range of solvent radii
varlist.append("rsolvrange")

if cfg.has_option("limits", "eumbrallimit"):  
    eumbrallimit = float(cfg.get("limits", "eumbrallimit"))
else:  
    eumbrallimit = 0.1
varlist.append("eumbrallimit")

if cfg.has_option("limits", "eumbralrange"):  
    eumbralrange = float(cfg.get("limits", "eumbralrange"))
else:  
    eumbralrange = 0.05
varlist.append("eumbralrange")

if cfg.has_option("limits", "krange"):  
    krange = float(cfg.get("limits", "krange"))
else:  
    krange = 0.5
varlist.append("krange")

if cfg.has_option("limits", "klimit"):  
    klimit = float(cfg.get("limits", "klimit"))
else:  
    klimit = 1.0
varlist.append("klimit")

#method
if cfg.has_option("method", "extrakeys"):  
    extrakeys =cfg.get("method", "extrakeys")
else:  
    extrakeys = "PM6 PRECISE 1SCF"
varlist.append("extrakeys")

if cfg.has_option("method", "extrakeyssolv"):  
    extrakeyssolv = cfg.get("method", "extrakeyssolv")
else:  
    extrakeyssolv = "EPS=78 COSWRT"
varlist.append("extrakeyssolv")

if cfg.has_option("method", "templateext"):  
    templateext = cfg.get("method", "templateext")
else:  
    templateext = DEFTEMPLATEEXT
varlist.append("templateext")


#system
if cfg.has_option("system", "calculationtype"):  
    calculationtype = cfg.get("system", "calculationtype")
else:  
    calculationtype="mc"
varlist.append("calculationtype")


# genetic algorithm parameters
if cfg.has_option("system", "numbermembers"):
    numbermembers = int(cfg.get("system", "numbermembers"))
else:
    numbermembers = 20
varlist.append("numbermembers")
    
if cfg.has_option("system", "maxgen"):
    maxgen = int(cfg.get("system", "maxgen"))
else:
    maxgen = 100
varlist.append("maxgen")



if cfg.has_option("system", "onlyneutral"):  
    onlyneutral = cfg.getboolean("system", "onlyneutral")
else:  
    onlyneutral=True
varlist.append("onlyneutral")

if cfg.has_option("system", "cyclicoptimization"):  
    cyclicoptimization = cfg.getboolean("system", "cyclicoptimization")
else:  
    cyclicoptimization = True #sequence of radii and gamma optimization
varlist.append("cyclicoptimization")

if cfg.has_option("system", "maxiter"):  
    maxiter = int(cfg.get("system", "maxiter"))
else:  
    maxiter = 20000 #maximum number of iterations in the minimization
varlist.append("maxiter")

if cfg.has_option("system", "nptype"):  
    nptype = cfg.get("system", "nptype")
else:  
    nptype = "claverie" #"electronmod" #"claverie"
varlist.append("nptype")

if cfg.has_option("system", "excludeatomlist"):  
    excludeatomlist = cfg.get("system", "excludeatomlist").split()
    varlist.append("excludeatomlist")

if cfg.has_option("system", "temperature"):  
    temperature = float(cfg.get("system", "temperature"))
else:  
    temperature = 0.10 #temperature of Monte Carlo simulation
varlist.append("temperature")

if cfg.has_option("system", "rangeslope"):  
    rangeslope = float(cfg.get("system", "rangeslope"))
else:  
    rangeslope = 0.995 #velocity of decrement of range that take values
varlist.append("rangeslope")

if cfg.has_option("system", "fixlist"):  
    fixlist = cfg.get("system", "fixlist").split()
    varlist.append("fixlist")

#set the formats of optimization, step1 p1 p2 p3 pn N, where p1...pn are paremeters, and N number of cycles.
stepslist = []
try:
    for n in range(1,1000): #max of 1000 optimization step formats
        if cfg.has_option("system", "step" + str(n)):
            stepall = cfg.get("system", "step"+str(n)).split()
            stepcycles = int(stepall.pop())
            stepslist.append([stepall,stepcycles]) #store: parameters, number of cycles
except: 
    if len(stepslist)==0:
        print "ERROR: You must define step settings in config.in, eg. step1 = p1 p2 p3 100 where p1, p2, p3 are parameters to optimize, and 100 the number of cycles per step."
        exit()
    else:
        varlist.append("stepslist")
    

#~ limitdic={"radii":radiilimit, "gamma":gammalimit, "rsolv":rsolvlimit} #limits determine fix extremes of values center in "initial" values.
rangesdic = {"radii":radiirange, "gamma":gammarange, "rsolv":rsolvrange, "eumbral":eumbralrange, "cosmoradii":cosmoradiirange, "k":krange} #moveable ranges, centers are "current" values in MC run



###########
# PROGRAM
###########

# NOTATION SUFFIX: 
#         0: initial
#      test: value in current step
#      best: the lowest totalerror step
#   current: store last trajectory values of MC run


# to debbug, copy source code to directory
os.system("cp " + __file__ + " source_code.py")

#initial time
start_time = time.time()

#initialize report file
#outfile = open(outfilename, "w",0)
outfile = open(outfilename, "w")


#read initial parameter values to adjust
param0dic = parameters_read(paramfile)

#read referential values to compare
datadic = exp_read(reffilename)
if onlyneutral==True:
    datadic = {your_key: datadic[your_key] for your_key in datadic if "anion" not in your_key }
    datadic = {your_key: datadic[your_key] for your_key in datadic if "protonated" not in your_key }
    

#erase compounds with a particular element of dgrefdic and template lists
try:
    excludeatomlist
except NameError:
    pass
else:
    newdatadic = {}
    for key, values in datadic.iteritems():
        if check_excludeatom(templatedir + "/" + key + templateext,excludeatomlist)==False:
            newdatadic[key]=datadic[key]
        else:
            print "Exclude:"  + key
    datadic = newdatadic

fixlimitdic={}
minlimitdic={}
maxlimitdic={}
for key, value in param0dic.iteritems():
    if "r@" in key:
        limitvalue = radiilimit
    elif "g@" in key:
        limitvalue = gammalimit
    elif "rsolv@" in key:
        limitvalue = rsolvlimit
    elif "e@" in key: #number of electrons umbral
        limitvalue = eumbrallimit
    elif "rc@" in key: #cosmo radius
        limitvalue = cosmoradiilimit
    elif "k@" in key: #constant multiplica radii
        limitvalue = klimit
    else:
        limitvalue = 0
    fixlimitdic[key]=[value-limitvalue, value+limitvalue]
    minlimitdic[key]=value-limitvalue
    maxlimitdic[key]=value+limitvalue

#make GAS phase calculation
for key, value in datadic.iteritems():
    datadic[key]["template"]=templatedir + "/" + key + templateext
    inputfile_make(value["template"],extrakeys)

#run GAS phase calculation
gasinputlist = list_files(gasdir)
commands = []
for gasinput in gasinputlist:
    commands.append([MOPACPATH, gasinput])
exec_commands(commands)

#extract GAS Heat of Formation
for key, values in datadic.iteritems():
    datadic[key]["hofgas"],datadic[key]["cosmoarea"] = mopacout_read(gasdir + "/" + key + templateext.replace(".mop","_gas.out")) #extract HOF

#head of report file
outfile.write(str(datetime.datetime.now())+"\n")

#write settings in report file
for name in varlist:
    try:
        outfile.write(name + "\t\t=\t" + str(eval(name)) + "\n")
    except:
        pass

outfile.write("\nKEYS:" + extrakeys + " SOLVENT:" + extrakeyssolv +"\n")
outfile.write("REFERENCES:\n")
for key, value in datadic.iteritems():
    outfile.write(key + ":" + str(value["dgexp"]) + "\n")
outfile.write(20*"=" + "\n")

if cyclicoptimization == True:
    period = 0
    maxcyclebystep = stepslist[period][1] # number of cycle by step #add to fix GA
    freeparamlist = stepslist[period][0] #add to fix GA
else:
    freeparamlist=[]
    for value,key in param0dic.iteritems():
        freeparamlist.append(value)


rsolvtest = param0dic["rsolv@"] #take rsolv parameter defined in input parameter file

if nptype=="claverie" or nptype=="claverietype":    
    yrel0 = 4*PI*NS*(rsolvtest*1.0E-10)*(rsolvtest*1.0E-10)*(rsolvtest*1.0E-10)/3.0 # ratio to use in claverie cavitation term
else:
    yrel0 = 0

#calculate initial value
mintotalerror, mae, rmse, bias, r2, slope, intercept, datadic = calc_error("0000", param0dic, datadic, extrakeys, extrakeyssolv + " RSOLV=% .3f" % (rsolvtest), outfile, nptype, yrel0)
outfile.write("\n")
# write summary file with initial values for first step
print_summary(mintotalerror, mae, rmse, bias, r2, slope, intercept,0,datadic,nptype)
paramtestdic = param0dic #copy initial parameter values to a new dictonary (test) to modified
beststep = [mintotalerror, paramtestdic, datadic] # first element: minimum error, second element: list of values of best step
currentstep = [mintotalerror, paramtestdic, datadic] # first element: minimum error, second element: list of values of current step


#Monte Carlo run
if calculationtype == 'mc':
    acumcycle = 0 # accumulate the number of cycles
    for ncycle in range(1,maxiter):
            mcmark = "\n"
            
            #number of step, to generate name of files
            numberstep = next(step)
            
            if cyclicoptimization == True:
                maxcyclebystep = stepslist[period][1] # number of cycle by step
                freeparamlist = stepslist[period][0]
                if acumcycle >= maxcyclebystep:
                    acumcycle = 0 # reset count
                    if period == len(stepslist)-1:
                        period = 0 #restart to first step
                    else:
                        period += 1 #next step
                    paramtestdic = beststep[1]
                    datadic = beststep[2]
                    print period
            
            
            #check restrictions, if not assign parameters again
            tag1sttry = True # the first time into while loop
            while (check_restrictions(paramtestdic,fixlimitdic)==0) or (tag1sttry == True):
                tag1sttry = False # first time into while loop
                # generate new values of parameters with a Gaussian distribution probability from current values
                paramtestdic = modified_values(currentstep[1], freeparamlist, rangesdic)
                
                
                # multiplied initial radii with a constant k, which is optimizing.
                try:
                    if paramtestdic["k@"]:
                        for key, value in paramtestdic.iteritems():
                            if "rc@" in key:
                                paramtestdic[key] = param0dic[key] * paramtestdic["k@"] 
                except:
                    pass
                
                #fix values
                try:
                    for key in fixlist:
                        paramtestdic[key]=param0dic[key]
                except:
                    pass

            rsolvtest = paramtestdic["rsolv@"] #take rsolv parameter 
            if nptype=="claverie" or nptype=="claverietype":    
                yrel = 4*PI*NS*(rsolvtest*1.0E-10)*(rsolvtest*1.0E-10)*(rsolvtest*1.0E-10)/3.0 # ratio to use in claverie cavitation term
            else:
                yrel = 0

            #=================================================================================================
            # check if the new value of error function is below respect to the best value store in beststep[0]

            totalerror, mae, rmse, bias, r2, slope, intercept, datadic = calc_error(numberstep, paramtestdic, datadic, extrakeys, extrakeyssolv + " RSOLV=% .3f" % (rsolvtest), outfile, nptype, yrel)
            
            if totalerror < mintotalerror:
                mintotalerror = totalerror
                # write summary file with the low totalerror step
                print_summary(mintotalerror, mae, rmse, bias, r2, slope, intercept, ncycle,datadic,nptype)
        
                # write parameter values
                paramout = open("list_param.out","w")
                for key, value in sorted(paramtestdic.iteritems()):
                    paramout.write("%-3s \t %.3f\n" % (key,value))
                paramout.close()
                            
                # ====================================================================
                
            #~ else:    
                #~ shutil.rmtree(numberstep)
            shutil.rmtree(numberstep) # erase current directory
                
            prob = np.random.rand()
            

            if totalerror < currentstep[0] or prob < math.exp(-(totalerror-currentstep[0])/temperature):
                mcmark = " @PROB\n"# mark that indicates this is a step selected by Monte Carlo because decent in energy
                # "current" vars store last trajectory values of MC run.
                currentstep = [totalerror,paramtestdic,datadic]
                            
                if totalerror < beststep[0]:
                    mcmark = " @DESC\n" # change mark that indicates this is a step selected by Monte Carlo because decent in energy
                    beststep = [totalerror,paramtestdic,datadic]

            outfile.write(mcmark) #end of line with a mark, @DESC: recenter because below value of function, @PROB: select from probability

            #scale the temperature as function of times of @PROB mcmark
            #mcmarklist store the n last mcmark
            temperature, mcmarklist = temperature_control(temperature,mcmarklist,mcmark) 
            print temperature 
            print mcmarklist
            
            print "PASO:" + str(ncycle)
            acumcycle += 1
            
#Genetic Algorithm run
elif calculationtype == 'ga':
    acumcycle = 0 # accumulate the number of cycles
    
    memberlist = [] #list that store members (list of paramtestdic)
    #building parents
    for n in range(0,numbermembers):

                        
        #check restrictions, if not assign parameters again
        tag1sttry = True # the first time into while loop
        while (check_restrictions(paramtestdic,fixlimitdic)==0) or (tag1sttry == True):
            tag1sttry = False # first time into while loop
            # generate new values of parameters with a distribution probability from current values
            paramtestdic = modified_values(currentstep[1], freeparamlist, rangesdic)
            
            # multiplied initial radii with a constant k, which is optimizing.
            try:
                if paramtestdic["k@"]:
                    for key, value in paramtestdic.iteritems():
                        if "rc@" in key:
                            paramtestdic[key] = param0dic[key] * paramtestdic["k@"] 
            except:
                pass
            
            #fix values
            try:
                for key in fixlist:
                    paramtestdic[key]=param0dic[key]
            except:
                pass
                
        memberlist.append(ga_interface.memberobj(paramtestdic))


     
    #generations
    for ncycle in range(1,maxgen):
        
        if cyclicoptimization == True:
            maxcyclebystep = stepslist[period][1] # number of cycle by step
            freeparamlist = stepslist[period][0]
            if acumcycle >= maxcyclebystep:
                acumcycle = 0 # reset count
                if period == len(stepslist)-1:
                    period = 0 #restart to first step
                else:
                    period += 1 #next step
                paramtestdic = beststep[1]
                datadic = beststep[2]
                print period
        
        
        #number of step, to generate name of files
        numberstep = next(step)
        
        #fitness calculation
        fitnesslist = [] #store fitness
        for member in memberlist:
            paramtestdic = member.paramdic
            rsolvtest = paramtestdic["rsolv@"] #take rsolv parameter 
            if nptype=="claverie" or nptype=="claverietype":    
                yrel = 4*PI*NS*(rsolvtest*1.0E-10)*(rsolvtest*1.0E-10)*(rsolvtest*1.0E-10)/3.0 # ratio to use in claverie cavitation term
            else:
                yrel = 0

            #=================================================================================================
            # check if the new value of error function is below respect to the best value store in beststep[0]

            if member.fitness == None:
                totalerror, mae, rmse, bias, r2, slope, intercept, datadic = calc_error(numberstep, paramtestdic, datadic, extrakeys, extrakeyssolv + " RSOLV=% .3f" % (rsolvtest), outfile, nptype, yrel)
            
                if totalerror < mintotalerror:
                    mintotalerror = totalerror
                    # write summary file with the low totalerror step
                    print_summary(mintotalerror, mae, rmse, bias, r2, slope, intercept, ncycle,datadic,nptype)
            
                    # write parameter values
                    paramout = open("list_param.out","w")
                    for key, value in sorted(paramtestdic.iteritems()):
                        paramout.write("%-3s \t %.3f\n" % (key,value))
                    paramout.close()
                                
                    # ====================================================================
                    
                #~ else:    
                    #~ shutil.rmtree(numberstep)
                shutil.rmtree(numberstep) # erase current directory
                fitness = 1.0/(totalerror*totalerror*totalerror*totalerror)
                fitnesslist.append(fitness)
                member.set_fitness(fitness)
            else:
                fitnesslist.append(member.fitness)
                #add to fix write all individues in report file
                prevline = "%-5s Err: %3.4f MAE: %3.4f RMSE: %3.4f BIAS: %3.4f R2: %1.5f " % (numberstep, totalerror, mae, rmse, bias, r2)
                outfile.write(prevline + print_param(paramtestdic))
        
        

        outfile.write("\n")
        genetic = ga_interface.SetsGeneration(memberlist, minlimitdic, maxlimitdic)
        memberlist = genetic.next()
  
        outfile.write("\n") #end of line with a mark

        acumcycle += 1
        print "PASO: %i %.3f"%(ncycle,np.mean(fitnesslist))
        salida = ""
        for fitn in sorted(fitnesslist):
            salida += " %.3f"%(fitn)
        print salida

# Minimize Algorithm run
elif calculationtype == 'minimize':
    from lmfit import minimize, Parameters, Parameter, report_fit, Minimizer
    acumcycle = 0
    limitcycle = 2
    cycleswithoutdown = 0

    
    def fcn2min(params, extrakeys, extrakeyssolv, outfile, nptype, datacompoundnamelist):
        global mintotalerror
        global acumcycle
        global datadic
        global cycleswithoutdown
        global limitcycle

        #rebuild paramdic
        paramdic={}
        for key, values in params.iteritems():
            key = key.replace("zzz","@") #change to original key formats
            key = key.replace("xxx",".")
            paramdic[key]=values.value
        
        rsolvtest = paramdic["rsolv@"]
                
        if nptype=="claverie" or nptype=="claverietype":
            yrel = 4*PI*NS*(rsolvtest*1.0E-10)*(rsolvtest*1.0E-10)*(rsolvtest*1.0E-10)/3.0 # ratio to use in claverie cavitation term
        else:
            yrel = 0
            
        totalerror, mae, rmse, bias, r2, slope, intercept, datadic = calc_error("0000", paramdic, datadic, extrakeys, extrakeyssolv + " RSOLV=% .3f" % (rsolvtest), outfile, nptype, yrel)
        if totalerror < mintotalerror:
            cycleswithoutdown = 0
            mintotalerror = totalerror
            # write summary file with the low totalerror step
            print_summary(mintotalerror, mae, rmse, bias, r2, slope, intercept, acumcycle,datadic,nptype)
            # write parameter values
            paramout = open("list_param.out","w")
            for key, value in sorted(paramdic.iteritems()):
                print "%-3s \t %.5f\n" % (key,value)
                paramout.write("%-3s \t %.5f\n" % (key,value))
            paramout.close()
        else:
            cycleswithoutdown += 1
            if cycleswithoutdown > limitcycle:
                return #exit function with error
            #~ shutil.rmtree(numberstep)
        shutil.rmtree("0000") # erase current directory
        print "PASO %i: %f"%(acumcycle,totalerror)
        acumcycle += 1
        errors = []
        for i in range(0,len(datacompoundnamelist)):
            errors.append(datadic[datacompoundnamelist[i]]['dgexp']-datadic[datacompoundnamelist[i]]['dgcalc'])
        #~ return errors
        
        return totalerror
        
    params = Parameters()
    for ikey, ivalue in param0dic.iteritems():
        maxlimit = maxlimitdic[ikey]
        minlimit = minlimitdic[ikey]
        ikey = ikey.replace("@","zzz") #replace @ with zzz because this character is not supported by lmfit library
        ikey = ikey.replace(".","xxx") #replace . with xxx because this character is not supported by lmfit library
        #~ params.add(ikey, value=ivalue)
        #~ if "rczzz" in ikey:
        #~ params.add(ikey, value=ivalue, min=minlimit, max=maxlimit)
        #~ else:
        #~ params.add(ikey, ivalue, False)
        params.add(ikey, value=ivalue, min=minlimit, max=maxlimit)
        
    print params
    #experimental data
    datacompoundnamelist = []
    for ikey, ivalue in datadic.iteritems():
        datacompoundnamelist.append(ikey) #to convert in error function to a dictionary
    
    #~ extra_kwargs={}
    #~ extra_kwargs['epsfcn'] = 0.5
    #~ result = minimize(fcn2min, params, args=(extrakeys, extrakeyssolv, rsolvtest, outfile, nptype, datacompoundnamelist), method='lbfgs', **extra_kwargs)
    #~ extra_kwargs={}
    #~ extra_kwargs['T'] = 300.0
    #~ extra_kwargs['stepsize']=0.1
    #~ result = minimize(fcn2min, params, args=(extrakeys, extrakeyssolv, rsolvtest, outfile, nptype, datacompoundnamelist), method='basinhopping', **extra_kwargs)
    
    options = {}
    options["maxiter"] = 3
    
    kws = {}
    kws["options"]=options
    print kws
    
    #~ myfit = Minimizer(fcn2min, params, fcn_args=(extrakeys, extrakeyssolv, outfile, nptype, datacompoundnamelist), maxfun=5, **options)
    #~ myfit.prepare_fit()
    #~ init = my_residual(p_fit, x)
    #~ pylab.plot(x, init, 'b--')

    #~ myfit.fmin()
    
    
    try:
        result = minimize(fcn2min, params, args=(extrakeys, extrakeyssolv, outfile, nptype, datacompoundnamelist), method='powell')
    except:
        pass
    #~ result = minimize(fcn2min, params, args=(extrakeys, extrakeyssolv, outfile, nptype, datacompoundnamelist), method='powell', **kws)
    #~ #calculate final result
    #~ final = data + result.residual
    # write error report
    report_fit(params)
    #~ # try to plot results
    #~ try:
    #~ import pylab
    #~ pylab.plot(x, data, 'k+')
    #~ pylab.plot(x, final, 'r')
    #~ pylab.show()
    #~ except:
    #~ pass

# genetic algorithm + minimization run
elif calculationtype == 'gamin':
    from lmfit import minimize, Parameters, Parameter, report_fit
    maxcycleminimization = 30 #set minimization period
    maxcyclega = 20 #set genetic algorithm period

    
    #function to minimize
    def fcn2min(params, extrakeys, extrakeyssolv, outfile, nptype, datacompoundnamelist):
        global mintotalerror #the most minimal error
        global acumcycle #general count of cycles
        global datadic #store experimental, calculated and contributions to solvation energies for all compounds
        global paramtestdic # params to test
        global cyclemin # cycle number of the period of minimization
        global mintotalerrormember # minimal error of a member
        global minparamtestdic #store params with mintotalerrormemeber
        #rebuild paramtestdic
        paramtestdic={}
        for key, values in params.iteritems():
            key = key.replace("zzz","@") #change to original key formats
            key = key.replace("xxx",".")
            paramtestdic[key]=values.value
        rsolvtest = paramtestdic["rsolv@"]
        if nptype=="claverie" or nptype=="claverietype":
            yrel = 4*PI*NS*(rsolvtest*1.0E-10)*(rsolvtest*1.0E-10)*(rsolvtest*1.0E-10)/3.0 # ratio to use in claverie cavitation term
        else:
            yrel = 0
        totalerror, mae, rmse, bias, r2, slope, intercept, datadic = calc_error("0000", paramtestdic, datadic, extrakeys, extrakeyssolv + " RSOLV=% .3f" % (rsolvtest), outfile, nptype, yrel)
        if totalerror < mintotalerror:
            mintotalerror = totalerror
            # write summary file with the low totalerror step
            print_summary(mintotalerror, mae, rmse, bias, r2, slope, intercept, acumcycle,datadic,nptype)
            # write parameter values
            paramout = open("list_param.out","w")
            for key, value in sorted(paramtestdic.iteritems()):
                paramout.write("%-3s \t %.5f\n" % (key,value))
            paramout.close()
        
        if totalerror < mintotalerrormember:
            minparamtestdic = deepcopy(paramtestdic)
        shutil.rmtree("0000") # erase current directory
        print "PASO %i: %f"%(cyclemin,totalerror)        
        
        if cyclemin >= maxcycleminimization:
            return
        cyclemin += 1
        return totalerror
        
    
    memberlist = [] #list that store members (list of paramtestdic)
    #building parents
    for n in range(0,numbermembers):
        #check restrictions, if not assign parameters again
        tag1sttry = True # the first time into while loop
        while (check_restrictions(paramtestdic,fixlimitdic)==0) or (tag1sttry == True):
            tag1sttry = False # first time into while loop
            # generate new values of parameters with a distribution probability from current values
            paramtestdic = modified_values(currentstep[1], freeparamlist, rangesdic)
            
            # multiplied initial radii with a constant k, which is optimizing.
            try:
                if paramtestdic["k@"]:
                    for key, value in paramtestdic.iteritems():
                        if "rc@" in key:
                            paramtestdic[key] = param0dic[key] * paramtestdic["k@"] 
            except:
                pass
            
            #fix values
            try:
                for key in fixlist:
                    paramtestdic[key]=param0dic[key]
            except:
                pass
                
        memberlist.append(ga_interface.memberobj(paramtestdic))

    acumcycle = 0 # accumulate the number of cycles
    cyclega = 0 # accumulate the number of cycles of GA period
    
    #generations
    for ncycle in range(1,maxgen):
        
        #number of step, to generate name of files
        numberstep = next(step)
        
        #fitness calculation
        fitnesslist = [] #store fitness
        for member in memberlist:
            paramtestdic = deepcopy(member.paramdic)
            minparamtestdic = deepcopy(member.paramdic)
            rsolvtest = paramtestdic["rsolv@"] #take rsolv parameter 
            if nptype=="claverie" or nptype=="claverietype":    
                yrel = 4*PI*NS*(rsolvtest*1.0E-10)*(rsolvtest*1.0E-10)*(rsolvtest*1.0E-10)/3.0 # ratio to use in claverie cavitation term
            else:
                yrel = 0

            if (cyclega == maxcyclega) or (acumcycle == maxgen - 1) :
                # do minimization
                print "minimization..."
                cyclemin = 0

                params = Parameters()
                for ikey, ivalue in paramtestdic.iteritems():
                    maxlimit = maxlimitdic[ikey]
                    minlimit = minlimitdic[ikey]
                    ikey = ikey.replace("@","zzz") #replace @ with zzz because this character is not supported by lmfit library
                    ikey = ikey.replace(".","xxx") #replace . with xxx because this character is not supported by lmfit library
                    #~ params.add(ikey, value=ivalue)
                    #~ if "rczzz" in ikey:
                    #~ params.add(ikey, value=ivalue, min=minlimit, max=maxlimit)
                    #~ else:
                    #~ params.add(ikey, ivalue, False)
                    params.add(ikey, value=ivalue, min=minlimit, max=maxlimit)
                    
                #experimental data
                datacompoundnamelist = []
                for ikey, ivalue in datadic.iteritems():
                    datacompoundnamelist.append(ikey) #to convert in error function to a dictionary
                
                try:
                    mintotalerrormember = 1 / (member.fitness * member.fitness * member.fitness * member.fitness)
                except:
                    mintotalerrormember = 1000.0
      
                try:
                    result = minimize(fcn2min, params, args=(extrakeys, extrakeyssolv, outfile, nptype, datacompoundnamelist), method='powell')
                except:
                    pass
                
               
                member.paramdic = deepcopy(minparamtestdic)
                member.fitness = None

            if member.fitness == None:
                totalerror, mae, rmse, bias, r2, slope, intercept, datadic = calc_error(numberstep, paramtestdic, datadic, extrakeys, extrakeyssolv + " RSOLV=% .3f" % (rsolvtest), outfile, nptype, yrel)
            
                if totalerror < mintotalerror:
                    mintotalerror = totalerror
                    # write summary file with the low totalerror step
                    print_summary(mintotalerror, mae, rmse, bias, r2, slope, intercept, ncycle,datadic,nptype)
                    
                    # write parameter values
                    paramout = open("list_param.out","w")
                    for key, value in sorted(paramtestdic.iteritems()):
                        paramout.write("%-3s \t %.3f\n" % (key,value))
                    paramout.close()
                                
                shutil.rmtree(numberstep) # erase current directory
                fitness = 1.0/(totalerror*totalerror*totalerror*totalerror)
                fitnesslist.append(fitness)
                member.set_fitness(fitness)
            else:
                fitnesslist.append(member.fitness)
                #add to fix write all individues in report file
                prevline = "%-5s Err: %3.4f MAE: %3.4f RMSE: %3.4f BIAS: %3.4f R2: %1.5f " % (numberstep, totalerror, mae, rmse, bias, r2)
                outfile.write(prevline + print_param(paramtestdic))
        
        if cyclega >= maxcyclega:
            cyclega = 0

        outfile.write("\n")
        genetic = ga_interface.SetsGeneration(memberlist, minlimitdic, maxlimitdic)
        memberlist = genetic.next()
  
        outfile.write("\n") #end of line with a mark
        cyclega = cyclega + 1
        acumcycle += 1
        print "PASO: %i %.3f"%(ncycle,np.mean(fitnesslist))
        salida = ""
        for fitn in sorted(fitnesslist):
            salida += " %.3f"%(fitn)
        print salida

# static to manual test  
elif calculationtype == 'statics':

    print "MAE: %2.3f RMSE: %2.3f BIAS: %2.3f R2: %1.5f Slope: %2.3e Intercept: %2.3e" % (mae, rmse, bias, r2, slope, intercept)
    print "symbol rc@ g@ MAE RMSE BIAS slope intercept R2 eslope eintercept eR2"
    
    for elementsymbol in ["C", "O", "N", "H", "F", "Cl", "Br", "P", "S"]:
        try:
            mae, rmse, bias, r2, slope, intercept, errorr2, errorslope, errorintercept = calc_staticsbyelement("0000", datadic, elementsymbol)
            radio = "rc@" + elementsymbol + ".a"
            gamma = "g@" + elementsymbol + ".a"
            print "<%s> %1.3f %1.3f %3.3f %3.3f %3.3f %3.3f %3.3f %1.4f %3.3f %3.3f %1.4f" % (elementsymbol, param0dic[radio], param0dic[gamma], mae, rmse, bias, slope, intercept, r2, errorslope, errorintercept, errorr2)

            #~ print "<%s> %1.3f %1.3f MAE: %3.3f RMSE: %3.3f BIAS: %3.3f  CALCvsEXP: %3.3f x + %3.3f, R2: %1.4f ERRORvsELECTRON: %3.3f x + %3.3f R2 %1.4f" % (elementsymbol, param0dic[radio], param0dic[gamma], mae, rmse, bias, slope, intercept, r2, errorslope, errorintercept, errorr2)
        except:
            print "<%s> sin compuestos" %(elementsymbol)
                
elif calculationtype == 'minstatics':
    def make_modification(value0,rangevalue,stepsize):
        currentvalue = value0 - rangevalue
        valuelist = []
        while currentvalue < value0 + rangevalue:
            valuelist.append([currentvalue,None])
            currentvalue += stepsize
        return valuelist
    
    
    from lmfit import minimize, Parameters, Parameter, report_fit
    maxcycleminimization = 30 #set minimization period

    
    #function to minimize
    def fcn2min(params, extrakeys, extrakeyssolv, outfile, nptype, datacompoundnamelist):
        global mintotalerror #the most minimal error
        global ncycle #general count of cycles
        global datadic #store experimental, calculated and contributions to solvation energies for all compounds
        global paramtestdic # params to test
        global cyclemin # cycle number of the period of minimization
        #rebuild paramtestdic
        paramtestdic={}
        for key, values in params.iteritems():
            key = key.replace("zzz","@") #change to original key formats
            key = key.replace("xxx",".")
            paramtestdic[key]=values.value
        rsolvtest = paramtestdic["rsolv@"]
        if nptype=="claverie" or nptype=="claverietype":
            yrel = 4*PI*NS*(rsolvtest*1.0E-10)*(rsolvtest*1.0E-10)*(rsolvtest*1.0E-10)/3.0 # ratio to use in claverie cavitation term
        else:
            yrel = 0
        totalerror, mae, rmse, bias, r2, slope, intercept, datadic = calc_error("0000", paramtestdic, datadic, extrakeys, extrakeyssolv + " RSOLV=% .3f" % (rsolvtest), outfile, nptype, yrel)
        if totalerror < mintotalerror:
            mintotalerror = totalerror
            # write summary file with the low totalerror step
            print_summary(mintotalerror, mae, rmse, bias, r2, slope, intercept, ncycle,datadic,nptype)
            # write parameter values
            paramout = open("list_param.out","w")
            for key, value in sorted(paramtestdic.iteritems()):
                paramout.write("%-3s \t %.5f\n" % (key,value))
            paramout.close()

        print "PASO %i: %f"%(cyclemin,totalerror)        
        
        if cyclemin >= maxcycleminimization:
            cyclemin = 0
            return
        cyclemin += 1
        return totalerror
        
           

    
    
    def run_minstatics(paramtestdic):
        global datadic
        staticlist = []
        print "\n symbol rc@ g@ MAE RMSE BIAS slope intercept R2 eslope eintercept eR2"
        for elementsymbol in ["C", "O", "N", "H", "F", "Cl", "Br", "P", "S"]:
            #~ try:
            mae, rmse, bias, r2, slope, intercept, errorr2, errorslope, errorintercept = calc_staticsbyelement("0000", datadic, elementsymbol)
            radio = "rc@" + elementsymbol + ".a"
            gamma = "g@" + elementsymbol + ".a"
            print "<%s> %1.4f %1.4f %3.3f %3.3f %3.3f %3.3f %3.3f %1.4f %3.3f %3.3f %1.4f" % (elementsymbol, paramtestdic[radio], paramtestdic[gamma], mae, rmse, bias, slope, intercept, r2, errorslope, errorintercept, errorr2)
            staticlist.append([mae,elementsymbol, paramtestdic[radio], paramtestdic[gamma],rmse, bias, slope, intercept, r2, errorslope, errorintercept, errorr2])
        #~ except:
            #~ print "<%s> sin compuestos" %(elementsymbol)
        newstaticlist = sorted(staticlist, key=lambda lista: lista[0], reverse=True)
        return newstaticlist
    
    
    # select worst element (with the highest individual MAE contribution)
    staticlist = run_minstatics(paramtestdic)
    print staticlist
    elementsymbol = staticlist[0][1]
    minmaebyelement = staticlist[0][0]
    
    initialradiistepsize = 0.08
    initialgammastepsize = 0.02
    
    radiistepsize = initialradiistepsize
    gammastepsize = initialgammastepsize
    
    initialradiirange = radiirange
    initialgammarange = gammarange    
    
    for ncycle in range(0,maxiter):
        
        radioselect = "rc@" + elementsymbol + ".a"
        gammaselect = "g@" + elementsymbol + ".a"
        radiilist = make_modification(paramtestdic[radioselect],radiirange,radiistepsize)
        gammalist = make_modification(paramtestdic[gammaselect],gammarange,gammastepsize)
        
        print "SELECTED %s : number test %i, radii step %f , gamma step %f" % (elementsymbol, len(radiilist)*len(gammalist), radiistepsize, gammastepsize)
        
        for radio in radiilist:
            for gamma in gammalist:
                paramtestdic[radioselect] = radio[0]
                paramtestdic[gammaselect] = gamma[0]
                rsolvtest = paramtestdic["rsolv@"] #take rsolv parameter 
                if nptype=="claverie" or nptype=="claverietype":    
                    yrel = 4*PI*NS*(rsolvtest*1.0E-10)*(rsolvtest*1.0E-10)*(rsolvtest*1.0E-10)/3.0 # ratio to use in claverie cavitation term
                else:
                    yrel = 0
                totalerror, mae, rmse, bias, r2, slope, intercept, datadic = calc_error("0000", paramtestdic, datadic, extrakeys, extrakeyssolv + " RSOLV=% .3f" % (rsolvtest), outfile, nptype, yrel)
                staticlist = run_minstatics(paramtestdic)
                
                if totalerror < mintotalerror:
                    mintotalerror = totalerror
                    # write summary file with the low totalerror step
                    print_summary(mintotalerror, mae, rmse, bias, r2, slope, intercept, ncycle,datadic,nptype)
                    
                    # write parameter values
                    paramout = open("list_param.out","w")
                    for key, value in sorted(paramtestdic.iteritems()):
                        paramout.write("%-3s \t %.4f\n" % (key,value))
                    paramout.close()
                
                for staticbyelement in staticlist:
                    if elementsymbol in staticbyelement[1]:
                        maebyelement = staticbyelement[0]
                        if maebyelement < minmaebyelement:
                            minmaebyelement = maebyelement
                            minradio = radio[0]
                            mingamma = gamma[0]
                            nextelementsymbol = staticlist[0][1]
                            nextmae = staticlist[0][0]
                            minstaticlist = staticlist
                            
                            
        paramtestdic[radioselect] = minradio
        paramtestdic[gammaselect] = mingamma
                            
        if minmaebyelement < nextelementsymbol:
            if nextelementsymbol == elementsymbol:
                
                #minimize
                params = Parameters()
                for ikey, ivalue in paramtestdic.iteritems():
                    maxlimit = maxlimitdic[ikey]
                    minlimit = minlimitdic[ikey]
                    ikey = ikey.replace("@","zzz") #replace @ with zzz because this character is not supported by lmfit library
                    ikey = ikey.replace(".","xxx") #replace . with xxx because this character is not supported by lmfit library
                    #~ params.add(ikey, value=ivalue)
                    #~ if "rczzz" in ikey:
                    #~ params.add(ikey, value=ivalue, min=minlimit, max=maxlimit)
                    #~ else:
                    #~ params.add(ikey, ivalue, False)
                    params.add(ikey, value=ivalue, min=minlimit, max=maxlimit)
                    
                print params
                #experimental data
                datacompoundnamelist = []
                for ikey, ivalue in datadic.iteritems():
                    datacompoundnamelist.append(ikey) #to convert in error function to a dictionary

                
                cyclemin = 0
                try:
                    result = minimize(fcn2min, params, args=(extrakeys, extrakeyssolv, outfile, nptype, datacompoundnamelist), method='powell')
                except:
                    pass
                
                staticlist = run_minstatics(paramtestdic)
                for staticbyelement in staticlist:
                    if elementsymbol in staticbyelement[1]:
                        maebyelement = staticbyelement[0]
                        if maebyelement < minmaebyelement:
                            minmaebyelement = maebyelement
                            minradio = radio[0]
                            mingamma = gamma[0]
                            nextelementsymbol = staticlist[0][1]
                            nextmae = staticlist[0][0]
                            minstaticlist = staticlist

                #~ gammastepsize = gammastepsize * 0.8
                #~ radiistepsize = radiistepsize * 0.8
                #~ radiirange = radiirange * 0.8
                #~ gammarange = gammarange * 0.8
            else:
                radiistepsize = initialradiistepsize
                gammastepsize = initialgammastepsize
                #~ radiirange = initialradiirange
                #~ gammarange = initialgammarange
                
                
            elementsymbol = nextelementsymbol
            minmaebyelement = nextmae
                    

    
    print "\n symbol rc@ g@ MAE RMSE BIAS slope intercept R2 eslope eintercept eR2"
    print minstaticlist
                        
    
    
#print best cycle
outfile.write("BEST CYCLE:\n")
outfile.write(str(beststep[0])+str(beststep[1])+"\n")
        
outfile.write("\n" + 20*"=" + "\n")
outfile.write(str(datetime.datetime.now()))
outfile.write(" TOTAL TIME:" + str(time.time() - start_time) + "seconds")

outfile.close() #close report file







