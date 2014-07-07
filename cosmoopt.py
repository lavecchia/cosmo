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
#parser.add_option("-k", "--keys", dest="extrakeys", default=defextrakeys,
#                  help="KEYS to run MOPAC software", metavar="KEYS")
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
    
#~ if cfg.has_option("system", "step1"):
    #~ stepall = cfg.get("system", "step1").split()
    #~ stepcycles = int(stepall.pop())
    #~ step1 = [stepall,stepcycles]
    #~ stepslist.append(step1)
#~ 
#~ if cfg.has_option("system", "step2"):
    #~ stepall = cfg.get("system", "step2").split()
    #~ stepcycles = int(stepall.pop())
    #~ step2 = [stepall,stepcycles]
    #~ stepslist.append(step2)
#~ 
#~ if cfg.has_option("system", "step3"):
    #~ stepall = cfg.get("system", "step3").split()
    #~ stepcycles = int(stepall.pop())
    #~ step3 = [stepall,stepcycles]
    #~ stepslist.append(step3)
#~ 
#~ if cfg.has_option("system", "step4"):
    #~ stepall = cfg.get("system", "step4").split()
    #~ stepcycles = int(stepall.pop())
    #~ step4 = [stepall,stepcycles]
    #~ stepslist.append(step4)
    #~ 
#~ varlist.append("stepslist")
#~ 
#~ if cfg.has_option("system", "step5"):
    #~ stepall = cfg.get("system", "step5").split()
    #~ stepcycles = int(stepall.pop())
    #~ step5 = [stepall,stepcycles]
    #~ stepslist.append(step5)
#~ varlist.append("stepslist")

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

            acumcycle += 1
            
#Genetic Algorithm run
elif calculationtype == 'ga':
    #initial generation
    memberlist = []
    for member in range(0,nmember):
        memberparamdic = param0dic #initial parameters
        #check restrictions, if not assign parameters again
        tag1sttry = True # the first time into while loop
        while (check_restrictions(memberparamdic,fixlimitdic)==0) or (tag1sttry == True):
            tag1sttry = False # first time into while loop
            # generate new values of parameters with a Gaussian distribution probability from current values
            memberparamdic = [modified_values(currentstep[1], freeparamlist, rangesdic),None] #last element for score function
            #fix values
            try:
                for key in fixlist:
                    memberparamdic[key]=param0dic[key]
            except:
                pass
            
            rsolvtest = paramtestdic["rsolv@"] #take rsolv parameter 
            if nptype=="claverie" or nptype=="claverietype":    
                yrel = 4*PI*NS*(rsolvtest*1.0E-10)*(rsolvtest*1.0E-10)*(rsolvtest*1.0E-10)/3.0 # ratio to use in claverie cavitation term
            else:
                yrel = 0
            
            errorvalue, mad, datadic = calc_error(numberstep, minMAD, memberparamdic[0], datadic, extrakeys, extrakeyssolv + " RSOLV=% .3f" % (rsolvtest), outfile, nptype, yrel)
            memberparamdic[1]=0
            memberlist.append(memberparamdic)

            if mad < minMAD:
                minMAD = mad
                # write summary file with the low MAD step
                print_summary(mad,ncycle,datadic,nptype)
        
                # write parameter values
                paramout = open("list_param.out","w")
                for key, value in sorted(paramtestdic.iteritems()):
                    paramout.write("%-3s \t %.3f\n" % (key,value))
                paramout.close()
            shutil.rmtree(numberstep) # erase current directory
            
            for generation in range(0,ngeneration):
                newgeneration = genetic()
           
            outfile.write(mcmark) #end of line with a mark, @DESC: recenter because below value of function, @PROB: select from probability

            mcmark = "\n"
            
            #number of step, to generate name of files
            numberstep = next(step)
            acumcycle += 1
            
            

#print best cycle
outfile.write("BEST CYCLE:\n")
outfile.write(str(beststep[0])+str(beststep[1])+"\n")

#bottom of report file
#~ outfile.write("\n" + 20*"=" + "\n")
#~ outfile.write("%-30s \t%s \t%s \t%s \t%s\n" % ("compoundname", "dgexp", "dgcalc", "error", "abserror"))

#~ reportlist =[]
#~ for nelement in range(0,len(dgrefdic)):
    #~ reportlist.append([dgrefdic[nelement][0],dgrefdic[nelement][1],dgcalcdic[nelement],(dgrefdic[nelement][1]-dgcalcdic[nelement]), abs(dgrefdic[nelement][1]-dgcalcdic[nelement])])
#~ reportlist.sort(reverse=True, key=lambda a: a[4])
#~ 
#~ tagfirst = True
#~ for compoundname, dgexp, dgcalc, error, abserror in reportlist:
    #~ outfile.write("%-30s \t%.2f \t%.2f \t%.2f \t%.2f\n" % (compoundname, dgexp, dgcalc, error, abserror))
    #~ if abserror < 3.0 and tagfirst==True: # write ============ to divide error values under 3.0 kcal/mol
        #~ outfile.write(10*"="+"\n") 
        #~ tagfirst = False
        
outfile.write("\n" + 20*"=" + "\n")
outfile.write(str(datetime.datetime.now()))
outfile.write(" TOTAL TIME:" + str(time.time() - start_time) + "seconds")

outfile.close() #close report file


print beststep





