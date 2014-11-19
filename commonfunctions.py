# common functions

from commonparameters import *
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
import genetic
from scipy import stats
from pylab import plot,show

# ===============================================================================================
# READ FILES FUNCTIONS
# ===============================================================================================

#read experimental data from a file and return dic with format [[moleculename:deltaGhyd]]
def exp_read(namefile):
    datadic = {}
    try:
        file = open(namefile,"r")
    except:
        print "FILE ERROR: %s has some problems or is not defined. Check it please" % (namefile)
        exit()
    for line in file:
        linediv = line.split()
        datadic[linediv[2]] = {"dgexp":float(linediv[0])} #format of exp data: deltaGhyd ! moleculename
    file.close()
    return datadic
        
        
#read parameters from a file and return dictonary with format [atsymbol: parameter value...]
def parameters_read(namefile):
    paramdic = {}
    try:
        file = open(namefile,"r")
    except:
        print "FILE ERROR: %s has some problems or is not defined. Check it please" %(namefile)
        exit()
    for line in file:
        print line
        linediv = line.split()
        paramdic[linediv[0]] = float(linediv[1]) 
    file.close()
    return paramdic

#read output of MOPAC
def mopacout_read(namefile):
    
    
    # open file
    mopacfile = open(namefile,"r")
    heatformation = None
    for line in mopacfile:
        if KEYFHOF in line:
            linediv = line.split()
            #extract FINAL HEAT OF FORMATION in KCAL/MOL
            heatformation = float(linediv[5])
        elif KEYCA in line:
            linediv = line.split()
            #extract COSMO AREA in SQUARE ANGSTROMS
            cosmoarea = float(linediv[3])
    mopacfile.close()
    if heatformation:
        return heatformation,cosmoarea
    else:
        # try again,open file (horrible esto)
        sleep(1.0)
        mopacfile = open(namefile,"r")
        heatformation = None
        for line in mopacfile:
            if KEYFHOF in line:
                linediv = line.split()
                #extract FINAL HEAT OF FORMATION in KCAL/MOL
                heatformation = float(linediv[5])
            elif KEYCA in line:
                linediv = line.split()
                #extract COSMO AREA in SQUARE ANGSTROMS
                cosmoarea = float(linediv[3])
        mopacfile.close()
        if heatformation:
            return heatformation,cosmoarea
        print "ERROR: A mistake produced to try read the file " + namefile
        exit()
    
    
    
#read COS (cosmo information) files of MOPAC
def cosmoout_read(namefile):
    # open file
    cosmofile = open(namefile,"r")
    cosmoparamlist = []
    passindicator = 0
    for line in cosmofile:
        if passindicator == 1:
            try:
                linediv=line.split()
                atomicnumber = int(linediv[1])
                atomiccosmoarea = float(linediv[7])
                cosmoparamlist.append([atomicnumber,atomiccosmoarea])
            except: #if find empty line
                cosmofile.close()
                if len(cosmoparamlist)==0:
                    print "ERROR: A mistake produced to try read COSMO output " + namefile
                    exit()
                return cosmoparamlist
        if KEYCOSMO in line:
            passindicator = 1

# read .out file and extract number of electron by atoms
def electron_read(namefile):
    electronlist = []
    inputfile = open(namefile, "r")
    tag = 0
    for line in inputfile:
        if "DIPOLE" in line:
            tag = 0
        elif tag == 1:
            linediv = line.split()
            electronlist.append(float(linediv[3])) #electron number
        if "ATOM NO.   TYPE          CHARGE      No. of ELECS.   s-Pop       p-Pop" in line:
            tag = 1
    inputfile.close()
    return electronlist
            
# read .out file and extract number of electron by atoms + atomic number
def atomtypeelectron_read(namefile):
    electronlist = []
    inputfile = open(namefile, "r")
    tag = 0
    for line in inputfile:
        if "DIPOLE" in line:
            tag = 0
        elif tag == 1:
            linediv = line.split()
            electronlist.append([linediv[1],float(linediv[3])]) #atom type - electron number
        if "ATOM NO.   TYPE          CHARGE      No. of ELECS.   s-Pop       p-Pop" in line:
            tag = 1
    inputfile.close()
    return electronlist

    
# ===============================================================================================
# PRINT AND CONVERSION FUNCTIONS
# ===============================================================================================
    
#transform symbol to atomicnumber
def symbol_to_atomicnumber(symbol):
    dictionary = {"ALL":0, "H":1, "C":6, "N":7, "O":8, "F":9, "P":15, "S":16, "Cl":17, "Br": 35, "I": 53}
    try:
        return dictionary[symbol]
    except:
        print "ERROR: the symbol " + symbol + " is not setting in function symbol_to_atomicnumber. Please check this"
        exit()


#transform atomicnumber to symbol
def atomicnumber_to_symbol(atomicnumber):
    dictionary = { 1:"H", 6:"C", 7:"N", 8:"O", 9:"F", 15:"P", 16:"S", 17:"Cl", 35:"Br", 53:"I"}
    try:
        return dictionary[atomicnumber]
    except:
        print "ERROR: the symbol " + str(atomicnumber) + " is not setting in function atomicnumber_to_symbol. Please check this"
        exit()


#transform a key to atomic symbol
#symbols have the follow format: type@Symbol.modified, example: r@S.b, symbol of radii of sulfur atom mod b.
def key_to_symbol(key):
    typeparam,symbolmod = key.split("@")
    try:
        symbol,mod = symbolmod.split(".")
    except:
        symbol = symbolmod
    return symbol
        
        
#transform atomic number and electron number
#enumber: number of electrons
#return symbol + mod, example: S.b
# X.a ----> e@X.a <--- X.b ---> e@X.b <---X.c
def atomicnumber_to_symbolmod(atomicnumber, enumber, paramdic):
    symbol = atomicnumber_to_symbol(atomicnumber)
    try:
        umbralb = paramdic["e@"+symbol + ".b"]
        umbrala =  paramdic["e@"+symbol + ".a"]
        if  enumber > umbralb:
            mod = "c"
        elif umbrala < enumber <= umbralb:
            mod = "b"
        else:
            mod = "a"
    except:
        try:
            umbral = paramdic["e@"+symbol + ".a"]
            if enumber > umbral:
                mod = "b"
            else:
                mod = "a"
        except:
            mod = "a"   #"a" is the default mod
    return symbol + "." + mod
        
        
#transform parameter dictionary to MOPAC format VDW(Cl=2.33;Br=2.50...)
def radii_to_text(paramdic):
    line = "VDW("
    for key,value in paramdic.iteritems():
        if "rc@" in key:
            symbol = key_to_symbol(key)
            line += str(symbol)+"="+str("{0:.4f}".format(value)+";")
    line = line.strip(";") + ")"
    line.strip()
    return line
    

#print parameters dictionary to text
def print_param(paramdic):
    line=""
    for key, value in sorted(paramdic.iteritems()):
        line += " %-2s % .5f " % (key,value)
    line = line.strip()
    return line
    
    
def print_summary(totalerror, mae, rmse, bias, r2, slope, intercept,ncycle,datadic,nptype):
    summaryfile = open(DEFSUMMARYFILE,"w")
    summaryfile.write("Cycle: %i TotError: %2.4f MAE: %2.4f RMSE: %2.4f BIAS: %2.4f R2: %1.5f Slope: %2.3e Intercept: %2.3e" % (ncycle, totalerror, mae, rmse, bias, r2, slope, intercept) + "\n")
    
    if nptype == "claverie":
        summaryfile.write("%-48s %s %s %s %s %s %s %s %s\n" % ("compoundname", "dgexp", "dgcalc", "error", "abserror","hof_gas","hof_cosmo","delta_hof","npterm(disp,cavitation)"))
        reportlist =[]
        for compoundname,value in datadic.iteritems():
            reportlist.append([compoundname,datadic[compoundname]["dgexp"],datadic[compoundname]["dgcalc"],(datadic[compoundname]["dgexp"]-datadic[compoundname]["dgcalc"]), abs(datadic[compoundname]["dgexp"]-datadic[compoundname]["dgcalc"]),datadic[compoundname]["hofgas"],datadic[compoundname]["hofcosmo"], datadic[compoundname]["hofcosmo"] - datadic[compoundname]["hofgas"],datadic[compoundname]["cavitation"], datadic[compoundname]["dispterm"]])
        reportlist.sort(reverse=True, key=lambda a: a[4])

        #~ tagfirst = True
        for compoundname, dgexp, dgcalc, error, abserror, hofgas, hofcosmo, deltahof, cavitation, dispterm in reportlist:
            summaryfile.write("%-48s % +.2f % +.2f % +.2f % +.2f % +.2f %+.2f % +.2f % +.2f % +.2f % +.2f\n" % (compoundname, dgexp, dgcalc, error, abserror, hofgas, hofcosmo, deltahof, cavitation + dispterm,cavitation,dispterm))
            #~ if abserror < 3.0 and tagfirst==True: # write ============ to divide error values under 3.0 kcal/mol
                #~ summaryfile.write(10*"="+"\n") 
                #~ tagfirst = False
        
    elif nptype == "gammasasa":
        summaryfile.write("%-48s %s %s %s %s %s %s %s %s\n" % ("compoundname", "dgexp", "dgcalc", "error", "abserror","hof_gas","hof_cosmo","delta_hof","npterm"))
        reportlist =[]
        for compoundname,value in datadic.iteritems():
            reportlist.append([compoundname,datadic[compoundname]["dgexp"],datadic[compoundname]["dgcalc"],(datadic[compoundname]["dgexp"]-datadic[compoundname]["dgcalc"]), abs(datadic[compoundname]["dgexp"]-datadic[compoundname]["dgcalc"]),datadic[compoundname]["hofgas"],datadic[compoundname]["hofcosmo"], datadic[compoundname]["hofcosmo"] - datadic[compoundname]["hofgas"],datadic[compoundname]["npterm"]])
        reportlist.sort(reverse=True, key=lambda a: a[4])

        #~ tagfirst = True
        for compoundname, dgexp, dgcalc, error, abserror, hofgas, hofcosmo, deltahof, npterm in reportlist:
            summaryfile.write("%-48s % +.2f % +.2f % +.2f % +.2f % +.2f %+.2f % +.2f % +.2f\n" % (compoundname, dgexp, dgcalc, error, abserror, hofgas, hofcosmo, deltahof, npterm))
            #~ if abserror < 3.0 and tagfirst==True: # write ============ to divide error values under 3.0 kcal/mol
                #~ summaryfile.write(10*"="+"\n") 
                #~ tagfirst = False
        
    summaryfile.close()


#generate list of files inside directory with extension ext
def list_files(directory,ext=DEFINPUTEXT):
    filelist = glob.glob(directory + "/*" + ext)
    return sorted(filelist)


#make MOPAC gas or cosmo input file taking a file as template
def inputfile_make(templatename, extrakeys, extrakeyssolv = "", paramdic = "", step = ""):
    template = open(templatename, "r")
    #if there is define a COSMO calculation
    if extrakeyssolv: 
        try:
            os.stat(step)
        except:
            os.makedirs(step) # create dir 00X
        inputfile = open(step+"/"+os.path.basename(templatename).replace(".mop","_" + step + ".mop"), "w") #save the file inside dir 001/archivo_001.mop
    #if there is not define a COSMO calculation, then is a GAS calculation
    else:
        try:
            os.stat("gas")
        except:
            os.makedirs("gas") # create dir gas
        inputfile = open("gas/"+os.path.basename(templatename).replace(".mop","_gas.mop"), "w") #save the file inside dir gas/archivo_gas.mop
    
    if paramdic != "":
        radiitext = radii_to_text(paramdic)
    else:
        radiitext = ""
    
    #read template
    for line in template:
        if KEYDEFAULT in line:
            line = line.replace(KEYDEFAULT, extrakeys + " " + extrakeyssolv + " " + radiitext)
        inputfile.write(line)
        
    template.close()
    inputfile.close()
    return 0
   
    
#run calculation in serial or parallel
def exec_commands(cmds, cores = DEFCORE):
    ''' Exec commands in parallel in multiple process 
    (as much as we have CPU)
    '''
    if not cmds: return # empty list

    def done(p):
        return p.poll() is not None
    def success(p):
        return p.returncode == 0
    def fail():
        sys.exit(1)

    max_task = cores
    processes = []
    while True:
        while cmds and len(processes) < max_task:
            task = cmds.pop()
            #print list2cmdline(task)
            processes.append(Popen(task))

        for p in processes:
            if done(p):
                if success(p):
                    processes.remove(p)
                else:
                    fail()

        if not processes and not cmds:
            break
        else:
            pass
            #~ time.sleep(1E-13)
            


# ===============================================================================================
# MAIN AND SOLVATION FUNCTIONS
# ===============================================================================================
    
#calculation of cavitation term with simplified Pierotti equation
def calc_pierotti(rsolv,rsphere,yrel):
    rsolv = rsolv * 1.0E-10 # angstrom to m
    rsphere = rsphere * 1.0E-10 # angstrom to m
    yratio = yrel/(1-yrel) 
    return RCONST*TEMP*(-math.log(1-yrel) + (3*yratio)*(rsphere/rsolv)+(3*yratio+4.5*(yratio)*(yratio))*(rsphere/rsolv)*(rsphere/rsolv)) * J2KCAL


#calculation of no-polar term as sum (gamma_i * area_i)
#input:
#   comosparamlist: [[atomicnumber_i, cosmoarea_i],...
#   gammadic: {atomicnumber_i:gamma_i, ...
#output:
#   no-polar term (float)
def np_calc(cosmoparamlist,paramdic,nptype,yrel,electronlist=None):
    if nptype == "claverie":  # if claverie-pierotti term is activated
        dispterm = 0
        cavitationterm = 0
        i = 0
        rsolv = paramdic["rsolv@"]
       
        
        for atomicnumber,cosmoarea in cosmoparamlist:
            symbolmod = atomicnumber_to_symbolmod(atomicnumber,electronlist[i],paramdic)
            gkey = "g@" + symbolmod
            gamma = paramdic[gkey]
            
            # check if there r@ param, else look for rc@, then else print warning and exit
            try:
                rkey = "r@" + symbolmod
                radii = paramdic[rkey]
            except:
                try:
                    rkey = "rc@" + symbolmod
                    radii = paramdic[rkey]
                except:
                    print "There is not radius defined (r@ or rc@) for " + symbolmod + ". Please check it"
                    exit()
            
            #~ dispterm += gamma*cosmoarea
            dispterm += gamma*cosmoarea*(rsolv+radii)*(rsolv+radii)/(radii*radii) #change to scale SES area to SAS: (r+Rsolv)^2/r^2
            
            #~ cavitationterm += calc_pierotti(rsolv,radii,yrel)*cosmoarea/((radii+rsolv)*(radii+rsolv)*PI*4/3) # sum of Claverie terms
            cavitationterm += calc_pierotti(rsolv,radii,yrel)*cosmoarea/((radii)*(radii)*PI*4/3) # sum of Claverie terms
            i += 1
            
        npterm = dispterm + cavitationterm
        return npterm, dispterm, cavitationterm
        
    elif nptype == "electronmod": # if electron modification term is activated (is a test)
        npterm = 0
        i = 0
        for atomicnumber,cosmoarea in cosmoparamlist:
            npterm += paramdic[key]*cosmoarea*electronlist[i]
            i += 1
        return npterm        
    elif nptype == "gammasasa":  # or simple no-polar term: sum_i (gamma_i x cosmoarea_i)
        dispterm = 0
        cavitationterm = 0
        i = 0
        rsolv = paramdic["rsolv@"]
        npterm = 0
           
        for atomicnumber,cosmoarea in cosmoparamlist:
            symbolmod = atomicnumber_to_symbolmod(atomicnumber,electronlist[i],paramdic)
            gkey = "g@" + symbolmod
            gamma = paramdic[gkey]

            #~ npterm += gamma*cosmoarea
            try:
                rkey = "rc@" + symbolmod
                radii = paramdic[rkey]
            except:
                print "There is not radius defined (r@ or rc@) for " + symbolmod + ". Please check it"
                exit()
            
            npterm += gamma*cosmoarea
            #~ npterm += gamma*cosmoarea*(rsolv+radii)*(rsolv+radii)/(radii*radii) #change to scale SES area to SAS: (r+Rsolv)^2/r^2
            
            i += 1
        return npterm
    

# calculation of no-polar term with only one gamma value
def np_calc_simple(cosmoarea,gammavalue):
    return gammavalue*cosmoarea
    
        
#error between experimental and calculated data
# r0list original radii
# paramtestlist list with new radii and gamma to test
# dgreflist experimental Gibbs energy
# dgcalc calculated Gibbs energy to test
# hfgaslist: list of heat of formations in gas phase
def calc_error(numberstep, paramdic, datadic, extrakeys, extrakeyssolv, outfile, nptype ,yrel=None):
        
    #generate MOPAC input files corresponding to a numberstep generation

    for key, value in datadic.iteritems():
        inputfile_make(value["template"],extrakeys,extrakeyssolv,paramdic,numberstep)

    #run COSMO phase calculation
    cosmoinputlist = list_files(numberstep)
    commands = []
    for cosmoinput in cosmoinputlist:
        commands.append([MOPACPATH, cosmoinput])
    exec_commands(commands)        
        
    gasindex = 0  
    errorlist = []
    xdata = []
    ydata = []
     
    for key, value in datadic.iteritems():
        cosmoinput =str(numberstep)+"/"+os.path.basename(datadic[key]["template"]).replace(".mop","_" + str(numberstep) + ".mop")
        
        hofcosmo, areacosmo = mopacout_read(cosmoinput.replace(".mop",".out")) #read MOPAC output from .out
        cosmoparamlist = cosmoout_read(cosmoinput.replace(".mop",".cos")) #read cosmo param from .cos
        electronlist = electron_read(cosmoinput.replace(".mop",".out")) #read number of electron by atom from MOPAC output 
        
        datadic[key]["hofcosmo"] = hofcosmo #extract HOF
        datadic[key]["areacosmo"] = areacosmo #extract AREA
        
        #no-polar term
           
        if nptype == "claverie":
            npterm, dispterm, cavitationterm = np_calc(cosmoparamlist,paramdic,nptype,yrel,electronlist)
            datadic[key]["dispterm"] = dispterm
            datadic[key]["cavitation"] = cavitationterm
        #~ elif nptype == "electronmod":
            #~ npterm = np_calc(cosmoparamlist,gammatestdic,rtestdic,rsolv,nptype,yrel,electronlist)
            #~ npdetaillist.append(npterm)
        elif nptype == "gammasasa":
            npterm = np_calc(cosmoparamlist,paramdic,nptype,yrel,electronlist)
            datadic[key]["npterm"] = npterm
        else:
            #mono gamma
            gamma = paramtestlist[len(rtestlistcomplete)]
            npterm = np_calc_simple(areacosmo,gamma)
            npdetaillist.append(npterm)
    
        #Hf_gas
        hofgas = value["hofgas"]
        gasindex += 1
        
        # \delta G_calc = Hf_cosmo - Hf_gas + no-polar_term
        dgcalc = hofcosmo - hofgas + npterm 
        #~ 
        #~ dgcalc = npterm #reemplazado para optimizar solo parte no polar
        #~ dgcalc = hofcosmo - hofgas #reemplazado para optimizar solo la parte electrostatica
        
        datadic[key]["dgcalc"]= dgcalc
        xdata.append(datadic[key]["dgexp"]) #experimental data as x-axis
        ydata.append(datadic[key]["dgcalc"]) #calc data as y-axis
        
        datadic[key]["error"] = datadic[key]["dgexp"]-datadic[key]["dgcalc"]
        errorlist.append(datadic[key]["error"])
        
    errorarray = np.array(errorlist)
    mae = np.mean(abs(errorarray))
    rmse = math.sqrt(np.mean((errorarray)**2))
    bias = np.mean(errorarray)
    slope, intercept, r2 = fit_lineal(xdata,ydata) #linear fit

    
    totalerror = mae #var to optimize
    prevline = "%-5s Err: %3.4f MAE: %3.4f RMSE: %3.4f BIAS: %3.4f R2: %1.5f " % (numberstep, totalerror, mae, rmse, bias, r2)
    outfile.write(prevline + print_param(paramdic))
     
     
    return totalerror, mae, rmse, bias, r2, slope, intercept, datadic


#calculation of metrics by element
def calc_staticsbyelement(numberstep, datadic, elementsymbol):
    errorlist = []
    electronlistbyelementlist =[]
    xdata = []
    ydata = []
    
    for key, value in datadic.iteritems():
        cosmoinput =str(numberstep)+"/"+os.path.basename(datadic[key]["template"]).replace(".mop","_" + str(numberstep) + ".mop")
        
        #~ hofcosmo, areacosmo = mopacout_read(cosmoinput.replace(".mop",".out")) #read MOPAC output from .out
        cosmoparamlist = cosmoout_read(cosmoinput.replace(".mop",".cos")) #read cosmo param from .cos
        atomtypeelectronlist = atomtypeelectron_read(cosmoinput.replace(".mop",".out")) #read atomic symbol - number of electron from MOPAC output 
        tag = False
        # atomtypeelectronlist = [[atomsymbol1,electronnumber1], [atomsymbol2,electronnumber2], ...]
        for atomtypeelectron in atomtypeelectronlist:
            if elementsymbol in atomtypeelectron[0]:
                errorlist.append(datadic[key]["error"])
                electronlistbyelementlist.append(atomtypeelectron[1])
                
                if tag == False: # store dgexp vs dgcalc point only once
                    tag = True
                    xdata.append(datadic[key]["dgexp"]) #experimental data as x-axis
                    ydata.append(datadic[key]["dgcalc"]) #calc data as y-axis
        
    errorarray = np.array(errorlist)
    try:
        mae = np.mean(abs(errorarray))
        rmse = math.sqrt(np.mean((errorarray)**2))
        bias = np.mean(errorarray)
    
        # calc vs experimental
        slope, intercept, r2 = fit_lineal(xdata,ydata) #linear fit
    
        # error vs electron number
        errorslope, errorintercept, errorr2 = fit_lineal(electronlistbyelementlist,errorlist) #linear fit
    except:
            pass
    return mae, rmse, bias, r2, slope, intercept, errorr2, errorslope, errorintercept




#
#  name: check_restrictions
#  Check if a serie of conditions are True or False
#
#  @param paramdic: parameters to tested
#  @type  paramtestlist: dictionary
#
#  @return If check is true or not
#  @rtype  binary
def check_restrictions(paramdic, fixlimitdic):
    checker=0
    #check radii F < Cl < Br < I
    if (paramdic["rc@F.a"]<=paramdic["rc@Cl.a"]) and (paramdic["rc@Cl.a"]<=paramdic["rc@Br.a"]) and (paramdic["rc@Br.a"]<=paramdic["rc@I.a"]):
        checker = 1
    else:
        #~ #~#print "fail order between F < Cl < Br < I"
        return 0
    
    #check radii O <= S 
    if paramdic["rc@O.a"] < paramdic["rc@S.a"]:
        checker = 1
    else:
        #~ #print "fail order between O < S"
        return 0
    
    #check radii N < P
    if paramdic["rc@N.a"] <= paramdic["rc@P.a"]:
        checker = 1
    else:
        #~ #~#print "fail order between N < P"
        return 0
      
    #check radii of second period C > N > O > F 
    if (paramdic["rc@C.a"]>=paramdic["rc@N.a"]) and (paramdic["rc@N.a"]>=paramdic["rc@O.a"]) and (paramdic["rc@O.a"]>=paramdic["rc@F.a"]):
        checker = 1
    else:
        #~ #~#print "fail order between first row"
        return 0
        
    for key, value in paramdic.iteritems():
        if fixlimitdic[key][0] <= value <= fixlimitdic[key][1]:
            checker = 1
        else:
            #~ print key + " limites:" + str(fixlimitdic[key][0]) + " - " + str(fixlimitdic[key][1])
            return 0
	        
    return checker
    

#  
#  name: make_gaussmodification
#  Return an aleatory number with a Gaussian distribution obtained from a number. This number is inside a range determined by rangevalue.
#
#  @param value0: Number where the Gaussian distribution is center 
#  @type value0: number
#
#  @param rangevalue: Specify the modification range.
#  @type rangevalue: number
#
#  @return Modified value.
#  @rtype  number
#
def make_gaussmodification(value0,rangevalue):
    #gauss
    #~ return value0 + rangevalue * np.random.normal(0,0.2)
    
    #uniform
    return value0 + rangevalue * np.random.uniform(-1,1)


#  
#  name: modified_values
#  Return a new dictionary with modified values
#
#  @param paramdic: contain name and value of parametes to optimize
#  @type paramdic: dictonary
#
#  @param freeparamlist: name of parameters free of modified
#  @type rangevalue: list
#
#  @return rrange, gammarange, rsolvrange
#  @rtype  number
#
def modified_values(paramdic, freeparamlist, rangesdic):
    newparamdic = {} #store new values of parameters
    for key, value in paramdic.iteritems():
        if key in freeparamlist: #if the parameter is free, select the range to modified
            if "r@" in key:     #radii range
                paramrange = rangesdic["radii"]
            elif "g@" in key:   #gamma range
                paramrange = rangesdic["gamma"]
            elif "rsolv@" in key:   #rsolv range
                paramrange = rangesdic["rsolv"]
            elif "rc@" in key: #cosmo radii range
                paramrange = rangesdic["cosmoradii"]
            elif "e@" in key:
                paramrange = rangesdic["eumbral"]
            elif "k@" in key:
                paramrange = rangesdic["k"]
            newparamdic[key] = make_gaussmodification(paramdic[key],paramrange)
        else:
            newparamdic[key] = paramdic[key] #if the parameter is fix, then take the old value
    return newparamdic
        
        
#
#  name: check_excludeatom
#  Check if an atom/atoms is/are present in a compound
#
#  @param filename: file to check
#  @type  filename: string
#
#  @param excludelist: list with symbols of atoms to exclude
#  @type  excludelist: list of string (eg. ["S","I"])
#
#  @return If there are atoms of excludelist in filename then is True, else is False
#  @rtype  binary
def check_excludeatom(filename, excludelist):
    filein = open(filename, "r")
    for line in filein:
        if any ((symbol+" ") in line[0:3] for symbol in excludelist) == True:
            filein.close()
            return True
    filein.close()
    return False
    

def fit_lineal(x,y):  
    slope, intercept, r_value, p_value, slope_std_error = stats.linregress(x, y)
    return slope, intercept, r_value*r_value
    
# control MC temperature
def temperature_control(temp,mcmarklist,lastmark):
    limitstore = 20
    numberstoremark = float(len(mcmarklist))
    mcmarklist.append(lastmark)
    
    if numberstoremark<limitstore:
        return temp, mcmarklist
    else:    
        mcmarklist.pop(0)
        numberdesc = mcmarklist.count(' @DESC\n')
        numberprob = mcmarklist.count(' @PROB\n')
        print "%s %f"%("prob",(float(numberprob + numberdesc)/numberstoremark))
        if numberprob/numberstoremark > 0.4:
            return temp * 0.95, mcmarklist
        elif (numberprob + numberdesc)/numberstoremark < 0.2:
            return temp * 1.05, mcmarklist
        else:
            return temp, mcmarklist
        
     



#Genetic Algorithm Core
#~ def runGA(ngeneration,

# ===============================================================================================
# CLASSES
# ===============================================================================================

class Datatest:
    def __init__(self,radiilist,gammalist,rsolv):
        self.radiilist = radiilist
        self.gammalist = gammalist
        self.rsolv = rsolv
        self.fixparam = fixparam
        



