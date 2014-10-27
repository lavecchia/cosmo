""" minimize_interface.py

"""

from commonfunctions import * # common function
from lmfit import minimize, Parameters, Parameter, report_fit
import numpy as np
#~ import cosmoopt.py
acumcycle = 0
limitcycle = 25
cycleswithoutdown = 0





def run_minimize(paramdic, maxlimitdic, minlimitdic, datadic, extrakeys, extrakeyssolv, rsolvtest, outfile, nptype, mintotalerror):
    
    mintotalerror = mintotalerror
    params = Parameters()
    for ikey, ivalue in paramdic.iteritems():
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


    result = minimize(fcn2min, params, args=(extrakeys, extrakeyssolv, outfile, nptype, datacompoundnamelist, limitcycle), method='powell')

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
        
    return datadic

cycleswithoutdown = 0
acumcycle = 0

def fcn2min(params, extrakeys, extrakeyssolv, outfile, nptype, datacompoundnamelist, limitcycle=1000):
    global mintotalerror
    global acumcycle
    global datadic
    global cycleswithoutdown
    #~ global limitcycle
    

    print cycleswithoutdown
    print acumcycle
    print mintotalerror
	
	


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
        #~ if cycleswithoutdown > limitcycle:
            #~ exit()
        #~ shutil.rmtree(numberstep)
    shutil.rmtree("0000") # erase current directory
    print "PASO %i: %f"%(acumcycle,totalerror)
    acumcycle += 1
    errors = []
    for i in range(0,len(datacompoundnamelist)):
        errors.append(datadic[datacompoundnamelist[i]]['dgexp']-datadic[datacompoundnamelist[i]]['dgcalc'])
    #~ return errors
    return totalerror

