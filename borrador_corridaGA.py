#Genetic Algorithm run
elif calculationtype == 'ga':
    acumcycle = 0 # accumulate the number of cycles
    
    memberlist = [] #list that store members
    #building parents
    for n in range(0,numbermembers):
        #~ if cyclicoptimization == True:
            #~ maxcyclebystep = stepslist[period][1] # number of cycle by step
            #~ freeparamlist = stepslist[period][0]
            #~ if acumcycle >= maxcyclebystep:
                #~ acumcycle = 0 # reset count
                #~ if period == len(stepslist)-1:
                    #~ period = 0 #restart to first step
                #~ else:
                    #~ period += 1 #next step
                #~ paramtestdic = beststep[1]
                #~ datadic = beststep[2]
                #~ print period
                        
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
                
        memberlist.append([paramtestdic])

        
        #number of step, to generate name of files
        numberstep = next(step)
        
     
    #generations
    for ncycle in range(1,maxgen):
        
        #fitness calculation
        fitnesslist = [] #store fitness
        for member in memberlist:
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
               
            fitness = 1/totalerror
            fitnesslist.append(fitness)
            
        genetic = SetsGeneration(memberlist, mindic, maxdix, fitnesslist)
        memberlist = genetic.next()
  
        outfile.write(mcmark) #end of line with a mark

        acumcycle += 1
