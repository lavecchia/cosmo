"""ga_interface.py
connect simulation with genetic library.
"""

from genetic import *
import numpy as np


class memberobj(object):
    def __init__(self, paramdic, fitness=None):
        self.paramdic = paramdic
        self.fitness = fitness
        
    def set_fitness(self,newfitness):
        self.fitness = newfitness
        
        

class SetsGeneration(object):
    """Use to connect with genetic library.
    Keyword arguments:
    memberlst -- list of member objects
    mindic -- list with minima values with same format of parameters chromosome
    maxdic -- list with maxima values with same form of parameters chromosome
    """
   
    def __init__(self, memberlst, mindic, maxdic):
        self.memberlst = memberlst
        minlst = []
        maxlst = []
        keylst = []
        paramdic = memberlst[0].paramdic # extrae key del primer miembro (objeto) para armar la lista de minimos y maximos
        for key, value in paramdic.iteritems(): 
            minlst.append(mindic[key])
            maxlst.append(maxdic[key])
            keylst.append(key)
        
        chromosomelst = []
        fitnesslst = []
        for member in self.memberlst:
            fitnesslst.append(member.fitness)
            paramlst = []
            for key in keylst:
                paramlst.append(member.paramdic[key])
            chromosomelst.append(paramlst)
            
        
        self.minlst = minlst 
        self.maxlst = maxlst
        self.keylst = keylst
        self.chromosomelst = chromosomelst
        self.fitnesslst = fitnesslst
        
        
    def next(self):
        """ Return a list of dictionaries with new parameter values
        """
        ga = Generation(self.chromosomelst, self.fitnesslst,
                self.minlst, self.maxlst)
        newgeneration = ga.next()
        newmemberlst = []
        
        for member in newgeneration:
            paramdic = {}
            chromosome = member.get_chrom()
            for n in range(0,len(self.keylst)):
                paramdic[self.keylst[n]] = chromosome[n]
            
            #temporal solution
            #check radii of second period C > N > O > F 
            if (paramdic["rc@C.a"]<paramdic["rc@N.a"]):
                paramdic["rc@N.a"] = paramdic["rc@C.a"]
            
            if (paramdic["rc@N.a"]<paramdic["rc@O.a"]):
                paramdic["rc@O.a"] = paramdic["rc@N.a"]
                
            if (paramdic["rc@O.a"]<paramdic["rc@F.a"]):
                paramdic["rc@F.a"] = paramdic["rc@O.a"]
                
            if (paramdic["rc@S.a"]<paramdic["rc@O.a"]):
                paramdic["rc@S.a"] = paramdic["rc@O.a"]
                
            if (paramdic["rc@P.a"]<paramdic["rc@N.a"]):
                paramdic["rc@P.a"] = paramdic["rc@N.a"]

            newmemberlst.append(memberobj(paramdic))

        #compare input list of members (memberlst) with the new list of members (newmemberlst), to avoid calculate againg the fitness of new old members
        for member2 in newmemberlst:
            for member1 in self.memberlst:
                if member2.paramdic == member1.paramdic: #if the lists of members are equal
                    member2.set_fitness(member1.fitness)


        

        return newmemberlst

