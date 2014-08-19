"""ga_interface.py
connect simulation with genetic library.
"""

from genetic import *



class SetsGeneration(object):
    """Use to connect with genetic library.
    Keyword arguments:
    memberlst -- list of member' parameters 
    mindic -- list with minima values with same format of parameters chromosome
    maxdic -- list with maxima values with same form of parameters chromosome
    fitnesslst -- list with the fitness
    """
   
    def __init__(self, memberlst, mindic, maxdic, fitnesslst):
        
        minlst = []
        maxlst = []
        keylst = []
        for key, value in memberlst[0].iteritems: #extrae del primer elemento para armar la lista de minimos y maximos
            minlst.append(mindic[key])
            maxlst.append(maxdic[key])
            keylst.append(key)
        
        chromosomelst = []
        
        for member in memberlst:
            paramlst = []
            for key in self.keylst:
                paramlst.append(member[key])
            chromosomelst.append(paramlst)
            
        
        self.minlst = minlst 
        self.maxlst = maxlst
        self.keylst = keylst
        self.chromosomelst = chromosomelst
        self.fitnesslst = fitnesslst
        
        
    def next(self):
        """ Return a list with new parameter values
        """
        ga = Generation(self, self.chromosomelst, self.fitnesslst,
                self.minlst, self.maxlst)
        newgeneration = ga.next()
        newmemberlst = []
        
        for member in newgeneration:
            paramdic = {}
            chromosome = member.get_chrom()
            for n in range(0,len(self.keylst)):
                paramdic[self.keylst[n]] = chromosome[n]
            newmemberlst.append(paramdic)
        return newmemberlst

