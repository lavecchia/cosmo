"""genetic.py 
Library for Genetic Algorithm 
"""
DEFLIMIT = 2.0

import random

class Generation(object):
    """Generation is used as easy interface to interact with the 
    library.
    
    Keyword arguments:
    vallstlst -- list of list with values of chromosomes
    fitnesslst -- list of fitness respectivily of vallstlst
    minlst -- list with minima values with same form of chromosome
    maxlst -- list with maxima values with same form of chromosome
    """
    def __init__(self, vallstlst=[], fitnesslst=[], minlst=None, maxlst=None):
        self.memberlst = []
        self.fitnesslst = fitnesslst
        self.minlst = minlst
        self.maxlst = maxlst
        for i in range(0,len(vallstlst)):
            self.build_member(vallstlst[i],fitnesslst[i],minlst,maxlst)
                
    def set_memberlst(self, memberlst):
        ''' Set '''
        self.memberlst = memberlst
    
    def add_memberlst(self, memberlst):
        self.memberlst += memberlst
        
    
    def add_member(self, member):
        self.memberlst.append(member)
    
    
    def build_member(self, vallst, fitness, minlst=None, maxlst=None):
        chromosome = []
        if minlst == None:
            minlst = self.minlst
        if maxlst == None:
            maxlst = self.maxlst
        for i in range(len(vallst)):
            val = vallst[i]
            try:
                minval = minlst[i]
            except:
                minval = None
            try:
                maxval = maxlst[i]
            except:
                maxval = None
            gene = Gene(val, minval, maxval)
            chromosome.append(gene)
        member = Member(chromosome,fitness)        
        self.memberlst.append(member)
    
    
    def next(self):
        ga = GeneticAlgorithm(self.memberlst)
        return ga.run_ga()


class GeneticAlgorithm(object):
    """Contain the main functions of the library.
    
    Keyword arguments:
    memberlst -- list of member object
    """
    def __init__(self, memberlst):
        self.MAXRULETTEWHEELCYCLES = 400
        self.ALTERCONST = 0.1
        self.memberlst = memberlst
        self.MUTPERCENT = 0.05
        self.NEIGHBORPERCENT = 0.05
        self.MAXTOURNAMENTCYCLES = 400
        self.TOURNAMENTWEAKWINPROB = 0.05
    
    def check_limits(self,val,minval,maxval):
        """Common function to check if value is between limits.
    
        Keyword arguments:
        val -- float number
        minval -- minimum float value 
        maxval -- maximum float value 
        """
        if val < minval:
            return minval
        elif val > maxval:
            return  maxval
        else:
            return val
    
    
    def run_ga(self):
        """Run the genetic algorithm.
        
        Output:
        list of members
        """
        # select the half of parents
        parents = self.selection_parent()
        # and complete with children
        children = self.build_randomchildren(parents)
        return parents + children
        

    """ 
    SELECTION FUNCTIONS 
    """
    
    def selection_parent(self):
        """Selection of half of parents in function of fitness values.

        (only is implemented roulette wheel selector)
        """
        #~ return self.selection_roulettewheel()
        return self.selection_half()
    
    
    def selection_half(self):
        """Parent selector alternative
        Implementation of half selector.
        Take the best half of members
        Output:
        Return selected parent
        """
        
        memberlist = sorted(self.memberlst, key=lambda x: x.fitness, reverse=True)
        parentlst = []
        
        for nmember in range(len(memberlist)/2):
            parentlst.append(memberlist[nmember])
        return parentlst
        
        
    def selection_roulettewheel(self):
        """Implementation of Roulette Wheel selector.
        
        Output:
        Return selected parent
        """
        sumfitness =  sum(member.fitness for member in self.memberlst)
        parentlst = []
        memberlst = sorted(self.memberlst, key=lambda x: x.fitness, reverse=True)
        
        #the best solution at the moment is add automatically
        parentlst.append(memberlst[0])
        cycles = 0
        while len(parentlst)< (len(self.memberlst)/2):
            current = 0
            pick = random.uniform(0, sumfitness)
            for member in memberlst:
                current += member.fitness
                if (current > pick) and (member not in parentlst) and len(parentlst)< (len(self.memberlst)/2):
                    parentlst.append(member)
            if cycles > self.MAXRULETTEWHEELCYCLES:
                print "Warning! Cycles superate MAXRULETTEWHEELCYCLES"
                return parentlst
            cycles += 1
        return parentlst
        
        


    def selection_tournament(self):
        """Parent selector alternative
        Implementation of tournament selector.
        
        Output:
        Return selected parent
        """
        parentlst = []
        memberlist = sorted(self.memberlst, key=lambda x: x.fitness, reverse=True)
        #the best solution at the moment is add automatically
        parentlst.append(memberlst[0])
        
        while len(parentlst)<(len(self.memberlst)/2):
            #select parent
            parent1 = random.choice(memberlst)
            parent2 = random.choice(memberlst)
            if (random.random < self.TOURNAMENTWEAKWINPROB) or (parent1.fitness >= parent2.fitness):
                if parent1 not in parenlst:
                    parentlst.append(parent1)
            else:
                if parent2 not in parentlst:
                    parentlst.append(parent2)

            if cycles > self.MAXTOURNAMENTCYCLES:
                print "Warning! Cycles superate MAXTOURNAMENTCYCLES"
                return parentlst
            cycles += 1
        return parentlst
        
        
    def build_randomchildren(self, parentlst):
        """Return number of children object same as half of parents.
        
        Keyword arguments:
        parentlst -- list of selected parents
        
        Output
        childrenlst -- list of new member objects
        """
        childrenlst=[]
        #copy chromosome list
        for nchild in range(len(parentlst)):
            #build a child near a parent (neighbor child)
            if random.random() < self.NEIGHBORPERCENT:
                #select the parent
                parent = random.choice(parentlst)
                childrenlst.append(self.neighbor_chromosome(parent.chromosome))

            #build a child from crossover its parent chromosomes
            else:
                #select the parent 1
                parent1 = random.choice(parentlst)
                #~ parent1 = parentlst[random.randint(0,len(parentlst)-1)]
                #~ #select the parent 2
                parent2 = random.choice(parentlst)
                #~ parent2 = cpmemberlst.pop(random.randint(0,len(cpmemberlst)-1))
                #~ parent2 = parentlst[random.randint(0,len(parentlst)-1)]
                childrenlst.append(self.crossover_chromosome(parent1.chromosome,
                    parent2.chromosome))
        return childrenlst
    
    
    """
    CROSSOVER FUNCTIONS
    """
    
    def crossover_chromosome(self, chromosome1, chromosome2):
        """Chrossover the chromosomes.
        
        Keywords arguments:
        chromosome1 -- list of gene objects
        chromosome2 -- list of gene objects
        
        Output:
        
        """
        childgenlst = []
        for igen in range(0,len(chromosome1)):
            
            crossgen = self.crossover_gene(chromosome1[igen],chromosome2[igen])
            # mutation of gen 
            crossgen = self.mutation_gene(crossgen)
            childgenlst.append(crossgen)
            
        return Member(childgenlst)

    
    def crossover_gene(self,gene1,gene2):
        """Interchange genes, redirect depend of type of gene.
        
        Keyword arguments:
        gene1 -- gene object
        gene2 -- gene object
        
        Output:
        new gene object
        """
        if type(gene1.val)!=type(gene2.val):
            if "float" in (str(type(gene1.val))) and ("float" in str(type(gene2.val))):
                pass
            else:
                print "ERROR: genes with different types"
                exit()
        if isinstance(gene1.val,float):
            return self.crossover_floatgene(gene1,gene2)
        elif isinstance(gene1.val,list):
            return self.crossover_lstgene(gene1,gene2)
        else:
            # print gene1.val
            print "ERROR: Type of gene is not defined"
            exit()

    
    def crossover_lstgene(self,lstgene1,lstgene2):
        """Crossover list type genes.
        
        Keyword arguments:
        lstgene1 -- gene object with list as value
        lstgene2 -- gene object with list as value
        
        Output:
        A gene gene object with list as value
        """
        rndselector = random.random()
        # gaussian center in 0 x ALTERCONST
        rndalter = random.gauss(0,0.3) * self.ALTERCONST
        newgenelst = []
        # random number select one of the two input list of genes
        if rndselector < 0.5:
            lstgeneselect = lstgene1
        else:
            lstgeneselect = lstgene2
        # iterate over all genes in the selected list of genes
        valuelst = []
        for i in range(len(lstgeneselect.val)):
            minval = lstgeneselect.minval[i]
            maxval = lstgeneselect.maxval[i]
            value = lstgeneselect.val[i]
            newval = self.check_limits(value + value * rndalter, minval, maxval)
            valuelst.append(newval)
        lstnewgene = Gene(valuelst, lstgeneselect.minval, lstgeneselect.maxval)
        return lstnewgene
            
        
    def crossover_floatgene(self,floatgene1,floatgene2):
        """Crossover genes with float values.
        
        Keyword arguments:
        floatgene1 -- gene object with float value
        floatgene2 -- gene object with float value

        Output:
        New gene
        """
        rndselector = random.random()
        # gaussian center in 0 x ALTERCONST
        rndalter = random.gauss(0,0.3) * self.ALTERCONST

        if rndselector < 0.5:
            newval = self.check_limits(floatgene1.val + floatgene1.val * rndalter,floatgene1.minval,floatgene1.maxval)
            return Gene(newval,floatgene1.minval,floatgene1.maxval)
        else:
            val = floatgene2.val + floatgene2.val * rndalter
            minval = floatgene2.minval
            maxval = floatgene2.maxval
            # check if value is between limits
            newval = self.check_limits(val,minval,maxval)
            return Gene(newval,floatgene2.minval,floatgene2.maxval)
           

            
    """
    MUTATION FUNCTIONS
    """
   
    def mutation_gene(self,gene):
        """Mutation of a gene, that depends of type of its value.
        """
        # do mutation on val
        rnd = random.random() 
        if rnd > self.MUTPERCENT:
            return gene
        else:
            if isinstance(gene.val,float):
                gene.val = (self.mutation_floatgene(gene))
                return gene
            elif isinstance(gene.val,list):
                gene.val = (self.mutation_lstgene(gene))
                return gene
            
    
    def mutation_lstgene(self,lstgene):
        """Mutation a gene with list value.
        
        Keyword arguments:
        lstgene -- gene object with list as value
        
        Output:
        List of new values restricted by minval a maxval
        """
        mutvallst = []
        for i in range(len(lstgene.val)):
            newval = random.uniform(lstgene.minval[i], lstgene.maxval[i])
            mutvallst.append(newval)
        return mutvallst
            
    
    def mutation_floatgene(self,floatgene):
        """Mutation a gene with float value.
        
        Keyword arguments:
        lstgene -- gene object with float as value
        
        Output:
        Float as new value restricted by minval a maxval
        """
        minval = floatgene.minval
        maxval = floatgene.maxval
        mutval = random.uniform(minval,maxval)
        return mutval
        
    
    """
    NEIGHBOR FUNCTIONS
    """
    
    def neighbor_chromosome(self, chromosome):
        """Buid a neighbor chromosome from a parent chromosome.
        Modified a little the input chromosome
        
        Keywords arguments:
        chromosome -- list of gene objects
        
        Output:
        List of new gene objects
        """
        childgenlst = []
        for igen in range(0,len(chromosome)):
            neighborgen = self.neighbor_gene(chromosome[igen])
            childgenlst.append(neighborgen)
            
        return Member(childgenlst)
        
        
    def neighbor_gene(self,gene):
        """Build a neighbor gene, redirect depend of type of gene.
        
        Keyword arguments:
        gene -- gene object

        Output:
        new gene object
        """
 
        if isinstance(gene.val,float):
            return self.neighbor_floatgene(gene)
        elif isinstance(gene.val,list):
            return self.crossover_lstgene(gene)
        else:
            # print gene.val
            print "ERROR: Type of gene is not defined"
            exit()
 

    def neighbor_floatgene(self,floatgene):
        """The value of float gene is modified a little.
        
        Keyword arguments:
        floatgene -- gene object with float as value
        
        Output:
        New gene
        """
        rndselector = random.random()
        # gaussian center in 0 x ALTERCONST
        rndalter = random.gauss(0,0.3) * self.ALTERCONST

        newval = self.check_limits(floatgene.val + floatgene.val * rndalter,floatgene.minval,floatgene.maxval)
        return Gene(newval,floatgene.minval,floatgene.maxval)
 
  
    def neighbor_lstgene(self,lstgene):
        """The value of lstgene is modified a little.
        
        Keyword arguments:
        lstgene -- gene object with list as value
        
        Output:
        New gene
        """
        rndselector = random.random()
        # gaussian center in 0 x ALTERCONST
        rndalter = random.gauss(0,0.3) * self.ALTERCONST
        newgenelst = []
        valuelst = []
        for i in range(len(lstgene.val)):
            value = lstgene.val[i]
            newval = self.check_limits(value + value * rndalter, minval, maxval)
            valuelst.append(newval)
        return Gene(valuelst, lstgene.minval, lstgene.maxval)
              
                
        
class Member(object):
    """Member is an individue.
    
    Keyword arguments:
    genelst -- list of genes
    fitness -- positive float value
    """
    def __init__(self, genelst, fitness=None):
        self.chromosome = genelst
        self.fitness = fitness
        
    def get_chrom(self):
        genevallst = []
        for gene in self.chromosome:
            genevallst.append(gene.val)
        return genevallst
        
        
class Gene(object):
    """The smaller object of algorithm.
    
    Keyword arguments:
    val -- float o list of float
    minval -- float o list of float mimimum limit 
    maxval -- float o list of float maximum limit 
    """
    def __init__(self, val, minval=None, maxval=None):
        self.val = val
        
        if minval == None:
            if isinstance(val,float):
                self.minval = (val-DEFLIMIT)
            elif isinstance(val,list):
                #buil a list of dimension of gene
                self.minval = ([val-DEFLIMIT]*len(val))
        else:
            self.minval = minval
            
        
        if maxval == None:
            if isinstance(val,float):
                self.maxval = (val+DEFLIMIT)
            elif isinstance(val,list):
                #build a lst of dimension of gene
                self.maxval = ([val+DEFLIMIT]*len(val))
        else:
            self.maxval = maxval
        
 
