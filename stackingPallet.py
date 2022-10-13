import dimod
import math
from dwave.system import EmbeddingComposite, DWaveSampler
import time
import argparse
import sys
import dwave.inspector

from neal.sampler import SimulatedAnnealingSampler
from qaUtils import saveSampleset

def iterN(items, n):
    """! Generator that iterates over a given collection in slices of size n.
    If len(items) is not divisible by n the last slice returned will contain
    the remaining elements.

    @params items The value to slice
    @params n Size of slices
    """
    while len(items) > n :
        yield items[:n]
        items = items[n:]
    if len(items) > 0:
        yield items

class PalletQUBOGenerator:
    def constructSequenceGraph(this):
        """!
          \brief Constructs the sequence Graph for the problem

          The sequence graph contains a directed edge between labels iff
          the label at the start of the edge must be opened before the label
          at the end of the edge can be closed.
        """
        this.sequenceGraph = set() 

        for sequence in this.sequences:
            earlierLabels = set()
            for label in sequence:
                for eLabel in earlierLabels:
                    if eLabel != label:
                        this.sequenceGraph.add((eLabel,label))
                earlierLabels.add(label)
        
        #Convert the sequenceGraph to list for conistent ordering
        this.sequenceGraph = [edge for edge in this.sequenceGraph] 

    def __init__(this, sequences, autoGenerate=True, penaltyMul=50):
        """!
          Constructs a generator for pallet-solution bqms
        
          \param sequences List of sequences to stack from. The sequences are ordered lists of labels.
          \param autoGenerate Whether to immediately generate the full bqm during construction
        """
        this.sequences = sequences

        labels = set()
        for sequence in this.sequences:
            for label in sequence:
                labels.add(label)
        this.numLabels = len(labels) #Number of different labels
 
        this.bqm = dimod.BinaryQuadraticModel(vartype='BINARY')

        this.auxSize = math.floor(math.log(this.numLabels,2))+1#Number of auxiliary variables to hold a number
       
        this.penaltyFactor = pow(2,this.auxSize)*penaltyMul #Penalty larger than maximum possible p

        if autoGenerate:
            this.generateBQM()
    
    def varName(this, i, j):
        """!
          \brief Returns the BQM-Variable name for the plan variable with the given indices
          
          \param i Index i
          \param j Index j
        """
        return 'x(' + str(i) + ',' + str(j) +')'

    def permutationConstraint(this):
        """! 
        \brief Models the constraint which ensures that each position
        is used by exactly one label and each label uses exactly one
        position"""
        for k in range(0, this.numLabels):
            for i in range(0, this.numLabels):
                iName = this.varName(k, i)
                invIName = this.varName(i,k)

                this.bqm.add_variable(iName, -this.penaltyFactor)
                this.bqm.add_variable(invIName, -this.penaltyFactor)
                for j in range(i+1, this.numLabels):
                    jName = this.varName(k,j)
                    invJName = this.varName(j,k)

                    this.bqm.add_interaction(iName, jName, 2*this.penaltyFactor)
                    this.bqm.add_interaction(invIName, invJName,2*this.penaltyFactor)

        this.bqm.offset = 2*this.penaltyFactor*this.numLabels
    
    def modelOr(this, left, right, auxName):
        """!
          \brief Models the boolean expression auxName = left OR right in the BQM
          
          \param left One of the variables of the expression
          \param right One of the variables of the expression
          \param auxName Name of the auxiliary variable which holds the result of the expression
        """
        this.bqm.add_variable(left, this.penaltyFactor)
        this.bqm.add_variable(right, this.penaltyFactor)
        this.bqm.add_variable(auxName, this.penaltyFactor)
        this.bqm.add_interaction(left,right,this.penaltyFactor)
        this.bqm.add_interaction(left, auxName, -2*this.penaltyFactor)
        this.bqm.add_interaction(right, auxName, -2*this.penaltyFactor)

    def yName(this, j, c):
        """!
          \brief Returns name for auxiliary variable holding value of Y(j,c)
        
          \param j Value of j
          \param c Value of c
        """
        return 'Y(' + str(j) + ',' + str(c) + ')'

    def y(this,j,c):
        """! 
            \brief Model Y(j,c) for the given j and c.
        
            \param j Value of j
            \param c Value of c
        """
        test  = 0
        #Y(j,c) = Y(j,c+1) OR (verodert alle Konjunktionen von i, i')
        j2 = c+1
        conjunctions = []
        for edge in this.sequenceGraph:
            left = this.varName(edge[1],j)
            right = this.varName(edge[0],j2)
            auxName = left + 'and' + right
            #AND Bedingung
            if not this.bqm.has_variable(auxName):
                this.bqm.add_interaction(left, right, this.penaltyFactor)
                this.bqm.add_interaction(left, auxName, -2*this.penaltyFactor)
                this.bqm.add_interaction(right, auxName, -2*this.penaltyFactor)
                this.bqm.add_variable(auxName, 3*this.penaltyFactor)

            conjunctions.append(auxName)
            test += 1

        while(len(conjunctions) > 1):
            pair = []
            left = conjunctions.pop(0)
            right = conjunctions.pop(0)

            auxName = '(' + left +')or('+right + ')'
            test += 1
            if len(conjunctions) == 0 and c==(this.numLabels-1):
                auxName = this.yName(j,c)

            if not this.bqm.has_variable(auxName):
                this.modelOr(left, right, auxName)
                
            conjunctions.append(auxName)
             
        #The expression can be modeled recursively 
        #since it grows longer with smaller c but the edges don't change
        if c < this.numLabels-2:
            left = conjunctions[0]
            right = this.yName(j, c+1)
            auxName = this.yName(j,c)
            test += 1

            this.modelOr(left,right,auxName)
        else:
            this.bqm.relabel_variables({conjunctions[0]:this.yName(j,c)})
        

    def yjc(this):
        """!
          \brief Generates all relevant expressions for Y(j,c)
        """
        #Reversed to avoid having to combine variables 
        #created by recursion
        #(relabeling with an exisiting name is not permitted)
        for c in reversed(range(0, this.numLabels-1)):
            for j in range(0, c+1):
                this.y(j,c)

    def squareAux(this, auxName, factor=1):
        """! 
          \brief Add expression to represent the square of an auxiliary variable,
        which is a natural number represented by multiple qubits
        in binary notation.

        \param auxName Name of the number to square
        \param factor Factor to multiply the squared variable by
        """
        auxName+='_'
        for i in range(0, this.auxSize):
            this.bqm.add_variable(auxName+str(i), (pow(2, i)**2)*factor)
            for j in range(i+1, this.auxSize):
                this.bqm.add_interaction(auxName+str(i), auxName+str(j), pow(2,i+j+1)*factor)

   
    def inequalityConstraints(this):
        """! Models the necessary inequalities to set w to the correct value"""
        for c in range(0, this.numLabels-1):
            for j in range(0, c+1): 
                this.bqm.add_variable(this.yName(j,c), this.penaltyFactor)
                for j2 in range(j+1, c+1):
                    this.bqm.add_interaction(this.yName(j,c), this.yName(j2,c), 2*this.penaltyFactor)

                for i in range(0, this.auxSize):
                    this.bqm.add_interaction(this.yName(j,c),'s'+str(c)+'_'+str(i), this.penaltyFactor*2*pow(2,i))
                    this.bqm.add_interaction(this.yName(j,c),'w_'+str(i), -this.penaltyFactor*2*pow(2,i))

            for i in range(0, this.auxSize):
                for j in range(0, this.auxSize):
                    this.bqm.add_interaction('s'+str(c)+'_'+str(i), 'w_'+str(j), -this.penaltyFactor*2*pow(2,i)*pow(2,j))

            this.squareAux('w', this.penaltyFactor)
            this.squareAux('s'+str(c), this.penaltyFactor)


    def generateBQM(this):
        """!
          \brief Performs all neccessary steps to fully model the problem
        """
        this.constructSequenceGraph()
        this.permutationConstraint()
        this.yjc()
        this.inequalityConstraints()
 
        for i in range(0, this.auxSize):
            this.bqm.add_variable('w_'+str(i), pow(2,i))
    
    def breakDownVariables(this):
        """!
          \brief Prints information about variable usage to console
        """
        remaining = len(this.bqm)
        plan = this.numLabels**2
        remaining -= plan
        print('Number of plan variables:', plan)
        numbers = this.auxSize*this.numLabels
        remaining -= numbers
        print('Number of variables that model numbers:', numbers)
        print('Number of variables that model boolean expressions:', remaining)
    
    def getMaxBias(this):
        """!
          \brief Returns the biggest coupler bias
        """
        maxKey = ''
        maxBias = 0
        for key1, key2 in this.bqm.iter_interactions():
            bias =  abs(this.bqm.get_quadratic(key1, key2))
            if bias > maxBias:
                    maxKey = key1+','+key2 
                    maxBias = bias

        for key in this.bqm.iter_variables():
            bias = abs(this.bqm.get_linear(key))
            if bias > maxBias:
                maxKey = key
                maxBias = bias

        print('Max bias is',maxBias,'at',maxKey)

        return maxBias

    def interpretSample(this, sample):
        """!
          Interprets the solution described by the given sample

          \param sample The sample to examine 
        """
        if sample.energy > (pow(2,this.auxSize)-1):
            print('WARNING: There appear to be violated constraints in the given sample,\
making the solution invalid!')

        print('The pallets are opened in this order:')
        for j in range(0, this.numLabels):
            for i in range(0, this.numLabels):
                if sample.sample[this.varName(i,j)] == 1:
                    print( str(j+1)+'.', i)
        print('The number of stacking places required is (according to the sample)', sample.energy+1)


def solveDWave(sequences, num_reads, penaltyMul=50, **args):
    """! 
    \brief Approximate a solutions of the Stacking Problem with the given sequences
    using a DWave Quantum Annealer

    \param sequences The sequences of the problem instance
    \param num_reads Number of samples to generate
    \param penaltyMul Value to mutiply the minimum possible penalty for violation of constraints by
    \param **args Additional keyword arguments are forwarded to DwaveSampler.sample()
    """

    test = PalletQUBOGenerator(sequences, penaltyMul = penaltyMul)
    print("Generated bqm")
    print("Number of Variables: ", len(test.bqm))
   
    sampler = EmbeddingComposite(DWaveSampler())

    # parameter auto_scale=true, ist default, skaliert alle Größen in das Intervall [-1, +1]
    # Parameter chain_strength=chain_strength_value könnte was helfen
    sampleset = sampler.sample(test.bqm, num_reads=num_reads,  return_embedding=True,warnings='save', **args)#PARAMETERS HERE
    dwave.inspector.show(sampleset) 
    sampleset.info['bqm'] = test.bqm
    sampleset.info['sequences'] = sequences
    sampleset.info['penaltyFactor'] = test.penaltyFactor
    saveSampleset(sampleset, "data/pallet/QA-")
    #print(sampleset)
    
    print('Lowest energy:', sampleset.first.energy)
    test.interpretSample(sampleset.first)
    test.breakDownVariables()
    return sampleset

def solveSimAnneal(sequences,num_reads, penaltyMul=50, **args):
    """! 

    \brief Approximate a solution of the Stacking Problem with the given sequences
        using Simulated Annealing with a QUBO-Formulation of the Energy Function

    \param sequences The sequences of the problem instance
    \param num_reads Number of samples to generate
    \param penaltyMul Value to mutiply the minimum possible penalty for violation of constraints by
    \param **args Additional keyword arguments are forwarded to SimulatedAnnealingSampler.sample()
    """

    test = PalletQUBOGenerator(sequences, penaltyMul = penaltyMul)
    print("Generated bqm")
    print("Number of variables: ", len(test.bqm))

    sampler = SimulatedAnnealingSampler()
    start = time.time()
    sampleset = sampler.sample(test.bqm, num_reads=num_reads, **args)
    end = time.time()
    sampleset.info['bqm'] = test.bqm
    sampleset.info['sequences'] = sequences
    sampleset.info['penaltyFactor'] = test.penaltyFactor
    saveSampleset(sampleset, "data/pallet/SA-")

    print('Lowest energy:', sampleset.first.energy)
    test.interpretSample(sampleset.first)
    test.breakDownVariables()

    return [end - start, sampleset, test]


def parseSequences(text):
    parts = text.split('-')
    return [[int(x) for x in part.split(',')] for part in parts]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Solve the stacking Problem using a Quantum Annealer or Simulated Annealing')
    
    requiredNamed = parser.add_argument_group('required arguments')

    requiredNamed.add_argument('-s', type=str, action='store', dest='seqs', 
            metavar='Sequences. Entries are separated by commas. Sequences are\
 separated by -.Labels are numbers', required = True)
    requiredNamed.add_argument('-m', type=str, action='store', dest='method', metavar='Method to use. Either SA or QA.', required = True)
    requiredNamed.add_argument('-nr', type=int, action='store', dest='num_reads', metavar='Number of samples to generate.', required = True)

    parser.add_argument('-p', type=int, action='store', dest='penalty', metavar='Factor to multiply lowest possible penalty A by', default = 50)

    args = parser.parse_args(sys.argv[1:])
    sequences = parseSequences(args.seqs)
    print("Solving instance " + str(sequences))
    
    if args.method == 'SA':
        solveSimAnneal(sequences, args.num_reads, args.penalty)
    elif args.method == 'QA':
        solveDWave(sequences, args.num_reads, args.penalty)
    else:
        print('Method (-m) must be either SA or QA!') 
