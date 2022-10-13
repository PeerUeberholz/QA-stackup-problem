from datetime import datetime
import pickle
import dimod
import networkx as nx
import matplotlib.pyplot as plt
import math

def varNameToLatex(name):
    """!
      \brief Converts the given BQM variable name to LaTeX notation
    """
    if name[0] == 'f' or name[0] == 'Y':
        firstLetter = name[0]
        firstIndex = str(int(name[name.find('(')+1:name.find(',')])+1)
        secondIndex = str(int(name[name.find(',')+1:name.find(')')])+1)
        return '$'+firstLetter+'('+firstIndex + ',' + secondIndex + ')$'
    elif name[0] != 'x':
        firstIndexStr = name[1:name.find('_')]
        firstIndex = ''
        if len(firstIndexStr) > 0:
            firstIndex = str(int(firstIndexStr)+1)

        secondIndex = str(int(name[name.find('_')+1:])+1)
        return '$'+name[0]+firstIndex+'_'+secondIndex+'$'

    res = '$'
    
    while len(name)>0:
        name = name[2:]
        firstIndex = str(int(name[:name.find(',')])+1)
        secondIndex = str(int(name[name.find(',')+1:name.find(')')])+1)

        res += 'x_{' + firstIndex+ '}^{' + secondIndex+ '}'
        name = name[name.find(')')+1:]
        
        if name[:2] == 'or':
            name = name[2:]
            res += '\\lor '

        if name[:3] == 'and':
            name = name[3:]
            res += '\\land '

    return res + '$'

def outputbqm(path, bqm, sep=',', latexMode=False):
    out = open(path, 'w')
    out.write(' ')
    for key in bqm:
        outKey = key
        if latexMode:
            outKey = varNameToLatex(key)
        out.write(sep+str(outKey))
    
    if latexMode:
        out.write('\\\\')

    out.write('\n')
    
    rowCounter = 0
    for rowKey in bqm:
        rowOutKey = rowKey
        if latexMode:
            rowOutKey = varNameToLatex(rowKey)

        out.write(rowOutKey)
        colCounter = 0
        for colKey in bqm:
            field = 0
            if colKey == rowKey:
                field = bqm.linear[rowKey]
            elif colCounter < rowCounter:
                field = ' '
            elif (rowKey, colKey) not in bqm.quadratic:
                field = 0
            else:
                field = bqm.adj[rowKey][colKey]
            
            if not isinstance(field, str) and field%1 == 0:
                field = math.floor(field)

            out.write(sep+str(field))
            colCounter += 1
            print(rowCounter, colCounter)
        if latexMode:
            out.write('\\\\')

        out.write('\n')
        rowCounter += 1

def saveSampleset(sampleset, prefix="", timestamp=True):
    """!Saves the sampleset with the given prefix and a timestamp
    @param sampleset The sampleset to save
    @param prefix The prefix of the path
    @param timestamp Whether to add a timestamp to the filename"""

    now = datetime.now()
    
    path = prefix+now.strftime("%Y-%m-%d-%H-%M-%S")+'.dat'

    pickle.dump(sampleset.to_serializable(), open(path, 'wb'))

def loadSampleset(path):
    """!Loads and returns a serialized sampleset at the given location"""
    return dimod.SampleSet.from_serializable(pickle.load(open(path, 'rb')))

def extractNeighborhood(bqm, var):
    """!Generates a new BQM that only contains the node var and neighboring nodes

    @param bqm The BinaryQuadraticModel to extract the neighborhood from
    @param var The variable to extract"""

    res = dimod.BinaryQuadraticModel(vartype='BINARY')
    knownKeys = [var]
    
    res.add_variable(var, bqm.linear[var])
    for key in bqm.adj[var]:
        knownKeys.append(key)
        res.add_variable(key, bqm.linear[key])
        res.add_interaction(var, key, bqm.quadratic[var,key])

        for innerKey in bqm.adj[key]:
            if innerKey in knownKeys:
                res.add_interaction(key, innerKey, bqm.quadratic[key, innerKey])

    return res

def drawNeighborhood(path, var):
    sampleS = loadSampleset(path)
    bqm = sampleS.info['bqm']
    part = extractNeighborhood(bqm, var)
    graph = part.to_networkx_graph()
    
    cm = []
    for key in part.variables:
        if key == var:
            cm.append('red')
        else:
            cm.append('blue')

    nx.draw(graph, with_labels=True, node_color=cm)
    plt.show()
