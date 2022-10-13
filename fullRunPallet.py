"""
    Convenience script to run through multiple plots and generate a plot and LaTeX table
    samplesets are still saved to data/pallet/QA-{TIMESTAMP}.dat
"""

import plotResultsPal as plotting
import stackingPallet

print("Running this script will run 10 instances using approx. 7.5 seconds of computation time(depending on parameters, num_reads etc")
input("Press Enter to continue")

additional_params = {} #Add additional parameters here
                       #e.g additional_params = {'anneal_schedule':somevalue}
instances = [
        [[0,1],[1,0]],
        [[0,1,1],[1,0,1]],
        [[0,2,1],[1,0,2]],
        [[0,1,0,1],[1,1,0,0]],
        [[0,2],[1,1],[2,0]],
        [[0,2,1],[1,0,2],[1,2]],
        [[0,2,1],[1,0,2],[1,2,0]],
        [[1,2,1,0],[1,0,2,0]],
        [[0,1,3,2],[3,1,0,2]]
]#Problem instances to solve

instanceIds = [1,2,3,4,5,6,7,8,9]#x-Axis labels and leftmost table column

num_reads = 10000
resultSamplesets = []

for instance in instances:
    print("Solving instance", instance)
    resultSamplesets.append(stackingPallet.solveDWave(instance, num_reads, **additional_params));

plotting.plotResults(resultSamplesets, instanceIds);
