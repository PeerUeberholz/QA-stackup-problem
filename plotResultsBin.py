import numpy as np
import qaUtils
from collectConstStats import calcConstraintStats
from matplotlib import pyplot as plt

stackedWidth = 0.6
groupedWidth = stackedWidth/4

files = [('data/2Lab2Bin1Dec.dat',1), ('data/3Seq2Lab3Bin1Dec.dat',1), ('data/3Lab2Bin1Dec.dat',1), ('data/3Lab2Bin2Dec.dat',2), ('data/2Seq2Lab4Bin1Dec.dat',1), ('data/3Seq3Lab2Bin1Dec.dat',1), ('data/3Seq3Lab2Bin2Dec.dat',2), ('data/3Seq3Lab8BinTot1Dec.dat',1), ('data/3Seq3Lab8BinTot2Dec.dat',2)]
#labels = list((np.arange(len(files))+1).astype(str))
labels = ['1,1','2,1','3,1','3,2','4,1','5,1','5,2','6,1','6,2']
x = np.arange(len(files))

correct = []
incorrect = []
permut = []
seq = []
ftc = []
count = []


for instance in files:
    ss = qaUtils.loadSampleset(instance[0])
    stats = calcConstraintStats(ss, instance[1])

    print('==========')
    print(ss.info['sequences'])
    print(stats)
    print(len(ss.info['bqm']))

    permut.append(stats['Permutation'])
    seq.append(stats['SequenceOrder'])
    ftc.append(stats['f(t,c)'])
    count.append(stats['Count'])

    correctCount = np.sum(ss.record[ss.record['energy'] < 10]['num_occurrences'])
    print(correctCount)
    correct.append(correctCount)
    incorrect.append(10000-correctCount)

fig, ax = plt.subplots()
ax.bar(labels, incorrect, stackedWidth+.1, label='Invalid', color='tab:red')
ax.bar(labels, correct, stackedWidth+.1, bottom=incorrect, label='Valid', color='tab:olive')

ax.bar(x-groupedWidth*1.5, permut, groupedWidth, label='PERMUTATION', color='aqua')
ax.bar(x-groupedWidth*0.5, seq, groupedWidth, label='SEQUENCE_ORDER', color='darkturquoise')
ax.bar(x+groupedWidth*0.5, ftc, groupedWidth, label='f(t,c)', color='cadetblue')
ax.bar(x+groupedWidth*1.5, count, groupedWidth, label='Inequality', color='magenta')

plt.xlabel("Instance ID, Value of k")
plt.ylabel("Number of Samples")

plt.legend(bbox_to_anchor=(.95, 1), loc='upper left', fontsize='small')

plt.show()
