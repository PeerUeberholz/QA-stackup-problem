import stackingPallet

cases = [[[0,1],[1,0]], [[0,2,1],[1,0,2]], [[0,1,0,1], [1,1,0,0]],[[0,2],[1,1],[2,0]],[[0,2,1],[1,0,2],[1,2]], [[0,1,1],[1,0,1]], [[0,1,3,2],[3,1,0,2]]]

for case in cases:
    stackingPallet.solveDWave(case, 10000)
