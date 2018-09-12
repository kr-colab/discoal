from io import StringIO
from Bio import Phylo
from msTools import msOutToTreeStringArrays
import sys
import numpy as np

trees, pos = msOutToTreeStringArrays(sys.argv[1])
obs = np.zeros(len(trees))
for i in range(len(trees)):
    times = []
    posSum = np.sum(pos[i])
    for j in range(len(trees[i])):
        tree = Phylo.read(StringIO(trees[i][j]), "newick")
        t_tot = tree.total_branch_length()
        times.append(t_tot*(pos[i][j]/posSum))
    times = np.array(times)
    obs[i] = np.sum(times)
print(np.mean(obs))
