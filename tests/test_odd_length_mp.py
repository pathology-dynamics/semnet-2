import sys
import pandas as pd

# sys.path.insert(0, 'path_to_semnet')
from semnet.offline import HetGraph
from semnet.offline_hetesim import find_all_metapaths, hetesim, mean_hetesim_scores

toy_graph_odd = pd.read_csv('toy_graph_odd.csv', sep=",", header=0)

edgelist = toy_graph_odd.to_dict(orient='records')
rel2inv = {'R1': 'R1_inv', 'R2':'R2_inv'}

hg = HetGraph(edgelist, rel2inv)

mps = find_all_metapaths(hg, ['S'], ['T'], 3)
print(mps)

hs = hetesim(hg, ['S'], ['T'], mps)
print(hs)

mhs = mean_hetesim_scores(hg, ['S'], 'T', 3)
print(mhs)

print('Testing done.')