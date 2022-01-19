"""Test of offline_hetesim.py."""

import sys
import pandas as pd

# sys.path.insert(0, 'path_to_semnet')
from semnet.offline import HetGraph
from semnet.offline_hetesim import hetesim, hetesim_all_metapaths, mean_hetesim_scores, approximate_mean_hetesim_scores

def test_hetesim(graph, mp, true_value):
    hs = hetesim(graph, ['s'], ['t'], [mp])[str(mp)]['s']['t']
    print("Computed hetesim: " + str(hs))
    print("True hetesim: " + str(true_value))

    assert( true_value - hs < 0.001 and hs - true_value < 0.001)


def test_hetesim_all_metapaths(graph, path_len, metapath, true_hs_value):

    assert(abs(hetesim_all_metapaths(graph, ['s'],['t'], path_len)[str(metapath)]['s']['t'] -true_hs_value)< 0.001)


def test_mean_hetesim_scores(graph, path_len, true_mean_hs_value):

    assert(abs( mean_hetesim_scores (graph, ['s'], 't', path_len)['s'] - true_mean_hs_value) < 0.01)


def test_approximate_mean_hetesim_scores(graph, path_len, true_mean_hs_value):
    approx_mean_hs = approximate_mean_hetesim_scores(graph, ['s'],'t',path_len, 0.05, 0.95)['s']
    print("approx mean hs: " + str(approx_mean_hs))
    print("True mean hs: " + str(true_mean_hs_value))


if __name__ == '__main__':

    toy_graph_1_df = pd.read_csv('toy_graph_1.tsv', sep="\t", header=0)
    toy_graph_1 = HetGraph(toy_graph_1_df.to_dict(orient='records'))

    toy_graph_2_df = pd.read_csv('toy_graph_2.tsv', sep="\t", header=0)
    toy_graph_2 = HetGraph(toy_graph_2_df.to_dict(orient='records'))

    toy_graph_3_df = pd.read_csv('toy_graph_3.tsv', sep="\t", header=0)
    toy_graph_3 = HetGraph(toy_graph_3_df.to_dict(orient='records'))

    toy_graph_4_df = pd.read_csv('toy_graph_4.tsv', sep="\t", header=0)
    toy_graph_4 = HetGraph(toy_graph_4_df.to_dict(orient='records'))

    mp1 = ['t1', 'r1', 't2', 'r2', 't3', 'r3', 't1', 'r1', 't4']
    mp2 = ['t1', 'r1', 't2', 'r1', 't3', 'r1', 't4', 'r1', 't5', 'r1', 't6', 'r1', 't7']
    mp3 = ['t1', 'r1', 't2', 'r1', 't1', 'r2', 't1', 'r3', 't1']
    
    test_hetesim(toy_graph_1, mp1, 0.5774)
    test_hetesim(toy_graph_2, mp2, 0.8437)
    test_hetesim(toy_graph_3, mp3, 0.8333)

    test_hetesim_all_metapaths(toy_graph_1, 4, mp1, 0.5774)

    test_mean_hetesim_scores(toy_graph_4, 4, 0.6007)
    test_approximate_mean_hetesim_scores(toy_graph_4, 4, 0.6007)

    print('Testing complete.')