"""Test of randomized_hetesim.py."""

import sys
import pandas as pd

# sys.path.insert(0, 'path_to_semnet')
from semnet.offline import HetGraph
from semnet.randomized_hetesim import randomized_pruned_hetesim, restricted_random_walk_on_metapath, randomized_pruned_hetesim_all_metapaths, approximate_mean_pruned_hetesim

def test_restricted_random_walk_on_metapath(tg1):
    mp1 = ['t1', 'r1', 't2', 'r2', 't3', 'r3', 't1', 'r1', 't4']
    mp1l = mp1[:5]
    bad_nodes = [set() for i in range(2)]
    for _ in range(100):
        print(restricted_random_walk_on_metapath(tg1, 's', mp1l , bad_nodes))
    

def test_randomized_pruned_hetesim(graph, mp, epsilon, k, r, true_value, filename, N, plot_title):
    results = []

    for i_ in range(N):
       results.append(randomized_pruned_hetesim(graph, ['s'], ['t'], [mp], k, epsilon, r)[str(mp)]['s']['t'])

    results_df = pd.DataFrame(results)
    results_df.columns = ["approximate pruned hetesim"]

    num_within_epsilon = len([x for x in results if true_value - epsilon <= x and true_value + epsilon >= x])
    percent_within_epsilon = num_within_epsilon / N * 100

    print("Of " + str(N) + " iterations, " + str(num_within_epsilon) + " (" + str(percent_within_epsilon) + "% ) had error less than epsilon.")

    
def test_randomized_pruned_hetesim_all_metapaths(graph, mp, path_len, epsilon, r, true_value, N):
    results = []

    for _ in range(N):
       results.append(randomized_pruned_hetesim_all_metapaths(graph, ['s'], ['t'], path_len, epsilon, r)[str(mp)]['s']['t'])

    num_within_epsilon = len([x for x in results if true_value - epsilon <= x and true_value + epsilon >= x])
    percent_within_epsilon = num_within_epsilon / N * 100

    print("Of " + str(N) + " iterations, " + str(num_within_epsilon) + " (" + str(percent_within_epsilon) + "% ) had error less than epsilon.")


def test_approximate_mean_pruned_hetesim(graph, path_len, epsilon, r, true_value, N):
    results = []

    for _ in range(N):
        approx_mean_hs =approximate_mean_pruned_hetesim(graph, ['s'], 't',  path_len, epsilon, r)['s']
        results.append(approx_mean_hs)
        print("Approximate mean hs is " + str(approx_mean_hs))

    num_within_epsilon = len([x for x in results if true_value - epsilon <= x and true_value + epsilon >= x])
    percent_within_epsilon = num_within_epsilon / N * 100

    print("Of " + str(N) + " iterations, " + str(num_within_epsilon) + " (" + str(percent_within_epsilon) + "% ) had error less than epsilon.")


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

    test_randomized_pruned_hetesim(toy_graph_1, mp1, 0.05, 4, 0.95, 0.5774, "toy_graph_1_test", 1, "Computed approximate pruned HeteSim values for toy graph 1")
    test_randomized_pruned_hetesim(toy_graph_2, mp2, 0.05, 3, 0.95, 0.8944, "toy_graph_2_test", 1, "Computed approximate pruned HeteSim values for toy graph 2")
    test_randomized_pruned_hetesim(toy_graph_3, mp3, 0.05, 3, 0.95, 0.8333, "toy_graph_3_test", 1, "Computed approximate pruned HeteSim values for toy graph 3")

    print('Testing complete.')