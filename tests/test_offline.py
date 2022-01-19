"""Test of offline.py."""

import sys
import pandas as pd 

# sys.path.insert(0, 'path_to_semnet')
from semnet.offline import HetGraph

# Create toy graph with 6 nodes, 2 node types, 2 relation types
toy_graph = pd.DataFrame([['a','t1','b','t1','r1',1],
                          ['a','t1','c','t2','r1',1],
                          ['c','t2','d','t1','r1',1],
                          ['b','t1','e','t2','r1',1],
                          ['a','t1','e','t2','r1',1],
                          ['c','t2','f','t2','r1',1],
                          ['d','t1','b','t1','r2',1],
                          ['e','t2','c','t2','r2',1],
                          ['a','t1','f','t2','r2',1],
                          ['f','t2','d','t1','r2',1],
                          ['b','t1','a','t1','r2',1],
                          ['a','t1','a','t1','r1',1]], 
                         columns=['start_node',
                                  'start_type',
                                  'end_node',
                                  'end_type',
                                  'relation',
                                  'weight'])

def test_fan_out_0(hg):
    node, path = next(hg._fan_out('a', depth=0))

    assert node == {'t1':{'a'}}
    assert path == []


def test_schema_fan_out_0(hg):
    node, path = next(hg._schema_fan_out('t1', depth=0))

    assert node == {'t1'}
    assert path == []


def test_fan_out_1(hg):
    nbhrs = [(hg._path_to_string(x[1]), x[0]) for x in hg._fan_out('a', depth=1)]
    nbhr_nodes = set([])

    for path, node_dict in nbhrs:
        for nodes in node_dict.values():
            nbhr_nodes |= nodes

    assert nbhr_nodes == {'a','b','c','e','f'}

def test_schema_fan_out_1(hg):
    nbhrs = [(hg._path_to_string(x[1]), x[0]) for x in hg._schema_fan_out('t1', depth=1)]

    nbhr_nodes = set([])
    for path, nodes in nbhrs:
        nbhr_nodes |= nodes

    assert nbhr_nodes == {'t1', 't2'}

def test_fan_out_2(hg):
    nbhrs = [(hg._path_to_string(x[1]), x[0]) for x in hg._fan_out('a', depth=2)]
    nbhr_nodes = set([])

    for path, node_dict in nbhrs:
        for nodes in node_dict.values():
            nbhr_nodes |= nodes

    assert nbhr_nodes == {'a','b','c','d','e','f'}

    
def test_schema_fan_out_2(hg):
    nbhrs = [(hg._path_to_string(x[1]), x[0]) for x in hg._schema_fan_out('t1', depth=2)]

    nbhr_nodes = set([])
    for path, node_set in nbhrs:
        nbhr_nodes |= node_set

    assert nbhr_nodes == {'t1','t2'}


def test_fan_in_0(hg):
    node, path = next(hg._fan_in('a', depth=0))

    assert node == {'t1':{'a'}}
    assert path == []


def test_fan_in_1(hg):
    nbhrs = [(x[0], hg._path_to_string(x[1])) for x in hg._fan_in('a', depth=1)]

    nbhr_nodes = set([])
    for node_dict, path in nbhrs:
        for nodes in node_dict.values():
            nbhr_nodes |= nodes

    assert nbhr_nodes == {'a','b','c','e','f'}


def test_fan_in_2(hg):
    nbhrs = [(x[0], hg._path_to_string(x[1])) for x in hg._fan_in('a', depth=2)]
    nbhr_nodes = set([])
    for node_dict, path in nbhrs:
        for nodes in node_dict.values():
            nbhr_nodes |= nodes

    assert nbhr_nodes == {'a','b','c','d','e','f'}


def test_merge_paths(hg):
    merged = hg._merge_paths('a',['b'], 'c')

    assert merged ==['a','b','c']


def test_fixed_length_paths(hg):
    paths = hg.compute_fixed_length_paths('a','c',length=2)
    path_strings = [hg._path_to_string(p) for p in paths]
    
    assert 'a->r1->e->r2->c' in path_strings
    assert 'a->r2->f->r1->c' in path_strings
    assert all([x.startswith('a') for x in path_strings])
    assert all([x.endswith('c') for x in path_strings])
    
def test_fixed_length_schema_walks(hg):
    paths = hg.compute_fixed_length_schema_walks('t1','t2',length=2)
    path_strings = [hg._path_to_string(p) for p in paths]

    assert 't1->r1->t1->r2->t2' in path_strings

def test_compute_metapath_reachable_nodes(hg):
    mp = ['t1', 'r1', 't1']
    reachable_nodes = hg.compute_metapath_reachable_nodes('a', mp)

    assert 'a' in reachable_nodes
    assert 'b' in reachable_nodes
    assert not 'd' in reachable_nodes
    assert not 'c' in reachable_nodes
    
    mp2 = ['t1', 'r1', 't2', 'r1', 't1'] 
    reachable_nodes_2 = hg.compute_metapath_reachable_nodes('a', mp2)

    assert 'a' in reachable_nodes_2
    assert 'b' in reachable_nodes_2
    assert 'd' in reachable_nodes_2
    assert not 'e' in reachable_nodes_2 


def test_compute_fixed_length_metapaths(hg):
    mps = [hg._path_to_string(mp) for mp in hg.compute_fixed_length_metapaths('a', 'b', length=2)]

    assert 't1->r1->t2->r1->t1' in mps
    assert 't1->r1->t1->r1->t2' in mps
    

if __name__ == '__main__':

    edgelist = toy_graph.to_dict(orient='records')
    rel2inv = {'r1': 'r1', 'r2':'r2_inv'}

    test_graph = HetGraph(edgelist, rel2inv)

    test_fan_out_0(test_graph)
    test_fan_out_1(test_graph)
    test_fan_out_2(test_graph)

    test_fan_in_0(test_graph)
    test_fan_in_1(test_graph)
    test_fan_in_2(test_graph)

    test_fixed_length_paths(test_graph)

    test_compute_metapath_reachable_nodes(test_graph)

    test_schema_fan_out_0(test_graph)
    test_schema_fan_out_1(test_graph)    
    test_schema_fan_out_2(test_graph)

    test_fixed_length_schema_walks(test_graph)

    print('Testing complete.')