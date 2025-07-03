"""Networkx common operations."""
#    Copyright (C) 2011-2012 by
#    Leonardo Bermeo <lbermeoc@unal.edu.co>
#    Joao Carlos Basilio <basilio@poli.ufrj.br>
#    GNU license.
import networkx as nx
from networkx import MultiDiGraph, DiGraph
from deslab.src.automatadefs import fsa


def strconncomps(G):  
    """This function returns the stronly connected components in the 
    graph of automaton G"""   
    if isinstance(G, MultiDiGraph): 
        sccs = nx.algorithms.components.strongly_connected.strongly_connected_components(G)
    if isinstance(G, DiGraph): 
        sccs = nx.algorithms.components.strongly_connected.strongly_connected_components(G)
    elif isinstance(G, fsa):
        sccs = nx.algorithms.components.strongly_connected.strongly_connected_components(G.Graph)
    else:
        raise invalidArgument('G must be an automaton or a graph')
    return sccs

def selfloopnodes(G):
    """ This function looks for nodes with self loops in the graph of 
    automaton G """    
    if isinstance(G, MultiDiGraph): 
        selfloopnodes = G.nodes_with_selfloops() 
    elif isinstance(G, DiGraph): 
        selfloopnodes = G.nodes_with_selfloops() 
    elif isinstance(G, fsa):  
        selfloopnodes = G.Graph.nodes_with_selfloops()   
    else:
        raise invalidArgument('G must be an automaton or a graph')
    return selfloopnodes

def condensation(G):
    """ Returns the condensation graph corresponding to the
    subjacent graph of automaton G """ 

    if isinstance(G, MultiDiGraph): 
        C = nx.condensation(G) 
    elif isinstance(G, DiGraph): 
        C = nx.condensation(G) 
    elif isinstance(G, fsa):  
        C =  nx.condensation(G.Graph) 
    return C



