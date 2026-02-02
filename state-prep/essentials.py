import fusion_blossom as fb
import numpy as np
from joblib import Parallel, delayed
from pymatching import Matching
from tqdm import tqdm 
import networkx as nx

import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['text.color'] = 'k'
plt.rcParams['axes.labelcolor'] = 'k'
plt.rcParams['text.usetex'] = True

pal = ["003049","d62828","f77f00","fcbf49","eae2b7"]
pal = ['#' + i for i in pal]

marks = ['o', 's', '^', 'D', 'v', '>', '<', 'p', '*', 'h', 'H', '+', 'x']


import seaborn as sns


def repetition_stabilisers(N):
    check_matrix = np.zeros((N-1, N), dtype=int)
    for i in range(N-1):
        check_matrix[i, i] = 1
        check_matrix[i, i+1] = 1
    return check_matrix

def x_logicals(L):
    return np.eye(L,dtype = int) 



class RepetitionCode:

    """A customize repetition code with non-i.i.d. noise model"""
    def __init__(self, d, L, p, q):
        
        self.d = d # number of qubits
        self.p = p # qubit error rate
        self.L = L # number of repetitions
        self.q = q # measurement error rate 
        
        self.space_vertex_num = (self.d - 1) + 2
        self.time_vertex_num = (self.L + 2)
        self.num_check = self.d - 1
        self.tot_nodes = self.space_vertex_num*self.L + 2*self.num_check
        self.mid_nodes = self.space_vertex_num * self.L 
        self.nodes = self._initialize_nodes()
        
    def _initialize_nodes(self):
        
        """Helper method to initialize the nodes"""

        nodes = np.concatenate([
            [0], 
            np.arange(0, self.num_check, 1), 
            [0], 
            np.arange(self.num_check, self.mid_nodes + self.num_check, 1), 
            [0], 
            np.arange(self.mid_nodes + self.num_check, self.mid_nodes + 2 * (self.num_check)), 
            [0]
        ])
        
        nodes = np.array(nodes, dtype=int)
        nodes = nodes.reshape((self.time_vertex_num, self.space_vertex_num))
        nodes = nodes.T # grid layout for the node labelling 
        
        return nodes

    def _edge_labels(self):
        
        spacelike_labels = np.arange(0,self.L*self.d)
        spacelike_labels = spacelike_labels.reshape(self.L,self.d)
        
        timelike_labels = self.L*self.d + np.arange(0,(self.d-1)*(self.L+1))
        timelike_labels = timelike_labels.reshape(self.d-1,self.L+1)        
        return spacelike_labels.T, timelike_labels
        
        
    def get_initializer(self, max_half_weight=500):
        
        nodes = self.nodes
        
        virtual_vertices = np.hstack((nodes[1:-1,0],nodes[0,1:-1],nodes[-1,1:-1],nodes[1:-1,-1]))
        
        spacelike_weights = np.log((1-self.p)/self.p)
        timelike_weights = np.log((1-self.q)/self.q) 
        
        max_weight = max([spacelike_weights,timelike_weights])
               
        spacelike_half_weights = np.round(max_half_weight*(spacelike_weights/max_weight))
        
        spacelike_weighted_edges = []

        for i in range(1,self.L+1,1): # choose column
            for j in range(self.space_vertex_num-1): # choose row
                spacelike_weighted_edges.append((nodes[j,i],nodes[j+1,i],int(2 * spacelike_half_weights)))


        timelike_half_weights = np.round(max_half_weight*(timelike_weights/max_weight)) 

        timelike_weighted_edges = []

        for i in range(1,self.d,1): # rows
            for j in range(self.time_vertex_num-1): # column
                timelike_weighted_edges.append((nodes[i,j],nodes[i,j+1],int(2*timelike_half_weights)))

        weighted_edges = spacelike_weighted_edges + timelike_weighted_edges
        
        return fb.SolverInitializer(self.tot_nodes, weighted_edges, virtual_vertices)

    def get_positions(self):
        
        nodes = self._initialize_nodes()
        
        pos = []
        
        pos += [fb.VisualizePosition(i,0, 0) for i in range(1,self.d)]
        pos += [fb.VisualizePosition(i,j+1, 0) for j in range(nodes.shape[1]-2) for i in range(nodes.shape[0])]
        pos += [fb.VisualizePosition(i,self.L+1,0) for i in range(1,self.d)]
        
        return pos


