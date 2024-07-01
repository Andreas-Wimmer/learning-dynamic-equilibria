#In this code we just want to play around with the packages graph, network, dynamic flow, rightConstant, 
#Piecewise linear and network loader to see how it works before we implement regularized fictitious play 
#for general networks.

from __future__ import annotations

from graph import *
from network import Network
from dynamic_flow import DynamicFlow 
from network_loader import NetworkLoader
from piecewise_linear import PiecewiseLinear
from right_constant import RightConstant

test_graph = DirectedGraph()
s = Node(1,test_graph)
v = Node(2, test_graph)
t = Node(3, test_graph)

test_graph.nodes = {1:s, 2:v, 3:t}
test_network = Network()
test_network.graph = test_graph

test_network.add_edge(1, 2, 1, 1)
test_network.add_edge(1, 2, 0, 3)
test_network.add_edge(2, 3, 0, 2)

net_inflow = RightConstant([0,2],[2,0],(0, float("inf")))
test_network.add_commodity({1: net_inflow}, 1, 1)

edge_1 = Edge(s, v, 1, test_graph)
edge_2 = Edge(s, v, 2, test_graph)
edge_3 = Edge(v, t, 3, test_graph)

s.outgoing_edges = [edge_1, edge_2]
v.outgoing_edges = [edge_3]
v.incoming_edges = [edge_1, edge_2]
t.incoming_edges = [edge_3]

path_1 = [edge_1, edge_3]
path_2 = [edge_2, edge_3]

inflow_1 = RightConstant([0,0.5,1],[1.5,0.5,0],(0,2))
inflow_2 = RightConstant([0,0.5,1],[0.5,1.5,2],(0,2))

loader = NetworkLoader(test_network, [(path_1, inflow_1),(path_2, inflow_2)])
result = loader.build_flow()
print("Hello World")






