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





