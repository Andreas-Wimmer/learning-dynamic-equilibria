#In this code we just want to play around with the packages graph, network, dynamic flow, rightConstant, 
#Piecewise linear and network loader to see how it works before we implement regularized fictitious play 
#for general networks.

import core.graph as graph
import core.network as network
import dynamic_flow as df
import network_loader as nl
import utilities.piecewise_linear as pl
import utilities.right_constant as rc

test_graph = graph.DirectedGraph()
s = graph.Node(1,test_graph)
v = graph.Node(2, test_graph)
t = graph.Node(3, test_graph)

test_graph.nodes = [s, v, t]
test_network = network.Network()
test_network.graph = test_graph

test_network.add_edge(1, 2, 1, 1)
test_network.add_edge(1, 2, 0, 3)
test_network.add_edge(2, 3, 0, 2)

net_inflow = rc.RightConstant([0,2],[2])
test_network.add_commodity((1, net_inflow), 1)



