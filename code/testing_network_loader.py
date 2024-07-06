#In this code we just want to play around with the packages graph, network, dynamic flow, rightConstant, 
#Piecewise linear and network loader to see how it works before we implement regularized fictitious play 
#for general networks.

from __future__ import annotations

from graph import *
from network import Network
from dynamic_flow import DynamicFlow 
from network_loader import NetworkLoader
from piecewise_linear import PiecewiseLinear, identity
from right_constant import RightConstant
from machine_precision import eps


test_graph = DirectedGraph()
s = Node(0, test_graph)
v = Node(1, test_graph)
t = Node(2, test_graph)

test_graph.nodes = {0:s, 1:v, 2:t}
test_network = Network()
test_network.graph = test_graph

test_network.add_edge(0, 1, 1, 1)
test_network.add_edge(0, 1, 0, 3)
test_network.add_edge(1, 2, 0, 2)

net_inflow = RightConstant([0,1,1.75],[2.5,1,3],(0, float("inf")))
test_network.add_commodity({1: net_inflow}, 0, 1)

edge_1 = Edge(s, v, 0, test_graph)
edge_2 = Edge(s, v, 1, test_graph)
edge_3 = Edge(v, t, 2, test_graph)

s.outgoing_edges = [edge_1, edge_2]
v.outgoing_edges = [edge_3]
v.incoming_edges = [edge_1, edge_2]
t.incoming_edges = [edge_3]

path_1 = [edge_1, edge_3]
path_2 = [edge_2, edge_3]

inflow_1 = RightConstant([0,0.5,1],[1.5,0.5,0],(0,2))
inflow_2 = RightConstant([0,0.5,1,1.75],[1,2,1,3],(0,2))
inflow_3 = RightConstant([0,1],[1,0],(0,2))
inflow_4 = RightConstant([0,1,1.75],[1.5,1,3],(0,2))

test_network.add_path(path_1)
test_network.add_path(path_2)

loader_1 = NetworkLoader(test_network, [(path_1, inflow_1),(path_2, inflow_2)])
loader_2 = NetworkLoader(test_network, [(path_1, inflow_3),(path_2, inflow_4)])
result_1 = loader_1.build_flow()
result_2 = loader_2.build_flow()
flow_1 = next(result_1)
flow_2 = next(result_2)
queues_1 = flow_1.get_queues()
arrivals_1 = loader_1.expected_arr()
arrivals_2 = loader_2.expected_arr()
delays_1 = loader_1.path_delay()
delays_2 = loader_2.path_delay()
delays = [delays_1, delays_2]
diff_delay = [delays_1[0] - delays_2[0], delays_1[1] - delays_2[1]]
diff_inflow = [inflow_1 - inflow_3, inflow_2 - inflow_4]

product_test = diff_inflow[1].multiply(diff_delay[1], 0, 0.5)

steps = []
for i in range(2):
    steps.append([])
for i in range(len(diff_delay[0].times)):
    steps[0].append(diff_delay[0].times[i])
for i in range(len(diff_delay[1].times)):
    steps[1].append(diff_delay[1].times[i])
for i in range(len(diff_inflow[0].times)):
    if diff_inflow[0].times[i] not in steps[0]:
        steps[0].append(diff_inflow[0].times[i])
for i in range(len(diff_inflow[1].times)):
    if diff_inflow[1].times[i] not in steps[1]:
        steps[1].append(diff_inflow[1].times[i])

steps[0].sort()
steps[1].sort()
sorted_steps = [steps[0], steps[1]]

integral_1 = 0
integral_2 = 0
for i in range(len(sorted_steps[0]) - 1):
    start = sorted_steps[0][i]
    end = sorted_steps[0][i+1] - 2*eps
    value = diff_inflow[0].multiply(diff_delay[0], start, end).integrate(start, end, True)
    integral_1 = integral_1 + value
for i in range(len(sorted_steps[1]) - 1):
    start = sorted_steps[1][i]
    end = sorted_steps[1][i+1] - 2*eps
    value = diff_inflow[1].multiply(diff_delay[1], start, end).integrate(start, end, True)
    integral_2 = integral_2 + value

scalar_product = integral_1 + integral_2


print("Hello World")






