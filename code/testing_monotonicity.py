#Here we implement a routine that checks for a given network, a given network inflow rate and 
#given path inflows, if we have monotonicity for the given parameters
from __future__ import annotations

from graph import *
from network import Network,Commodity
from dynamic_flow import DynamicFlow 
from network_loader import NetworkLoader, Path
from piecewise_linear import PiecewiseLinear, identity
from right_constant import RightConstant
from machine_precision import eps
import gap_function

def monotonicity_check(graph: DirectedGraph, capacities: List[float], travel_times: List[float],
                       network_inflow: RightConstant, T: float, paths = List[Path], 
                        inflows_1 = List[RightConstant], inflows_2 = List[RightConstant]) -> float:
    network = Network()
    network.graph = graph
    network.capacity = capacities
    network.travel_time = travel_times
    s = graph.nodes[0]
    t = graph.nodes[len(graph.nodes) - 1]
    net_inflow = Commodity({s: network_inflow},t,1)
    network.commodities = [net_inflow]
    network.paths = paths 

    inflow_dict_1 = []
    inflow_dict_2 = []
    for i in range(len(paths)):
        inflow_dict_1.append((paths[i], inflows_1[i]))
        inflow_dict_2.append((paths[i], inflows_2[i]))

    loader_1 = NetworkLoader(network, inflow_dict_1)
    loader_2 = NetworkLoader(network, inflow_dict_2)

    result_1 = loader_1.build_flow()
    result_2 = loader_2.build_flow()

    flow_1 = next(result_1)
    flow_2 = next(result_2)

    delays_1 = loader_1.path_delay(T)
    delays_2 = loader_2.path_delay(T)

    diff_delays = []
    diff_inflows = []
    for i in range(len(paths)):
        diff_delays.append(delays_1[i] - delays_2[i])
        diff_inflows.append(inflows_1[i] - inflows_2[i])
    
    steps = []
    for i in range(len(paths)):
        steps.append([])
        for j in range(len(diff_delays[i].times)):
            if diff_delays[i].times[j] <= T:
                steps[i].append(diff_delays[i].times[j])
        for k in range(len(diff_inflows[i].times)):
            if diff_inflows[i].times[k] not in steps[i] and diff_inflows[i].times[k] <= diff_delays[i].times[-1] :
                steps[i].append(diff_inflows[i].times[k])
    
    for i in range(len(paths)):
        steps[i].sort()

    integrals = []
    for i in range(len(paths)):
        integrals.append(0)

    for i in range(len(paths)):
        for j in range(len(steps[i]) - 1):
            if steps[i][j+1] - steps[i][j] >= 10*eps:
                start = steps[i][j]
                end = steps[i][j+1] - 2*eps
                value = diff_inflows[i].multiply(diff_delays[i], start, end).integrate(start, end, True)
                integrals[i] = integrals[i] + value

    scalar_product = sum(integrals)
    #term = gap_function.gap_function(inflows_1,delays_1,inflows_2,T)
    if abs(scalar_product) < eps:
        scalar_product = 0
    print(str(scalar_product))
    return scalar_product

test_graph = DirectedGraph()
s = Node(0,test_graph)
v = Node(1,test_graph)
w = Node(2,test_graph)
x = Node(3,test_graph)
t = Node(4,test_graph)
edge_1 = Edge(s,v,0,test_graph)
edge_2 = Edge(s,x,1,test_graph)
edge_3 = Edge(v,w,2,test_graph)
edge_4 = Edge(x,v,3,test_graph)
edge_5 = Edge(w,x,4,test_graph)
edge_6 = Edge(w,t,5,test_graph)
edge_7 = Edge(x,t,6,test_graph)


test_graph.nodes = {0:s,1:v,2:w,3:x,4:t}
test_graph.edges = [edge_1,edge_2,edge_3,edge_4,edge_5,edge_6,edge_7]

capacities = [2,3,1,1,1,2,1]
travel_times = [0,1,1,0,1,0,0]
net_inflow = RightConstant([0,1,2,3,4],[6,2,2,4,0],(0,4))

path_1 = [edge_1,edge_3,edge_6]
path_2 = [edge_2,edge_4,edge_3,edge_6]
path_3 = [edge_1,edge_3,edge_5,edge_7]
path_4 = [edge_2,edge_7]


inflow_1 = RightConstant([0,1,2,3,4],[1.5,0.5,1,1,0],(0,4))
inflow_2 = RightConstant([0,1,2,3,4],[1.5,0.5,1,1,0],(0,4))
inflow_3 = RightConstant([0,1,2,3,4],[1.5,0.5,0,1,0],(0,4))
inflow_4 = RightConstant([0,1,2,3,4],[1.5,0.5,0,1,0],(0,4))
inflow_5 = RightConstant([0,1,2,3,4],[2,1,0.5,2,0],(0,4))
inflow_6 = RightConstant([0,1,2,3,4],[1,0,0.5,0,0],(0,4))
inflow_7 = RightConstant([0,1,2,3,4],[2,1,0.5,2,0],(0,4))
inflow_8 = RightConstant([0,1,2,3,4],[1,0,0.5,0,0],(0,4))

inflow_f = [inflow_1,inflow_2,inflow_3,inflow_4]
inflow_g = [inflow_5,inflow_6,inflow_7,inflow_8]

monotonicity_check(test_graph, capacities, travel_times, net_inflow, 4, [path_1,path_2,path_3,path_4], inflow_f, inflow_g)