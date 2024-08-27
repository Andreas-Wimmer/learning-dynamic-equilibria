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

graph = DirectedGraph()
s = Node(0,graph)
v = Node(1,graph)
w = Node(2,graph)
x = Node(3,graph)
y = Node(4,graph)
z = Node(5,graph)
u = Node(6,graph)
t = Node(7,graph)
e_1 = Edge(s,v,0,graph)
e_2 = Edge(s,w,1,graph)
e_3 = Edge(s,x,2,graph)
e_4 = Edge(s,y,3,graph)
e_5 = Edge(v,z,4,graph)
e_6 = Edge(w,z,5,graph)
e_7 = Edge(x,u,6,graph)
e_8 = Edge(y,u,7,graph)
e_9 = Edge(z,t,8,graph)
e_10 = Edge(u,t,9,graph)


graph.nodes = {0:s,1:v,2:w,3:x,4:y,5:z,6:u,7:t}
graph.edges = [e_1,e_2,e_3,e_4,e_5,e_6,e_7,e_8,e_9,e_10]

capacities = [1,3,1,3,1,3,1,3,2,2]
travel_times = [1,1,1,1,1,1,1,1,1,1]
net_inflow = RightConstant([0,1,1.75,2],[5,2,6,0],(0,2))

path_1 = [e_1,e_5,e_9]
path_2 = [e_2,e_6,e_9]
path_3 = [e_3,e_7,e_10]
path_4 = [e_4,e_8,e_10]


inflow_1 = RightConstant([0,0.5,1,2],[1.5,0.5,0,0],(0,2))
inflow_2 = RightConstant([0,0.5,1,1.75,2],[1,2,1,3,0],(0,2))
inflow_3 = RightConstant([0,0.5,1,2],[1.5,0.5,0,0],(0,2))
inflow_4 = RightConstant([0,0.5,1,1.75,2],[1,2,1,3,0],(0,2))
inflow_5 = RightConstant([0,1,2],[1,0,0],(0,2))
inflow_6 = RightConstant([0,1,1.75,2],[1.5,1,3,0],(0,2))
inflow_7 = RightConstant([0,1,2],[1,0,0],(0,2))
inflow_8 = RightConstant([0,1,1.75,2],[1.5,1,3,0],(0,2))

inflow_f = [inflow_1,inflow_2,inflow_3,inflow_4]
inflow_g = [inflow_5,inflow_6,inflow_7,inflow_8]

monotonicity_check(graph, capacities, travel_times, net_inflow, 2, [path_1,path_2,path_3,path_4], inflow_f, inflow_g)