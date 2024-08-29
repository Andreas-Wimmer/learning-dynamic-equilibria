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
#x = Node(2,graph)
#y = Node(3,graph)
t = Node(2,graph)
e1 = Edge(s,v,0,graph)
e2 = Edge(s,v,1,graph)
e3 = Edge(v,t,2,graph)
#e4 = Edge(x,y,3,graph)
#e5 = Edge(x,y,4,graph)
#e6 = Edge(y,t,5,graph)
#e7 = Edge(x,t,6,graph)

graph.nodes = {0:s,1:v,2:t}
graph.edges = [e1,e2,e3]

capacities = [1,2,2]
travel_times = [1,1,1]
net_inflow = RightConstant([0,1,1.75,2],[2,1,2,0],(0,2))

p1 = [e1,e3]
p2 = [e2,e3]
#p3 = [e1,e3,e4,e6]
#p4 = [e1,e3,e5,e6]



i1 = RightConstant([0,0.5,1,2],[1.5,0.5,0,0],(0,2))
i2 = RightConstant([0,0.5,1,1.75,2],[0.5,1.5,1,2,0],(0,2))
#i3 = RightConstant([0,0.5,1,1.5,3],[0.75,1.25,0.5,1.5,0],(0,3))
#i4 = RightConstant([0,0.5,1,1.5,3],[0.75,1.25,0.5,1.5,0],(0,3))
#i5 = RightConstant([0,1,1.5,2],[0,1,0,0],(0,2))
#i6 = RightConstant([0,1,1.5,2],[0,1,0,0],(0,2))

i7 = RightConstant([0,1,2],[1,0,0],(0,2))
i8 = RightConstant([0,1,1.75,2],[1,1,2,0],(0,2))
#i9 = RightConstant([0,1,1.5,3],[1,0.5,1.5,0],(0,3))
#i10 = RightConstant([0,1,1.5,3],[1,0.5,1.5,0],(0,3))
#i11 = RightConstant([0,1,1.5,2],[0,1,0,0],(0,2))
#i12 = RightConstant([0,1,1.5,2],[0,1,0,0],(0,2))

in_f = [i1,i2]
in_g = [i7,i8]

monotonicity_check(graph, capacities, travel_times, net_inflow,2, [p1,p2], in_f, in_g)