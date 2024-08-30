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
u = Node(2,graph)
w = Node(3,graph)
#x = Node(4,graph)
#y = Node(5,graph)
#z = Node(6,graph)
#r = Node(7,graph)
t = Node(4,graph)
e1 = Edge(s,v,0,graph)
e2 = Edge(s,w,1,graph)
e3 = Edge(v,u,2,graph)
e4 = Edge(u,s,3,graph)
e5 = Edge(w,u,4,graph)
e6 = Edge(u,t,5,graph)
#e7 = Edge(w,y,6,graph)
#e8 = Edge(w,z,7,graph)
#e9 = Edge(y,x,8,graph)
#e10 = Edge(y,r,9,graph)
#e11 = Edge(y,z,10,graph)
#e12 = Edge(x,t,11,graph)
#e13 = Edge(r,t,12,graph)
#e14 = Edge(z,t,13,graph)

graph.nodes = {0:s,1:v,2:u,3:w,4:t}
graph.edges = [e1,e2,e3,e4,e5,e6]

capacities = [1,3,1,2,3,2]
travel_times = [1,1,1,1,1,1]
net_inflow = RightConstant([0,1,1.75,3],[2.5,1,3,0],(0,3))

p1 = [e1,e3,e6]
p2 = [e2,e5,e6]
#p3 = [e1,e5,e10,e13]
#p4 = [e1,e5,e11,e14]
#p5 = [e2,e6,e9,e12]
#p6 = [e2,e6,e10,e13]
#p7 = [e2,e6,e11,e14]
#p8 = [e3,e7,e9,e12]
#p9 = [e3,e7,e10,e13]
#p10 = [e3,e7,e11,e14]
#p11 = [e3,e8,e14]


i1 = RightConstant([0,0.5,1,3],[1.5,0.5,0,0],(0,3))
i2 = RightConstant([0,0.5,1,1.75,3],[1,2,1,3,0],(0,3))
#i3 = RightConstant([0,4],[0,0],(0,4))
#i4 = RightConstant([0,4],[0,0],(0,4))
#i5 = RightConstant([0,4],[0,0],(0,4))
#i6 = RightConstant([0,0.5,1,4],[1.5,0.5,0,0],(0,4))
#i7 = RightConstant([0,4],[0,0],(0,4))
#i8 = RightConstant([0,4],[0,0],(0,4))
#i9 = RightConstant([0,0.5,1,1.75,4],[1,2,1,3,0],(0,4))
#i10 = RightConstant([0,4],[0,0],(0,4))
#i11 = RightConstant([0,4],[0,0],(0,4))

i12 = RightConstant([0,1,3],[1,0,0],(0,3))
i13 = RightConstant([0,1,1.75,3],[1.5,1,3,0],(0,3))
#i14 = RightConstant([0,4],[0,0],(0,4))
#i15 = RightConstant([0,4],[0,0],(0,4))
#i16 = RightConstant([0,4],[0,0],(0,4))
#i17 = RightConstant([0,1,4],[1,0,0],(0,4))
#i18 = RightConstant([0,4],[0,0],(0,4))
#i19 = RightConstant([0,4],[0,0],(0,4))
#i20 = RightConstant([0,1,1.75,4],[1.5,1,3,0],(0,4))
#i21 = RightConstant([0,4],[0,0],(0,4))
#i22 = RightConstant([0,4],[0,0],(0,4))

in_f = [i1,i2]
in_g = [i12,i13]

monotonicity_check(graph, capacities, travel_times, net_inflow,3,[p1,p2], in_f, in_g)