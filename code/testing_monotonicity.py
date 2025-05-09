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
import nguyen_network
import sioux_falls_network

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

graph = sioux_falls_network.sioux_graph
capacities = sioux_falls_network.capacities
travel_times = sioux_falls_network.travel_times
net_inflow = RightConstant([0,200],[8,0],(0,200))
paths = sioux_falls_network.paths
edges = sioux_falls_network.edges
horizon = 200

in_f = []
in_g = []
for path in paths:
    if path.edges == [edges[0],edges[3],edges[14],edges[12],edges[23],edges[19],edges[17],edges[55]]:
        in_f.append(RightConstant([0,50,100,200],[8,1,4,0],(0,200)))
        in_g.append(RightConstant([0,50,100,200],[7,4,4,0],(0,200)))
    elif path.edges == [edges[1],edges[5],edges[8],edges[12],edges[23],edges[19],edges[17],edges[55]]:
        in_f.append(RightConstant([0,50,100,200],[0,4,4,0],(0,200)))
        in_g.append(RightConstant([0,50,100,200],[1,2,4,0],(0,200)))
    elif path.edges == [edges[1],edges[6],edges[36],edges[38],edges[74],edges[63]]:
        in_f.append(RightConstant([0,50,100,200],[0,3,0,0],(0,200)))
        in_g.append(RightConstant([0,50,100,200],[0,2,0,0],(0,200)))
    else:
        in_f.append(RightConstant([0,200],[0,0],(0,200)))
        in_g.append(RightConstant([0,200],[0,0],(0,200)))

#i1 = RightConstant([0,300],[0,0],(0,300))
#i2 = RightConstant([0,300],[0,0],(0,300))
#i3 = RightConstant([0,300],[0,0],(0,300))
#i4 = RightConstant([0,100,200,300],[3,1,2,0],(0,300))
#i5 = RightConstant([0,300],[0,0],(0,300))
#i6 = RightConstant([0,300],[0,0],(0,300))
#i7 = RightConstant([0,100,200,300],[1,3,2,0],(0,300))
#i8 = RightConstant([0,300],[0,0],(0,300))
#i9 = RightConstant([0,1,2,4,5],[0.5,1,0.5,1.5,0],(0,5))
#i10 = RightConstant([0,4],[0,0],(0,4))
#i11 = RightConstant([0,4],[0,0],(0,4))
#i12 = RightConstant([0,4],[0,0],(0,4))
#i13 = RightConstant([0,4],[0,0],(0,4))
#i14 = RightConstant([0,4],[0,0],(0,4))
#i15 = RightConstant([0,4],[0,0],(0,4))
#i16 = RightConstant([0,4],[0,0],(0,4))
#i17 = RightConstant([0,1,1.75,4],[1,2,1,0],(0,4))

#i18 = RightConstant([0,300],[0,0],(0,300))
#i19 = RightConstant([0,300],[0,0],(0,300))
#i20 = RightConstant([0,300],[0,0],(0,300))
#i21 = RightConstant([0,100,200,300],[1,3,2,0],(0,300))
#i22 = RightConstant([0,300],[0,0],(0,300))
#i23 = RightConstant([0,300],[0,0],(0,300))
#i24 = RightConstant([0,100,200,300],[3,1,2,0],(0,300))
#i25 = RightConstant([0,300],[0,0],(0,300))
#i26 = RightConstant([0,300],[0,0],(0,300))
#i27 = RightConstant([0,300],[0,0],(0,300))
#i28 = RightConstant([0,300],[0,0],(0,300))
#i29 = RightConstant([0,4],[0,0],(0,4))
#i30 = RightConstant([0,4],[0,0],(0,4))
#i31 = RightConstant([0,4],[0,0],(0,4))
#i32 = RightConstant([0,4],[0,0],(0,4))
#i33 = RightConstant([0,4],[0,0],(0,4))
#i34 = RightConstant([0,1,1.75,4],[1,2,1,0],(0,4))

#in_f = [i1,i2,i3,i4,i5,i6,i7,i8]
#in_g = [i18,i19,i20,i21,i22,i23,i24,i25]

monotonicity_check(graph, capacities, travel_times, net_inflow, horizon, paths, in_f, in_g)