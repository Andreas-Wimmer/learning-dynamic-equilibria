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

def monotonicity_check(graph: DirectedGraph, capacities: List[float], travel_times: List[float],
                       network_inflow: RightConstant, T: float, paths = List[Path], 
                       inflows_1 = List[RightConstant], inflows_2 = List[RightConstant]):
    net_inflow = Commodity({0: network_inflow}, len(graph.edges) - 1, 1)
    network = Network()
    network.graph = graph
    network.capacity = capacities
    network.travel_time = travel_times
    network.commodities = [net_inflow]
    network.paths = paths 

    inflow_dict_1 = []
    for i in range(len(paths)):
        inflow_dict_1.append((paths[i], inflows_1[i]))

    inflow_dict_2 = []
    for i in range(len(paths)):
        inflow_dict_2.append((paths[i], inflows_2[i]))

    loader_1 = NetworkLoader(network, inflow_dict_1)
    loader_2 = NetworkLoader(network, inflow_dict_2)

    result_1 = loader_1.build_flow()
    result_2 = loader_2.build_flow()

    flow_1 = next(result_1)
    flow_2 = next(result_2)

    delays_1 = loader_1.path_delay()
    delays_2 = loader_2.path_delay()

    diff_delays = []
    for i in range(len(paths)):
        diff_delays.append(delays_1[i] - delays_2[i])

    diff_inflows = []
    for i in range(len(paths)):
        diff_inflows.append(inflows_1[i] - inflows_2[i])
    
    steps = []
    for i in range(len(paths)):
        steps.append([])
        for j in range(len(diff_delays[i].times)):
            steps[i].append(diff_delays[i].times[j])
        for k in range(len(diff_inflows[i].times)):
            if diff_inflows[i].times[j] not in steps[i]:
                steps[i].append(diff_inflows[i].times[k])
    
    for i in range(len(paths)):
        steps[i].sort()

    integrals = []
    for i in range(len(paths)):
        integrals.append(0)

    for i in range(len(paths)):
        for j in range(len(steps[j]) - 1):
            start = steps[i][j]
            end = steps[i][j+1] - 2*eps
            value = diff_inflows[i].multiply(diff_delays[i], start, end).integrate(start, end, True)
            integrals[i] = integrals[i] + value

    scalar_product = sum(integrals)
    return scalar_product
    