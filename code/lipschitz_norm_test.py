#A little routine that computes for a given network, a given time horizon, a given network inflow rate and two
#given feasible path flows the path delay functions and the norm difference of the path flows and the path delay functions
#and also their ratio

from __future__ import annotations
import math
from graph import DirectedGraph, Node, Edge
from network import Network, Commodity, Path
from dynamic_flow import DynamicFlow
from piecewise_linear import PiecewiseLinear 
from right_constant import RightConstant 
import scipy
from network_loader import NetworkLoader
from machine_precision import eps
from typing import List
from arrays import *

def lipschitz_norm_test(graph: DirectedGraph, capacities: List[float], travel_times: List[float], network_inflow: RightConstant,
                        horizon: float, paths: List[Path], times: List[float], values_1: List[List[float]], values_2: List[List[float]]):
    network = Network()
    network.graph = graph
    network.capacity = capacities
    network.travel_time = travel_times
    s = graph.nodes[0]
    t = graph.nodes[len(graph.nodes) - 1]
    network_inflow = Commodity({s, network_inflow}, t, 1)
    network.commodities = [network_inflow]
    network.paths = paths
    
    inflows_1 = []
    inflows_2 = []
    for i in range(len(network.paths)):
        inflows_1.append(RightConstant(times, values_1[i], (0,horizon)))
        inflows_2.append(RightConstant(times, values_2[i],(0,horizon)))

    inflow_dict_1 = []
    inflow_dict_2 = []
    for i in range(len(network.paths)):
        inflow_dict_1.append((network.paths[i],inflows_1[i]))
        inflow_dict_2.append((network.pahts[i],inflows_2[i]))


    loader_1 = NetworkLoader(network, inflow_dict_1)
    loader_2 = NetworkLoader(network, inflow_dict_2)
    result_1 = loader_1.build_flow()
    result_2 = loader_2.buiild_flow()
    flow_1 = next(result_1)
    flow_2 = next(result_2)
    delays_1 = loader_1.path_delay(horizon)
    delays_2 = loader_2.path_delay(horizon)

    sum_1 = 0
    for i in range(len(network.paths)):
        for j in range(len(times)):
            sum_1 = sum_1 + math.abs(inflows_1[i].values[j] - inflows_2[i].values[j])**2

    norm_difference_inflows =  math.sqrt(sum_1)

    time_delays = []
    for i in range(len(network.paths)):
        for j in range(len(delays_1[i].times)):
            time_delays.append(delays_1[i].times[j])
    
    for i in range(len(network.paths)):
        for j in range(len(delays_2[i].times)):
            if delays_2[i].times[j] not in time_delays:
                time_delays.append(delays_2[i].times[j])

    time_delays.sort()

    #sum_2 = 0
    #for i in range(len(network.paths)):
    #    for j in range(len(time_delays) - 1):
    #        value = 