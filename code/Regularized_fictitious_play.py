#Here we want to implement the learning dynamics on general networks
from __future__ import annotations

from graph import DirectedGraph, Node, Edge
from network import Network, Commodity
from dynamic_flow import DynamicFlow
from piecewise_linear import PiecewiseLinear 
from right_constant import RightConstant 
import numpy 
from network_loader import NetworkLoader, Path
from machine_precision import eps
from typing import List
import math
from arrays import *


def reg_fictitious_play(graph: DirectedGraph, cap: List[float], travel: List[float],
                         net_inflow: RightConstant, paths: List[Path], horizon: float,
                         delta: float, epsilon: float, numSteps: int, lamb: float):
    #Steps that need to be taken:
    #1. Initialization:
    network = Network()
    network.graph = graph
    network.capacity = cap
    network.travel_time = travel
    s = graph.nodes[0]
    t = graph.nodes[-1]
    network_inflow = Commodity({s, net_inflow}, t, 1)
    network.commodities = [network_inflow]
    network.paths = paths

    if numSteps == None:
        numSteps = math.inf

    if lamb == None:
        lamb = 0

    min_cap = []
    for i in range(len(network.paths)):
        min_cap.append(0)
        caps = []
        for j in range(len(network.paths[i])):
            index = network.paths[i][j].id
            caps.append(network.capacity[index])
        min_cap.append(min(caps))

    sum_min_caps = sum(min_cap)
    portions = []
    for i in range(len(network.paths)):
        portions.append(min_cap[i]/sum_min_caps)

    breaks_net_inflow = net_inflow.times
    values = []
    for i in range(len(network.paths)):
        values.append([])
        for j in range(len(breaks_net_inflow)):
            values[i].append(portions[i]*net_inflow.eval(breaks_net_inflow[j]))

    inflows = []
    inflow_dict = []
    for i in range(len(network.paths)):
        inflows.append(RightConstant(breaks_net_inflow, values[i], (0, horizon)))
        inflow_dict.append((network.paths[i], inflows[i]))

    loader_beg = NetworkLoader(network, inflow_dict)
    result_beg = loader_beg.build_flow()
    flow_beg = next(result_beg)
    delays_beg = loader_beg.path_delay()
    flow_avg = flow_beg
    delays_avg = delays_beg
    inflow_avg = inflows.copy()
    
    counter_steps = 1
    while counter_steps < numSteps:
        stop = 0
        steps = []
        while stop*delta <= horizon - delta:
            for i in range(len(network.paths)):
                steps.append(stop*delta)
            stop = stop + 1

        breaks = []
        for i in range(len(network.paths)):
            breaks.append(steps.copy())
            for j in range(len(delays_avg[i].times[j])):
                if delays_avg[i].times[j] not in breaks[i]:
                    breaks[i].append(delays_avg[i].times[j])
            breaks[i].sort()
    #2. Best-response problem: 
    def obj(h):
        sum = 0
        for i in range(len(network.paths)):
            for j in range(len(breaks[i]) - 1):
                value_1 = delays_avg[i].eval(breaks[i][j])
                value_2 = delays_avg[i].eval(breaks[i][j+1])
                if breaks[i][j] in steps:
                    index = steps.index(breaks[i][j])
                    var_index = (len(steps) - 1)*i + index
                else: 
                    index = elem_rank(steps, breaks[i][j])
                    var_index = (len(steps) - 1)*i + index
                value_3 = ((value_1 + value_2)/2)*h[var_index]
                value_4 = epsilon*(h[var_index] - inflows[i].eval(breaks[i][j]))
                sum = sum + value_3 + value_4
    #3. Update the flow and run the edge-loading procedure
    #4. Update the population average and run the edge-loading procedure again
    #5. Calculate the difference of the new flow and the last flow for checking convergence
    #6. Calculate the (regularized) gap for checking, if we get close to a dynamic equilibrium
