#Here we want to implement the learning dynamics on general networks
from __future__ import annotations

import scipy.optimize

from graph import DirectedGraph, Node, Edge
from network import Network, Commodity
from dynamic_flow import DynamicFlow
from piecewise_linear import PiecewiseLinear 
from right_constant import RightConstant 
import numpy 
import scipy
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
    accuracy_reached = False
    equilibrium_reached = False
    while counter_steps < numSteps and not accuracy_reached and not equilibrium_reached:
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
                    value_4 = epsilon*(h[var_index] - inflow_avg[i].eval(breaks[i][j]))
                    sum = sum + value_3 + value_4
        
        A = []
        for j in range(len(steps) - 1):
            A.append([])
            for k in range(len(network.paths)):
                for g in range(len(steps)):
                    if j == g:
                        A[j].append(1)
                    else:
                        A[j].append(0)

        B = []
        for j in range(len(network.paths)):
            B.append([])
            for k in range(len(network.paths)):
                for g in range(len(steps)):
                    if g == len(steps) - 1 and k == j:
                        B[j].append(1)
                    else:
                        B[j].append(0)

        net_bound = []
        for i in range(len(steps)):
            net_bound.append(net_inflow.eval(steps[i]))
        constraint_1 = scipy.optimize.LinearConstraint(A, net_bound, net_bound)
        constraint_2 = scipy.optimize.LinearConstraint(B, 0, 0)
        bounds = []
        for j in range(len(network.paths)):
            for k in range(len(steps)):
                bounds.append((0, None))

        start = []
        for i in range(len(network.paths)):
            for j in range(len(steps)):
                start.append(inflows[i].eval(steps[j]))
        
        sol = scipy.optimize.minimize(obj, start, bounds = bounds, constraints = [constraint_1, constraint_2])
    
        #3. Update the path inflows and run the edge-loading procedure
        old_inflows = inflows.copy()
        old_inflows_dict = inflow_dict.copy()

        inflows = []
        values = []
        for i in range(len(network.paths)):
            values.append([])
            for j in range(len(steps)):
                values[i].append(sol.x[len(steps)*i + j])

        for i in range(len(network.paths)):
            inflows.append(RightConstant(steps, values[i], (0, horizon)))

        inflow_dict = []
        for i in range(len(network.paths)):
            inflow_dict.append((network.paths[i], inflows[i]))

        loader_new = NetworkLoader(network, inflow_dict)
        result_new = loader_new.build_flow()
        flow_new = next(result_new)
        delays_new = loader_new.path_delay()
        counter_steps = counter_steps + 1
        #4. Update the population average and run the edge-loading procedure again
        old_avg = inflow_avg.copy()
        old_avg_delays = delays_avg.copy()

        values_avg = []
        for i in range(len(network.paths)):
            values_avg.append([])
            for j in range(len(steps)):
                new_value_1 = (1/counter_steps)*inflows[i].values[j]
                new_value_2 = ((counter_steps - 1)/counter_steps)*old_avg[i].values[j]
                values_avg[i].append(new_value_1 + new_value_2)

        inflow_avg = []
        for i in range(len(network.paths)):
            inflow_avg.append(RightConstant(steps, values_avg[i], (0, horizon)))
    
        inflow_dict_avg = []
        for i in range(len(network.paths)):
            inflow_dict_avg.append((network.paths[i], inflow_avg[i]))

        loader_avg = NetworkLoader(network, inflow_dict_avg)
        result_avg = loader_avg.build_flow()
        flow_avg = next(result_avg)
    
        #5. Calculate the difference of the new flow and the last flow for checking convergence
        diff_inflows = []
        for i in range(len(network.paths)):
            curr_diff = inflows[i] - old_inflows[i]
            diff_inflows.append(curr_diff.abs())
    
        diff = 0
        for i in range(len(network.paths)):
            diff = diff + diff_inflows[i].integral()
    
        if round(diff, 12) <= lamb:
            accuracy_reached = True

        print("Norm difference to the last flow" + str(diff))
        #6. Calculate the (regularized) gap for checking, if we get close to a dynamic equilibrium
        gap_steps = []
        for i in range(len(network.paths)):
            gap_steps.append([])
            for j in range(len(steps)):
                gap_steps[i].append(steps[j])
            for k in range(len(delays_new[i].times)):
                if delays_new[i].times[k] not in gap_steps[i]:
                    gap_steps.append(delays_new[i].times[k])
            gap_steps[i].sort()

        def g(h):
            sum_1 = 0
            for i in range(len(network.paths)):
                for j in range(len(gap_steps[i]) - 1):
                    if gap_steps[i][j] in steps:
                        varindex = steps.index(gap_steps[i][j])
                    else:
                        varindex = elem_lrank(steps, gap_steps[i][j])
                    val_1 = delays_new[i].eval(gap_steps[i][j+1])
                    val_2 = delays_new[i].eval(gap_steps[i][j])
                    val_3 = 2*epsilon*(-h[len(steps)*i + varindex] + inflows[i].eval(gap_steps[i][j]))
                    val_4 = h[len(steps)*i + varindex] - inflows[i].eval(gap_steps[i][j])
                    sum_1 = sum_1 + (((val_1 + val_2)/2) - val_3)*val_4
            return sum_1

        A = []
        for j in range(len(steps) - 1):
            A.append([])
            for k in range(len(network.paths)):
                for g in range(len(steps)):
                    if j == g:
                        A[j].append(1)
                    else:
                        A[j].append(0)

        B = []
        for j in range(len(network.paths)):
            B.append([])
            for k in range(len(network.paths)):
                for g in range(len(steps)):
                    if g == len(steps) - 1 and k == j:
                        B[j].append(1)
                    else:
                        B[j].append(0)

        net_bound = []
        for i in range(len(steps)):
            net_bound.append(net_inflow.eval(steps[i]))
        constraint_1 = scipy.optimize.LinearConstraint(A, net_bound, net_bound)
        constraint_2 = scipy.optimize.LinearConstraint(B, 0, 0)
        bounds = []
        for j in range(len(network.paths)):
            for k in range(len(steps)):
                bounds.append((0, None))

        start = []
        for i in range(len(network.paths)):
            for j in range(len(steps)):
                start.append(inflows[i].eval(steps[j]))
        constraints = [constraint_1, constraint_2]
        sol_gap = scipy.optimize.minimize(g, start, bounds=bounds, constraints=constraints)
        
        if round(sol.gap, 5) == 0:
            equilibrium_reached = True
            print("Regularized equilibrium reached")

        print("Gap at " + str(sol.gap))

    if accuracy_reached == False or equilibrium_reached == False:
        print("The learning dynamics did not converge to a equilibrium with the given accuracy")
    else:
        print("The learning dynamics did converge to a regularized equilibrium with the given accuracy")