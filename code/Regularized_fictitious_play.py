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
    t = graph.nodes[len(graph.nodes) - 1]
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
        min_cap[i] = min(caps)

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
    delays_beg = loader_beg.path_delay(horizon)
    delays_avg = delays_beg
    inflow_avg = inflows.copy()
    
    counter_steps = 1
    accuracy_reached = False
    equilibrium_reached = False
    while counter_steps < numSteps and not accuracy_reached and not equilibrium_reached:
        step = 0
        steps = []
        while step*delta <= horizon - delta:
            steps.append(step*delta)
            step = step + 1
        steps.append(horizon)

        breaks = steps.copy()
        for k in range(len(breaks_net_inflow)):
            if breaks_net_inflow[k] not in breaks:
                breaks.append(breaks_net_inflow)

        for i in range(len(network.paths)):
            for j in range(len(delays_avg[i].times)):
                if delays_avg[i].times[j] not in breaks and delays_avg[i].times[j] <= horizon:
                    breaks.append(delays_avg[i].times[j])
        
        breaks.sort()

        time_steps = []
        time_steps.append(0)
        for i in range(len(breaks)):
            if breaks[i] - time_steps[-1] >= delta/5:
                time_steps.append(breaks[i])
        
        time_steps.sort()
        #2. Best-response problem: 
        def obj(h):
            sums = 0
            for i in range(len(network.paths)):
                for j in range(len(breaks) - 1):
                    value_1 = delays_avg[i].eval(breaks[j])
                    value_2 = delays_avg[i].eval(breaks[j+1])
                    if breaks[j] in time_steps:
                        varindex = time_steps.index(breaks[j])
                    else:
                        varindex = elem_lrank(time_steps, breaks[j])
                    value_3 = ((value_1 + value_2)/2)*h[len(time_steps)*i + varindex]
                    value_4 = epsilon*pow((h[len(time_steps)*i + varindex] - inflow_avg[i].eval(breaks[j])),2)
                    sums = sums + value_3 + value_4

            return sums 
        
        A = []
        for j in range(len(time_steps)):
            A.append([])
            for k in range(len(network.paths)):
                for g in range(len(time_steps)):
                    if j == g:
                        A[j].append(1)
                    else:
                        A[j].append(0)

        net_bound = []
        for i in range(len(time_steps)):
            net_bound.append(net_inflow.eval(time_steps[i]))
        constraint_1 = scipy.optimize.LinearConstraint(A, net_bound, net_bound)
        bounds = []
        for j in range(len(network.paths)):
            for k in range(len(time_steps)):
                bounds.append((0, float("inf")))

        start = []
        for i in range(len(network.paths)):
            for j in range(len(time_steps)):
                start.append(inflows[i].eval(time_steps[j]))
        
        sol = scipy.optimize.minimize(obj, start, bounds = bounds, constraints = constraint_1)
    
        #3. Update the path inflows and run the edge-loading procedure
        old_inflows = inflows.copy()

        inflows = []
        values = []
        for i in range(len(network.paths)):
            values.append([])
            for j in range(len(time_steps)):
                values[i].append(sol.x[len(time_steps)*i + j])

        for i in range(len(network.paths)):
            inflows.append(RightConstant(time_steps, values[i], (0, horizon)))

        inflow_dict = []
        for i in range(len(network.paths)):
            inflow_dict.append((network.paths[i], inflows[i]))

        loader_new = NetworkLoader(network, inflow_dict)
        result_new = loader_new.build_flow()
        flow_new = next(result_new)
        delays_new = loader_new.path_delay(horizon)
        counter_steps = counter_steps + 1
        #4. Update the population average and run the edge-loading procedure again
        old_avg = inflow_avg.copy()

        values_avg = []
        for i in range(len(network.paths)):
            values_avg.append([])
            for j in range(len(breaks)):
                new_value_1 = (1/counter_steps)*inflows[i].eval(breaks[j])
                new_value_2 = ((counter_steps - 1)/counter_steps)*old_avg[i].eval(breaks[j])
                values_avg[i].append(new_value_1 + new_value_2)

        inflow_avg = []
        for i in range(len(network.paths)):
            inflow_avg.append(RightConstant(breaks, values_avg[i], (0, horizon)))
    
        inflow_dict_avg = []
        for i in range(len(network.paths)):
            inflow_dict_avg.append((network.paths[i], inflow_avg[i]))

        loader_avg = NetworkLoader(network, inflow_dict_avg)
        result_avg = loader_avg.build_flow()
        flow_avg = next(result_avg)
        delays_avg = loader_avg.path_delay(horizon)
        #5. Calculate the difference of the new flow and the last flow for checking convergence
        diff_inflows = []
        for i in range(len(network.paths)):
            curr_diff = inflow_avg[i] - old_avg[i]
            diff_inflows.append(curr_diff.abs_val())
    
        diff = 0
        for i in range(len(network.paths)):
            diff = diff + diff_inflows[i].integral().eval(horizon)
    
        if round(diff, 12) <= lamb:
            accuracy_reached = True

        print("Norm difference to the last flow: " + str(diff))
        #6. Calculate the (regularized) gap for checking, if we get close to a dynamic equilibrium
        gap_breaks = []
        gap_breaks = steps.copy()
        for k in range(len(breaks_net_inflow)):
            if breaks_net_inflow[k] not in gap_breaks:
                gap_breaks.append(breaks_net_inflow[k])
        
        for i in range(len(network.paths)):
            for j in range(len(delays_avg[i].times)):
                if delays_avg[i].times[j] not in gap_breaks and delays_avg[i].times[j] <= horizon:
                    gap_breaks.append(delays_avg[i].times[j])
        
        gap_breaks.sort()

        gap_steps = []
        gap_steps.append(0)
        for i in range(len(gap_breaks)):
            if gap_breaks[i] - gap_steps[-1] >= delta/5:
                gap_steps.append(gap_breaks[i])

        gap_steps.sort() 

        def obj_gap(h):
            sum_1 = 0
            for i in range(len(network.paths)):
                for j in range(len(gap_breaks) - 1):
                    val_1 = delays_avg[i].eval(gap_breaks[j+1])
                    val_2 = delays_avg[i].eval(gap_breaks[j])
                    if gap_breaks[j] in gap_steps:
                        varindex = gap_steps.index(gap_breaks[j])
                    else:
                        varindex = elem_lrank(gap_steps, gap_breaks[j])
                    #val_3 = 2*epsilon*(-h[len(gap_steps)*i + varindex] + inflow_avg[i].eval(gap_breaks[j]))
                    val_4 = h[len(gap_steps)*i + varindex] - inflow_avg[i].eval(gap_breaks[j])
                    sum_1 = sum_1 + ((val_1 + val_2)/2)*val_4
            return sum_1

        A = []
        for j in range(len(gap_steps)):
            A.append([])
            for k in range(len(network.paths)):
                for g in range(len(gap_steps)):
                    if j == g:
                        A[j].append(1)
                    else:
                        A[j].append(0)

        net_bound = []
        for i in range(len(gap_steps)):
            net_bound.append(net_inflow.eval(gap_steps[i]))
        constraint_1 = scipy.optimize.LinearConstraint(A, net_bound, net_bound)
        bounds = []
        for j in range(len(network.paths)):
            for k in range(len(gap_steps)):
                bounds.append((0, float("inf")))

        start = []
        for i in range(len(network.paths)):
            for j in range(len(gap_steps)):
                start.append(inflow_avg[i].eval(gap_steps[j]))
        #sol_gap = scipy.optimize.minimize(obj_gap, start, bounds=bounds, constraints=constraint_1)
        
        #if round(sol_gap.fun, 5) == 0:
        #    equilibrium_reached = True
        #    print("Regularized equilibrium reached")

        #print("Gap at " + str(sol_gap.fun))

    if equilibrium_reached:
        print("The learning dynamics reached a regularized equilibrium")
    elif accuracy_reached:
        print("The learning dynamics converged")
    else:
        print("The learning dynamics neither reached a regularized equilibrium nor converged wtih the given accuracy in the given number of steps")

graph = DirectedGraph
s = Node(0,graph)
v = Node(1,graph)
w = Node(2,graph)
x = Node(3,graph)
y = Node(4,graph)
t = Node(5,graph)
e_1 = Edge(s,v,0,graph)
e_2 = Edge(s,w,1,graph)
e_3 = Edge(v,x,2,graph)
e_4 = Edge(v,y,3,graph)
e_5 = Edge(w,x,4,graph)
e_6 = Edge(w,y,5,graph)
e_7 = Edge(x,t,6,graph)
e_8 = Edge(y,t,7,graph)

graph.nodes = {0:s,1:v,2:w,3:x,4:y,5:t}
graph.edges = [e_1,e_2,e_3,e_4,e_5,e_6,e_7,e_8]
graph.reversed = False

capacities = [2,1,2,2,1,1,3,3]
travel_times = [1,1,1,2,2,1,1,1]

p_1 = [e_1,e_3,e_7]
p_2 = [e_1,e_4,e_8]
p_3 = [e_2,e_5,e_7]
p_4 = [e_2,e_6,e_8]

paths = [p_1, p_2,p_3,p_4]
net_inflow = RightConstant([0,4],[8,0], (0, 4))
horizon = 4
delta = 0.5
epsilon = 0.5
numSteps = 200
lamb = 0.00001


reg_fictitious_play(graph, capacities, travel_times,
                    net_inflow, paths, horizon, delta, epsilon, numSteps, lamb)