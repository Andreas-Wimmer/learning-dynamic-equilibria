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
    gap_values = []
    storage_values = []
    norm_differences = []

    breaks_net_inflow = net_inflow.times
    values = []

    caps = []
    for i in range (len(network.paths)):
        max_cap = 0
        for j in range(len(network.paths[i])):
            if capacities[network.paths[i][j].id] > max_cap:
                max_cap = capacities[network.paths[i][j].id]
        caps.append(max_cap)

    max_path = paths[caps.index(max(caps))]
    for i in range(len(network.paths)):
        values.append([])
        for j in range(len(breaks_net_inflow)):
            if network.paths[i] == max_path:
                values[i].append(net_inflow.values[j])
            else:
                values[i].append(0)
            

    inflows = []
    inflow_dict = []
    for i in range(len(network.paths)):
        inflows.append(RightConstant(net_inflow.times, values[i], (0, horizon)))
        inflow_dict.append((network.paths[i], inflows[i]))

    loader_beg = NetworkLoader(network, inflow_dict)
    result_beg = loader_beg.build_flow()
    flow_beg = next(result_beg)
    delays_beg = loader_beg.path_delay(horizon)
    delays_avg = delays_beg.copy()
    inflow_avg = inflows.copy()
    
    counter_steps = 1
    accuracy_reached = False
    equilibrium_reached = False
    stop = 0
    steps = []
    while stop*delta <= horizon - delta:
        steps.append(stop*delta)
        stop = stop + 1
    steps.append(horizon)

    for i in range(len(net_inflow.times)):
        if net_inflow.times[i] not in steps:
            steps.append(net_inflow.times[i])
    steps.sort()
    while counter_steps < numSteps and not accuracy_reached and not equilibrium_reached:
        #2. Best-response problem: 
        def obj(h):
            sums = 0
            for i in range(len(network.paths)):
                for j in range(len(steps) - 1):
                    sum_delays = 0
                    count_delays = 0
                    steps_in = []
                    for k in range(len(delays_avg[i].times) - 1):
                        if delays_avg[i].times[k] > steps[j] and delays_avg[i].times[k] < steps[j+1]:
                            steps_in.append(delays_avg[i].times[k])
                            count_delays = count_delays + 1
                    if count_delays != 0:
                        weight = ((steps_in[0] - steps[j])/(steps[j + 1] - steps[j]))
                        value = ((delays_avg[i].eval(steps_in[0]) + delays_avg[i].eval(steps[j]))/2)
                        sum_delays = sum_delays + weight*value
                        weight = ((steps[j + 1] - steps_in[-1])/(steps[j + 1] - steps[j]))
                        value = ((delays_avg[i].eval(steps_in[-1]) + delays_avg[i].eval(steps[j + 1]))/2)
                        sum_delays = sum_delays + weight*value
                        for l in range(len(steps_in) - 1):
                            weight = ((steps_in[l + 1] - steps_in[l])/(steps[j + 1] - steps[j]))
                            value = ((delays_avg[i].eval(steps_in[l]) + delays_avg[i].eval(steps_in[l + 1]))/2)
                            sum_delays = sum_delays + weight*value
                    if count_delays == 0:
                        sum_delays = ((delays_avg[i].eval(steps[j]) + delays_avg[i].eval(steps[j + 1]))/2)
                        count_delays = 1
                    value_1 = sum_delays*h[len(steps)*i + j]
                    value_2 = epsilon*(h[len(steps)*i + j])**2
                    sums  = sums + value_1 + value_2
            return sums 
        
        A = []
        for j in range(len(steps)):
            A.append([])
            for k in range(len(network.paths)):
                for g in range(len(steps)):
                    if j == g:
                        A[j].append(1)
                    else:
                        A[j].append(0)

        net_bound = []
        for i in range(len(steps)):
            net_bound.append(net_inflow.eval(steps[i]))
        constraint_1 = scipy.optimize.LinearConstraint(A, net_bound, net_bound)
        
        bounds = []
        start = []
        for i in range(len(network.paths)):
            for j in range(len(steps)):
                bounds.append((0, float("inf")))
                start.append(inflow_avg[i].eval(steps[j]))

        sol = scipy.optimize.minimize(obj, start, bounds = bounds, constraints = constraint_1)
    
        #3. Update the path inflows and run the edge-loading procedure
        inflows = []
        values = []
        inflow_dict = []
        for i in range(len(network.paths)):
            values.append([])
            for j in range(len(steps)):
                values[i].append(sol.x[len(steps)*i + j])

        for i in range(len(network.paths)):
            inflows.append(RightConstant(steps, values[i], (0, horizon)))
            inflow_dict.append((network.paths[i], inflows[i]))

        loader_new = NetworkLoader(network, inflow_dict)
        result_new = loader_new.build_flow()
        flow_new = next(result_new)
        counter_steps = counter_steps + 1
        #4. Update the population average and run the edge-loading procedure again
        old_avg = inflow_avg.copy()

        inflow_avg = []
        inflow_dict_avg = []
        for i in range(len(network.paths)):
            func_1 = inflows[i].mult_scalar((1/counter_steps))
            func_2 = old_avg[i].mult_scalar(((counter_steps - 1)/counter_steps))
            new_avg = func_1.__add__(func_2)
            inflow_avg.append(new_avg)
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

        print("Norm difference: " + str(diff))
        norm_differences.append(diff) 
        #6. Calculate the (regularized) gap for checking, if we get close to a dynamic equilibrium
        gap_steps = steps.copy()

        def obj_gap(h):
            sum_1 = 0
            for i in range(len(network.paths)):
                for j in range(len(gap_steps) - 1):
                    sum_delays = 0
                    count_delays = 0
                    steps_in = []
                    for k in range(len(delays_avg[i].times) - 1):
                        if delays_avg[i].times[k] > gap_steps[j] and delays_avg[i].times[k] < gap_steps[j+1]:
                            steps_in.append(delays_avg[i].times[k])
                            count_delays = count_delays + 1
                    if count_delays != 0:
                        weight = ((steps_in[0] - gap_steps[j])/(gap_steps[j + 1] - gap_steps[j]))
                        value = ((delays_avg[i].eval(steps_in[0]) + delays_avg[i].eval(gap_steps[j]))/2)
                        sum_delays = sum_delays + weight*value
                        weight = ((gap_steps[j + 1] - steps_in[-1])/(gap_steps[j + 1] - gap_steps[j]))
                        value = ((delays_avg[i].eval(steps_in[-1]) + delays_avg[i].eval(gap_steps[j + 1]))/2)
                        sum_delays = sum_delays + weight*value
                        for l in range(len(steps_in) - 1):
                            weight = ((steps_in[l + 1] - steps_in[l])/(gap_steps[j + 1] - gap_steps[j]))
                            value = ((delays_avg[i].eval(steps_in[l]) + delays_avg[i].eval(steps_in[l + 1]))/2)
                            sum_delays = sum_delays + weight*value
                    if count_delays == 0:
                        sum_delays = ((delays_avg[i].eval(gap_steps[j]) + delays_avg[i].eval(gap_steps[j + 1]))/2)
                        count_delays = 1
                    value_1 = sum_delays
                    value_2 = 2*epsilon*(inflow_avg[i].eval(gap_steps[j]))
                    value_3 = h[len(gap_steps)*i + j] - inflow_avg[i].eval(gap_steps[j])
                    sum_1 = sum_1 + (value_1 + value_2)*value_3
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
        start = []
        for i in range(len(network.paths)):
            for j in range(len(gap_steps)):
                bounds.append((0, float("inf")))
                start.append(inflow_avg[i].eval(gap_steps[j]))

        sol_gap = scipy.optimize.minimize(obj_gap, start, bounds=bounds, constraints=constraint_1)
        print("Gap: " + str((-1)*sol_gap.fun))
        gap_values.append((-1)*sol_gap.fun)
        
        if round(sol_gap.fun, 4) == 0:
            equilibrium_reached = True
            print("The empirical frequency has reached a regularized equilbrium")

        #7. We introduce another measure for closeness to a dynamic equilibrium, namely the storage function based on the very definition of
        #dynamic equilibrium (will have to be changed for regularized FP)
        theta = []
        for i in range(len(network.paths)):
            theta.append([])
            for j in range(len(network.paths)):
                theta[i].append([])
                for k in range(len(inflow_avg[i].times)):
                    theta[i][j].append(inflow_avg[i].times[k])
                for l in range(len(delays_avg[i].times)):
                    if delays_avg[i].times[l] not in theta[i][j]:
                        theta[i][j].append(delays_avg[i].times[l])
                for m in range(len(delays_avg[j].times)):
                    if delays_avg[j].times[m] not in theta[i][j]:
                        theta[i][j].append(delays_avg[j].times[m])
                for n in range(len(inflow_avg[j].times)):
                    if inflow_avg[j].times[n] not in theta[i][j]:
                        theta[i][j].append(inflow_avg[j].times[n])
                theta[i][j].sort()

        storage = 0
        for i in range(len(network.paths)):
            for j in range(len(network.paths)):
                if j != i:
                    for k in range(len(theta[i][j]) - 1):
                        start = theta[i][j][k] + eps
                        end = theta[i][j][k+1] - eps
                        value = 0
                        if end - start >= 100*eps and end <= horizon and start <= horizon:
                            diff_inf = inflow_avg[i] - inflow_avg[j]
                            diff_inf_1 = diff_inf.mult_scalar(2*epsilon)
                            diff_un = delays_avg[i] - delays_avg[j]
                            new_times = []
                            for o in range(len(diff_un.times)):
                                new_times.append(diff_un.times[o])
                            for p in range(len(diff_inf_1.times)):
                                if diff_inf_1.times[p] not in new_times:
                                    new_times.append(diff_inf_1.times[p])
        
                            new_times.sort()
                            new_values = []
                            for q in range(len(new_times)):
                                if new_times[q] <=  diff_inf_1.domain[1]:
                                    new_values.append(diff_un.eval(new_times[q]) + diff_inf_1.eval(new_times[q]))
                                else: 
                                    new_values.append(diff_un.eval(new_times[q]))

                            diff_un_1 = PiecewiseLinear(new_times,new_values,diff_un.first_slope,diff_un.last_slope,diff_un.domain)
                            diff = diff_un_1.restrict((start, end))
                            value1 = round(diff.eval(start),12)
                            value2 = round(diff.eval(end),12)
                            if value1 > 0 and value2 > 0:
                                if end not in diff.times:
                                    diff.times.append(end)
                                    diff.values.append(value2)
                                if start not in diff.times:
                                    diff.times.insert(0,start)
                                    diff.values.insert(0,value1)
                                value = inflow_avg[i].multiply(diff,start,end).integrate(start,end,False)
                            elif value1 <= 0 and value2 <= 0:
                                value = 0
                            elif value1 > 0 and value2 <= 0:
                                gradient = (value2 - value1)/(end - start)
                                point = start - value1/gradient
                                if point > end:
                                    point = end
                                if point not in diff.times:
                                    diff.times.append(point)
                                    diff.times.sort()
                                    diff.values.insert(diff.times.index(point),0)
                                if start not in diff.times:
                                    diff.times.insert(0,start)
                                    diff.values.insert(0,value1)
                                value = inflow_avg[i].multiply(diff,start,point).integrate(start,point,False)
                            else:
                                gradient = (value2 - value1)/(end - start)
                                point = start - value1/gradient
                                if point < start:
                                    point = start
                                if point not in diff.times:
                                    diff.times.append(point)
                                    diff.times.sort()
                                    diff.values.insert(diff.times.index(point),0)
                                if end not in diff.times:
                                    diff.times.append(end)
                                    diff.values.append(value2)
                                value = inflow_avg[i].multiply(diff,point,end).integrate(point,end,False)
                        storage = storage + value
        
        print("")
        print("Storage :" + str(storage))
        print("")
        storage_values.append(storage)

        if round(storage,4) == 0:
            equilibrium_reached = True
            print("the dynamics reached a dynamic equilibrium")
                    


    diff_delays_avg = []
    for i in range(len(network.paths)):
        for j in range(len(network.paths)):
            if i != j:
                diff_delays_avg.append(delays_avg[i] - delays_avg[j])

    if equilibrium_reached and accuracy_reached:
        print("The sequence of empirical frequencies converged in the given precision to a regularized equilibrium with epsilon = " + str(epsilon))
    elif equilibrium_reached:
        print("The sequence of empirical frequencies reached a regularized equilibrium with epsilon = " + str(epsilon) + " but not in the given accuracy")
    elif accuracy_reached:
        print("The sequence of empirical frequencies converged in the given precision, but not to a regularized equilibrium")
    else:
        print("The sequence of empirical frequencies neither converged nor reached a regularized equilibrium")

graph = DirectedGraph()
s = Node(0,graph)
v = Node(1,graph)
w = Node(2,graph)
x = Node(3,graph)
y = Node(4,graph)
t = Node(5,graph)
e1 = Edge(s,v,0,graph)
e2 = Edge(s,w,1,graph)
e3 = Edge(v,x,2,graph)
e4 = Edge(v,y,3,graph)
e5 = Edge(w,x,4,graph)
e6 = Edge(w,y,5,graph)
e7 = Edge(x,t,6,graph)
e8 = Edge(y,t,7,graph)

graph.nodes = {0:s,1:v,2:w,3:x,4:y,5:t}
graph.edges = [e1,e2,e3,e4,e5,e6,e7,e8]

capacities =   [1,1,1,1,1,1,2,2]
travel_times = [1,1,1,1,1,1,2,1]
net_inflow = RightConstant([0,5],[8,0],(0,5))

p1 = [e1,e3,e7]
p2 = [e1,e4,e8]
p3 = [e2,e5,e7]
p4 = [e2,e6,e8]

paths = [p1,p2,p3,p4]
horizon = 5
delta = 0.5
epsilon = 0
numSteps = 500
lamb = 0.0001

reg_fictitious_play(graph, capacities, travel_times,
                    net_inflow, paths, horizon, delta, epsilon, numSteps, lamb)