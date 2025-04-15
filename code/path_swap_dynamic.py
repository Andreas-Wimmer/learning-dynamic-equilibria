#Here we implement the path swap dynamic by Richard Mounce
from __future__ import annotations

from graph import DirectedGraph, Node, Edge
from network import Network, Commodity
from dynamic_flow import DynamicFlow
from piecewise_linear import PiecewiseLinear 
from right_constant import RightConstant 
from network_loader import NetworkLoader, Path
from machine_precision import eps
from typing import List
import math
from arrays import *

def path_swap(graph: DirectedGraph, cap: List[float], travel: List[float], paths: List[Path], horizon: float, net_inflow: RightConstant,
               delta: float, numSteps: int, lamb: float) -> List[RightConstant]:
    #Initialization
    network = Network()
    network.graph = graph
    network.capacity = cap
    network.travel_time = travel
    s = graph.nodes[0]
    t = graph.nodes[len(graph.nodes) - 1]
    network_inflow = Commodity({s, net_inflow}, t, 1)
    network.commodities = [network_inflow]
    network.paths = paths
    storMou_values = []
    norm_differences = []

    breaks_net_inflow = net_inflow.times
    values = []

    caps = []
    for i in range (len(network.paths)):
        max_cap = 0
        for j in range(len(network.paths[i])):
            if cap[network.paths[i][j].id] > max_cap:
                max_cap = cap[network.paths[i][j].id]
        caps.append(max_cap)

    max_path = paths[caps.index(max(caps))]
    
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

    for i in range(len(network.paths)):
        values.append([])
        for j in range(len(steps)):
            if network.paths[i] == max_path:
                values[i].append(net_inflow.eval(steps[j]))
            else:
                values[i].append(0)
    
    inflows = []
    inflow_dict = []
    for i in range(len(network.paths)):
        inflows.append(RightConstant(steps, values[i], (0, horizon)))
        inflow_dict.append((network.paths[i], inflows[i]))

    inflow = inflows.copy()
    loader = NetworkLoader(network, inflow_dict)
    result = loader.build_flow()
    flow = next(result)
    delays = loader.path_delay(horizon)
    
    counter_steps = 1
    accuracy_reached = False
    equilibrium_reached = False
    
    while counter_steps < numSteps and not accuracy_reached and not equilibrium_reached:
        #Compute Update
        update = []
        for j in range(len(steps) - 1):
            current = []
            for z in range(len(network.paths)):
                current.append(0)
            for p in range(len(network.paths)):
                for q in range(len(network.paths)):
                    sum_delays_1 = 0
                    sum_delays_2 = 0
                    count_delays_1 = 0
                    count_delays_2 = 0
                    steps_in_1 = []
                    steps_in_2 = []
                    for k in range(len(delays[p].times) - 1):
                        if delays[p].times[k] > steps[j] and delays[p].times[k] < steps[j+1]:
                            steps_in_1.append(delays[p].times[k])
                            count_delays_1 = count_delays_1 + 1
                    if count_delays_1 != 0:
                        weight = ((steps_in_1[0] - steps[j])/(steps[j + 1] - steps[j]))
                        value = ((delays[p].eval(steps_in_1[0]) + delays[p].eval(steps[j]))/2)
                        sum_delays_1 = sum_delays_1 + weight*value
                        weight = ((steps[j + 1] - steps_in_1[-1])/(steps[j + 1] - steps[j]))
                        value = ((delays[p].eval(steps_in_1[-1]) + delays[p].eval(steps[j + 1]))/2)
                        sum_delays_1 = sum_delays_1 + weight*value
                        for l in range(len(steps_in_1) - 1):
                            weight = ((steps_in_1[l + 1] - steps_in_1[l])/(steps[j + 1] - steps[j]))
                            value = ((delays[p].eval(steps_in_1[l]) + delays[p].eval(steps_in_1[l + 1]))/2)
                            sum_delays_1 = sum_delays_1 + weight*value
                    if count_delays_1 == 0:
                        sum_delays_1 = ((delays[p].eval(steps[j]) + delays[p].eval(steps[j + 1]))/2)
                        count_delays_1 = 1
                    value_1 = sum_delays_1
                    for k in range(len(delays[q].times) - 1):
                        if delays[q].times[k] > steps[j] and delays[q].times[k] < steps[j+1]:
                            steps_in_2.append(delays[q].times[k])
                            count_delays_2 = count_delays_2 + 1
                    if count_delays_2 != 0:
                        weight = ((steps_in_2[0] - steps[j])/(steps[j + 1] - steps[j]))
                        value = ((delays[q].eval(steps_in_2[0]) + delays[q].eval(steps[j]))/2)
                        sum_delays_2 = sum_delays_2 + weight*value
                        weight = ((steps[j + 1] - steps_in_2[-1])/(steps[j + 1] - steps[j]))
                        value = ((delays[q].eval(steps_in_2[-1]) + delays[q].eval(steps[j + 1]))/2)
                        sum_delays_2 = sum_delays_2 + weight*value
                        for l in range(len(steps_in_2) - 1):
                            weight = ((steps_in_2[l + 1] - steps_in_2[l])/(steps[j + 1] - steps[j]))
                            value = ((delays[q].eval(steps_in_2[l]) + delays[q].eval(steps_in_2[l + 1]))/2)
                            sum_delays_2 = sum_delays_2 + weight*value
                    if count_delays_2 == 0:
                        sum_delays_2 = ((delays[q].eval(steps[j]) + delays[q].eval(steps[j + 1]))/2)
                        count_delays_2 = 1
                    value_2 = sum_delays_2
                    value_4 = (value_1 - value_2)
                    value_5 = 0
                    if value_4 < 0:
                            value_5 = 0
                    else:
                            value_5 = value_4
                    value6 = inflow[p].eval(steps[j])*value_5
                    for y in range(len(network.paths)):
                        if y == p:
                            current[y] = current[y] - value6
                        elif y == q: 
                            current[y] = current[y] + value6
                        else:
                            current[y] = current[y] + 0
            update.append(current)
        update.append([0 for i in range(len(network.paths))])

        #Update flow
        inflow_old = inflow.copy()
        inflow_dict_old = []
        for i in range(len(network.paths)):
            inflow_dict_old.append((network.paths[i], inflow_old[i]))
        loader_old = NetworkLoader(network, inflow_dict_old)
        result_old = loader_old.build_flow()
        flow_old = next(result_old)
        delays_old = loader_old.path_delay(horizon)

        values_new = []
        for i in range(len(network.paths)):
            values_new.append([])
            for j in range(len(steps)):
                values_new[i].append(inflow[i].values[j] + update[j][i])

        inflows_new = []
        for i in range(len(network.paths)):
            inflows_new.append(RightConstant(steps, values_new[i],(0,horizon)))

        inflow = inflows_new.copy()
        new_dict = []
        for i in range(len(network.paths)):
            new_dict.append((network.paths[i],inflow[i]))
        loader = NetworkLoader(network, new_dict)
        result = loader.build_flow()
        flow = next(result)
        delays = loader.path_delay(horizon)
        #Check norm difference
        current_difference = 0
        for i in range(len(network.paths)):
            curr_path = 0
            for j in range(len(steps) - 1):
                diff = inflow[i].values[j] - inflow_old[i].values[j]
                if diff >= 0:
                    diff = diff
                else:
                    diff = (-1)*diff
                squared = diff*diff
                integral = squared*(steps[j+1] - steps[j])
                curr_path = curr_path + integral
            curr_path = math.sqrt(curr_path)
            current_difference = current_difference + curr_path

        norm_differences.append(current_difference)
        print("Current norm difference: " + str(current_difference))
        if current_difference <= lamb:
            accuracy_reached = True
        #Check value of Lyapunov function
        storage = 0
        for p in range(len(network.paths)):
            for q in range(len(network.paths)):
                for t in range(len(steps) - 1):
                    sum_delays_1 = 0
                    sum_delays_2 = 0
                    count_delays_1 = 0
                    count_delays_2 = 0
                    steps_in_1 = []
                    steps_in_2 = []
                    for k in range(len(delays[p].times) - 1):
                        if delays[p].times[k] > steps[t] and delays[p].times[k] < steps[t+1]:
                            steps_in_1.append(delays[p].times[k])
                            count_delays_1 = count_delays_1 + 1
                    if count_delays_1 != 0:
                        weight = ((steps_in_1[0] - steps[t])/(steps[t + 1] - steps[t]))
                        value = ((delays[p].eval(steps_in_1[0]) + delays[p].eval(steps[t]))/2)
                        sum_delays_1 = sum_delays_1 + weight*value
                        weight = ((steps[t + 1] - steps_in_1[-1])/(steps[t + 1] - steps[t]))
                        value = ((delays[p].eval(steps_in_1[-1]) + delays[p].eval(steps[t + 1]))/2)
                        sum_delays_1 = sum_delays_1 + weight*value
                        for l in range(len(steps_in_1) - 1):
                            weight = ((steps_in_1[l + 1] - steps_in_1[l])/(steps[t + 1] - steps[t]))
                            value = ((delays[p].eval(steps_in_1[l]) + delays[p].eval(steps_in_1[l + 1]))/2)
                            sum_delays_1 = sum_delays_1 + weight*value
                    if count_delays_1 == 0:
                        sum_delays_1 = ((delays[p].eval(steps[t]) + delays[p].eval(steps[t + 1]))/2)
                        count_delays_1 = 1
                    value_1 = sum_delays_1
                    for k in range(len(delays[q].times) - 1):
                        if delays[q].times[k] > steps[t] and delays[q].times[k] < steps[t+1]:
                            steps_in_2.append(delays[q].times[k])
                            count_delays_2 = count_delays_2 + 1
                    if count_delays_2 != 0:
                        weight = ((steps_in_2[0] - steps[t])/(steps[t + 1] - steps[t]))
                        value = ((delays[q].eval(steps_in_2[0]) + delays[q].eval(steps[t]))/2)
                        sum_delays_2 = sum_delays_2 + weight*value
                        weight = ((steps[t + 1] - steps_in_2[-1])/(steps[t + 1] - steps[t]))
                        value = ((delays[q].eval(steps_in_2[-1]) + delays[q].eval(steps[t + 1]))/2)
                        sum_delays_2 = sum_delays_2 + weight*value
                        for l in range(len(steps_in_2) - 1):
                            weight = ((steps_in_2[l + 1] - steps_in_2[l])/(steps[t + 1] - steps[t]))
                            value = ((delays[q].eval(steps_in_2[l]) + delays[q].eval(steps_in_2[l + 1]))/2)
                            sum_delays_2 = sum_delays_2 + weight*value
                    if count_delays_2 == 0:
                        sum_delays_2 = ((delays[q].eval(steps[t]) + delays[q].eval(steps[t + 1]))/2)
                        count_delays_2 = 1
                    value_2 = sum_delays_2
                    value_4 = (value_1 - value_2)
                    value_5 = 0
                    if value_4 < 0:
                        value_5 = 0
                    else:
                        value_5 = value_4*value_4
                    value6 = inflow[p].eval(steps[t])*value_5
                    storage = storage + value6
        
        storMou_values.append(storage)
        print("Current storage value: " + str(storage))
        if storage <= lamb:
            equilibrium_reached = True
        counter_steps = counter_steps + 1

    #Check convergence and output flow
    if counter_steps > numSteps:
        print("Number of learning steps reached")
    if equilibrium_reached or accuracy_reached:
        print("Convergence to an equilibrium")
    else: 
        print("No convergence in the given number of steps")
    
    return inflow

graph = DirectedGraph()
s = Node(0,graph)
u = Node(1,graph)
t = Node(2,graph)
e1 = Edge(s,u,0,graph)
e2 = Edge(s,u,1,graph)
e3 = Edge(u,t,2,graph)

graph.nodes = {0:s,1:u,2:t}
graph.edges = [e1,e2,e3]

capacities = [1,3,2]
travel_times = [1,1,1]
net_inflow = RightConstant([0,1,1.75,2],[2.5,1,4,0],(0,2))

p1 = [e1,e3]
p2 = [e2,e3]

paths = [p1,p2]
horizon = 2
delta = 0.25
numSteps = 100
lamb = 0.00001

path_swap(graph, capacities, travel_times, paths, horizon, net_inflow, delta, numSteps, lamb)
