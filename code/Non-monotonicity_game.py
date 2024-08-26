#Here we try to implement a two player game that should lead to the computation of two path inflows
#such that the monotonicity condition is violated.

from graph import Edge,Node,DirectedGraph
from right_constant import RightConstant
from piecewise_linear import PiecewiseLinear
from machine_precision import eps
from network_loader import NetworkLoader, Path
from network import Network, Commodity
from dynamic_flow import DynamicFlow
from typing import List
from arrays import *
import scipy.optimize
from testing_monotonicity import monotonicity_check


def non_monotonicity_game(graph: DirectedGraph, cap: List[float], travel: List[float], net_inflow: RightConstant,
                          paths: List[Path], delta: float, horizon: float, epsilon: float):
    #Steps that need to be taken:
    #1. Initialize the network
    network = Network()
    network.graph = graph
    network.capacity = cap
    network.travel_time = travel
    s = graph.nodes[0]
    t = graph.nodes[len(graph.nodes) - 1]
    network_inflow = Commodity({s: net_inflow}, t, 1)
    network.commodities = [network_inflow]
    network.paths = paths

    #2. Initialize the two flows
    min_cap = []
    for i in range(len(network.paths)):
        caps = []
        for j in range(len(network.paths[i])):
            index = network.paths[i][j].id
            caps.append(network.capacity[index])
        min_cap.append(min(caps))

    steps = []
    stop = 0
    while stop*delta <= horizon - delta:
        steps.append(stop*delta)
        stop = stop + 1
    steps.append(horizon)

    for i in range(len(net_inflow.times)):
        if net_inflow.times[i] not in steps:
            steps.append(net_inflow[i])
    steps.sort()

    values_1 = []
    values_2 = []
    for i in range(len(network.paths)):
        values_1.append([])
        values_2.append([])
        for j in range(len(steps)):
            if i == min_cap.index(min(min_cap)):
                values_1[i].append(net_inflow.eval(steps[j]))
            else:
                values_1[i].append(0)
        for k in range(len(steps)):
            if i == min_cap.index(max(min_cap)):
                values_2[i].append(net_inflow.eval(steps[j]))
            else:
                values_2[i].append(0)

    inflows_min = []
    inflows_max = []
    inflow_dict_min = []
    inflow_dict_max = []
    for i in range(len(network.paths)):
        inflows_min.append(RightConstant(steps,values_1[i],(0,horizon)))
        inflows_max.append(RightConstant(steps,values_2[i],(0,horizon)))
        inflow_dict_min.append((network.paths[i],inflows_min[i]))
        inflow_dict_max.append((network.paths[i],inflows_max[i]))

    loader_min = NetworkLoader(network,inflow_dict_min)
    loader_max = NetworkLoader(network,inflow_dict_max)
    result_min = loader_min.build_flow()
    result_max = loader_max.build_flow()
    flow_min = next(result_min)
    flow_max = next(result_max)
    delays_min = loader_min.path_delay(horizon)
    delays_max = loader_max.path_delay(horizon)
    non_mono_reached = False
    
    while not non_mono_reached:
        #3. Best response of first flow
        def obj_min(h):
            sums = 0
            for i in range(len(network.paths)):
                for j in range(len(steps) - 1):
                    sum_delays = 0
                    count_delays = 0
                    for k in range(len(delays_max[i].times) - 1):
                        length = 0
                        if delays_max[i].times[k] >= steps[j] and delays_max[i].times[k] <= steps[j+1]:
                            if delays_max[i].times[k + 1] > steps[j + 1]:
                                end = steps[j + 1]
                            else:
                                end = delays_max[i].times[k + 1]
                            start = delays_max[i].times[k]
                            length = end - start 
                            average = ((delays_max[i].eval(end) + delays_max[i].eval(start))/2)
                            sum_delays = sum_delays + (length/(steps[j + 1] - steps[j]))*average
                            count_delays = count_delays + 1
                    if count_delays == 0:
                        sum_delays = ((delays_max[i].eval(steps[j]) + delays_max[i].eval(steps[j + 1]))/2)
                        count_delays = 1
                    value_1 = (sum_delays/count_delays)*(inflows_max[i].eval(steps[j]) - h[len(steps)*i +j])
                    value_2 = epsilon*(h[len(steps)*i + j] - inflows_min[i].eval(steps[j]))**2
                    sums  = sums + value_1 - value_2
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
                start.append(inflows_min[i].eval(steps[j]))

        sol_min = scipy.optimize.minimize(obj_min, start, bounds = bounds, constraints = constraint_1)
        
        #4. Update the first flow to respond
        values_1 = []
        inflows_min = []
        inflow_dict_min = []
        for i in range(len(network.paths)):
            values_1.append([])
            for j in range(len(steps)):
                values_1[i].append(sol_min.x[len(steps)*i + j])

        for i in range(len(network.paths)):
            inflows_min.append(RightConstant(steps,values_1[i],(0,horizon)))
            inflow_dict_min.append((network.paths[i],inflows_min[i]))

        loader_min = NetworkLoader(network,inflow_dict_min)
        result_min = loader_min.build_flow()
        flow_min = next(result_min)
        delays_min = loader_min.path_delay(horizon)

        #5. Best response of the second flow
        def obj_max(h):
            sums = 0
            for i in range(len(network.paths)):
                for j in range(len(steps) - 1):
                    sum_delays = 0
                    count_delays = 0
                    for k in range(len(delays_min[i].times) - 1):
                        length = 0
                        if delays_min[i].times[k] >= steps[j] and delays_min[i].times[k] <= steps[j+1]:
                            if delays_min[i].times[k + 1] > steps[j + 1]:
                                end = steps[j + 1]
                            else:
                                end = delays_min[i].times[k + 1]
                            start = delays_min[i].times[k]
                            length = end - start 
                            average = ((delays_min[i].eval(end) + delays_min[i].eval(start))/2)
                            sum_delays = sum_delays + (length/(steps[j + 1] - steps[j]))*average
                            count_delays = count_delays + 1
                    if count_delays == 0:
                        sum_delays = ((delays_min[i].eval(steps[j]) + delays_min[i].eval(steps[j + 1]))/2)
                        count_delays = 1
                    value_1 = (sum_delays/count_delays)*(h[len(steps)*i + j] - inflows_min[i].eval(steps[j]))
                    value_2 = epsilon*(h[len(steps)*i + j] - inflows_max[i].eval(steps[j]))**2
                    sums  = sums + value_1 - value_2
            return -sums 
        
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
                start.append(inflows_max[i].eval(steps[j]))

        sol_max = scipy.optimize.minimize(obj_max, start, bounds = bounds, constraints = constraint_1)

        #6. Update the second flow to respond
        values_2 = []
        inflows_max = []
        inflow_dict_max = []
        for i in range(len(network.paths)):
            values_2.append([])
            for j in range(len(steps)):
                values_2[i].append(sol_max.x[len(steps)*i + j])

        for i in range(len(network.paths)):
            inflows_max.append(RightConstant(steps,values_2[i],(0,horizon)))
            inflow_dict_max.append((network.paths[i],inflows_max[i]))

        loader_max = NetworkLoader(network,inflow_dict_max)
        result_max = loader_max.build_flow()
        flow_max = next(result_max)
        delays_max = loader_max.path_delay(horizon)

        #7. Evaluation of the variational inequality
        monotonicity_gap = monotonicity_check(graph, cap, travel, net_inflow, horizon, paths, inflows_max, inflows_min)
        print(str(monotonicity_gap))

        if monotonicity_gap < -1000*eps:
            non_mono_reached = True

    print(str(values_1))
    print(str(values_2))
    return [values_1, values_2]

graph = DirectedGraph()
s = Node(0,graph)
v = Node(1,graph)
t = Node(2,graph)
e_1 = Edge(s,v,0,graph)
e_2 = Edge(s,v,1,graph)
e_3 = Edge(v,t,2,graph)

graph.nodes = {0:s,1:v,2:t}
graph.edges = [e_1,e_2,e_3]
graph.reversed = False

capacities = [1,3,2]
travel_times = [1,0,0]
net_inflow = RightConstant([0,1,1.75,2],[2.5,1,3,0],(0,2))

path_1 = [e_1,e_3]
path_2 = [e_2,e_3]
paths = [path_1, path_2]

non_monotonicity_game(graph,capacities,travel_times,net_inflow,paths,0.25,2,0.2)