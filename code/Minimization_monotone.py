#Here I want to try to implement the following: for a given network with given capcities 
#given travel times and a given (constant) network inflow rate and a given time horizon
#I want to minimize the scalar product of the difference of two delay operator vector and the 
#difference of two path inflow vectors in the given situation with the discretization as
#in the best response problem

import scipy.optimize
from graph import *
from network import Network,Commodity
from dynamic_flow import DynamicFlow 
from network_loader import NetworkLoader, Path
from piecewise_linear import PiecewiseLinear, identity
from right_constant import RightConstant
from machine_precision import eps
import scipy
import math
import arrays

def minimize_monotone(graph: DirectedGraph, capacities: List[float], travel_times: List[float],
                      net_inflow: RightConstant, paths: List[Path], 
                      horizon: float, delta: float) -> float:
    network = Network()
    network.graph = graph
    network.capacity = capacities
    network.travel_time = travel_times
    network_inflow = Commodity({s: net_inflow}, t, 1)
    network.commodities = [network_inflow]
    network.paths = paths
    
    break_points = []
    counter = 0
    while counter*delta <= horizon - delta:
        break_points.append(counter*delta)
        counter = counter + 1
    break_points.append(horizon)

    def obj(h):
        h_1 = []
        for i in range(len(network.paths)):
            for j in range(len(break_points)):
                h_1.append(h[(len(break_points) - 1)*i+j])
        h_2 = []
        for i in range(len(network.paths)):
            for j in range(len(break_points)):
                h_2.append(h[len(network.paths)*len(break_points) + (len(break_points) - 1)*i + j])
        
        
        values_1 = []
        for i in range(len(network.paths)):
            values_1.append([])
            for j in range(len(break_points)):
                values_1[i].append(h[(len(break_points) - 1)*i+j])
        values_2 = []
        for i in range(len(network.paths)):
            values_2.append([])
            for j in range(len(break_points)):
                values_2[i].append(h[len(network.paths)*len(break_points) + (len(break_points) - 1)*i + j])
        
        
        inflow_funcs_1 = []
        for i in range(len(network.paths)):
            inflow_funcs_1.append(RightConstant(break_points,values_1[i],(0,horizon)))
        inflow_funcs_2 = []
        for i in range(len(network.paths)):
            inflow_funcs_2.append(RightConstant(break_points, values_2[i],(0,horizon)))
        
        
        inflow_dict_1 = []
        for i in range(len(network.paths)):
            inflow_dict_1.append((network.paths[i], inflow_funcs_1[i]))
        inflow_dict_2 = []
        for i in range(len(network.paths)):
            inflow_dict_2.append((network.paths[i], inflow_funcs_2[i]))


        loader_1 = NetworkLoader(network, inflow_dict_1)
        loader_2 = NetworkLoader(network, inflow_dict_2)

        result_1 = loader_1.build_flow()
        result_2 = loader_2.build_flow()

        flow_1 = next(result_1)
        flow_2 = next(result_2)

        delays_1 = loader_1.path_delay(horizon)
        delays_2 = loader_2.path_delay(horizon)

        diff_delays = []
        for i in range(len(network.paths)):
            diff_delays.append(delays_1[i] - delays_2[i])

        diff_inflows = []
        for i in range(len(network.paths)):
            diff_inflows.append(inflow_dict_1[i] - inflow_dict_2[i])

        steps = []
        for i in range(len(network.paths)):
            steps.append([])
            for j in range(len(diff_delays[i].times)):
                steps[i].append(diff_delays[i].times[j])
            for k in range(len(break_points)):
                if diff_inflows[i].times[k] not in steps[i]:
                    steps[i].append(diff_inflows[i].times[k])
        
        for i in range(len(network.paths)):
            steps[i].sort()

        sums = []
        for i in range(len(network.paths)):
            sums.append(0)

        for i in range(len(network.paths)):
            for j in range(len(steps[i]) - 1):
                value_s = diff_delays[i].eval(steps[i][j])
                value_e = diff_delays[i].eval(steps[i][j+1])
                value_p = (value_s + value_e)/2
                index_step = arrays.elem_lrank(break_points, steps[i][j])
                value_i = diff_inflows[i].eval[steps[i[j]]]
                sums[i] = sums[i] + value_p*value_i
        
        return sum(sums)
    
    A = []
    for i in range(2):
        for j in range(len(break_points) - 1):
            A.append([])
            for h in range(2):
                for k in range(len(network.paths)):
                    for g in range(len(break_points)):
                        if j == g and h == i:
                            A[2*i+j].append(1)
                        else:
                            A[2*i+j].append(0)

    B = []
    for i in range(2):
        for j in range(len(network.paths)):
            B.append([])
            for h in range(2):
                for k in range(len(network.paths)):
                    for g in range(len(break_points)):
                        if g == len(break_points) - 1 and h == i and k == j:
                            B[2*i + j].append(1)
                        else:
                            B[2*i + j].append(0)

    constraint_1 = scipy.optimize.LinearConstraint(A, net_inflow.values[0], net_inflow.values[0])
    constraint_2 = scipy.optimize.LinearConstraint(B, 0, 0)
    bounds = []
    for i in range(2):
        for j in range(len(network.paths)):
            for k in range(len(break_points)):
                bounds.append((0, None))
        
    min_cap = []
    for i in range(len(network.paths)):
        caps = []
        for j in range(len(network.paths[i])):
            index = network.graph.edges.index(network.paths[i][j])
            caps.append(network.capacity[index])
        min_cap.append(min(caps))

    path_min = min_cap.index(min(min_cap))
    path_max = min_cap.index(max(min_cap))
    start = []
    for i in range(len(network.paths)):
        for j in range(len(break_points)):
            if i == path_min and j == 0:
                start.append(net_inflow.values[0])
            else:
                start.append(0)
    for i in range(len(network.paths)):
        for j in range(len(break_points)):
            if i == path_max and j == 0:
                start.append(net_inflow.values[0])
            else: 
                start.append(0)
        
    constraints = [constraint_1, constraint_2]
    sol = scipy.optimize.minimize(obj, x0 = start, bounds = bounds, constraints = constraints )


    print(str(sol.obj))
    return sol.obj

graph = DirectedGraph()
s = Node(0, graph)
v = Node(1, graph)
t = Node(2, graph)

e_1 = Edge(s, v, 0, graph)
e_2 = Edge(s, v, 1, graph)
e_3 = Edge(v, t, 2, graph)

graph.nodes = {s:0,v:1,t:2}
graph.edges = [e_1,e_2,e_3]

capacities = [1,3,2]
travel_times = [1,0,0]

horizon = 2
network_inflow = RightConstant([0,2],[2,0],(0,2))
delta = 1

p_1 = [e_1, e_3]
p_2 = [e_2, e_3]

paths = [p_1, p_2]
minimize_monotone(graph, capacities, travel_times,network_inflow, paths, horizon, delta)
