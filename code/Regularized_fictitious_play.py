#Here we want to implement the learning dynamics on general networks
from __future__ import annotations

import scipy.optimize

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
import nguyen_network
import sioux_falls_network
import time

#Graph contains the directed graph, cap the capacitites of the edges in the order of the edges in the edge list of the graph, similarly for travel containing the travel times,
#net_inflow is the right-constatn network inflow rate, paths is the list of s-t-paths in the network, horizon is the time horizon, delta is the size of the interval the path 
#inflows should be constant at, epsilon is the regularization weight, numSteps is the opitonal number of learning steps before termination, lamb is the optional accuracy to be 
#reached for termination, size the step size (after all we are dealing with the discretization of an ODE) for the Euler scheme (other schemes doable; choosing this not to be 1
# will probably slow down the convergence process, but will maybe enhance stability)
def reg_fictitious_play(graph: DirectedGraph, cap: List[float], travel: List[float],
                         net_inflow: RightConstant, paths: List[Path], horizon: float,
                         delta: float, epsilon: float, numSteps: int, lamb: float, size: float) -> List[RightConstant]:
    #Initialize the network, the commodity and various lists for saving values
    start_time = time.time()
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
    norm_differences = []
    values = []
    step_size = size

    #Initialization of the path flow (can vary); for example uniform or everything into the path with highest capcacity
    #Optional: search for path with highest capactiy
    caps = []
    for i in range (len(network.paths)):
        max_cap = 0
        for j in range(len(network.paths[i].edges)):
            if capacities[network.paths[i].edges[j].id] > max_cap:
                max_cap = capacities[network.paths[i].edges[j].id]
        caps.append(max_cap)

    max_path = paths[caps.index(max(caps))]

    #Here: initialization by path with highest capacity
    for i in range(len(network.paths)):
        values.append([])
        for j in range(len(net_inflow.times)):
            if network.paths[i] == max_path:
                values[i].append(net_inflow.values[j])
            else:
                values[i].append(0)

    #Here: uniform initialization
    #for i in range(len(network.paths)):
    #    values.append([])
    #    for j in range(len(net_inflow.times)):
    #        values[i].append((1/len(network.paths))*net_inflow.values[j])
            
    #Initialize inflow dictionary
    inflows = []
    inflow_dict = []
    for i in range(len(network.paths)):
        inflows.append(RightConstant(net_inflow.times, values[i], (0, horizon)))
        inflow_dict.append((network.paths[i], inflows[i]))

    #Initialize population average and the path delays according to the average path inflows
    loader_beg = NetworkLoader(network, inflow_dict)
    result_beg = loader_beg.build_flow()
    flow_beg = next(result_beg)
    delays_beg = loader_beg.path_delay(horizon)
    delays_avg = delays_beg.copy()
    inflow_avg = inflows.copy()
    
    #Set a counter for the number of steps (can be infinity) and flags, if the accuracy regarding the norm difference or the 
    #Lyapunov function has been reached (accuracy can be 0)
    counter_steps = 1
    accuracy_reached = False
    equilibrium_reached = False

    #Create set of intervals equidistantly according to delta (set of intervals can also be constructed differently)
    #and accordingly the set of time steps
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
        #Define the objective of the best response problem
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
        
        #Define the side-constraint securing the solution to be feasible w.r.t. the network inflow rate
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
        
        #Set the bounds such that the regularized best response is non-negative
        bounds = []
        start = []
        for i in range(len(network.paths)):
            for j in range(len(steps)):
                bounds.append((0, float("inf")))
                start.append(inflow_avg[i].eval(steps[j]))

        #Compute the regularized best response
        sol = scipy.optimize.minimize(obj, start, bounds = bounds, constraints = constraint_1)
    
        #Extract the solution of above best response problem into a dictionary of path inflows and run network loading
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
        #Update the population average and recompute the path delay according to it
        old_avg = inflow_avg.copy()

        inflow_avg = []
        inflow_dict_avg = []
        for i in range(len(network.paths)):
            func_1 = inflows[i].mult_scalar((size/counter_steps))
            func_2 = old_avg[i].mult_scalar(((counter_steps - size)/counter_steps))
            new_avg = func_1.__add__(func_2)
            inflow_avg.append(new_avg)
            inflow_dict_avg.append((network.paths[i], inflow_avg[i]))

        loader_avg = NetworkLoader(network, inflow_dict_avg)
        result_avg = loader_avg.build_flow()
        flow_avg = next(result_avg)
        delays_avg = loader_avg.path_delay(horizon)
        #Calculate the norm difference between the old population average and the new population average
        diff_inflows = []
        for i in range(len(network.paths)):
            curr_diff = inflow_avg[i] - old_avg[i]
            diff_inflows.append(curr_diff.abs_val())
    
        diff = 0
        for i in range(len(network.paths)):
            diff = diff + diff_inflows[i].integral().eval(horizon)
    
        #Check convergence w.r.t. the accuracy for the norm difference and save norm difference
        if round(diff, 12) <= lamb:
            accuracy_reached = True

        print("Norm difference: " + str(diff))
        norm_differences.append(diff) 
        #Compute the value of the Lyapunov function
        gap_steps = steps.copy()

        #Define the objective of the Lyapunov function
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

        #Define the feasibility set to be the set of feasible path inflows according to the given network inflow rate
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
        
        #Give the bounds such that the values of the minimizer of the above objective are non-negative and set the 
        #initial value for the computation of the Lyapunov function to be the population average (can be changed)
        bounds = []
        start = []
        for i in range(len(network.paths)):
            for j in range(len(gap_steps)):
                bounds.append((0, float("inf")))
                start.append(inflow_avg[i].eval(gap_steps[j]))

        #Save the value of the Lyapunov function
        sol_gap = scipy.optimize.minimize(obj_gap, start, bounds=bounds, constraints=constraint_1)
        print("Gap: " + str((-1)*sol_gap.fun))
        gap_values.append((-1)*sol_gap.fun)
        
        #Check convergence w.r.t. to the accuracy for the Lyapunov function
        if (-1)*sol_gap.fun <= lamb:
            equilibrium_reached = True
            print("The empirical frequency has reached a regularized equilbrium")

    end_time = time.time()
    print("Time take: " + str(end_time - start_time))

    size_support = 0
    for i in range(len(network.paths)):
        in_support = False
        for j in range(len(inflow_avg[i].values)):
            if round(inflow_avg[i].values[j], 6) > 0:
                in_support = True
        if in_support:
            size_support = size_support + 1
    print("Size of support: " + str(size_support))



    #Print the overall outcome of the computation
    if equilibrium_reached and accuracy_reached:
        print("The sequence of empirical frequencies converged in the given precision to a regularized equilibrium with epsilon = " + str(epsilon))
    elif equilibrium_reached:
        print("The sequence of empirical frequencies reached a regularized equilibrium with epsilon = " + str(epsilon) + " but not in the given accuracy")
    elif accuracy_reached:
        print("The sequence of empirical frequencies converged in the given precision, but not to a regularized equilibrium")
    else:
        print("The sequence of empirical frequencies neither converged nor reached a regularized equilibrium")

    #Return the computed population average
    return inflow_avg


#Initialize any network instance here
graph = DirectedGraph
s = Node(0,graph)
b = Node(1,graph)
c = Node(2,graph)
d = Node(3,graph)
e = Node(4,graph)
t = Node(5,graph)

e_1 = Edge(s,b,0,graph)
e_2 = Edge(s,c,1,graph)
e_3 = Edge(b,d,2,graph)
e_4 = Edge(b,e,3,graph)
e_5 = Edge(c,d,4,graph)
e_6 = Edge(c,e,5,graph)
e_7 = Edge(d,t,7,graph)
e_8 = Edge(e,t,7,graph)

graph.nodes = {0:s,1:b,2:c,3:d,4:e,5:t}
graph.edges = [e_1,e_2,e_3,e_4,e_5,e_6,e_7,e_8]

capacities = [1,1,1,1,1,1,2,2]
travel_times = [1,1,1,1,1,1,2,1]

net_inflow = RightConstant([0,10],[5,0],(0,10))
p_1 = Path([e_1,e_3,e_7])
p_2 = Path([e_1,e_4,e_8])
p_3 = Path([e_2,e_5,e_7])
p_4 = Path([e_2,e_6,e_8])
paths_in = [p_1,p_2,p_3,p_4]

horizon = 10
delta = 0.5
numSteps = 100000
lamb = 0.00001
epsilon = 0.05
size = 0.1

reg_fictitious_play(graph, capacities, travel_times,
                    net_inflow, paths_in, horizon, delta, epsilon, numSteps, lamb, size)