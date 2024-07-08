#Here I want to try to implement the following: for a given network with given capcities 
#given travel times and a given (constant) network inflow rate and a given time horizon
#I want to minimize the scalar product of the difference of two delay operator vector and the 
#difference of two path inflow vectors in the given situation with the discretization as
#in the best response problem

from graph import *
from network import Network,Commodity
from dynamic_flow import DynamicFlow 
from network_loader import NetworkLoader, Path
from piecewise_linear import PiecewiseLinear, identity
from right_constant import RightConstant
from machine_precision import eps
import scipy
import math

def minimize_monotone(graph: DirectedGraph, capacities: List[float], travel_times: List[float],
                      net_inflow: RightConstant, paths: List[Path], 
                      horizon: float, delta: float) -> float:
    network = Network()
    network.graph = graph
    network.capacity = capacities
    network.travel_time = travel_times
    network_inflow = Commodity({graph.nodes[0]: net_inflow}, graph.nodes[len(graph.nodes) - 1], 1)
    network.commodities = [network_inflow]
    network.paths = paths
    
    break_points = []
    counter = 0
    while counter*delta < horizon - delta:
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

        
        
    return 0