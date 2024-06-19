#Here we want to implement the learning dynamics on general networks: for that we make use of the already implemented code: we need a network (and a graph?) we need to somehow get all the paths from the origin to the destination, we need an edge_loading procedure and thats basically it

import graph
import network as net
import piecewise_linear as pl
import right_constant as rc
import numpy as np

#graph is the directed graph that forms the network, we only assume one commodity and we assume its source is always indexed by 0 and its sink always has the highest index, cap contains the capacities of the edges and travel the travel times, the order on cap and travel has to be the same as the indices in the graph, u is either just a constant or a right constant function, T is the time horizon, delta is the length of an interval, epsilon is the weight of the regularizer, numSteps is the number of learning steps that we want to have and lamb is the accuracy we want to achieve in the convergence
def reg_fictitious_play(graph, cap, travel, u, T, delta, epsilon, numSteps, lamb):
    #Steps that need to be taken:
    #1. Initialization: maybe random or maybe also distribute the inflow according to the minimum capacitiy on a path
    #2. Best-response problem: this should be not so different from that that we got for the parallel link networks
    #3. Update the flow and run the edge-loading procedure
    #4. Update the population average and run the edge-loading procedure again
    #5. Calculate the difference of the new flow and the last flow for checking convergence
    #6. Calculate the (regularized) gap for checking, if we get close to a dynamic equilibrium
