#Here we just want to implement the gap function (regularized and not)

from piecewise_linear import PiecewiseLinear
from right_constant import RightConstant
from typing import List,Tuple 
from machine_precision import eps
from network_loader import Path

#inflows is the path inflow vector of the flow, where we want to compute the gap function, path_delays is the vector of its path
#delay operator and epsilon is the regularization weight
def gap_function(inflows: List[RightConstant], path_delays: List[PiecewiseLinear], inflows_other: List[RightConstant], horizon: float, espilon: float = 0):
    diff_inflows = []
    for i in range(len(inflows)):
        diff_inflows.append(inflows[i] - inflows_other[i])
    
    steps = []
    for i in range(len(inflows)):
        steps.append([])
        for j in range(len(diff_inflows[i].times)):
            if diff_inflows[i].times[j] <= horizon and diff_inflows[i].times[j] not in steps[i] :
                steps[i].append(diff_inflows[i].times[j])
        for k in range(len(path_delays[i].times)):
            if path_delays[i].times[k] not in steps[i] and path_delays[i].times[k] <= horizon :
                steps[i].append(path_delays[i].times[k])
    
    for i in range(len(inflows)):
        steps[i].sort()

    integrals = []
    for i in range(len(inflows)):
        integrals.append(0)

    for i in range(len(inflows)):
        for j in range(len(steps[i]) - 1):
            start = steps[i][j]
            end = steps[i][j+1] - 2*eps
            value = diff_inflows[i].multiply(path_delays[i], start, end).integrate(start, end, True)
            integrals[i] = integrals[i] + value

    scalar_product = sum(integrals)
    if abs(scalar_product) < eps:
        scalar_product = 0
    print(str(scalar_product))
    return scalar_product