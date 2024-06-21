#Here we implement the learning dynamics of fictitious play in a dynamic network to see, if it
#converges and if it converges against a dynamic equilibrium.

import utilities_c as u_c
import math
import sympy as sy
import scipy as sci
import numpy as np
import pdb
import matplotlib.pyplot as plt
import edge_loading_pln as elp

pdb.set_trace()

#We assume a parallel link network: m is the number of edges, cap contains the capacities of the
#edges, travel contains the travel times of the edges, T gives the time horizon, u gives the
#(constant) network inflow rate, delta could give the size of an interval in the discretization
#of the best-response problem, epsilon gives the proximity of convergence and numsteps could give
#the maximal amount of steps we make and lambda the precision of convergence we want to have
def learning_pln(m, cap, travel, T, u, delta, epsilon, numsteps, lamb):
    #1. Initialization
    E = []
    for i in range(m):
        E.append(["s","t"])

    #if we do not give a number of steps we want to reach the given accuracy for sure and if we do not give a
    #value for lambda we want exact convergence or want so see how close we get in the given number of steps
    if numsteps == None:
        numsteps = math.inf

    if lamb == None:
        lamb = 0

    #for the initialization we just distribute the flow over the edges according to their
    #capacity for the whole interval
    sum_cap = sum(cap)
    values = []
    for i in range(m):
        values.append((cap[i]/sum_cap)*u)

    [flow, queues, pathop] = elp.edge_loading_pln(m, cap, travel, [0,T], values,T)

    #breakpoint()
    #We need to store the time average flow; in the beginning it is just the same as our
    #initial flow
    flow_avg = flow.copy()
    queues_avg = queues.copy()
    pathop_avg = pathop.copy()
    counter_steps = 1
    while counter_steps < numsteps:
        stop = 0
        step_points = [];
        while stop < T - delta:
            step_points.append(stop)
            stop = stop + delta

        step_points.append(T)

        breaks = []
        for i in range(m):
            for j in range(len(pathop[i].segmentBorders)):
                breaks.append(pathop[i].segmentBorders[j])

        breaks.sort()
        for i in range(len(breaks)):
            if breaks[i] not in step_points:
                counter = 0
                right = 0
                for i in range(len(step_points)):
                    if step_points[counter] > breaks[i]:
                        right = counter
                        break
                    counter = counter + 1
                step_points.insert(right-1,breaks[i])

        #breakpoint()
        #2. Best-response problem
        #At first we try it out with scipy and if that does not work, then we change to sympy
        #Objective function
        def f(h):
            sum_1 = 0
            sum_2 = 0
            for i in range(m):
                for j in range(len(step_points) - 1):
                    val_1 = pathop_avg[i].getValueAt(step_points[j+1])
                    val_2 = pathop_avg[i].getValueAt(step_points[j])
                    sum_1 = sum_1 + ((val_1+val_2)/2)*h[(len(step_points)-1)*i+j]
                    sum_2 = sum_2 + epsilon*(h[(len(step_points)-1)*i+j] - flow_avg[i].getValueAt(step_points[j]))**2
            return sum_1 + sum_2

        A = []
        for i in range(len(step_points) - 1):
            A.append([])
            for j in range(m):
                for k in range(len(step_points) - 1):
                    if k == i:
                        A[i].append(1)
                    else:
                        A[i].append(0)

        constraints = sci.optimize.LinearConstraint(A, lb=u, ub=u)
        bounds = []
        for i in range(m):
            for j in range(len(step_points) - 1):
                bounds.append((0,None))

        start = []
        for i in range(m):
            for j in range(len(step_points) - 1):
                start.append(flow[i].getValueAt(step_points[j]))
        h0 = start.copy()
        sol = sci.optimize.minimize(f, h0, bounds=bounds, constraints=constraints)
        #breakpoint()
        #3. Updating the flow
        old_flow = flow.copy()
        old_queues = queues.copy()
        old_pathop = pathop.copy()
        [flow, queues, pathop] = elp.edge_loading_pln(m, cap, travel, step_points, sol.x, T)

        counter_steps = counter_steps + 1
        #breakpoint()
        #4. Updating the time average
        old_avg = flow_avg.copy()
        flow_avg = []
        for i in range(m):
            flow_avg.append(u_c.PWConst([0],[],0,True))

        breakpoints = []
        for i in range(m):
            for j in range(old_avg[i].noOfSegments + 1):
                if old_avg[i].segmentBorders[j] not in breakpoints:
                    breakpoints.append(old_avg[i].segmentBorders[j])

        for i in range(m):
            for j in range(flow[i].noOfSegments + 1):
                if flow[i].segmentBorders[j] not in breakpoints:
                    breakpoints.append(flow[i].segmentBorders[j])

        list_points = breakpoints.copy()
        list_points.sort()

        value_list = []
        for i in range(m):
            for j in range(len(list_points) - 1):
                value_1 = (1/counter_steps)*flow[i].getValueAt(list_points[j])
                value_2 = ((counter_steps-1)/counter_steps)*old_avg[i].getValueAt(list_points[j])
                value_list.append(value_1 + value_2)

        [flow_avg, queues_avg, pathop_avg] = elp.edge_loading_pln(m, cap, travel, list_points, value_list, T)

        breakpoint()
        #5. Check convergence
        diff_func = []
        diff = 0
        for i in range(m):
            neg_old = old_flow[i].smul(-1)
            diff_func.append(flow[i].__add__(neg_old))
            diff = diff + diff_func[i].norm()

        #to not have to worry with rounding errors and stuff like that
        if lamb != 0:
            if diff < lamb:
                break


        print("Difference to the last flow: " + str(diff))

        breakpoint()

        #6. Gap function: Gap = zero iff we have a dynamic equilibrium
        #Objective function
        def g(h):
            sum_1 = 0
            for i in range(m):
                for j in range(len(step_points) - 1):
                    val_1 = pathop[i].getValueAt(step_points[j+1])
                    val_2 = pathop[i].getValueAt(step_points[j])
                    sum_1 = sum_1 + (((val_1+val_2)/2) - 2*epsilon*(-h[(len(step_points)-1)*i+j] + flow[i].getValueAt(step_points[j])))*(h[(len(step_points)-1)*i+j] -flow[i].getValueAt(step_points[j]))
            return sum_1

        A = []
        for i in range(len(step_points) - 1):
            A.append([])
            for j in range(m):
                for k in range(len(step_points) - 1):
                    if k == i:
                        A[i].append(1)
                    else:
                        A[i].append(0)

        constraints = sci.optimize.LinearConstraint(A, lb=u, ub=u)
        bounds = []
        for i in range(m):
            for j in range(len(step_points) - 1):
                bounds.append((0,None))

        start = []
        for i in range(m):
            for j in range(len(step_points) - 1):
                start.append(flow[i].getValueAt(step_points[j]))
        h0 = start.copy()
        sol_gap = sci.optimize.minimize(g, h0, bounds=bounds, constraints=constraints)
        print("Value of gap function: " + str(sol_gap.fun))

        breakpoint()

    breakpoint()
    if counter_steps >= numsteps:
        print("The learning dynamics did not converge within the given number of steps and within the given precision")
    else:
        print("The learning dynamics did converge")


learning_pln(2, [2,1],[0.5,1],2,3,0.1,0.3,100,None)