#Here we implement the learning dynamics of fictitious play in a dynamic network to see, if it
#converges and if it converges against a dynamic equilibrium.

import utilities_c as u_c
import math
import sympy as sy
import scipy as sci
import numpy as np
import pdb
import matplotlib.pyplot as plt

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

    if numsteps == None:
        numsteps = math.inf

    #for the initialization we just distribute the flow over the edges according to their
    #capacity for the whole interval
    sum_cap = sum(cap)
    values = []
    for i in range(m):
        values.append((cap[i]/sum_cap)*u)

    flow = []
    for i in range(m):
        flow.append(u_c.PWConst([0,T],[values[i]],0,True))


    queues = []
    for i in range(m):
        if values[i] >= cap[i]:
            queues.append(u_c.PWLin([0,T],[values[i] - cap[i]],[0],True))
        else:
            queues.append(u_c.PWLin([0,T],[0],[0],True))


    pathop = []
    for i in range(m):
        pathop.append(u_c.PWLin([0,T],[queues[i].segmentMvalues[0]/cap[i]],[travel[i]],True))

    breakpoint()
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

        breakpoint()
        #2. Best-response problem
        #Optimization with scipy is somehow too slow, so here we do it with sympy and based on the KKT-conditions
        #At first we need a variable for every time interval of every path-flow
        h_v = []
        for i in range(m):
            for j in range(len(step_points) - 1):
                ij = str(i) + str(j)
                h_ij = sy.symbols("h_" + str(i) + str(j))
                h_v.append(h_ij)


        lamb_0 = sy.symbols('lamb_0')
        lamb = []
        mu = []
        for i in range(m):
            for j in range(len(step_points) - 1):
                ij = str(i) + str(j)
                lamb_ij = sy.symbols('lamb_' + str(i) + str(j))
                lamb.append(lamb_ij)

        lamb_vec = sy.Matrix(lamb)

        for j in range(len(step_points) - 1):
            ij = str(j)
            mu_ij = sy.symbols('mu_' + str(j))
            mu.append(mu_ij)

        mu_vec = sy.Matrix(mu)


        equalities = []
        for i in range(m):
            for j in range(len(step_points) - 1):
                val_1 = pathop_avg[i].getValueAt(step_points[j])
                val_2 = pathop_avg[i].getValueAt(step_points[j+1])
                val_3 = 2*epsilon*(h_v[(len(step_points) - 1)*i + j] - flow_avg[i].getValueAt(step_points[j]))
                val_4 = lamb[(len(step_points) - 1)*i + j]
                val_5 = mu[j]
                val_6 = lamb_0*((val_1 + val_2)/2 + val_3) - val_4 + val_5
                equalities.append(sy.Eq(val_6,0))

        for j in range(len(step_points) - 1):
            sum_interval = 0
            for i in range(m):
                sum_interval = sum_interval + h_v[(len(step_points) - 1)*i + j]
            sum_interval = sum_interval - u
            #equalities.append(sy.Eq(sum_interval,0))

        complem = 0
        for i in range(m):
            for j in range(len(step_points) - 1):
                complem = complem - lamb[(len(step_points) - 1)*i + j]*h_v[(len(step_points) - 1)*i + j]
        #equalities.append(sy.Eq(complem,0))
        eq_tuple = tuple(equalities)

        variables = []
        variables.append(lamb_0)
        for i in range(m):
            for j in range(len(step_points) - 1):
                variables.append(h_v[(len(step_points) - 1)*i + j])

        for i in range(m):
            for j in range(len(step_points) - 1):
                variables.append(lamb[(len(step_points) - 1)*i + j])

        for j in range(len(step_points) - 1):
            variables.append(mu[j])

        var_tuple = tuple(variables)
        breakpoint()



        sol = sy.solve(eq_tuple,var_tuple)



        breakpoint()
        #3. Updating the flow
        old_flow = flow.copy()
        old_queues = queues.copy()
        old_pathop = pathop.copy()
        flow = []
        for i in range(m):
            flow.append(u_c.PWConst([0],[],0,True))

        for i in range(m):
            for j in range(len(step_points) - 1):
                flow[i].addSegment(step_points[j+1], sol.x[(len(step_points)-1)*i+j])

        queues = []
        for i in range(m):
            queues.append(u_c.PWLin([0,0],[0],[0],True))

        for i in range(m):
            for j in range(flow[i].noOfSegments):
                if flow[i].getValueAt(flow[i].segmentBorders[j]) >= cap[i]:
                    queues[i].addSegment(flow[i].segmentBorders[j+1], flow[i].getValueAt(flow[i].segmentBorders[j]) - cap[i])
                else:
                    if queues[i].getValueAt(flow[i].segmentBorders[j]) == 0:
                        queues[i].addSegment(flow[i].segmentBorders[j+1],0)
                    else:
                        empty = flow[i].segmentBorders[j] + ((queues[i].getValueAt(flow[i].segmentBorders[j]))/(cap[i] - flow[i].getValueAt(flow[i].segmentBorders[j])))
                        if T == min(T,flow[i].segmentBorders[j+1],empty):
                            queues[i].addSegment(T,cap[i] - flow[i].getValueAt(flow[i].segmentBorders[j]))
                        elif flow[i].segmentBorders[j+1] == min(T, flow[i].segmentBorders[j+1],empty):
                            queues[i].addSegment(flow[i].segmentBorders[j+1],cap[i] - flow[i].getValueAt(flow[i].segmentBorders[j]))
                        else:
                            queues[i].addSegment(empty, flow[i].getValueAt(flow[i].segmentBorders[j]) - cap[i])
                            queues[i].addSegment(flow[i].segmentBorders[j+1],0)

        pathop = []
        for i in range(m):
            pathop.append(u_c.PWLin([0,0],[0],[travel[i]],True))

        for i in range(m):
            for j in range(queues[i].get_noOfSegments()):
                pathop[i].addSegment(queues[i].segmentBorders[j+1],(queues[i].segmentMvalues[j]/cap[i]))

        counter_steps = counter_steps + 1
        breakpoint()
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
        for i in range(m):
            for j in range(len(list_points) - 1):
                value_1 = (1/counter_steps)*flow[i].getValueAt(list_points[j])
                value_2 = ((counter_steps-1)/counter_steps)*old_avg[i].getValueAt(list_points[j])
                flow_avg[i].addSegment(list_points[j+1], value_1 + value_2)

        queues_avg = []
        for i in range(m):
            queues_avg.append(u_c.PWLin([0,0],[0],[0],True))

        for i in range(m):
            for j in range(flow_avg[i].get_noOfSegments()):
                if flow_avg[i].getValueAt(flow_avg[i].segmentBorders[j]) >= cap[i]:
                    queues_avg[i].addSegment(flow_avg[i].segmentBorders[j+1], flow_avg[i].getValueAt(flow_avg[i].segmentBorders[j]) - cap[i])
                else:
                    if queues_avg[i].getValueAt(flow_avg[i].segmentBorders[j]) == 0:
                        queues_avg[i].addSegment(flow_avg[i].segmentBorders[j+1],0)
                    else:
                        empty = flow_avg[i].segmentBorders[j] + ((queues_avg[i].getValueAt(flow_avg[i].segmentBorders[j]))/(cap[i] - flow_avg[i].getValueAt(flow_avg[i].segmentBorders[j])))
                        if T == min(T,flow_avg[i].segmentBorders[j+1],empty):
                            queues_avg[i].addSegment(T,cap[i] - flow_avg[i].getValueAt(flow_avg[i].segmentBorders[j]))
                        elif flow_avg[i].segmentBorders[j+1] == min(T, flow_avg[i].segmentBorders[j+1],empty):
                            queues_avg[i].addSegment(flow_avg[i].segmentBorders[j+1],cap[i] - flow_avg[i].getValueAt(flow_avg[i].segmentBorders[j]))
                        else:
                            queues_avg[i].addSegment(empty,  flow_avg[i].getValueAt(flow_avg[i].segmentBorders[j])- cap[i])
                            queues_avg[i].addSegment(flow_avg[i].segmentBorders[j+1],0)

        pathop_avg = []
        for i in range(m):
            pathop_avg.append(u_c.PWLin([0,0],[0],[travel[i]],True))

        for i in range(m):
            for j in range(queues_avg[i].get_noOfSegments()):
                pathop_avg[i].addSegment(queues_avg[i].segmentBorders[j+1],(queues_avg[i].segmentMvalues[j]/cap[i]))

        breakpoint()
        #5. Check convergence
        diff_func = []
        diff = 0
        for i in range(m):
            neg_old = old_flow[i].smul(-1)
            diff_func.append(flow[i].__add__(neg_old))
            diff = diff + diff_func[i].norm()

        if diff < lamb:
            break

        print("Difference to the last flow: " + str(diff))

        breakpoint()

        #6. Gap function: Gap = zero iff we have a dynamic equilibrium, also calculated based on the KKT-conditions; but it is not exactly the gap-function, because we can not expect to get a true equilibrium, only an epsilon-regularized equilibrium
        #again we need for every path and every time interval a variable
        m_v = []
        for i in range(m):
            for j in range(len(step_points) - 1):
                ij = str(i) + str(j)
                m_ij = sy.symbols('m_' + str(i) + str(j))
                m_v.append(m_ij)

        m_vec = sy.Matrix(m_v)

        #Objective function
        obj_gap = 0
        for i in range(m):
            for j in range(len(step_points) - 1):
                val_1 = pathop[i].getValueAt(step_points[j])
                val_2 = pathop[i].getValueAt(step_points[j+1])
                val_3 = flow[i].getValueAt(step_points[j]) - m_ij
                val_4 = ((val_1 + val_2)/2 + 2*epsilon*val_3)*val_3
                obj_gap = obj_gap + val_4

        #Gradient of the objective function
        grad_obj_gap = sy.Matrix([obj_gap]).jacobian(sy.Matrix(list(obj_gap.free_symbols)))

        #Equality constraint
        result_eq_gap = []
        for j in range(len(step_points) - 1):
            sum_p = 0
            for i in range(m):
                ij = str(i) + str(j)
                sum_p = sum_p + m_ij
            sum_p = sum_p - u
            result_eq_gap.append(sum_p)
        eq_const_gap = sy.Matrix(result_eq_gap)


        #Jacobian of equality constraint
        grad_eq_gap = eq_const_gap.jacobian(sy.Matrix(list(obj_gap.free_symbols)))

        #Inequality constraint
        result_ineq_gap = []
        for i in range(m):
            for j in range(len(step_points) - 1):
                ij = str(i) + str(j)
                result_ineq_gap.append(-h_ij)
        ineq_const_gap = sy.Matrix(result_ineq_gap)

        #Jacobian of inequality constraint
        grad_ineq_gap = ineq_const_gap.jacobian(sy.Matrix(list(obj_gap.free_symbols)))

        lambd_0 = sy.symbols('lambd_0')
        for i in range(m):
            for j in range(len(step_points) - 1):
                ij = str(i) + str(j)
                lambd_ij = sy.symbols('lambd_' + str(i) + str(j))

        for j in range(len(step_points) - 1):
            ij = str(j)
            muu_ij = sy.symbols('muu_' + str(j))

        #Lagrangian function
        lambd = []
        for i in range(m):
            for j in range(len(step_points) - 1):
                ij = str(i) + str(j)
                lambd.append(lambd_ij)
        lambd_vec = sy.Matrix(lambd)

        muu = []
        for j in range(len(step_points) - 1):
            ij = str(j)
            muu.append(muu_ij)
        muu_vec = sy.Matrix(muu)

        Lag_gap = sy.Matrix([lambd_0 * obj_gap]) + lambd_vec.T * ineq_const_gap + muu_vec.T * eq_const_gap
        Grad_lag_gap = lambd_0 * grad_obj_gap + lambd_vec.T * grad_ineq_gap + muu_vec.T * grad_eq_gap

        zero_1 = sy.zeros(1, m*(len(step_points) - 1))
        zero_2 = sy.zeros(1, (len(step_points) - 1))
        zero_3 = sy.zeros(1,(m + 1)*(len(step_points) - 1))

        stat_gap = sy.Eq(Grad_lag_gap,zero_1)
        pfeas_eq_gap = sy.Eq(eq_const_gap,zero_2)
        comp_gap = sy.Eq(lambd_vec.T * ineq_const_gap, 0)

        sol_gap = sy.solveset((stat_gap,pfeas_eq_gap,ineq_const_gap <= zero_1,lambd_vec >= zero_1, muu_vec >= zero_2, comp_gap),(m_vec,lambd_vec,muu_vec), sy.Reals)


    #breakpoint()

    breakpoint()
    if counter_steps >= numsteps:
        print("The learning dynamics did not converge within the given number of steps and within the given precision")
    else:
        print("The learning dynamics did converge")


learning_pln(2, [1,1], [0.5,0.25], 1, 8, 0.1, 0.1, None, 0.000001)