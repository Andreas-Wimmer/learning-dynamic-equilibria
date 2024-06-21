import utilities_c as u_c
import math

#Here we implement a small routine that automatically does the edge-loading procedure for us in parallel link networks; only for piecewise constant functions; given are the number of edges, the capacities, the travel times, the step points and the constant values on the intervals; in our case steps is the same for every edge, only for the values we need a list of values for every edge
def edge_loading_pln(m, cap, travel, steps, values, T):
    E = []
    for i in range(m):
        E.append(["s","t"])

        flow = []
        for i in range(m):
            flow.append(u_c.PWConst([0],[],0,True))

        for i in range(m):
            for j in range(len(steps) - 1):
                flow[i].addSegment(steps[j+1], values[(len(steps)-1)*i+j])

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

        return [flow, queues, pathop]