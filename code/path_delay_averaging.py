#Here I extra implement the averaging process for the path delays operators to see what goes wrong
#in the regularized fictitious play routine

from typing import List 
from piecewise_linear import PiecewiseLinear

def averaging(horizon: float, times: List[List[float]], values: List[List[float]], f_slopes: List[float], l_slopes: List[float], 
              steps: List[float], numb_paths: float):
    
    delays_avg = []
    for i in range(numb_paths):
        delays_avg.append(PiecewiseLinear(times[i],values[i],f_slopes[i],l_slopes[i],(0,horizon)))
    
    avg = []
    for i in range(numb_paths):
        avg.append([])
        for j in range(len(steps) - 1):
            sum_delays = 0
            count_delays = 0
            for k in range(len(delays_avg[i].times) - 1):
                length = 0
                if delays_avg[i].times[k] >= steps[j] and delays_avg[i].times[k] <= steps[j+1]:
                    if delays_avg[i].times[k + 1] > steps[j + 1]:
                        end = steps[j + 1]
                    else:
                        end = delays_avg[i].times[k + 1]
                    start = delays_avg[i].times[k]
                    length = end - start 
                    average = ((delays_avg[i].eval(end) + delays_avg[i].eval(start))/2)
                    sum_delays = sum_delays + (length/(steps[j + 1] - steps[j]))*average
                    count_delays = count_delays + 1
            if count_delays == 0:
                sum_delays = ((delays_avg[i].eval(steps[j]) + delays_avg[i].eval(steps[j + 1]))/2)
                count_delays = 1
            avg[i].append((sum_delays/count_delays))

    return avg

