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

horizon = 2
steps = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2]
numb_paths = 2
times = [[0.0, 0.8571428571428572, 1.0, 1.25, 1.75, 1.9375000000000002, 2.0, 2.0833333333333335, 2.6250000000000004],[0.0, 1.0, 1.25, 1.75, 2.0, 2.0833333333333335, 2.625]]
values = [[1.0, 1.8095238095238093, 1.8055555555555554, 1.625, 1.375, 1.6458333333333333, 1.625, 1.5416666666666665, 0.9999999999999996],[1.0, 1.6666666666666665, 1.625, 1.375, 1.5833333333333335, 1.5416666666666665, 1.0]]
f_slopes = [0,0]
l_slopes = [0.0]

averaging(horizon, times, values, f_slopes, l_slopes, steps, numb_paths)