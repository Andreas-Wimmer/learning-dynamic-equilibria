from __future__ import annotations
from typing import List, Tuple
import pdb
import extendedRational

#import matplotlib.pyplot as plt

from extendedRational import *

pdb.set_trace();

# A right-constant function with finitely many steps
class PWConst:
    noOfSegments: int
    segmentBorders: List[float]
    segmentValues: List[float]
    defaultValue: float
    autoSimplify: bool

    def __init__(self, borders: List[float], values: List[float],
                 defaultValue: float = None, autoSimplyfy: bool = True):
        # autoSimplify=True means that adjacent segments with the same value are automatically unified
        # to a single segment
        # If a defaultValue is given this is the value of the function outside of the given borders
        # If none is given the function is undefined there
        self.defaultValue = defaultValue
        assert (len(borders) == len(values) + 1)
        self.noOfSegments = len(values)
        # TODO: check that borders are non-decreasing
        self.segmentBorders = borders
        self.segmentValues = values

        self.autoSimplify = autoSimplyfy

    def info(self):
        print("The function has the segment borders: " + str(self.segmentBorders) + " and the segment values " + str(self.segmentValues) +
              " and a default value of " + str(self.defaultValue))

    def addSegment(self, border: float, value: float):
        # Adds a new constant segment at the right side
        assert (self.segmentBorders[-1] <= border + 10*numPrecision)

        # Only add new segment if it is larger than the given precision
        if self.segmentBorders[-1] - numPrecision < border < self.segmentBorders[-1]+numPrecision:
            return

        # If autoSimplify is active and the new intervals value is the same (up to the given precision) as the one of
        # the last interval, we extend the last interval instead of creating a new one
        if self.autoSimplify and len(self.segmentValues) > 0 and value - numPrecision <= self.segmentValues[-1] <= value + numPrecision:
            self.segmentBorders[-1] = border
        else:
            self.segmentBorders.append(border)
            self.segmentValues.append(value)
            self.noOfSegments += 1

    #TODO: make changes to this function analogously to the changes of the variant for the PWLin functions
    def changeSegment(self, old_borders : [float, float], new_step : float, new_value : float):
        #Only change segment, if the changed segment is bigger than the given precision
        if new_borders[1] - new_borders < numPrecision:
            return

        if self.autoSimplify and len(self.segmentValues) > 0 and new_value - numPrecision <= self.segmentValues[-1] <= new_value + numPrecision:
            self.segmentBorders[-1] = new_step

        self.segmentBorders[-1] = old_borders[0]
        self.segmentBorders[0] = new_step
        self.segmentBorders.insert(1,old_borders[1])
        last_value = self.segmentValues[-1]
        old_value = self.segmentValues[0]
        self.segmentValues.insert(0,new_value)
        self.segmentValues[1] = old_value
        self.noOfSegments = self.noOfSegments + 1;

    def get_noOfSegments(self):
        return self.noOfSegments

    def get_borders(self):
        return self.segmentBorders

    def getValueAt(self, x: float) -> float:
        # Returns the value of the function at x
        if x < self.segmentBorders[0] or x >= self.segmentBorders[-1]:
            # x is outside the range of the function
            return self.defaultValue
        else:
            for i in range(0, self.noOfSegments):
                if x < self.segmentBorders[i + 1]:
                    return self.segmentValues[i]

    def getNextStepFrom(self, x: float) -> float:
        # get the next step of the function strictly after x
        if x >= self.segmentBorders[-1]:
            if self.defaultValue is None:
                # TODO: Raise an error
                pass
            else:
                return infinity
        else:
            for i in range(0, self.noOfSegments + 1):
                if x < self.segmentBorders[i]:
                    return self.segmentBorders[i]

    def __add__(self, other: PWConst) -> PWConst:
        # Add two piecewise constant functions

        # If at least one of the functions is undefined outside its borders the sum of the two function can only be
        # defined within the boundaries of that function (the intersection of the two boundaries if both functions
        # are undefined outside their boundaries)
        if self.defaultValue is None and other.defaultValue is None:
            default = None
            leftMost = max(self.segmentBorders[0], other.segmentBorders[0])
            rightMost = min(self.segmentBorders[-1], other.segmentBorders[-1])
        elif self.defaultValue is None and not (other.defaultValue is None):
            default = None
            leftMost = self.segmentBorders[0]
            rightMost = self.segmentBorders[-1]
        elif not(self.defaultValue is None) and other.defaultValue is None:
            default = None
            leftMost = other.segmentBorders[0]
            rightMost = other.segmentBorders[-1]
        else:
            default = self.defaultValue + other.defaultValue
            leftMost = min(self.segmentBorders[0], other.segmentBorders[0])
            rightMost = max(self.segmentBorders[-1], other.segmentBorders[-1])

        sum = PWConst([leftMost], [], default, self.autoSimplify and other.autoSimplify)

        x = leftMost
        while x < rightMost:
            val = self.getValueAt(x) + other.getValueAt(x)
            x = min(self.getNextStepFrom(x), other.getNextStepFrom(x))
            sum.addSegment(x, val)

        return sum

    def smul(self, mu: float) -> PWConst:
        # Creates a new piecewise constant function by scaling the current one by mu

        if self.defaultValue is None:
            default = None
        else:
            default = mu * self.defaultValue
        scaled = PWConst([self.segmentBorders[0]], [], default, self.autoSimplify)
        for i in range(len(self.segmentValues)):
            scaled.addSegment(self.segmentBorders[i + 1], mu * self.segmentValues[i])

        return scaled

    def restrictTo(self, a: float, b: float, default: float = None) -> PWConst:
        # Creates a new piecewise constant function by restricting the current one to the interval [a,b)
        # and setting it to default outside [a,b)

        x = max(a,self.segmentBorders[0])
        restrictedF = PWConst([x], [], defaultValue=default)

        while x <= self.segmentBorders[-1] and x < b:
            val = self.getValueAt(x)
            x = min(self.getNextStepFrom(x), b)
            restrictedF.addSegment(x, val)

        if x < b and (not self.getValueAt(x) is None):
            restrictedF.addSegment(b, self.getValueAt(x))

        return restrictedF

    def isZero(self) -> bool:
        # Checks whether the function is zero wherever it is defined
        if self.defaultValue is not None and -infinity < self.segmentBorders[0] \
                and self.segmentBorders[-1] < infinity and self.defaultValue != zero:
            # If the default value is not zero, the function is not zero
            return False
        for y in self.segmentValues:
            if y != zero:
                # If there is one segment where the function is non-zero, the function is not zero
                # (this assume that there are no zero-length intervals!)
                return False
        return True

    def __abs__(self) -> PWConst:
        # Creates a new piecewise constant functions |f|
        if self.defaultValue is None:
            default = None
        else:
            default = self.defaultValue.__abs__()
        absf = PWConst([self.segmentBorders[0]], [], default, self.autoSimplify)
        for i in range(len(self.segmentValues)):
            absf.addSegment(self.segmentBorders[i + 1], self.segmentValues[i].__abs__())

        return absf

    def integrate(self, a: float, b: float) -> float:
        # Determines the value of the integral of the given piecewise function from a to b
        assert (self.defaultValue is not None or (a >= self.segmentBorders[0] and b <= self.segmentBorders[-1]))

        integral = zero
        x = a
        while x < b:
            y = min(self.getNextStepFrom(x), b)
            integral += (y - x) * self.getValueAt(x)
            x = y

        return integral

    def norm(self) -> float:
        # Computes the L1-norm of the function
        # requires the function to either be undefined or zero outside its borders
        # (otherwise the L1-norm would be undefined/+-infty)
        assert (self.defaultValue is None or self.defaultValue == 0)
        return self.__abs__().integrate(self.segmentBorders[0], self.segmentBorders[-1])

    def drawGraph(self, start: float, end: float):
        # Draws a graph of the function between start and end
        current = start
        x = []
        y = []
        while self.getNextStepFrom(current) < end:
            x.append(current)
            x.append(self.getNextStepFrom(current))
            y.append(self.getValueAt(current))
            y.append(self.getValueAt(current))
            current = self.getNextStepFrom(current)
        x.append(current)
        x.append(end)
        y.append(self.getValueAt(current))
        y.append(self.getValueAt(current))
        plt.plot(x, y)
        return plt

    def getXandY(self, start: float, end: float) -> Tuple[List[float],List[float]]:
        # Returns two vectors x and y representing the function between start and end in the following form:
        # x = [a_0,a_1,a_1,a_2,a_2,...,a_n], y = [b_0,b_0,b_1,b_1,...,b_{n-1}]
        # such that [a_i,a_{i+1}) form a partition of [start,nextStep(end))
        # into maximal (if autoSimplify=True) intervals of constant value b_i of the function
        # i.e. for even i x[i] is the left boundary of such an interval and y[i] the value in it
        #      for odd i  x[i] is the right boundary of such an interval and y[i] the value in it
        current = start
        x = []
        y = []
        while self.getNextStepFrom(current) < end:
            x.append(current)
            x.append(self.getNextStepFrom(current))
            y.append(self.getValueAt(current))
            y.append(self.getValueAt(current))
            current = self.getNextStepFrom(current)
        x.append(current)
        x.append(end)
        y.append(self.getValueAt(current))
        y.append(self.getValueAt(current))
        return x,y

    def __str__(self):
        f = "|" + str(round(float(self.segmentBorders[0]),2)) + "|"
        for i in range(len(self.segmentValues)):
            # f += "-" + str(self.segmentValues[i]) + "-|" + str(self.segmentBorders[i + 1]) + "|"
            f += " " + str(round(float(self.segmentValues[i]),2)) + " |"
            if self.segmentBorders[i+1] < infinity:
                f += str(round(float(self.segmentBorders[i + 1]),2)) + "|"
            else:
                f += str(self.segmentBorders[i + 1]) + "|"
        return f

#A small function that creates a list of piecewise constant fucntions out of a list of step values and step times and default values
def create_functions_c(values: List[List[float]], breaks: List[List[float]], default: List[float]) -> List[PWConst]:
    assert len(values) == len(breaks)
    assert len(default) == len(values)
    for i in range(len(values)):
        assert len(values[i]) + 1 == len(breaks[i])

    functions = []
    for i in range(len(values)):
        functions.append(PWConst([0],[], default[i], True))

    for i in range(len(values)):
        for j in range(len(values[i])):
            functions[i].addSegment(breaks[i][j+1], values[i][j])

    return functions
            


# A piecewise linear function with finitely many break points
# I.e. each piece of the function is of the form f(x) = mx + t for all x in some interval [a,b) 
class PWLin:
    noOfSegments: int
    autoSimplify: bool
    segmentBorders: List[float]
    segmentTvalues: List[float]
    segmentMvalues: List[float]

    def __init__(self, borders: List[float], mvalues: List[float],
                 tvalues: List[float], autoSimplify: bool = True):
        # autoSimplify=True means that adjacent segments are automatically unified whenever possible
        self.autoSimplify = autoSimplify

        self.noOfSegments = len(mvalues)
        assert (len(tvalues) == len(mvalues))
        assert (len(borders) == self.noOfSegments + 1)

        # TODO: check that borders are non-decreasing
        self.segmentBorders = borders
        self.segmentMvalues = mvalues
        self.segmentTvalues = tvalues

    def info(self):
        print("The function has the segment borders: " + str(self.segmentBorders) + " and the segment m-values " + str(self.segmentMvalues) + 
              " and the segment T-values " + str(self.segmentTvalues))

    def get_noOfSegments(self):
        return self.noOfSegments

    def get_borders(self):
        return self.segmentBorders

    def addSegment(self, border: float, m: float, t: float = None):
        # Adds a new segment on the right side
        # If no t value is provided the function is extended continuously
        if t is None:
            assert (self.noOfSegments > 0)
            t = self.segmentTvalues[-1] + (self.segmentBorders[-1] - self.segmentBorders[-2]) * self.segmentMvalues[-1]

        if self.autoSimplify and self.noOfSegments > 0 and self.segmentMvalues[-1] == m and self.getValueAt(
                self.segmentBorders[-1]) == t:
            self.segmentBorders[-1] = border
        else:
            self.segmentBorders.append(border)
            self.segmentMvalues.append(m)
            self.segmentTvalues.append(t)
            self.noOfSegments += 1

    #borders contains the borders of the interval we want to change the function on; TODO: add the depleting to the function
    def changeSegment(self, borders : [float, float], new_mvalue : float, deplete: bool=False):
        position = self.segmentBorders.index(self.getNextStepFrom(borders[0]))
        old_mvalue = self.segmentMvalues[position - 1]
        old_borders = [self.segmentBorders[position-1],self.segmentBorders[position]]
        if position < len(self.segmentMvalues):
            next_mvalue = self.segmentMvalues[position]
            

        if borders[0] - old_borders[0] < numPrecision and old_borders[1] - borders[1] < numPrecision:
            self.segmentBorders[position - 1] = borders[0]
            self.segmentBorders[position] = borders[1]
            self.segmentMvalues[position - 1] = new_mvalue
        elif borders[0] - old_borders[0] < numPrecision:
            self.segmentBorders[position - 1] = borders[0]
            self.segmentBorders.insert(position, borders[1])
            self.segmentMvalues.insert(position - 1, new_mvalue)
            self.noOfSegments = self.noOfSegments + 1
        elif old_borders[1] - borders[1] < numPrecision:
            self.segmentBorders.insert(position, borders[0])
            self.segmentBorders[position + 1] = borders[1]
            self.segmentMvalues.insert(position, new_mvalue)
            self.noOfSegments = self.noOfSegments + 1
        else:
            self.segmentBorders.insert(position - 1, borders[0])
            self.segmentBorders.insert(position, borders[1])
            self.segmentMvalues.insert(position - 2,old_mvalue)
            self.segmentMvalues.insert(position - 1, new_mvalue)
            self.noOfSegments = self.noOfSegments + 1

        self.segmentTvalues.insert(position, 0)                                 
             
        for i in range(self.noOfSegments):
            self.segmentTvalues[i] = self.getValueAt(self.segmentBorders[i])
            


    def getValueAt(self, x: float) -> float:
        # Returns the value of the function at x
        if x < self.segmentBorders[0] or x > self.segmentBorders[-1]:
            # x is outside the range of the function
            pass # TODO: Raise error
        else:
            for i in range(0, self.noOfSegments):
                if x <= self.segmentBorders[i + 1]:
                    return self.segmentTvalues[i] + (x - self.segmentBorders[i]) * self.segmentMvalues[i]

    def getNextStepFrom(self, x: float) -> float:
        # Returns the next break point strictly after x
        if x >= self.segmentBorders[-1]:
            # TODO: Implement default value and/or better error handling
            pass
        elif self.noOfSegments == 1:
            if x < self.segmentBorders[1]:
                return self.segmentBorders[1]
        else:
            for i in range(0, self.noOfSegments + 1):
                if x < self.segmentBorders[i]:
                    return self.segmentBorders[i]


    def drawGraph(self, start: float, end: float):
        # Draws a graph of the function between start and end
        x = [start]
        y = [self.getValueAt(start)]
        while self.getNextStepFrom(x[-1]) < end:
            x.append(self.getNextStepFrom(x[-1]))
            y.append(self.getValueAt(x[-1]))
        x.append(end)
        y.append(self.getValueAt(end))
        plt.plot(x, y)
        return plt

    def segment_as_str(self,i:int,ommitStart:bool=True) -> str:
        # Creates a string of the form |2|3 4|4| for the i-th segment
        # |2| is omitted if ommitStart=True (standard)
        assert(i<self.noOfSegments)
        s = ""
        if not ommitStart:
            s += "|" + str(round(float(self.segmentBorders[i]),2)) + "|"
        s += str(round(float(self.segmentTvalues[i]),2)) + " "
        if self.segmentBorders[i+1] < infinity:
            s += str(round(float(self.segmentTvalues[i] + (self.segmentBorders[i + 1] - self.segmentBorders[i]) * self.segmentMvalues[i]),2))
        else:
            if self.segmentMvalues[i] == 0:
                s += "0"
            elif self.segmentMvalues[i] > 0:
                s += "infty"
            else:
                s += "-infty"
        s += "|" + str(round(float(self.segmentBorders[i+1]),2)) + "|"
        return s

    def __str__(self) -> str:
        f = "|" + str(round(float(self.segmentBorders[0]),2)) + "|"
        for i in range(len(self.segmentMvalues)):
            f += self.segment_as_str(i)

        return f

    #We only add piecewise linear functions that are defined on the same domain (could be changed later)
    def __add__(self, other: PWLin) -> PWLin:
        assert(self.get_borders()[0] == other.get_borders()[0] and
               self.get_borders()[self.get_noOfSegments()] == other.get_borders()[other.get_noOfSegments()])
        sum_func = PWLin([0,0], [0], [0], False)
        x = 0
        while x < self.get_borders()[self.get_noOfSegments()]:
            val = self.segmentMvalues[self.segmentBorders.index(self.getNextStepFrom(x)) - 1] + other.segmentMvalues[other.segmentBorders.index(other.getNextStepFrom(x)) - 1]
            x = min(self.getNextStepFrom(x), other.getNextStepFrom(x))
            sum_func.addSegment(x, val)

        return sum_func



    #We only substract two piecewise linear functions, if they are defined on the same domain (could be changed later)
    def __sub__(self, other: PWLin) -> PWLin:
        assert(self.get_borders()[0] == other.get_borders()[0] and
               self.get_borders()[self.get_noOfSegments()] == other.get_borders()[other.get_noOfSegments()])
        sub_func = PWLin([0,0], [0], [0], False)
        x = 0
        while x < self.get_borders()[self.get_noOfSegments()]:
            val = self.segmentMvalues[self.segmentBorders.index(self.getNextStepFrom(x)) - 1] - other.segmentMvalues[other.segmentBorders.index(other.getNextStepFrom(x)) - 1]
            x = min(self.getNextStepFrom(x), other.getNextStepFrom(x))
            sub_func.addSegment(x, val)

        return sub_func

    #multiplies a piecewise linear and a piecewise constant functions; result: piecewise linear; the product of a continous piecewise linear and a
    #piecewise constant function is not continous anymore
    def multiply(self, other: PWConst) -> PWLin:
        if other.defaultValue is None:
            assert self.segmentBorders[0] == other.segmentBorders[0] and self.segmentBorders[-1] == other.segmentBorders[-1]

        mul_func = PWLin([0,0], [0], [0], False)
        x = 0
        while x < self.segmentBorders[-1]:
            valM = self.segmentMvalues[self.segmentBorders.index(self.getNextStepFrom(x)) - 1] * other.getValueAt(x)
            valT = self.segmentTvalues[self.segmentBorders.index(self.getNextStepFrom(x)) - 1] * other.getValueAt(x)
            x = min(self.getNextStepFrom(x), other.getNextStepFrom(x))
            mul_func.addSegment(x, valM,valT)

        return mul_func

    #integrates a piecewise linear function in a given interval [a, b]
    def integrate(self, start: float, end: float) -> float:
        assert start <= end
        assert self.segmentBorders[0] <= start and self.segmentBorders[-1] >= end
        if start == self.segmentBorders[0] and end != self.segmentBorders[-1]:
            first_b = 1
            last_b = self.segmentBorders.index(self.getNextStepFrom(end - numPrecision))
        elif start != self.segmentBorders[0] and end == self.segmentBorders[-1]:
            first_b = self.segmentBorders.index(self.getNextStepFrom(start))
            last_b = -1
        elif start == self.segmentBorders[0] and end == self.segmentBorders[-1]:
            first_b = 1
            last_b = -1
        else:
            first_b = self.segmentBorders.index(self.getNextStepFrom(start))
            last_b = self.segmentBorders.index(self.getNextStepFrom(end - numPrecision))

        integral = 0
        diff = self.segmentBorders.index(self.segmentBorders[last_b]) - first_b
        first_q = (1/2)*self.segmentMvalues[first_b-1]*((self.segmentBorders[first_b])**2 - (start)**2)
        first_l = (self.segmentTvalues[first_b-1] - self.segmentMvalues[first_b-1]*self.segmentBorders[first_b-1])*(self.segmentBorders[first_b] - start)
        integral = integral + first_q + first_l
        last_q = (1/2)*self.segmentMvalues[last_b-1]*((end)**2 - (self.segmentBorders[last_b])**2)
        last_l = (self.segmentTvalues[last_b-1] - self.segmentMvalues[last_b-1]*self.segmentBorders[last_b-1])*(end - self.segmentBorders[last_b])
        integral = integral + last_q + last_l
        for i in range(diff):
            curr_q = (1/2)*self.segmentMvalues[first_b+i]*((self.segmentBorders[first_b+i+1])**2 - (self.segmentBorders[first_b+i])**2)
            curr_l = (self.segmentTvalues[first_b+i] - self.segmentMvalues[first_b+i]*self.segmentBorders[first_b+i])*(self.segmentBorders[first_b+i+1] - self.segmentBorders[first_b+i])
            integral = integral + curr_q + curr_l

        return integral

    #We compute the composition of two piecewise linear functions g * f (first f, then g)
    def compose(self, other: PWLin) -> PWLin:
        f = other
        g = self
        assert f.image()[0] <= g.segmentBorders[0] and f.image()[1] >= g.segmentBorders[-1]
        comp = PWLin([0,0],[0],[0],True)
        theta = 0
        next_f = 0
        next_g = 0
        breakpoint()
        while theta <= g.segmentBorders[-1]:
            if theta > f.image()[1]:
                break
            next_f = f.getNextStepFrom(theta)
            next_g = g.getNextStepFrom(f.getValueAt(theta))
            next_step = min(next_f,next_g)
            comp.addSegment(next_step, g.segmentMvalues[g.segmentBorders.index(next_g) - 1]*f.segmentMvalues[f.segmentBorders.index(next_f) - 1])
            theta = next_step
            breakpoint()

        return comp
        
        
    #Computes the value and the position of a minimum of a piecewise linear functions; if there exist several, one is returned
    def minimum(self) -> [float,float]:
        #A minimum is always gonna be at one of the border points
        minimum = math.inf
        pos_min = 0
        for i in range(self.noOfSegments + 1):
            if self.getValueAt(self.segmentBorders[i]) <= minimum:
                minimum = self.getValueAt(self.segmentBorders[i])
                pos_min = self.segmentBorders[i]

        return [minimum, pos_min]
        

    #Computes the value and the position of a maximum of a piecewise linear function; if there exist several, one is returned
    def maximum(self) -> [float,float]:
        #A maximum is always gonna be at one of the border points
        maximum = 0
        pos_max = 0
        for i in range(self.noOfSegments + 1):
            if self.getValueAt(self.segmentBorders[i]) >= maximum:
                maximum = self.getValueAt(self.segmentBorders[i])
                pos_max = self.segmentBorders[i]

        return [maximum, pos_max]

    #Returns the borders of the image of a piecewise linear function
    def image(self) -> [float,float]:
        return [self.minimum()[0], self.maximum()[0]]
        
        
    # def __str__(self):
        # f = "|" + str(self.segmentBorders[0]) + "|"
        # for i in range(len(self.segmentMvalues)):
            # f += str(self.segmentTvalues[i]) + "-" \
                 # + str(
                # self.segmentTvalues[i] + (self.segmentBorders[i + 1] - self.segmentBorders[i]) * self.segmentMvalues[i]) \
                 # + "|" + str(self.segmentBorders[i + 1]) + "|"

        # return f

#A small function that creates a list of piecewise linear fucntions out of a list of slope values and step times and a initial t-values: the functions
#are created in such a way that they are extended continously at every break point
def create_functions_l(Mvalues: List[List[float]], breaks: List[List[float]], Tvalues: List[float]) -> List[utilities_c.PWLin]:
    assert len(Mvalues) == len(breaks)
    assert len(Tvalues) == len(Mvalues)
    for i in range(len(Mvalues)):
        assert len(Mvalues[i]) + 1 == len(breaks[i])

    functions = []
    for i in range(len(Mvalues)):
        functions.append(PWLin([0,0],[0], [Tvalues[i]], True))

    for i in range(len(Mvalues)):
        for j in range(len(Mvalues[i])):
            functions[i].addSegment(breaks[i][j+1], Mvalues[i][j])

    return functions
