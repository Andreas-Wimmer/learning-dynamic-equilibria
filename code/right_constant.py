from __future__ import annotations

from typing import List, Tuple

from machine_precision import eps
from arrays import elem_lrank, merge_sorted, merge_sorted_many
from piecewise_linear import PiecewiseLinear


class RightConstant:
    """
    This class defines right-continuous functions with finitely many break points (xᵢ, yᵢ).
    The breakpoints are encoded in the two lists times = [x₀, ..., xₙ] and values = [y₀,..., yₙ]
    where we assume that times is strictly increasing.
    This encodes the function f(x) = yᵢ with i maximal s.t. xᵢ <= x (or i = 0 if x₀ > x).
    """

    times: List[float]
    values: List[float]
    domain: Tuple[float, float] = (float("-inf"), float("inf"))

    def __json__(self):
        return {
            "times": self.times,
            "values": self.values,
            "domain": [
                "-Infinity" if self.domain[0] == float("-inf") else self.domain[0],
                "Infinity" if self.domain[1] == float("inf") else self.domain[1],
            ],
        }

    def __init__(
        self,
        times: List[float],
        values: List[float],
        domain: Tuple[float, float] = (float("-inf"), float("inf")),
    ):
        self.times = times
        self.values = values
        self.domain = domain
        assert len(self.values) == len(self.times)
        assert all(
            float("-inf") < self.values[i] < float("inf")
            for i in range(len(self.times))
        )
        assert all(
            self.domain[0] <= self.times[i] <= self.domain[1]
            for i in range(len(self.times))
        )

    def __call__(self, at: float) -> float:
        return self.eval(at)

    def eval_from_end(self, at: float) -> float:
        """
        Searches the lower rank of the element x in arr by going backwards from the last entry.
        The lower rank is the minimal number i in -1, ..., len(arr)-1,
        such that arr[i] <= x < arr[i+1] (with the interpretation arr[-1] = -inf and arr[len(arr)] = inf)
        """
        assert self.domain[0] <= at <= self.domain[1], f"Function not defined at {at}."
        rnk = len(self.times) - 1
        while rnk >= 0 and self.times[rnk] > at:
            rnk -= 1
        return self._eval_with_lrank(at, rnk)

    def eval(self, at: float) -> float:
        assert self.domain[0] <= at <= self.domain[1], f"Function not defined at {at}."
        rnk = elem_lrank(self.times, at)
        return self._eval_with_lrank(at, rnk)

    def _eval_with_lrank(self, at: float, rnk: int):
        assert self.domain[0] <= at <= self.domain[1], f"Function not defined at {at}."
        assert -1 <= rnk <= len(self.times)
        assert rnk == elem_lrank(self.times, at)

        if rnk == -1:
            return self.values[0]
        else:
            return self.values[rnk]

    def extend(self, start_time: float, value: float):
        assert start_time >= self.times[-1] - eps
        if abs(self.values[-1] - value) <= eps:
            return
        if abs(start_time - self.times[-1]) <= eps:
            #  Simply replace the last value
            self.values[-1] = value
        else:
            self.times.append(start_time)
            self.values.append(value)

    def equals(self, other):
        if not isinstance(other, RightConstant):
            return False
        return (
            self.values == other.values
            and self.times == other.times
            and self.domain == other.domain
        )

    def __radd__(self, other):
        if other == 0:
            return self
        if not isinstance(other, RightConstant):
            raise TypeError("Can only add a RightConstantFunction.")
        assert self.domain == other.domain

        new_times = merge_sorted(self.times, other.times)

        new_values = [0.0] * len(new_times)

        lptr = 0
        rptr = 0
        for i, time in enumerate(new_times):
            while lptr < len(self.times) - 1 and self.times[lptr + 1] <= time:
                lptr += 1
            while rptr < len(other.times) - 1 and other.times[rptr + 1] <= time:
                rptr += 1
            new_values[i] = self.values[lptr] + other.values[rptr]

        return RightConstant(new_times, new_values, self.domain)

    @staticmethod
    def sum(functions: List[RightConstant], domain=(0, float("inf"))) -> RightConstant:
        if len(functions) == 0:
            return RightConstant([0.0], [0.0], domain)
        new_times = merge_sorted_many([f.times for f in functions])
        new_values = [0.0] * len(new_times)
        ptrs = [0 for _ in functions]

        for i, time in enumerate(new_times):
            for j in range(len(ptrs)):
                while (
                    ptrs[j] < len(functions[j].times) - 1
                    and functions[j].times[ptrs[j] + 1] <= time
                ):
                    ptrs[j] += 1
            new_values[i] = sum(
                functions[j].values[ptrs[j]] for j in range(len(functions))
            )
        return RightConstant(new_times, new_values, domain)

    def __add__(self, other):
        return self.__radd__(other)

    def __neg__(self):
        return RightConstant(self.times, [-v for v in self.values], self.domain)

    def __sub__(self, other):
        if not isinstance(other, RightConstant):
            raise TypeError("Can only subtract a RightConstantFunction.")
        return self + (-other)

    def simplify(self) -> RightConstant:
        """
        This removes unnecessary timesteps
        """
        new_times = [self.times[0]]
        new_values = [self.values[0]]
        for i in range(0, len(self.times) - 1):
            # Add i+1, if it's necessary.
            if abs(self.values[i] - self.values[i + 1]) >= 1000 * eps:
                new_times.append(self.times[i + 1])
                new_values.append(self.values[i + 1])
        return RightConstant(new_times, new_values, self.domain)

    def integral(self) -> PiecewiseLinear:
        """
        Returns the integral starting from self.times[0] to x of self.
        """
        assert self.times[0] == self.domain[0] and self.domain[1] >= self.times[-1]
        times = self.times
        values = [0.0] * len(times)
        for i in range(len(times) - 1):
            values[i + 1] = values[i] + self.values[i] * (times[i + 1] - times[i])
        return PiecewiseLinear(
            times, values, self.values[0], self.values[-1], self.domain
        )
    
    #Computes the product of a piecewise linear function and a piecewise constant function,
    #but only on an interval, where both functions do not have a jump (to avoid discontinuities)
    def multiply(self, other: PiecewiseLinear, start: float, end: float) -> PiecewiseLinear:
        assert (self.domain[0] <= start and self.domain[1] >= end)
        assert (other.domain[0] <= start and other.domain[1] >= end) or (round(other.domain[0] - start, 12) >= 0 and round(other.domain[1] - end,12) >= 0)
        if start < self.times[-1]:
            self_index_s = 0
            while self.times[self_index_s] <= start:
                self_index_s = self_index_s + 1
            self_index_s = self_index_s - 1
            self_index_e = 0
            while self.times[self_index_e] < end:
                self_index_e = self_index_e + 1
            assert (self_index_s == self_index_e - 1)
        other_index_s = 0
        while other.times[other_index_s] <= start:
            other_index_s = other_index_s + 1
        other_index_s = other_index_s - 1
        other_index_e = 0
        while other.times[other_index_e] < end:
            other_index_e = other_index_e + 1
        
        assert (other_index_s == other_index_e - 1)
        times = [start, end]
        value_1 = self.eval(start)*other.eval(start)
        value_2 = self.eval(end)*other.eval(end)
        values = [value_1, value_2]
        first_slope = 0
        last_slope = 0
        product = PiecewiseLinear(times, values, first_slope, last_slope)
        return product

    #multiplies the given function with a arbitrary scalar    
    def mult_scalar(self, scalar: float) -> RightConstant:
        new_vals = []
        for i in range(len(self.values)):
            new_vals.append(self.values[i]*scalar)
        scalared_self = RightConstant(self.times, new_vals, self.domain)
        return scalared_self
    
    #returns the function in absolute value
    def abs_val(self) -> RightConstant:
        times = self.times
        domain = self.domain
        values = []
        for i in range(len(self.values)):
            if self.values[i] < 0:
                values.append(-self.values[i])
            else:
                values.append(self.values[i])
        abs_val_func = RightConstant(times, values, domain)
        return abs_val_func
