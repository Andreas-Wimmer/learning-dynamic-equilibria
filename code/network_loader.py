from __future__ import annotations

from typing import Dict, Generator, List, Set, Tuple

from dijkstra import dynamic_dijkstra, reverse_dijkstra
from dynamic_flow import DynamicFlow, FlowRatesCollection
from graph import Edge, Node
from machine_precision import eps
from network import Network
from queues import PriorityQueue
from right_constant import RightConstant
from piecewise_linear import PiecewiseLinear,identity

Path = List[Edge]


class NetworkLoader:
    network: Network
    path_inflows: List[Tuple[Path, RightConstant]]  # (path, inflow)
    _network_inflow_changes: PriorityQueue[
        Tuple[int, Node, float]
    ]  # (path_id, source, time)
    _net_inflow_by_node: Dict[Node, FlowRatesCollection]
    _path_indexes_by_sink: Dict[Node, Set[int]]

    def __init__(
        self, network: Network, path_inflows: List[Tuple[Path, RightConstant]]
    ):
        self.network = network
        self.path_inflows = path_inflows
        self._built = False
        self._flow = DynamicFlow(network)
        self._net_inflow_by_node = {
            v: FlowRatesCollection(
                {
                    i: path_inflow[1]
                    for i, path_inflow in enumerate(self.path_inflows)
                    if v == path_inflow[0][0].node_from
                }
            )
            for v in network.graph.nodes.values()
        }
        self._network_inflow_changes = PriorityQueue(
            [
                ((i, path[0].node_from, time), time)
                for i, (path, inflow) in enumerate(path_inflows)
                for time in inflow.times
            ]
        )

        self._path_indexes_by_sink = {
            v: set(
                i
                for i, path_inflow in enumerate(self.path_inflows)
                if path_inflow[0][-1].node_to == v
            )
            for v in network.graph.nodes.values()
        }
        self._handle_nodes = set()

    def build_flow(self) -> Generator[DynamicFlow,None,None]:
        if self._built:
            raise RuntimeError("Flow was already built. Initialize a new builder.")
        self._built = True

        #yield self._flow
        while self._flow.phi < float("inf"):
            while self._flow.phi >= self._network_inflow_changes.min_key():
                _, s, _ = self._network_inflow_changes.pop()
                self._handle_nodes.add(s)

            new_inflow = self._determine_new_inflow()
            max_ext_time = self._network_inflow_changes.min_key()
            edges_with_outflow_change = self._flow.extend(new_inflow, max_ext_time)
            self._handle_nodes = set(
                self.network.graph.edges[e].node_to for e in edges_with_outflow_change
            )

        yield self._flow
    
    
    def path_delay(self, T: float) -> List[PiecewiseLinear]:
        arr_funcs = self.expected_arr()
        minimal_delay = []
        for i in range(len(self.network.paths)):
            minimal_delay.append(0)
            for j in range(len(self.network.paths[i])):
                index = self.network.paths[i][j].id
                minimal_delay[i] = minimal_delay[i] + self.network.travel_time[index]
        path_delays = []
        for path in self.network.paths:
            delay_op = identity.restrict((0, float("inf")))
            for edge in path:
                index = edge.id
                delay_op = arr_funcs[index].compose(delay_op)
            
            for i in range(len(delay_op.times)):
                assert delay_op.values[i] >= minimal_delay[self.network.paths.index(path)]

            if delay_op.times[-1] >= T:
                delay_op = delay_op - identity.restrict((0,delay_op.times[-1]))
                delay_op.last_slope = 0.0
                path_delays.append(delay_op.restrict((0,delay_op.times[-1])))
            else:
                delay_op = delay_op - identity.restrict((0,T))
                delay_op.last_slope = 0.0
                path_delays.append(delay_op.restrict((0,T)).simplify())
                

        return path_delays

    def expected_arr(self) -> List[PiecewiseLinear]:
        queues = self._flow.get_queues()
        arr_funcs = []
        for i in range(len(self.network.graph.edges)):
            times = queues[i].times
            values = []
            for j in range(len(queues[i].values)):
                delay = queues[i].values[j]/self.network.capacity[i]
                curr_time = queues[i].times[j]
                values.append(self.network.travel_time[i] + curr_time + delay)
            first_slope = queues[i].first_slope/self.network.capacity[i] + 1
            last_slope = queues[i].last_slope/self.network.capacity[i] + 1
            new_func = PiecewiseLinear(times, values, first_slope, last_slope, (0, float("inf")))
            arr_funcs.append(new_func)

        return arr_funcs




    def _get_active_edges(self, i: int, s: Node) -> List[Edge]:
        path = self.path_inflows[i][0]
        edge = None
        for e in path:
            if e.node_from == s:
                edge = e
                break
        if edge is None:
            raise RuntimeError(
                "Node does not appear in path p (or is the last node of p)."
            )
        return [edge]

    def _determine_new_inflow(self) -> Dict[int, Dict[int, float]]:
        new_inflow = {}
        for v in self._handle_nodes:
            new_inflow.update({e.id: {} for e in v.outgoing_edges})

            outflows = {
                e.id: self._flow.outflow[e.id].get_values_at_time(self._flow.phi)
                for e in v.incoming_edges
            }

            net_inflow_by_com = self._net_inflow_by_node[v].get_values_at_time(
                self._flow.phi
            )

            used_commodities = (
                set(key for outflow in outflows.values() for key in outflow)
                .union(net_inflow_by_com.keys())
                .difference(self._path_indexes_by_sink[v])
            )

            for i in used_commodities:
                inflow = sum(
                    outflow[i] for outflow in outflows.values() if i in outflow
                )
                if i in net_inflow_by_com:
                    inflow += net_inflow_by_com[i]
                if inflow <= eps:
                    continue

                active_edges = self._get_active_edges(i, v)
                distribution = inflow / len(active_edges)
                for e in active_edges:
                    new_inflow[e.id][i] = distribution
        return new_inflow
