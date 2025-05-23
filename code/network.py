from __future__ import annotations

import pickle
from typing import Dict, List

import numpy as np

from graph import DirectedGraph, Edge, Node
from predictor_type import PredictorType
from right_constant import RightConstant
import math
from collections import deque

class Path:
    edges: List[Edge]

    def __init__(self, e_list: List[Edge]):
        self.edges = e_list

    def isNodeInPath(self, x: Node) -> bool:
        # print("nodes in path ", str(self.getNodesInPath()))
        if x in self.getNodesInPath():
            # print("False")
            return True
        else:
            # print("True")
            return False
        
    def getNodesInPath(self) -> List[Node]:
        nodeList = [self.edges[0]._node_from]
        for e in self.edges:
            nodeList.append(e._node_to)
        return nodeList
    
class Commodity:
    sources: Dict[Node, RightConstant]
    sink: Node
    predictor_type: PredictorType

    def __init__(
        self,
        sources: Dict[Node, RightConstant],
        sink: Node,
        predictor_type: PredictorType,
    ):
        self.sources = sources
        self.sink = sink
        self.predictor_type = predictor_type


class Network:
    graph: DirectedGraph
    paths: List[Path]
    capacity: np.ndarray[float]
    travel_time: np.ndarray[float]
    commodities: List[Commodity]

    def __init__(self):
        self.graph = DirectedGraph()
        self.capacity = np.array([])
        self.travel_time = np.array([])
        self.commodities = []
        self.paths = []

    def __getstate__(self):
        return {
            "graph": self.graph,
            "capacity": self.capacity,
            "travel_time": self.travel_time,
            "commodities": [
                {
                    "sources": {s.id: value for s, value in c.sources.items()},
                    "sink": c.sink.id,
                    "predictor_type": c.predictor_type,
                }
                for c in self.commodities
            ],
        }

    def __setstate__(self, state):
        self.graph = state["graph"]
        self.capacity = state["capacity"]
        self.travel_time = state["travel_time"]
        self.commodities = []
        for c in state["commodities"]:
            self.add_commodity(c["sources"], c["sink"], c["predictor_type"])

    def add_edge(
        self, node_from: int, node_to: int, travel_time: float, capacity: float
    ):
        self.graph.add_edge(node_from, node_to)
        self.travel_time = np.append(self.travel_time, travel_time)
        self.capacity = np.append(self.capacity, capacity)

    def add_commodity(
        self,
        sources: Dict[int, RightConstant],
        sink: int,
        predictor_type: PredictorType,
    ):
        nodes = self.graph.nodes
        assert all(
            s in nodes.keys() for s in sources
        ), f"No node with id#{sink} in the graph!"
        assert sink in nodes.keys(), f"No node with id#{sink} in the graph!"
        self.commodities.append(
            Commodity(
                {nodes[s]: v for s, v in sources.items()}, nodes[sink], predictor_type
            )
        )

    def add_path(self, edges: List[Edge]):
        new_path = edges
        self.paths.append(new_path)

    def _remove_edge(self, edge: Edge):
        edge.node_to.incoming_edges.remove(edge)
        edge.node_from.outgoing_edges.remove(edge)
        del self.graph.edges[edge.id]
        self.capacity = np.delete(self.capacity, edge.id)
        self.travel_time = np.delete(self.travel_time, edge.id)
        for i in range(len(self.graph.edges)):
            self.graph.edges[i].id = i

    def compress_lonely_nodes(self):
        """
        A node is lonely, if it is no source or sink and if it has a single incoming and a single outgoing edge.
        This function removes these useless nodes to speed up computation
        """
        remove_nodes = []
        for v in self.graph.nodes.values():
            if len(v.outgoing_edges) == 1 == len(v.incoming_edges) and all(
                c.source != v != c.sink for c in self.commodities
            ):
                edge1 = v.incoming_edges[0]
                edge2 = v.outgoing_edges[0]
                new_travel_time = (
                    self.travel_time[edge1.id] + self.travel_time[edge2.id]
                )
                new_capacity = min(self.capacity[edge1.id], self.capacity[edge2.id])
                self._remove_edge(edge1)
                self._remove_edge(edge2)
                if edge1.node_from != edge2.node_to:
                    self.add_edge(
                        edge1.node_from.id,
                        edge2.node_to.id,
                        new_travel_time,
                        new_capacity,
                    )
                remove_nodes.append(v)
        for v in remove_nodes:
            self.graph.nodes.pop(v.id)

    def remove_unnecessary_nodes(self):
        """
        This functions checks whether a node is necessary for any commodity.
        A node v is necessary for a commodity, if there is a path from its source to its sink passing through v.
        """
        important_nodes = set()
        for commodity in self.commodities:
            reaching_t = self.graph.get_nodes_reaching(commodity.sink)
            reachable_from_s = self.graph.get_reachable_nodes(commodity.source)
            important_nodes = important_nodes.union(
                reaching_t.intersection(reachable_from_s)
            )

        remove_nodes = set(self.graph.nodes.values()).difference(important_nodes)

        print(f"\rRemoving {len(remove_nodes)} unnecessary nodes.", end="\r")
        for v in remove_nodes:
            to_remove = v.outgoing_edges.copy()
            for edge in to_remove:
                self._remove_edge(edge)
            to_remove = v.incoming_edges.copy()
            for edge in to_remove:
                self._remove_edge(edge)
            self.graph.nodes.pop(v.id)
        print(f"\rRemoved {len(remove_nodes)} unnecessary nodes.")

    def remove_unnecessary_commodities(self, selected_commodity: int) -> int:
        """
        Remove all commodities that cannot interfere with commodity i.
        """
        selected_comm = self.commodities[selected_commodity]
        important_nodes = [
            self.graph.get_reachable_nodes(i.source).intersection(
                self.graph.get_nodes_reaching(i.sink)
            )
            for i in self.commodities
        ]
        remove_commodities = []
        for i, comm in enumerate(self.commodities):
            if (
                len(
                    important_nodes[i].intersection(important_nodes[selected_commodity])
                )
                == 0
            ):
                remove_commodities.append(comm)

        print(
            f"\rRemoving {len(remove_commodities)} non-interfering commodities.",
            end="\r",
        )
        for comm in remove_commodities:
            self.commodities.remove(comm)
        print(f"Removed {len(remove_commodities)} non-interfering commodities.")

        return self.commodities.index(selected_comm)

    def print_info(self):
        print(
            f"The network contains {len(self.graph.nodes)} nodes and {len(self.graph.edges)} edges."
        )
        print(f"Moreover, there are {len(self.commodities)} commodities.")
        print(
            f"Minimum/Average/Maximum capacity: {np.min(self.capacity)}/{np.average(self.capacity)}/{np.max(self.capacity)}"
        )
        print(
            f"Minimum/Average/Maximum transit time: {np.min(self.travel_time)}/{np.average(self.travel_time)}/{np.max(self.travel_time)}"
        )
        max_in_degree = 0
        max_out_degree = 0
        max_degree = 0
        for node in self.graph.nodes.values():
            max_degree = max(
                max_degree, len(node.incoming_edges) + len(node.outgoing_edges)
            )
            max_in_degree = max(max_in_degree, len(node.incoming_edges))
            max_out_degree = max(max_out_degree, len(node.outgoing_edges))
        print(f"Maximum indgree: {max_in_degree}")
        print(f"Maximum outdegree: {max_out_degree}")
        print(f"Maximum degree: {max_degree}")
        avg_demand = np.average(
            [
                inflow.values[0]
                for c in self.commodities
                for inflow in c.sources.values()
            ]
        )
        print(f"Average demand: {avg_demand}")

    def printPathInNetwork(self, p: Path) -> str:
        s = str()
        for e in p:
            if len([i for i in e._node_from.outgoing_edges\
                if i._node_to == e._node_to]) > 1:
                    s += str(self.graph.edges.index(e))
            s += str(e)
        return s

    def findPaths(self, src, dest,excludeSelfLoopNodes: bool=False, verbose: bool=False) -> List[Path]:
        if excludeSelfLoopNodes:
            # Find nodes with self loops
            selfLoopNodes = [e._node_from for e in self.graph.edges if e._node_from == e._node_to]
        else:
            selfLoopNodes = []
        # print('Nodes with self loop: ', *(n for n in selfLoopNodes))

        # Queue to store (partial) paths
        q = deque()

        # Add edges going out of the source to q
        for e in self.graph.edges:
            if e._node_from == src:
                q.append(Path([e]))
        # List to store the final paths
        pathList = []
        count = 0

        while q:
            count += 1
            if verbose: print("\ncount:%d"%count)

            # Get the (earliest generated) partial path
            if verbose: print("q before pop")
            for p in q:
                if verbose: print(printPathInNetwork(self, p))
            path = q.popleft()
            # print("after pop ", printPathInNetwork(path, self), q)
            for p in q:
                if verbose: printPathInNetwork(self, p)

            # Get the last node in the partial path
            last = path.edges[-1]._node_to
            if verbose: print("last ", last)

            # If the last node is the destination node then store the path
            if last == dest:
                if verbose: print("Found s-t Path:", printPathInNetwork(self,path))
                pathList.append(path)

            # Traverse all the nodes connected to the current node and push new partial
            # path to queue
            edgeListCurrNode = [e for e in self.graph.edges if (e._node_from == last and
                e._node_to not in selfLoopNodes)]
            if verbose: print("edgeListCurrNode ", len(edgeListCurrNode), edgeListCurrNode)
            for e in edgeListCurrNode:
                if verbose: print("edge %d" %self.graph.edges.index(e), e,
                        printPathInNetwork(self,path), path.isNodeInPath(e._node_to))
                if not path.isNodeInPath(e._node_to):
                    newpath = Path(path.edges.copy())
                    if verbose: print("newpath before append ", printPathInNetwork(self,newpath))
                    newpath.edges.append(e)
                    q.append(newpath)
                    if verbose: print("newpath after append ", printPathInNetwork(self,newpath))

        # Print pathList
        if verbose:
            print("\nTotal %d paths found from node %s to node %s:"%(len(pathList),src,dest))
            for i,p in enumerate(pathList):
                print(i, len(p.edges), printPathInNetwork(self,p))
        return pathList
    
    def to_file(self, file_path: str):
        with open(file_path, "wb") as file:
            pickle.dump(self, file)

    @staticmethod
    def from_file(file_path: str) -> Network:
        with open(file_path, "rb") as file:
            return pickle.load(file)
        
def printPathInNetwork(G: Network, p: Path) -> str:
        s = str()
        for e in p.edges:
            if len([i for i in e._node_from._outgoing_edges\
                if i._node_to == e._node_to]) > 1:
                    s += str(G.graph.edges.index(e))
            s += str(e)
        return s
