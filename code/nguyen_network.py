#Here we implement the ngyuen network

from graph import DirectedGraph, Edge, Node
from network import Network, Commodity, Path
from right_constant import RightConstant

nguyen_graph = DirectedGraph()
#Set nodes
a = Node(0, nguyen_graph)
b = Node(1, nguyen_graph)
c = Node(2, nguyen_graph)
d = Node(3, nguyen_graph)
e = Node(4, nguyen_graph)
f = Node(5, nguyen_graph)
g = Node(6, nguyen_graph)
h = Node(7, nguyen_graph)
i = Node(8, nguyen_graph)
j = Node(9, nguyen_graph)
k = Node(10, nguyen_graph)
l = Node(11, nguyen_graph)
m = Node(12, nguyen_graph)
#Set edges 
e_1 = Edge(a,l,0,nguyen_graph)
e_2 = Edge(a,e,1,nguyen_graph)
e_3 = Edge(l,f,2,nguyen_graph)
e_4 = Edge(l,h,3,nguyen_graph)
e_5 = Edge(d,e,4,nguyen_graph)
e_6 = Edge(e,f,5,nguyen_graph)
e_7 = Edge(f,g,6,nguyen_graph)
e_8 = Edge(g,h,7,nguyen_graph)
e_9 = Edge(d,i,8,nguyen_graph)
e_10 = Edge(e,i,9,nguyen_graph)
e_11 = Edge(f,j,10,nguyen_graph)
e_12 = Edge(g,k,11,nguyen_graph)
e_13 = Edge(h,b,12,nguyen_graph)
e_14 = Edge(i,j,13,nguyen_graph)
e_15 = Edge(j,k,14,nguyen_graph)
e_16 = Edge(k,b,15,nguyen_graph)
e_17 = Edge(i,m,16,nguyen_graph)
e_18 = Edge(k,c,17,nguyen_graph)
e_19 = Edge(m,c,18,nguyen_graph)
#Set nodes and edges in graoh object
nguyen_graph.nodes = {0:a,1:b,2:c,3:d,4:e,5:f,6:g,7:h,8:i,9:j,10:k,11:l,12:m}
nguyen_graph.edges = [e_1,e_2,e_3,e_4,e_5,e_6,e_7,e_8,e_9,e_10,e_11,e_12,e_13,e_14,e_15,e_16,e_17,e_18,e_19]
#Set capacities
capacities = [5/6,5/6,5/6,5/6,5/6,5/6,5/6,5/6,5/6,5/6,5/6,5/6,5/6,5/6,5/6,5/6,5/6,5/6,5/6]
#Set travel times
travel_times = [150,75,75,150,150,150,150,150,225,75,75,75,75,150,150,150,225,75,150]

#Source node is a and sink node is b
p_1 = Path([e_1,e_4,e_13])
p_2 = Path([e_1,e_3,e_7,e_8,e_13])
p_3 = Path([e_1,e_3,e_7,e_12,e_16])
p_4 = Path([e_1,e_3,e_11,e_15,e_16])
p_5 = Path([e_2,e_6,e_7,e_8,e_13])
p_6 = Path([e_2,e_6,e_7,e_12,e_16])
p_7 = Path([e_2,e_6,e_11,e_15,e_16])
p_8 = Path([e_2,e_10,e_14,e_15,e_16])
paths = [p_1,p_2,p_3,p_4,p_5,p_6,p_7,p_8]

#Network inflow rate
u = RightConstant([0,300],[4,0],(0,300))
