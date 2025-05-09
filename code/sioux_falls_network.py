#Here we implement the Sioux falls network

from graph import DirectedGraph,Node,Edge
from collections import deque
from network import Network,Commodity,Path

sioux_graph = DirectedGraph

#Set nodes
a = Node(0,sioux_graph)
b = Node(1,sioux_graph)
c = Node(2,sioux_graph)
d = Node(3,sioux_graph)
e = Node(4,sioux_graph)
f = Node(5,sioux_graph)
g = Node(6,sioux_graph)
h = Node(7,sioux_graph)
i = Node(8,sioux_graph)
j = Node(9,sioux_graph)
k = Node(10,sioux_graph)
l = Node(11,sioux_graph)
m = Node(12,sioux_graph)
n = Node(13,sioux_graph)
o = Node(14,sioux_graph)
p = Node(15,sioux_graph)
q = Node(16,sioux_graph)
r = Node(17,sioux_graph)
s = Node(18,sioux_graph)
t = Node(19,sioux_graph)
u = Node(20,sioux_graph)
v = Node(21,sioux_graph)
w = Node(22,sioux_graph)
x = Node(23,sioux_graph)

#Set edges 
e_1 = Edge(a,b,0,sioux_graph)
e_2 = Edge(a,c,1,sioux_graph)
e_3 = Edge(b,a,2,sioux_graph)
e_4 = Edge(b,f,3,sioux_graph)
e_5 = Edge(c,a,4,sioux_graph)
e_6 = Edge(c,d,5,sioux_graph)
e_7 = Edge(c,l,6,sioux_graph)
e_8 = Edge(d,c,7,sioux_graph)
e_9 = Edge(d,e,8,sioux_graph)
e_10 = Edge(d,k,9,sioux_graph)
e_11 = Edge(e,d,10,sioux_graph)
e_12 = Edge(e,f,11,sioux_graph)
e_13 = Edge(e,i,12,sioux_graph)
e_14 = Edge(f,b,13,sioux_graph)
e_15 = Edge(f,e,14,sioux_graph)
e_16 = Edge(f,h,15,sioux_graph)
e_17 = Edge(g,h,16,sioux_graph)
e_18 = Edge(g,r,17,sioux_graph)
e_19 = Edge(h,f,18,sioux_graph)
e_20 = Edge(h,g,19,sioux_graph)
e_21 = Edge(h,i,20,sioux_graph)
e_22 = Edge(h,p,21,sioux_graph)
e_23 = Edge(i,e,22,sioux_graph)
e_24 = Edge(i,h,23,sioux_graph)
e_25 = Edge(i,j,24,sioux_graph)
e_26 = Edge(j,i,25,sioux_graph)
e_27 = Edge(j,k,26,sioux_graph)
e_28 = Edge(j,o,27,sioux_graph)
e_29 = Edge(j,p,28,sioux_graph)
e_30 = Edge(j,q,29,sioux_graph)
e_31 = Edge(k,d,30,sioux_graph)
e_32 = Edge(k,j,31,sioux_graph)
e_33 = Edge(k,l,32,sioux_graph)
e_34 = Edge(k,n,33,sioux_graph)
e_35 = Edge(l,c,34,sioux_graph)
e_36 = Edge(l,k,35,sioux_graph)
e_37 = Edge(l,m,36,sioux_graph)
e_38 = Edge(m,l,37,sioux_graph)
e_39 = Edge(m,x,38,sioux_graph)
e_40 = Edge(n,k,39,sioux_graph)
e_41 = Edge(n,o,40,sioux_graph)
e_42 = Edge(n,w,41,sioux_graph)
e_43 = Edge(o,j,42,sioux_graph)
e_44 = Edge(o,n,43,sioux_graph)
e_45 = Edge(o,s,44,sioux_graph)
e_46 = Edge(o,v,45,sioux_graph)
e_47 = Edge(p,h,46,sioux_graph)
e_48 = Edge(p,j,47,sioux_graph)
e_49 = Edge(p,q,48,sioux_graph)
e_50 = Edge(p,r,49,sioux_graph)
e_51 = Edge(q,j,50,sioux_graph)
e_52 = Edge(q,p,51,sioux_graph)
e_53 = Edge(q,s,52,sioux_graph)
e_54 = Edge(r,g,53,sioux_graph)
e_55 = Edge(r,p,54,sioux_graph)
e_56 = Edge(r,t,55,sioux_graph)
e_57 = Edge(s,o,56,sioux_graph)
e_58 = Edge(s,q,57,sioux_graph)
e_59 = Edge(s,t,58,sioux_graph)
e_60 = Edge(t,r,59,sioux_graph)
e_61 = Edge(t,s,60,sioux_graph)
e_62 = Edge(t,u,61,sioux_graph)
e_63 = Edge(t,v,62,sioux_graph)
e_64 = Edge(u,t,63,sioux_graph)
e_65 = Edge(u,v,64,sioux_graph)
e_66 = Edge(u,x,65,sioux_graph)
e_67 = Edge(v,o,66,sioux_graph)
e_68 = Edge(v,t,67,sioux_graph)
e_69 = Edge(v,u,68,sioux_graph)
e_70 = Edge(v,w,69,sioux_graph)
e_71 = Edge(w,n,70,sioux_graph)
e_72 = Edge(w,v,71,sioux_graph)
e_73 = Edge(w,x,72,sioux_graph)
e_74 = Edge(x,m,73,sioux_graph)
e_75 = Edge(x,u,74,sioux_graph)
e_76 = Edge(w,x,75,sioux_graph)

#Set nodes and edges in the graph object
sioux_graph.nodes = {0:a,1:b,2:c,3:d,4:e,5:f,6:g,7:h,8:i,9:j,10:k,11:l,12:m,13:n,14:o,15:p,16:q,17:r,18:s,19:t,20:u,21:v,22:w,23:x}
sioux_graph.edges = [e_1,e_2,e_3,e_4,e_5,e_6,e_7,e_8,e_9,e_10,e_11,e_12,e_13,e_14,e_15,e_16,e_17,e_18,e_19,e_20,e_21,e_22,e_23,e_24,e_25,e_26,e_27,e_28,
                     e_29,e_30,e_31,e_32,e_33,e_34,e_35,e_36,e_37,e_38,e_39,e_40,e_41,e_42,e_43,e_44,e_45,e_46,e_47,e_48,e_49,e_50,e_51,e_52,e_53,e_54,e_55,e_56,
                     e_57,e_58,e_59,e_60,e_61,e_62,e_63,e_64,e_65,e_66,e_67,e_68,e_69,e_70,e_71,e_72,e_73,e_74,e_75,e_76]

edges = sioux_graph.edges.copy()

#Set capacities
capacities = [7.19450017,6.500964775,7.19450017,1.37727248,6.50096775,4.75292325,4.939665027,1.3635629805,4.939665027,1.37444318583,
              2.77,1.37727248,1.37444318583,1.3607187905,2.178209194,6.500964775,1.3607187905,2.1782809194,1.402831432,1.40161738417,
              2.77,1.402831432,3.865496783,3.865496783,2.77,3.75333376389,1.3485825472,1.38708630389,1.3635629805,2.77,1.3635629805,1.35458563527,
              6.500964775,1.3635629805,1.35458563527778,6.500964775,1.3635629805,7.194500177,7.194500177,1.41423782,1.354585635277,1.42431281083,
              1.367997390277,3.75333376388,1.42431281083,4.04576476388,2.66643904583,1.401617384166,1.348588254722,1.452752795277,5.466637975,
              1.38708630388,1.452752795277,1.339986341944,6.500964775,5.466637975,6.500964775,4.04576476388,1.33998634194,1.389613211944,
              6.500964775,1.38961321194,1.4055312055,1.40991588694,1.4055312055,1.452752795277,1.3570437677,2.666439045833,1.40991588694,
              1.452752795277,1.388,1.367997390277,1.3888,1.4106967877,1.41423782,1.3570437677,1.4106967877]

#Set travel times
travel_times = [360,240,360,300,240,240,240,240,120,360,120,240,300,300,240,120,180,120,120,180,600,300,300,600,
                180,180,300,360,240,480,360,300,360,240,240,360,180,180,240,240,300,240,360,300,180,180,300,240,120,180,480,120,
                120,120,180,240,180,120,240,240,240,360,300,360,120,180,180,300,120,240,240,240,120,240,180,120]

#Initiate network
sioux_graph.reversed = False
sioux_network = Network()
sioux_network.graph = sioux_graph
sioux_network.capacity = capacities
sioux_network.travel_time = travel_times

#Compute all paths from a to t
paths = sioux_network.findPaths(a, t)

#Optional (if total demand is also low enough): exclude the more expensive paths to make computations easier for RFP
#for i in range(0):
#    average_minium_delay = 0
#    minimal_delay = []
#    for i in range (len(paths)):
#        minimal_delay.append(0)
#        for j in range(len(paths[i].edges)):
#                index = paths[i].edges[j].id - 1
#                minimal_delay[i] = minimal_delay[i] + travel_times[index]
#        average_minium_delay = average_minium_delay + minimal_delay[i]
#
#    average_minium_delay = average_minium_delay/len(paths)
#    print(average_minium_delay)
#
#    new_paths = []
#    for i in range(len(paths)):
#           if minimal_delay[i] <= average_minium_delay:
#            new_paths.append(paths[i])
#    print(len(new_paths))
#    paths = new_paths.copy()

