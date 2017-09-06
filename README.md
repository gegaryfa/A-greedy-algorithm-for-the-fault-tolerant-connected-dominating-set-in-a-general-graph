# Implementation of: A greedy algorithm for the fault-tolerant connected dominating set in a general graph,(Jiao Zhou, Zhao Zhang,Weili Wu ,Kai Xing)

Using a connected dominating set (CDS) to serve as the virtual backbone of a wireless network
is an effective way to save energy and alleviate broadcasting storm. Since nodes may fail due to accidental damage
or energy depletion, it is desirable that the virtual backbone is fault tolerant. A node set C is an m-fold connected
dominating set (m-fold CDS) of graph G if every node in V(G)\C has at least m-neighbors in C and the subgraph of G
induced by C is connected. A greedy algorithm is presented to compute an m-fold CDS in a general graph, which has
size at most 2+ln(Delta + m - 2) times that of a minimum m-fold CDS, where Delta is the maximum degree of the graph.
This result improves on the previous best known performance ratio of 2*H(Delta + m - 1) for this problem, where H is
the harmonic number

Algorithm:
Input: A connected graph G = (V,E) and an integer m
Output: A(1,m)-CDS of G.
```
1. Set C <- 0
2. while there existsa node x in V\C such that ( -Delta(x)f(C) > 0  ) do
3.	Select x which maximizes -Delta(x)f(C).
4.	C <- C U {x}.
5. end while
6. Output: C(G) <- C
```