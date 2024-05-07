Thanks for your attention. This is the source code of our TKDE 2023 “Size-constrained Community Search on Large Graphs: An Effective and Efficient Solution”.

#### Experimental Environment            
The code can be run on a Linux environment. Our code is compiled on a CentOS Linux (release 7.9.2009) server with GCC 7.5.0.

#### Summary of Program         
TD: run Truss Decomposition algorithm and output the trussness of all edges in the graph.   
    usage: TD {input_orig} {output_truss}        
    
STExa: run ST-Exa algorithm and output the subgraph that satisfies the definition of STCS (Size-constrained Truss Community Search).     
    usage: STExa {input_truss} {query} {lower bound of size} {upper bound of size} {time limit}
    
SCBRB: run SC-BRB algorithm and output the subgraph taht satisfies the definition of SCS (Size-Bounded Community Search).
    usage: STExa {input_orig} {query} {lower bound of size} {upper bound of size} {time limit}

khash.h: a hash table from https://github.com/attractivechaos/klib

SCBRB is the state-of-the-art exact algorithm for SCS. And STExa algorithm proposed by us is the exact algorithm for STCS.

#### Parameters      
{input_orig}: input graph, the first line includes the total number of vertices and edges in the graph. After the first line each line is the endpoints of an edge       
{output_truss} output graph, each line is the endpoints of an edge with trussness.           
{input_truss} input graph, it is the same as "output_truss".                    
{query} the query vertex.     
{lower bound of size} the lower bound of size.     
{upper bound of size} the upper bound of size.  
{time limit} the limit of running time. If running time is over the time limit, terminate the program and output the subgraph with the maxinum min-support/min-degree among each computed subgraph satisfying size constraint and including query vertex.

#### Running Instance   
1. Running makefile to produce "STExa" and "TD", "make";
2. Running "TD" to produce "truss-file", "./TD ./dblp.txt ./tru-dblp.txt";
3. Running "STExa" to return resulting subgraph satisfying the limit of the STCS problem, "./STExa ./tru-dblp.txt 0 11 20 1000".

#### Note
If you need the code of SCBRB, please contact the authors of SCBRB [47]. And I attach the "makefile" of SCBRB to easily compile and run it after obtaining the code.

SCBRB* is SCBRB optimized with our pruning approaches, by modifying the code of "SCBRB.cpp" as follows. 
1. Function "CSSC_heu" in "SCBRB.cpp", is replaced by our heuristic algorithm (Algorithm 2, Functions "Findku" and "Heuristic" in "STExa.cpp"). The computation of Algorithm 2 is based on degree rather than support, e.g., for the computations of "k*" and "score function". 
2. Add a check function in "BB_dom_ustar" of "SCBRB.cpp". The function terminates current branch if there is a vertex with its coreness in G no larger than "kl" in current branch, where "kl" in "SCBRB.cpp" denotes the maximum min-degree among each computed subgraph satisfying size constraint and including query vertex.
