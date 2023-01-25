from turtle import back
import numpy as np
#import BLOSUM62
from Bio.Align import substitution_matrices
sub_mat = substitution_matrices.load("BLOSUM62")


SIGMA = 5

def read(file_name):
    
    with open(file_name, 'r+') as file:
        #k = int(file.readline().strip())
        all = [line.strip() for line in file.readlines()]
        return all[0], all[1]

def middle_edge(v,w):
    """ 
    Taking in sequences and returns the middle edge and node

    """
    m = len(w)
    middle = m//2
    n = len(v)
    linear_source = np.zeros(n+1, dtype=int) #creating 1D array to hold 
    linear_source[0] = 0 
    linear_sink = np.zeros(n+1, dtype=int) #creating 1D array to hold 
    linear_sink[0] = 0
    backtrack = np.zeros(n+1, dtype=int)
    sink_w = w[::-1]  #reversing strings
    sink_v = v[::-1] 
    midd = len(sink_w)//2

    for i in range(1, n+1):
        linear_source[i] = linear_source[i-1] - SIGMA
    
    for i, letter in enumerate(w):
        temp = linear_source[0]
        linear_source[0] += -SIGMA
        for k in range(1, n+1):
            diag = temp + sub_mat[letter, v[k-1]]   # 0
            vert = linear_source[k-1]  - SIGMA             # gap vertical # 1
            horz = linear_source[k]  - SIGMA               # gap horizontal # 2
            incoming = [diag, vert, horz]
            max_edge = np.argmax(incoming) #max score edge
            temp = linear_source[k]

            linear_source[k] = incoming[max_edge]
        
        if i + 1== middle:
            break
    
    for i in range(1, n+1):
        linear_sink[i] = linear_sink[i-1] - SIGMA

    #going from sink to middle is that same as reversing strings and going to middle
    for i, letter in enumerate(sink_w):
        temp2 = linear_sink[0]
        linear_sink[0] += -SIGMA
        for k in range(1, n+1):
            diag = temp2 + sub_mat[letter, sink_v[k-1]]     # 0
            vert = linear_sink[k-1]  - SIGMA                # 1 gap vertical
            horz = linear_sink[k]  - SIGMA                  # 2 gap horizontal
            incoming = [diag, vert, horz]
            max_edge = np.argmax(incoming) #max score edge
            backtrack[k] = max_edge 

            temp2 = linear_sink[k]

            linear_sink[k] = incoming[max_edge]
        
        if (i + 1) == (len(w) - middle):
            break
    
    #summing up the array
    flipped = linear_sink[::-1]
    ahh = linear_source + flipped
    #print(linear_source, flipped, ahh)
    #print(ahh)
    mid_node_i = np.argmax(ahh)

    from_node = (mid_node_i, middle)
    edge = backtrack[-(mid_node_i+1)] # backtrack is in reverse order, so count from end
    score = ahh[mid_node_i]
        
    #print(from_node + " " + str(edge))
    return from_node, edge, score

def LSA(v,w):
    """
    LinearSpaceAlignment(v, w, top, bottom, left, right)
        if left = right
            output path formed by bottom − top vertical edges
        if top = bottom
            output path formed by right − left horizontal edges
        middle ← ⌊ (left + right)/2⌋
        midEdge ← MiddleEdge(v, w, top, bottom, left, right)
        midNode ← vertical coordinate of the initial node of midEdge
        LinearSpaceAlignment(v, w, top, midNode, left, middle)
        output midEdge
        if midEdge = "→" or midEdge = "↘"
            middle ← middle + 1
        if midEdge = "↓" or midEdge ="↘"
            midNode ← midNode + 1
        LinearSpaceAlignment(v, w, midNode, bottom, middle, right)
  """
    if len(w) == 0:
        return v,"_" * len(v), -SIGMA * len(v)
    if len(v) == 0:
        return "_" * len(w), w, -SIGMA * len(w)
    #middle = (left + right)//2
    (midNodeI, midNodeJ), midEdge, score = middle_edge(v,w)
    #print(midNodeI,midNodeJ,midEdge)
    #output middle edge
    alignedV, alignedW, _ = LSA(v[0:midNodeI], w[0:midNodeJ])
    #midNodeI -= 1
    #midNodeJ -= 1
    # #, _ = LSA(v[midNodeI:bot],w[middle:right])
    
    if midEdge == 0:
        alignedV += v[midNodeI]
        alignedW += w[midNodeJ]
    elif midEdge == 1:
        alignedV += v[midNodeI]
        alignedW += "_"
    elif midEdge == 2:
        alignedV += "_"
        alignedW += w[midNodeI]
    

    if midEdge == 2 or midEdge == 0:
        midNodeJ += 1
    if midEdge == 1 or midEdge == 0:
        midNodeI += 1
    secV, secW, _ = LSA(v[midNodeI:], w[midNodeJ:])
    alignedV += secV
    alignedW += secW

    return alignedV, alignedW, score



def main():
    v, w = read('rosalind_ba5l.txt')
    # v = "PLEASANTLY"
    # w = "MEANLY"
    newV, newW, score = LSA(v, w)
    print(score)
    print(newV)
    print(newW)
    #find longest path (i, middle) from source to middle


if __name__ == '__main__':
    main()