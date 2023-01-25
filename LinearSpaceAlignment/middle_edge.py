from turtle import back
import numpy as np
#import BLOSUM62
from Bio.Align import substitution_matrices
sub_mat = substitution_matrices.load("BLOSUM62")


def read(file_name):
    
    with open(file_name, 'r+') as file:
        #k = int(file.readline().strip())
        all = [line.strip() for line in file.readlines()]
        return all[0], all[1]

def middle_edge(v,w):
    """ 
    PLEASANTLY
    MEASNLY

    """

    m = len(w)
    middle = m//2
    sigma = 5
    n = len(v)
    linear_source = np.zeros(n+1, dtype=int) #creating 1D array to hold 
    linear_source[0] = 0 
    for i in range(1, n+1):
        linear_source[i] = linear_source[i-1] - sigma
        
    print(w , " middle is " , middle)
    #print(linear)


    backtrack = np.zeros(n+1, dtype=int)

    for i, letter in enumerate(w):
        temp = linear_source[0]
        linear_source[0] += -sigma
        for k in range(1, n+1):
            diag = temp + sub_mat[letter, v[k-1]]   # 0
            vert = linear_source[k-1]  - sigma             # gap vertical # 1
            horz = linear_source[k]  - sigma               # gap horizontal # 2
            incoming = [diag, vert, horz]
            max_edge = np.argmax(incoming) #max score edge
            backtrack[k] = max_edge 
            temp = linear_source[k]

            linear_source[k] = incoming[max_edge]
        
        #print(linear)
        if i + 1== middle:
            #print(letter)
            break
    
    linear_sink = np.zeros(n+1, dtype=int) #creating 1D array to hold 
    linear_sink[0] = 0
    for i in range(1, n+1):
        linear_sink[i] = linear_sink[i-1] - sigma

    sink_w = w[::-1]  #reversing strings
    sink_v = v[::-1] 

    #print("sink")
    #print(linear_sink)
    midd = len(sink_w)//2

    #going from sink to middle is that same as reversing strings and going to middle
    for i, letter in enumerate(sink_w):
        temp2 = linear_sink[0]
        linear_sink[0] += -sigma
        for k in range(1, n+1):
            diag = temp2 + sub_mat[letter, sink_v[k-1]]# 0
            vert = linear_sink[k-1]  - sigma             # gap vertical # 1
            horz = linear_sink[k]  - sigma               # gap horizontal # 2
            incoming = [diag, vert, horz]
            max_edge = np.argmax(incoming) #max score edge

            temp2 = linear_sink[k]

            #
            
            linear_sink[k] = incoming[max_edge]
            #,missing backtrack arrays
        
        #print(linear_sink)
        if i == midd:
            #print(letter)
            break
    
    #summing up the arrays
    flipped = linear_sink[::-1]
    #print("linear source", linear)
    #print("linear sink ", flipped)
    ahh = linear_source + flipped
    #print(ahh)
    mid_node_j = np.argmax(ahh) 

    

    #print(mid_node_j)

    """
    from_source, to_sink, middle_column
    [-15 -11  -7  -3   6   1  -4  -9 -14 -19 -24] 
    [ -9  -4   1   6  11  13  12   7   1  -8 -20] 
    [-24 -15  -6   3  17  14   8  -2 -13 -27 -44]

    """
    
    from_node = ""
    #print(backtrack[mid_node_j])

    if(backtrack[mid_node_j]==0):
        from_node = "("+ str(mid_node_j+1) + ", " + str(middle+1) +")"
    
    elif(backtrack[mid_node_j]==1):
        from_node = "("+ str(mid_node_j) + ", " + str(middle+1) +")"
    else:
        from_node = "("+ str(mid_node_j+1) + ", " + str(middle) +")"
        
    print("("+str(mid_node_j) +", "+str(middle)+")" + " " + from_node)



def main():
    v, w = read('rosalind_ba5k.txt')
    middle_edge(v, w)
    #find longest path (i, middle) from source to middle


if __name__ == '__main__':
    main()
    