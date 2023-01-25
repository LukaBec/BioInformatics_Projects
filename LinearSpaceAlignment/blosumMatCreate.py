from Bio.Align import substitution_matrices
sub_mat = substitution_matrices.load("BLOSUM62")

def printMat():
    f = open("Blosum62_mat.txt", "w")
    print(type(sub_mat))
    f.writelines(str(sub_mat))
    f.close()

printMat()