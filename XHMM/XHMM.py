"""
XHMM
By Luka Becerra
"""
import numpy as np
from scipy.stats import norm
from math import log10
P = 1e-8
Q = 1/6
STATES = ["DEL", "DIP", "DUP"]

def genEmat(row):
    em = np.zeros((3,len(row)), np.longdouble)
    em[0,:] = norm.pdf(row, -3, 1)
    em[1,:] = norm.pdf(row, 0, 1)
    em[2,:] = norm.pdf(row, 3, 1)
    return em

def viterbi(trans, em):
    path = []
    score = np.zeros(em.shape,np.longdouble)
    back = np.zeros(em.shape,int)
    score[:,0] = (np.array([P,1-2*P,P]) * em[:,0])
    for i in range(1,em.shape[1]):
        for k in range(len(STATES)):
            edges = score[:,i-1] * trans[:,k] * em[k,i]
            score[k,i] = np.max(edges)
            back[k,i] = np.argmax(edges)
    idx = np.argmax(score[:,-1])
    for i in range(em.shape[1] - 1, -1, -1):
        path = [idx] + path
        idx = back[idx,i]
    return path

def XHMM(samples, targets, data):
    trans = np.array(  [[1-Q,   1,      0],
                        [P,     1-2*P,  P],
                        [0,     1,      1-Q]])
    for i,row in enumerate(data):
        em = genEmat(row)
        samplePath = viterbi(trans, em)
        t1 = 0
        t2 = -1
        for t1 in range(len(samplePath)):
            if samplePath[t1] != 1:
                for t2 in range(t1+1,len(samplePath)):
                    if samplePath[t2] != samplePath[t1]:
                        break
                t2 -= 1
                break
        
        #print non diploid
        if t2 != -1:
            #calculate qExact
            qExact = 0
            weight = float(1)
            forward = np.zeros((len(STATES),len(samplePath)),np.longdouble)
            back = np.zeros((len(STATES),len(samplePath)),np.longdouble)
            #forward mat
            forward[:,0] = (np.array([P,1-2*P,P]) * em[:,0])
            for j in range(1,len(samplePath)):
                for k in range(len(STATES)):
                    edges = forward[:,j-1] * trans[:,k] * em[k,j]
                    forward[k,j] = np.sum(edges)
            #weight
            for j in range(t1+1,t2+1):
                weight *= trans[samplePath[j-1], samplePath[j]] * em[samplePath[j],j] 
            #backward mat
            back[:,-1] = 1
            for j in range(len(samplePath)-2,t2-1,-1):
                for k in range(len(STATES)):
                    edges = back[:,j+1] * trans[k,:] * em[:,j+1]
                    back[k,j] = np.sum(edges)
            sink = np.sum(forward[:,-1])
            qExact = (forward[samplePath[t1],t1] * weight * back[samplePath[t1],t2]) / sink
            qExact = -10*log10(1-qExact)
            print(samples[i],STATES[samplePath[t1]],targets[t1], targets[t2], qExact)

def readData(file):
    
    with open(file, "r") as f:
        targets = f.readline().strip().split()
        targets = targets[1:]
        samples = np.genfromtxt(file, skip_header=1, usecols=0, dtype=str)
        data = np.genfromtxt(file, skip_header=1, usecols=range(1,1+len(targets)), dtype=np.longdouble)
        return samples, targets, data
    

def main():
    samples,targets,data = readData("XHMM.in.txt")
    XHMM(samples,targets,data)

if __name__ == '__main__':
    main()