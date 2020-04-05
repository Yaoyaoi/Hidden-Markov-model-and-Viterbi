import math as math

# P number of states; X:number of x; xs: set of x
nOP = nOX = xSet= None
# A : transition matrix; E : emission matrix; init: initial distribution
A = []
E = []
init = []
# genome
genome=''

#  reading files
def loadData(HMMFileName, FAFileName): 
    global nOP,nOX,xSet,init,A,E,genome
# get hmm parameters
    fid = open(HMMFileName)
    lines = fid.readlines()
    for i in range(len(lines)):
        line = lines[i]
        if i == 0:
            lineSplit = line.strip('\n').split()
            nOP = int(lineSplit[0])   #2
            nOX = int(lineSplit[1])   #4
            xSet = lineSplit[2]
            continue
        if i == 1:
            lineSplit = list(map(float,line.strip('\n').split()))
            for i in range(nOP):
                init.append(lineSplit[i])
            continue
# get A matrix and E matrix
        lineSplit = list(map(float,line.strip('\n').split()))
        A.append(lineSplit[:nOP])
        E.append(lineSplit[nOP:])

# get genome data from fasta
    with open(FAFileName, 'r') as f:
        for line in f:
            if not line[0] == '>':
                genome = genome + line.rstrip().upper()

# viterbi
def viterbi():
    global nOP,nOX,xSet,init,A,E,genome
    Q = ['state A','state B']
    length = len(genome)
# table initialization
    sKi = [[0 for i in range(length)] for j in range(nOP)]
    pathLog = [[0 for i in range(length)] for j in range(nOP)]

# table filling
    for i in range(nOP):
        sk0 = init[i] * E[i][xSet.index(genome[0])] 
        sKi[i][0] = math.log2(sk0) if sk0 else float('-inf')
    for i in range(1,length):
        for j in range(nOP):
            temp = [sKi[k][i-1] + math.log2(A[k][j]) if A[k][j] else float('-inf') for k in range(nOP)]
            maxI = 0
            for v in range(1,nOP):
                if temp[v] > temp[maxI]:
                    maxI = v
            TsProb = E[j][xSet.index(genome[i])]
            logTsProb = math.log2(TsProb) if TsProb else float('-inf') 
            sKi[j][i] = temp[maxI] + logTsProb
            pathLog[j][i] = maxI

# traceback
    start = end = length
    outputList = []
    maxState = 0
    for j in range(1,nOP):
        if sKi[j][length-1] > sKi[maxState][length-1]:
            maxState = j
    path = [maxState]
    countList = [0 for i in range(nOP)]
    countList[maxState] += 1
    current = maxState
    for i in range(length-1,0,-1):
        path.append(pathLog[path[-1]][i])
        if path[-1] != current:
            outputList.append((start, end, Q[current]))
            current = path[-1]
            end = start - 1
            countList[current] += 1
        start -= 1
    outputList.append((start,end,Q[current]))

    # output/write to file'output.txt'
    outputFile = open('output.txt', mode='w')
    for i in range(len(outputList)-1,-1,-1):
        line = ' '.join(list(map(str,outputList[i])))
        print(line)
        outputFile.write(line+'\n')
    
    for i in range(nOP):
        line = "{} segments of the genome are in ".format(countList[i]) + Q[i] 
        print(line)
        outputFile.write(line+'\n')

if __name__ == "__main__":
    #HMMFileName = "example.hmm"
    #FAFileName = "example.fa"
    HMMFileName = input("hmm file name:")
    FAFileName = input("fasta file name:")
    loadData(HMMFileName,FAFileName)
    viterbi()