import numpy as np
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
from dimod import BinaryQuadraticModel
import dwave.inspector

# Reads file, returns clauses, num qubits, num clauses and solution

def read_file(file):
    control=list(map(int, file.readline().split()))
    solution=list(map(str, file.readline().split()))
    
    qubits=control[0]
    numcl=control[1]
    
    sol_bin = 0
    for i in range(len(solution)):
        sol_bin+=int(solution[i])*2**(len(solution)-i-1)
        
    clauses=[]
    for i in range(numcl):
        clauses.append(list(map(int, file.readline().split())))   
        
    return qubits, numcl, sol_bin, clauses

# It counts the times each qubit appears in all the clauses

def times(qubits, clauses): 
    t=np.zeros(qubits)
    for clause in clauses:
        for i in clause:
            t[i-1]+=1
    return t

# Generates a list with the weight assigned to each clause: w = W*num clauses/sum over clauses W ; W = di+dj+dk-3 for i,j,k in clause

def weights(times, qubits, clauses, numcl): 
    temp=np.zeros(numcl)
    for clause in clauses:
        i=clauses.index(clause)
        temp[i]=sum(times[j-1] for j in clause)-3
    norm=sum(temp[i] for i in range(numcl))/numcl
    w=[k/norm for k in temp]
    return w
# WEIGHTS NORMALIZED SUCH THAT MAXIMUM COUPLING REMAINS THE SAME

def weights_max_coupling(times, qubits, clauses, numcl):
    temp=np.zeros(numcl)
# Genera 2 matrius de zeros    
    m=np.zeros([qubits,qubits])
    m1=np.zeros([qubits,qubits])
# (A matriu m: Compta número de cops que apareix una parella; m1: Compta número de cops que apareix parella amb el pes - matriu del coupling)
#  Recorre array que conté clàusules i per a cada clàusula calcula el pes wc = di+dj+dk-3, i suma +1 a m per a cada parella que apareix a la clàusula
# Ens queda matriu upper triangular amb num cops que apareixen parelles a les clàusules
# suma +pes clàusula per a cada parella que apareix a la clàusula a m1. m1 conté couplings (a H en veritat anirà multiplicat x2)      
    for clause in clauses:
        i=clauses.index(clause)
        temp[i]=sum(times[j-1] for j in clause)-3
        
        m[clause[0]-1,clause[1]-1]+=1
        m[clause[1]-1,clause[2]-1]+=1
        m[clause[0]-1,clause[2]-1]+=1
        
        m1[clause[0]-1,clause[1]-1]+=temp[i]
        m1[clause[1]-1,clause[2]-1]+=temp[i]
        m1[clause[0]-1,clause[2]-1]+=temp[i]

# Normalitza de forma que coupling max sigui el mateix a H weighted i non-weighted. Mira quin és el coupling max, busca el valor que tenia a non-weighted
# i defineix la normalització com coupling non-weighted/coupling weighted.        
    value=np.max(m1)
    index=np.where(m1 == value)
    norm=np.max(m1)/m[index[0][0]][index[1][0]]
    w=[k/norm for k in temp]
    return w
    
# Generates problem QUBO matrix for non-weighted problem Hamiltonian as a dictionary

def h_problem(qubits, clauses):
    m=np.zeros([qubits,qubits])

    for clause in clauses:
        m[clause[0]-1,clause[1]-1]+=2
        m[clause[1]-1,clause[2]-1]+=2
        m[clause[0]-1,clause[2]-1]+=2
    
        m[clause[0]-1,clause[0]-1]+=-1
        m[clause[1]-1,clause[1]-1]+=-1
        m[clause[2]-1,clause[2]-1]+=-1

    Q={}

    for i in range(qubits):
        for j in range(qubits):
            value=m[i][j]
            Q[(i,j)]=value

    return Q

# Generates problem QUBO matrix for weighted problem Hamiltonian as a dictionary

def h_weighted(qubits, clauses, w):
    m=np.zeros([qubits,qubits])

    for clause in clauses:
        k=clauses.index(clause)
        m[clause[0]-1,clause[1]-1]+=2*w[k]
        m[clause[1]-1,clause[2]-1]+=2*w[k]
        m[clause[0]-1,clause[2]-1]+=2*w[k]
    
        m[clause[0]-1,clause[0]-1]+=-w[k]
        m[clause[1]-1,clause[1]-1]+=-w[k]
        m[clause[2]-1,clause[2]-1]+=-w[k]

    Q={}

    for i in range(qubits):
        for j in range(qubits):
            value=m[i][j]
            Q[(i,j)]=value

    return Q


# Main program

file=open("n8i102.txt")
file2=open("results.txt", 'w+')

qubits, numcl, sol_bin, clauses = read_file(file)
t=times(qubits, clauses)
w=weights(t, qubits, clauses, numcl)

Q_w=h_weighted(qubits, clauses, w)

bqm_w=BinaryQuadraticModel.from_qubo(Q_w, offset = 0.0)

sampler=EmbeddingComposite(DWaveSampler())
sampleset=sampler.sample(bqm_w, num_reads=1000)
print(sampleset)
file.close()
file2.write("Instance: n8i102.txt")
file2.write("\n")
file2.write("num_reads = 1000")
file2.write("\n")
file2.write("weighted max energy")
file2.write("\n")
file2.write(str(sampleset))
file2.write("\n")
file2.write("\n")

Q=h_problem(qubits, clauses)
bqm=BinaryQuadraticModel.from_qubo(Q, offset = 0.0)
sampler=EmbeddingComposite(DWaveSampler())
sampleset=sampler.sample(bqm_w, num_reads=1000)
print(sampleset)

file2.write("num_reads = 1000")
file2.write("\n")
file2.write("non-weighted")
file2.write("\n")
file2.write(str(sampleset))
file2.write("\n")
file2.write("\n")

w=weights_max_coupling(t, qubits, clauses, numcl)
Q_w=h_weighted(qubits, clauses, w)

bqm_w=BinaryQuadraticModel.from_qubo(Q_w, offset = 0.0)

sampler=EmbeddingComposite(DWaveSampler())
sampleset=sampler.sample(bqm_w, num_reads=1000)
print(sampleset)

file2.write("num_reads = 1000")
file2.write("\n")
file2.write("weighted max coupling")
file2.write("\n")
file2.write(str(sampleset))


file2.close()

dwave.inspector.show(sampleset)
