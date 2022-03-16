from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
import dimod 
import dwave.inspector

# NOT GATE formulation qubo minor embedding unespecified sampler 

sampler = EmbeddingComposite(DWaveSampler())
Q = {('x','x'):-1, ('x','z'):2, ('z','z'):-1}

sampleset = sampler.sample_qubo(Q, num_reads=6000)
print(sampleset)


# AND gate formulation qubo minor embedding unespecified sampler 
# BQM formulated as QUBO, penality function: xy-2(x+y)z+3z

sampler = EmbeddingComposite(DWaveSampler())
Q = {('z','z'):3, ('x','z'):-2, ('y','z'):-2, ('x','y'):1}

sampleset = sampler.sample_qubo(Q, num_reads=6000)
print(sampleset)
dwave.inspector.show(sampleset)

# AND gate 

Q = {('z','z'):3, ('x','z'):-2, ('y','z'):-2, ('x','y'):1}
bqm = dimod.BinaryQuadraticModel.from_qubo(Q, offset=0.0)

# AND gate

from dimod.generators import and_gate
bqm=and_gate('in1','in2','out')
sampler=EmbeddingComposite(DWaveSampler())
sampleset=sampler.sample(bqm, num_reads=1000)
