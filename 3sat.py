from dwave.system import DWaveSampler, EmbeddingComposite
sampler = EmbeddingComposite(DWaveSampler())

linear = {('x','x'):-1, ('z','z'):-1, ('y','y'):-1}
quadratic = {('y','z'):2, ('x','y'):2, ('x','z'):2 }
Q = {**linear,**quadratic}

sampleset = sampler.sample_qubo(Q, num_reads=6000)
print(sampleset)
