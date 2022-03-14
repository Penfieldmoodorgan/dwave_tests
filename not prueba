from dwave.system import DWaveSampler, EmbeddingComposite
sampler = EmbeddingComposite(DWaveSampler())

Q = {('x','x'):-1, ('x','z'):2, ('z','x'):0, ('z','z'):-1}

sampleset = sampler.sample_qubo(Q, num_reads=6000)
print(sampleset)

sampleset.first.sample["x"] =! sampleset.first.sample["z"]
