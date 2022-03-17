from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
from dimod import BinaryQuadraticModel
import dwave.inspector

# Define problem QUBO as dictionary

linear = {('x','x'):-1, ('z','z'):-1, ('y','y'):-1}
quadratic = {('y','z'):2, ('x','y'):2, ('x','z'):2 }
Q = {**linear,**quadratic}

# Convert to BQM to submit to QPU

bqm=BinaryQuadraticModel.from_qubo(Q)

# Define sampler, auto minor-embedding 

sampler = EmbeddingComposite(DWaveSampler())

# Run the problem and print the results

sampleset = sampler.sample(bqm, num_reads=1000)
print(sampleset)

# Show problem details and solver details

dwave.inspector.show(sampleset)
