from dimod import ConstrainedQuadraticModel, Integer

a=Integer('a', upper_bound=5)
b=Integer('b', upper_bound=5)

cqm=ConstrainedQuadraticModel()
# ponemos el - porque minimiza y queremos maximizar
cqm.set_objective(-a*b)
cqm.add_constraint(2a+2b <= 10, "max perimeter")

from dwave.system import LeapHybridCQMSampler
sampler=CQMSampler()
sampleset=sampler.sample_cqm(cqm)

print(sampleset.first)
