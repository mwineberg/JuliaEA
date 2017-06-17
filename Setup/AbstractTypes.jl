typealias Chromosome Vector

# Problem Parameters
abstract Problem
abstract EvalSetup
abstract Encoding
abstract OptimumInfo

# Population and Components
abstract Population
abstract PopulationComponent
abstract GenotypeMembers  <: PopulationComponent  # single member called a "Chromosome"
abstract PhenotypeMembers <: PopulationComponent  # single member called a "Body"
abstract FitnessValues    <: PopulationComponent
abstract OptimalMembers   <: PopulationComponent

abstract Adjacent				# either an adjacency matrix or list - for future expansion

# Selection
abstract SelectionType
abstract RankSelection <: SelectionType
