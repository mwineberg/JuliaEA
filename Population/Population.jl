
####################################################################################################
# Population Types Definitions
####################################################################################################
#  Population
#	-  on the use of an inner constructor
#			-  phenotype and fitness and optimal are computed during decoding and fitness evaluation,
#			   which occurs afterthe creation of a populaiton.
# 			-  the inner constructor uses new() which allows fields to be left undefined
#			   i.e. the phenotype and fitness fields can be left undefined until decoding/evalutation
####################################################################################################

####################################################################################################
# GeneralPopn Types Definition and Constructuors
####################################################################################################
# General populations are panmictic (anyone can mate with anyone) as opposed to spatial
#		- if no crossover is used spatial is meaningly so General population is used as well
#			although population is technically not panmictic
####################################################################################################

type GeneralPopn <: Population
	chromosomes::Chromosomes
	phenotype::PhenotypeMembers
	fitness::FitnessValues
	optimal::OptimalMembers

	# Inner Constructor used for initial population creation and population creation after reproduction
	#	- i.e.  just chromosomes are supplied and phenotype, fitness and isoptimum are produced later during decode/evaluation steps
	function GeneralPopn(chr::Chromosomes)
		new(chr)
	end

	# Inner Constructor used for creating a population through selection withour reproduction
	#	or by combining multiple populations into a single population
	function GeneralPopn(chr::Chromosomes, phtype::PhenotypeMembers, fit::FitnessValues, opt::OptimalMembers)
		new(chr, phtype, fit, opt)
	end
end

function GeneralPopn(pi::PopulationInfo, encoding::Encoding)
	GeneralPopn(create_chromosomes(pi.popn_size, encoding))
end

function GeneralPopn(pi::PopulationInfo, values::Matrix)
	GeneralPopn(FixedLengthChromosomes(values))
end

function GeneralPopn(pi::PopulationInfo, values::Vector)
	GeneralPopn(VarLengthChromosomes(values))
end


####################################################################################################
# SpatialPopn Types Definition and Constructuors
####################################################################################################
# Spatial populations only allow mating within a neighborhood
#		- a neighborhood can be either direct neighbors, or neighbors of neighbors, etc.
#		- for future expansion
####################################################################################################

type SpatialPopn <: Population
	chromosomes::Chromosomes
	phenotype::PhenotypeMembers
	fitness::FitnessValues
	optimal::OptimalMembers
	structure::Adjacent				# either an adjacency matrix or list - Type not immutable (perhpas can change structure over generations)

	# Inner Constructor used for initial population creation and population creation after reproduction
	#	- i.e.  just chromosomes are supplied and phenotype, fitness and isoptimum are produced later during decode/evaluation steps
	function SpatialPopn(struct::Adjacent, chr::Chromosomes)
		new(struct, chr)
	end

	# Inner Constructor used for creating a population through selection withour reproduction
	#	or by combining multiple populations into a single population
	function SpatialPopn(struct::Adjacent, chr::Chromosomes, phtype::PhenotypeMembers, fit::FitnessValues, opt::OptimalMembers)
		new(struct, chr, phtype, fit, opt)
	end
end

function SpatialPopn(pi::PopulationInfo, struct::Adjacent, encoding::Encoding)
	SpatialPopn(struct, create_chromosomes(pi.popn_size, encoding))
end

function SpatialPopn(pi::PopulationInfo, struct::Adjacent, values::Matrix)
	SpatialPopn(struct, FixedLengthChromosomes(values))
end

function SpatialPopn(pi::PopulationInfo,struct::Adjacent, values::Vector)
	SpatialPopn(struct, VarLengthChromosomes(values))
end


####################################################################################################
# create_population Constructuors
####################################################################################################
#	- uses the PopulationInfo to determine which type of population to create
####################################################################################################

function create_population(pi::PopulationInfo, encoding::Encoding)
	PT = pi.PopnType
	PT(pi, encoding)
end

function create_population(pi::PopulationInfo, encoding::Encoding, chromosomes::Array)
	PT = pi.PopnType
	PT(pi, encoding, chromosomes)
end



####################################################################################################
# Population properties
####################################################################################################

# note: does not use popn.popn_info.popn_size
#			- because the size of the population can change during a run
#			- e.g. in (mu, lambda) Evolutionary Strategies the initial popn size == mu == popn_size,
#				   but the actual popn size becomes mu * lambda during the run
function popnsize(popn::Population)
	popnsize(popn.chromosomes)
end


# only use after isoptimal!() or evaluate!(popn, problem) has been run
#	otherwise an error is thrown
function goodenough(popn::Population)
	goodenough(popn.optimal)
end


####################################################################################################
# Evaluation methods
####################################################################################################

function decode!(popn::Population, encoding::Encoding)
	popn.phenotype = decode(popn.chromosomes, encoding)
end

function evaluate!(p::Population, setup::EvalSetup)
	p.fitness = (usesphenotype(p.phenotype) ? evaluate(p.phenotype, setup)
											: evaluate(p.chromosomes, setup))
end

function isoptimal!(popn::Population, optinfo::OptimumInfo)
	popn.optimal = isoptimal(popn.fitness, optinfo)
end

function evaluate!(popn::Population, problem::Problem)
	decode!(popn, problem.encoding)
	evaluate!(popn, problem.setup)
	isoptimal!(popn, problem.optinfo)
end

####################################################################################################
# getindex method
####################################################################################################
# 	- utility function - e.g. allows popn[1:5] to be called
#	- returns a tuble of (chromosomes, phenotype, fitness values) but stripped of all type information
# 	- therefore, not to be used to create new populations from selected members of an old population
# 		- use select_population() instead
####################################################################################################

# only use after evaluate!(popn, problem) or each of its compoent functions has have run;
#	otherwise an error is thrown
function getindex(p::Population, selection)
	chr    = p.chromosomes[selection]
	fit    = p.fitness[selection]
	phtype = usesphenotype(p.phenotype) ? p.phenotype[selection] : chr
	opt    = optunknown(p.optimal)     ? p.fitness[selection]   : UnknownOptimum()
	return (chr, phtype, fit, opt)
end

####################################################################################################
# duplicate methods:
####################################################################################################
#	-  selecting members from the population (includes chromosomes, phenotype and fitness)
#   -  creating the next generation's population from the reproduced chromosomes
#   -  evaluating a populaiton
####################################################################################################
function duplicate(chromosomes::Chromosomes, phenotype::PhenotypeMembers, fitness::FitnessValues, optimal::OptimalMembers, selection)
	chr 	= duplicate(chromosomes, selection)
	fit 	= duplicate(fitness, selection)
	phtype 	= usesphenotype(phenotype) ? duplicate(phenotype, selection) : NoPhenotype()
	opt 	= optknown(optimal)        ? duplicate(optimal, selection)   : UnknownOptimum()
	return (chr, phtype, fit, opt)
end

function duplicate{P <: Population}(p::P, selection; chromosome_only = false)
	if chromosome_only
		chr = duplicate(p.chromosomes)
		return P(chr)
	else
		chr, phtype, fit, opt = duplicate(p.chromosomes, p.phenotype, p.fitness, p.optimal, selection)
		return P(chr, phtype, fit, opt)
	end
end

function duplicate(popn::SpatialPopn, selection; chromosome_only = false)
	struct  = p.structure
	if chromosome_only
		chr = duplicate(p.chromosomes)
		return P(struct, chr)
	else
		chr, phtype, fit, opt = duplicate(p.chromosomes, p.phenotype, p.fitness, p.optimal, selection)
		return P(struct, chr, phtype, fit, opt)
	end
end


####################################################################################################
# combine methods:
####################################################################################################
#	-  selecting members from the population (includes chromosomes, phenotype and fitness)
#   -  creating the next generation's population from the reproduced chromosomes
#   -  evaluating a populaiton
####################################################################################################

function combine{P <: Population}(p1::P, p_rest...; chromosome_only = false)
	if chromosome_only
		combined_popn = chromosomes_only_combine(p1, p_rest)
	else
		combined_popn = all_fields_combine(p1, p_rest)
	end
	return combined_popn
end

function combine{P <: SpatialPopn}(p1::P, p_rest...; chromosome_only = false)
	if chromosome_only
		combined_popn = chromosomes_only_combine(p1, p_rest)
	else
		combined_popn = all_fields_combine(p1, p_rest)
	end
	combined_popn.structure = p1.structure
	return combined_popn
end

function all_fields_combine{P <: Population}(p1::P, p_rest)
	components = names(p1)
	popns = Array(PopulationComponent, length(p_rest) + 1, length(components))
	combined_pop = Array(PopulationComponent, length(components))

	for j = 1:length(components)
		component = components[j]
		popns[1,j] = p1.component
	end

	for i = 2:(length(p_rest) + 1), j = 1:length(components)
		component = components[j]
		popns[i,j] = p_rest[i].component
	end

	for j = 1:length(components)
		all_members = slice(popns, j)
		combined_pop[j] = combine(all_members...)
	end

	P(combined_pop...)
end

function chromosomes_only_combine{P <: Population}(p1::P, p_rest)
	all_members = Array(Chromosomes, length(p_rest) + 1)

	all_members[1] = p1.chromosomes

	for i = 2:(length(p_rest) + 1)
		all_members[i] = p_rest[i].chromosomes
	end

	combined_chromosomes = combine(all_members...)

	P(combined_chromosomes)
end
