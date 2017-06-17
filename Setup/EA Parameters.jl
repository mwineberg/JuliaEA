##########################################################################
#	EA Parameters
##########################################################################
#	popn_size						{initial size of the population and the size the population always returns to by the start of each generation}
#	ρ										{ratio of chromosome count after vs before reproduction}
#	ξ										{number of elites}
# λ 									{chromosome length; must be greater than 1 or, if variable length, be set to -1}
#	pμ_chr 							{percent/prob that a chromosome will be mutated}
#	pμ_gene							{percent/prob that a given gene within a chromosome will be mutated}
#	pχ 									{percent/prob of the max number of possible crossovers (AKA the probability of crossover)}
# parents_per_xover 	{number of parents used per individual crossover application}
#	offpring_per_xover	{number of offpring produced per individual crossover application}
#	include_parents 		{used to prepare for the reduction step}
##########################################################################
#	notes
#		True probability of gene mutation:		 Pμ,gene = pμ_chr * pμ_gene
#		True probability of chromosome mutation: Pμ,chr  = pμ_chr * (1 - (1 - pμ_gene)^λ)   where   λ = chromosome length
#		All variation (xover and mutation) ratio parameters {p_ratio = {pχ, pμ_xchr, pμ_nxchr, pμ_gene}
#			can be either probabilities or percentages
#			- as percentages, the counts are determined using p_ratio * largest_count
#				- count is floored if the result is not an integer
#				- the counts they produce never change throughout the entire run
#			- as probabilities, the counts are determined using Binom(largest_count, p_ratio)
#				- counts can change from generation to generation
#			- for pμ_gene, counts change from chromosome to chromosome
#				- consequently the paramenter is passed to the mutation method
#				  so the mutation count can be computed for each chromosome
##########################################################################
#	In practice:
#		- ρ
#			- it is >= 1 for all current methods
#			- GA: ρ = 1
#			- ES: ρ = λ/μ
#		- ξ
#			- GA: user provided, usually 1 or k where k is small relative to the population size
#			- ES: does not directly use elitism to my knowledge (must look into this further)
#					- relies on (μ+λ) with trunctation selection to automatically produce "elites"
#					- although if these "elites" have Gaussian mutation applied, they are not really the same
#		- pμ_chr and pμ_gene
#			- GA: some systems have pμ_chr = 1 and pμ_gene as a probability 		{Goldberg and other textbooks}
#				  others have pμ_chr as a probability, and pμ_gene as a percentage 	{many/most implementations}
#				  where pμ_gene is either fixed or randomly generated using a uniform distribution (often underspecified)
#				  this underspecified approach is not yet implemented here
#			- ES: pμ_chr = 1 and pμ_gene = 1, i.e. both are percentages, but not really meaningful parameters
#					- because Gaussian mutation is used instead of binary or symbolic/catagorical mutation
#				  	- the amount of mutation is determined by the variance of the Gaussian
#					- usually dynamically determined by the "evoluation strategy" used
#		- pχ
#			- GA: treats it as a probabily
#			- ES: treats it as a percentage, commonly set to 1
# 		- parents_per_crossover
#			- GA: always set to 2
#			- ES: user defined
#			- this implementation: user defined; can be any integer >= 2
# 		- offspring_per_crossover
#			- GA: usually set to 2, sometimes set to 1 (especially with spatial GAs)
#			- ES: alwasy set to 1
#			- this implementation: user defined; an integer >= 1, up to the number of parents per crossover
#		- include_parents::Bool
# 			- used to prepare for the reduction step (selection in ES)
#			- a form of implicit elitism
#			- ES: a (μ+λ) system; user specified
#			- GA: rarely performed except with NSGA2 (and later versions)
##########################################################################

abstract VariationParm

immutable  PercentageParm <: VariationParm
	value::Real
end

function VariatnParm(value; name = "")
	if length(value) >= 2
		return (value[2] == :prob || value[2] == :probability 	? ProbabilityParm(value[1]) :
				value[2] == :percent || value[2] == :percentage ? PercentageParm(value[1])  :
				error(":$(value[2]) is not recognized in $name = $value. The parameter must be set to either :prob or :percent"))
	else
		return PercentageParm(value)
	end
end

immutable ProbabilityParm <: VariationParm
	value::Real
end

function modificationcount(percent::PercentageParm, total::Integer)
	ifloor(percent.value * total)
end

function modificationcount(prob::ProbabilityParm, max_count::Integer)
	rand(Binomial(max_count, prob.value))
end

function update(percent::PercentageParm)
	false
end

function update(prob::ProbabilityParm)
	true
end

immutable EA_Parameters
	popn_size::Integer
	ρ::Real							# ratio of chr count after vs before reprodn: i.e size of the interim population before reduction
	ξ::Integer 						# number of elites
	λ::Integer 						# chromosome length; must be > 1 or -1 if variable length
	pμ_xchr::VariationParm 			# percent/prob that an offspring from the crossover operation will be mutated
	pμ_nxchr::VariationParm 		# percent/prob that a chromosome selected for reproduction, but is not a result of crossover, will be mutated
	pμ_gene::VariationParm 			# percent/prob that a given gene within a chromosome will be mutated
	pχ::VariationParm 				# percent/prob of the max number of possible crossovers (AKA probability of crossover)
	parents_per_xover::Integer
	offpring_per_xover::Integer
	include_parents::Bool 			# used to prepare for the reduction step (selection in ES); ES(μ+λ); a form of implicit elitism
	max_gen::Integer
	break_when_found::Bool
end

function chr_mutation_parms(pμ_chr, pμ_xchr, pμ_nxchr)
	if isnan(pμ_xchr) && isnan(pμ_nxchr)
		return (VariatnParm(pμ_chr,   name = "pμ_chr/prob_chr_mutated"),
						VariatnParm(pμ_chr,   name = "pμ_chr/prob_chr_mutated"))
	elseif isnan(pμ_xchr)
		warn("pμ_xchr is not provided, pμ_chr = $pμ_chr is used instead")
		return (VariatnParm(pμ_chr,   name = "pμ_chr/prob_chr_mutated"),
						VariatnParm(pμ_nxchr, name = "pμ_nxchr/prob_uncrossed_chr_mutated"))
	elseif isnan(pμ_nxchr)
		warn("pμ_nxchr is not provided, pμ_chr = $pμ_chr is used instead")
		return (VariatnParm(pμ_xchr,  name = "pμ_xchr/prob_xover_offspring_mutated"),
						VariatnParm(pμ_chr,   name = "pμ_chr/prob_chr_mutated"))
	else
		return (VariatnParm(pμ_xchr,  name = "pμ_xchr/prob_xover_offspring_mutated"),
						VariatnParm(pμ_nxchr, name = "pμ_nxchr/prob_uncrossed_chr_mutated"))
	end
end

function EA_Parameters(;popn_size = 100,
						include_parents = false,			# should not be set true if elite_size > 0 (allows elite even more representation)
						parents_per_xover = 2,
						offpring_per_xover = 1,
						max_gen = 200, 						# default is arbitrary; would not include it, except that keyword needs a default
						break_when_found = true,
						ρ = 1.0, 				reproduction_ratio = NaN,
						ξ = 1, 					elite_size = NaN,
						λ = 32, 				chr_length = NaN,
						pχ = 0.75, 			prob_xover = NaN,
						pμ_chr = 0.75, 	prob_chr_mutated = NaN,
						pμ_xchr = NaN,	prob_xover_offspring_mutated = NaN,
						pμ_nxchr = NaN,	prob_uncrossed_chr_mutated = NaN,
						pμ_gene = 0.05, prob_gene_mutated = NaN)

	ρ					= isnan(reproduction_ratio) 		   		?  ρ 	   		: reproduction_ratio
	ξ					= isnan(elite_size) 				   				?  ξ 	   		: elite_size
	λ 		 		= isnan(chr_length) 				   				?  λ 	   		: chr_length
	pμ_chr		= isnan(prob_chr_mutated) 			   		?  pμ_chr   : prob_chr_mutated
	pμ_xchr		= isnan(prob_xover_offspring_mutated) ?  pμ_xchr  : prob_xover_offspring_mutated
	pμ_nxchr	= isnan(prob_uncrossed_chr_mutated)   ?  pμ_nxchr : prob_uncrossed_chr_mutated
	pμ_gene		= isnan(prob_gene_mutated) 		   			?  pμ_gene  : prob_gene_mutated
	pχ				= isnan(prob_xover) 				   				?  pχ 	   	: prob_xover
	pμ_xchr, pμ_nxchr = chr_mutation_parms(pμ_chr, pμ_xchr, pμ_nxchr)
	pμ_gene		= VariatnParm(pμ_gene, name = "pμ_gene/prob_gene_mutated")
	pχ				= VariatnParm(pχ, 	   name = "pχ/prob_xover")
	EA_Parameters(popn_size, ρ, ξ, λ, pμ_xchr, pμ_nxchr, pμ_gene, pχ,
				  parents_per_xover, offpring_per_xover, include_parents, max_gen, break_when_found)
end
