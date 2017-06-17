# Uses package Distributions

# need Chromosomes to be preloaded for popnsize() and chrlength()

# 1) Generator locations
#	- fixed number
#		- number supplied by user
#		- usually a percentage of the chromosome length
#	- random number
#		- binomial distribution
#			- mutation rate (prob of mutation) supplied by user
#			- usually a percentage of the chromosome length
# 2) based on value of the gene at the location determine alternate
#	- for binary,
#		- use 1- gene
#	- for symbolic (catagorical) values
#		- use rand(Set(alphabet) - Set(gene))
#	- for integer finite
#	- for integer infinite
#		- creap mutation
#	- for real values
#		- could use Gaussian mutation it is gene[i + randn(0,sigma[i])
#			- sigma[i] could be a
#				- user specified constant (e.g. 1)
#				- dynamic constant based on proporties of the population
#				- evolved constant
#				- evolved chromosome
#					- extends the chromosome
#					- coevolving population
#		- could use a different distribution
#			- must be symmetric
#			- e.g. Cauchy, T, quadratic, etc.
# 3) replace gene with alternative

abstract MutationPerChr

type CountMutn <: MutationPerChr
	count::Integer
	chr_length::Integer
end

type BinomialMutn <: MutationPerChr
	binom::Binomial
	chr_length::Integer
end

type PercentageMutn <: MutationPerChr
	percent::Real
end

type ProbabiltyMutn <: MutationPerChr
	prob::Real
end

function MutnPerChr(p::Real, chr_length::Integer; ptype = :percentage)
	(ptype == :probability	? (0 <= p <= 1 ? BinomialMutn(Binomial(chr_length, p), chr_length)
					  					  	: BinomialMutn(Binomial(chr_length, p/chr_length), chr_length)) :
	(ptype == :percentage  	? (0 <= p <= 1 ? CountMutn(ifloor(p * chr_length), chr_length)
					  					  	: error("Mutation percentage p = $p isn't between 0 and 1")) :
	error("ptype = $ptype; must be either :percentage or :probability")))
end

function MutnPerChr(count::Integer, chr_length::Integer; ctype = :count)
	if !(0 <= count <= chr_length)
		error("count = $count is greater than the chromosome length = $chr_length")
	end

	(ctype == :count 	   		? CountMutn(count, chr_length) :
	(ctype == :probability 	? BinomialMutn(Binomial(chr_length, count / chr_length), chr_length) :
	error("ctype = $ctype; must be either :count or :probability")))
end

function MutnPerChr(p::Real; ptype = :percentage)
	if !(0 <= p <= 1)
		error("Mutation probability/percentage p = $p must be greater than 0")
	end

	(ptype == :probability ? ProbabiltyMutn(p) :
	(ptype == :percentage  ? PercentageMutn(p) :
	error("ptype = $ptype; must be either :percentage or :probability")))
end

# note: for chrs::Vector{Vector} using count::Integer does not make sense since each chromosome has a different length
#	- only a percentage of the chromosomes and or probability of selecting a gene given the specific chromosome size makes sense
#	- therefor not defined

#############################################

mutatecount(info::CountMutn)														= info.count
mutatecount(info::BinomialMutn)													= rand(info.binom)[1]
mutatecount(info::PercentageMutn, chr_length::Integer)	= ifloor(info.percent * chr_length)
mutatecount(info::ProbabiltyMutn, chr_length::Integer)	= rand(Binomial(chr_length, info.prob))[1]

randomloc(mpc::MutationPerChr)					    						= sample(1:mpc.chr_length, mutatecount(mpc), replace = false)
randomloc(mpc::MutationPerChr, chr_length::Integer) 		= sample(1:chr_length, mutatecount(mpc, chr_length), replace = false)

################################################################################################


abstract Mutation

type BinaryMutation <: Mutation
	per_chr::MutationPerChr
end

type SymbolicMutation <: Mutation
	per_chr::MutationPerChr
	alphabet						# is a set. Can be implemented as a range, integer vector, symbol vector
end

typealias CategoricalMutation SymbolicMutation

type GaussianMutation <: Mutation
	per_chr::MutationPerChr
	sigma
end

function alleles(mut_type::BinaryMutation, genes::Vector)
	1 - genes
end

function alleles(info::SymbolicMutation, genes::Vector)
	alt = similar(genes)
	for i = 1:length(genes)
		alt[i] = sample(setdiff(info.alphabet, genes[i]))
	end
	alt
end

function alleles(info::GaussianMutation, genes::Vector)
	randn(length(genes)) * info.sigma
end

function alleles(info::GaussianMutation, genes::Vector, sigma)
	randn(length(genes)) * sigma
end

function sigma!(mut::GaussianMutation, sigma)
	mut.sigma = sigma
end

function mutate(mutn::Mutation, chrs::Array)
	mchrs = deepcopy(chrs)
	loc, genes, alt = mutate!(mutn, mchrs)

	(mchrs, loc, genes, alt)
end

function mutate(mutn::Mutation, chrs::Array, members)
	mchrs = chrs[members]								# <- this line must change for Julia 0.4 - slices won't copy
	loc, genes, alt = mutate!(mutn, mchrs)

	(mchrs, loc, genes, alt)
end

function mutate!(mutn::Mutation, chrs::Matrix)
	popn_size = popnsize(chrs)
	chr_length = chrlength(chrs)

	loc = Array(Vector, popn_size)
	genes = Array(Vector, popn_size)
	alt = Array(Vector, popn_size)

	for i = 1:popn_size
#		println("i = $i")
		loc[i] = randomloc(mutn.per_chr)
#		println("loc[i] = $(loc[i])")
		genes[i] = vec(chrs[i,loc[i]])
#		println("genes[i] = $(genes[i])")
		alt[i] = alleles(mutn, genes[i])
#		println("alt[i] = $(alt[i])")
		chrs[i, loc[i]] = alt[i]
	end

	(loc, genes, alt)
end

function mutate!(mutn::Mutation, chrs::Matrix, members)
	chr_length = chrlength(chrs)

	loc = Array(Vector, length(members))
	genes = Array(Vector, length(members))
	alt = Array(Vector, length(members))

	for (i, m) in enumerate(members)
		loc[i] = randomloc(mutn.per_chr)
		genes[i] = vec(chrs[m,loc[i]])
		alt[i] = alleles(mutn, genes[i])
		chrs[m, loc[i]] = alt[i]
	end

	(loc, genes, alt)
end


function mutate!(mutn::Mutation, chrs::Vector)
	popn_size = popnsize(chrs)
	chr_lengths = chrlength(chrs)			# since each chromosome can be a different size, chrlength() returns a length for each chromosome

	loc = Array(Vector, popn_size)
	genes = Array(Vector, popn_size)
	alt = Array(Vector, popn_size)

	for i = 1:popn_size
		loc[i] = randomloc(mutn.per_chr, chr_lengths[i])
		genes[i] = chrs[i][loc[i]]
		alt[i] = alleles(mutn, genes[i])
		chrs[i][loc[i]] = alt[i]
	end

	(loc, genes, alt)
end

function mutate!(mutn::Mutation, chrs::Vector, members)
	chr_lengths = chrlength(chrs)			# since each chromosome can be a different size, chrlength() returns a length for each chromosome

	loc = Array(Vector, length(members))
	genes = Array(Vector, length(members))
	alt = Array(Vector, length(members))

	for (i, m) in enumerate(members)
		loc[i] = randomloc(mutn.per_chr, chr_lengths[i])
		genes[i] = chrs[m][loc[i]]
		alt[i] = alleles(mutn, genes[i])
		chrs[m][loc[i]] = alt[i]
	end

	(loc, genes, alt)
end
