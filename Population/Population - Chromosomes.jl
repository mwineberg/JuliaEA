####################################################################################################
# GenotypeMembers Type Definitions
####################################################################################################
#	-  genotype members (population of chromosomes) can have a fixed length, i.e. are homologous
#	-  or they can have varying lengths, i.e. non-homologous
#	-  ideally would want to use an iterator pattern here to have generalized code,
#		  but stepping through each member in the population is a time sensitive bottleneck in the EA
#		  and so needs to be optimized for different encodings
####################################################################################################


####################################################################################################
# Fixed Length GenotypeMembers (FL_GenotypeMembers) Type Definition and Constructuors
####################################################################################################

type FL_GenotypeMembers <: GenotypeMembers
	members::Matrix{Real}
end

function GenotypeMembers(popn_size::Integer, e::FixedLengthEncoding)
	members = Array(Real, popn_size, chr_length(e))

	for i = 1:popn_size
		chr =  create_chromosome(e)

		for locus = 1:e.chr_length
			members[i, locus] =  chr[locus]
		end
	end

	return FL_GenotypeMembers(members)
end

function GenotypeMemebrs(popn_size::Integer, e::FixedLengthEncoding)
	FL_GenotypeMembers(popn_size, e)
end

####################################################################################################
# Variable Length GenotypeMembers (VL_GenotypeMemebers) Type Definition and Constructuors
####################################################################################################

type VL_GenotypeMembers <: GenotypeMembers
	members::Vector{Chromosome}
end

function VL_GenotypeMembers(popn_size::Integer, e::Encoding)
	members = Array(Chromosome, popn_size)

	for i = 1:popn_size
		chr =  e.create_chromosome(e)
		members[i] =  chr
	end

	VL_GenotypeMembers(members)
end

function GenotypeMembers(popn_size::Integer, e::Encoding)
	VL_GenotypeMembers(popn_size, e)
end

####################################################################################################
# Chromosomes Constructuor
#	- depending on the encoding type (FixedLengthEncoding or VarLengthEncoding)
#		either General or FixedLength constructor is called and the chromosomes are created
####################################################################################################

function create_chromosomes(popn_size::Integer, e::Encoding)
	chromosomes_constructor(e)(popn_size, e)
end

function create_chromosomes(pi::PopulationInfo, e::Encoding)
	chromosomes_constructor(e)(pi.popn_size, e)
end

####################################################################################################
# Chromosomes methods  (works on all Chromosomes subtypes)
####################################################################################################
#	- popnsize is determined by the number of chromosome members
#		- unlike phenotype etc., these are always defined
# 	- duplicate select all chromosome(s) i, where i can be a range; returns a Chromosomes object
# 		- calls the appropriate getindex() function (i.e implements the [] operator)
#		- then wraps the result the appropriate Chromosome subtype C
#	- combine combines two chromosome matrices of the same type together using a vertical concatenation
#		and constructs a Chromosome instance of the same type that was passed in (either FixLength or General)
#	- show / print a single chromosome
####################################################################################################

function popnsize(genotypeMembers::Array)
	size(genotypeMembers, 1)
end

function popnsize(popn::GenotypeMembers)
	popnsize(popn.members)
end

function duplicate{G <: GenotypeMembers}(popn::G, i)
	G(popn[i])
end

function combine{G <: GenotypeMembers}(popn1::G, popn2::G)
	G(vcat(popn1.members, popn2.members))
end

import Base.print
function print(popn::GenotypeMembers)
	for i = 1:popnsize(popn)
		print("[", i, "] ")
		print(popn[i])
		print("\n")
	end
end

function print(popn::GenotypeMembers, i::Integer)
	print(popn[i])
end

#import Base.show
#function show(chromosomes::Chromosomes)
#	print(chromosomes)
#end

####################################################################################################
# FL_GenotypeMembers methods
####################################################################################################
#	chrlength(::FL_GenotypeMembers)
#		- computes the chromosome length uses by all chromosomes
#	getindex(::FL_GenotypeMembers, i)
#		- allowing indexing directly using [] to select chromosome members
#		- i can be a range e.g. 10:20 or 1:3:27, or it can be a vector of indicies
#		- no external type information is retained (i.e. not wrapped in the Chromosomes type)
#	getindex(::FL_GenotypeMembers, i, j)
#		- selects a range of genes across a selection of chromosomes
#		- no external type information is retained (i.e. not wrapped in the Chromosomes type)
#		- does not exist for GeneralChromosome because nothing is known about
#			the internal organsization of a chromosome
####################################################################################################

function chrlength(genotypeMembers::Matrix)
	size(genotypeMembers, 2)
end

function chrlength(popn::FL_GenotypeMembers)
	size(popn.members, 2)
end

import Base.getindex
function getindex(popn::FL_GenotypeMembers, i)
	popn.members[i,:]
end

function getindex(popn::FL_GenotypeMembers, i, j)
	popn.members[i,j]
end

####################################################################################################
# VL_GenotypeMembers methods
####################################################################################################
#	- getindex(::VL_GenotypeMembers, i)
#		- allows indexing directly using [] (by declare getindex())
# 		- selects the appropriate chromosome vectors, where i can be a range
#		- removes external type information (i.e. not wrapped in any Chromosomes type)
####################################################################################################

function chrlength(chr::Vector)
	popn_size = popnsize(chr)
	chr_length = Array(Integer, popn_size)
	for i = 1:popn_size
		chr_length[i] = length(chr[i])
	end
	chr_length
end


function getindex(chr::VL_GenotypeMembers, i)
	chr.members[i]
end

function getindex(chr::VL_GenotypeMembers, i, j)
	chr.members[i][j]
end
