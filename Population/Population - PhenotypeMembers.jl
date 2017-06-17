####################################################################################################
# PhenotypeMembers  Type Definitions - uses default Constructuors
####################################################################################################
abstract PhenotypeUsed <: PhenotypeMembers

type NoPhenotype <: PhenotypeMembers
end

type GeneralPhenotype <: PhenotypeUsed
	members::Vector
end

# converts a 2d matrix into a vector of vectors, i.e. from fixed length chromosomes 
# 	note: NoPhenotype would be faster then converting the chromosomes to a GeneralPhenotype
function GeneralPhenotype(phenotype_matrix::Matrix)
	popnsize, chrlngth = size(phenotype_matrix)
	phenotype_vector = Array(Vector, popnsize)
	for i = 1:popnsize
		phenotype_vector[i] = reshape(phenotype_matrix[i,:], chrlngth)
	end
	return GeneralPhenotype(phenotype_vector)
end

function popnsize(phtype::GeneralPhenotype)
	length(phtype.members)
end

####################################################################################################
# PhenotypeUsed methods  (works on all subtypes of PhenotypeMembersc except NoPhenotype)
####################################################################################################
#	methods inherrited unchanged from PopulationComponents
#		getindex(::PhenotypeMembers, i)
#			- allowing indexing directly using [] to select phenotype members
#			- i can be a range e.g. 10:20 or 1:3:27, or it can be a vector of indicies
#			- no external type information is retained (i.e. not wrapped in the Phenotype type)
# 		duplicate() select all phenotype(s) i, where i can be a range; returns a Phenotype object
# 			- calls the appropriate getindex() function (i.e implements the [] operator) 
#			- then wraps the result the appropriate Phenotype subtype PT
#		combine() combines two phenotype matrices together using a vertical concatenation 
#		  and constructs a PhenotypeMember instance of the same type that was passed in 
####################################################################################################

function usesphenotype(phtype::PhenotypeUsed)
	true
end

function nophenotype(phtype::PhenotypeUsed)
	false
end

import Base.print
function print(phtype::PhenotypeUsed)
	for i = 1:size(phtype.members, 1)
		print("[", i, "] ")
		print(phtype[i])
		print("\n")
	end
end

function print(phtype::PhenotypeUsed, i::Integer)
	print(phtype[i])
end

####################################################################################################
# NoPhenotype methods
####################################################################################################
#	methods overiding methods from PopulationComponents
#		getindex(::NoPhenotype, i)
#			- since there is no phenotype, subsetting it will give you an empty array
# 		duplicate(::NoPhenotype, i) 
# 			- the phentoype field always contains the NoPhenotype() singleton when there is no phenotype 
#		combine() combines two phenotype matrices together using a vertical concatenation 
# 			- the phentoype field always contains the NoPhenotype() singleton when there is no phenotype 
#	also overrides print method since there is no phenotype
####################################################################################################

function usesphenotype(phtype::NoPhenotype)
	false
end

function nophenotype(phtype::NoPhenotype)
	true
end

function getindex(phenotype::NoPhenotype, i)
	Array(Any,0)
end

function duplicate(phtype::NoPhenotype, i)
	NoPhenotype()
end 

function combine(pv1::NoPhenotype, pv_rest...)
	NoPhenotype()
end

function print(pheno::NoPhenotype)
end

function print(pheno::NoPhenotype, i::Integer)
end

####################################################################################################
# GeneralPhenotype methods
####################################################################################################
#	- none not included in Phenotype methods
#	- type used to distinguish from NoPhenotype, but methods kept general to allow for future subtyping 
####################################################################################################




