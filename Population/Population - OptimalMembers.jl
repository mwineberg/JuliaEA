####################################################################################################
# Type Definitions - uses default Constructuors
####################################################################################################
# stores whether or not each member from the population is optimal (is an optimum),
# i.e. has a fitness within epsilon of the optimal fitness
####################################################################################################

type KnownOptimum <: OptimalMembers
	optimal::Vector{Bool}
end

type UnknownOptimum <: OptimalMembers
end

####################################################################################################
# Optimum methods  (works on all Optimum subtypes unless overridden)
####################################################################################################
#	methods inherrited unchanged from PopulationComponents
#		getindex(::OptimalMembers, i)
#			- allowing indexing directly using [] to select phenotype members
#			- i can be a range e.g. 10:20 or 1:3:27, or it can be a vector of indicies
#			- no external type information is retained (i.e. not wrapped in the OptimalMembers subtype)
# 		duplicate() select all optimal indicators i, where i can be a range
# 			- calls the appropriate getindex() function
#			- then wraps the result the appropriate OptimalMembers subtype as determined during runtime
#		combine() combines two optimal members matrices together using a vertical concatenation
#		  and constructs a OptimalMembers instance of the same type that was passed in
####################################################################################################

function goodenough(opt::OptimalMembers)
	findfirst(opt.optimal, true) > 0
end

function print(opt::OptimalMembers)
	for i = 1:size(opt.optimal, 1)
		print("[", i, "] ")
		print(opt[i])
		print("\n")
	end
end

function print(opt::OptimalMembers, i::Integer)
	print(opt[i])
end


####################################################################################################
# KnownOptimum methods
####################################################################################################

function optknown(opt::KnownOptimum)
	true
end

function optunknown(opt::KnownOptimum)
	false
end

####################################################################################################
# UnknownOptimum methods
####################################################################################################
#	getindex(::UnknownOptimum, i)
#		- since there is no known optimum, subsetting it will give you
#			the tuple (empty array, empty array) for (isclose, isoptimum)
# 	duplicate(::UnknownOptimum, i)
# 		- the optimum field always contains the UnknownOptimum() singleton when there is no known optimum
#	combine() combines two optimum objects together
# 		- the optimum field always contains the UnknownOptimum() singleton when there is no known optimum
#	overrides print method since there is no known optimum
####################################################################################################

function optknown(opt::UnknownOptimum)
	false
end

function optunknown(opt::UnknownOptimum)
	true
end

function goodenough(opt::UnknownOptimum)
	false
end

function getindex(opt::UnknownOptimum, i)
	Array(Any,0)
end

function duplicate(opt::UnknownOptimum, i)
	UnknownOptimum()
end

function combine(opt::UnknownOptimum, pc_rest...)
	UnknownOptimum()
end

function print(pheno::UnknownOptimum)
end

function print(pheno::UnknownOptimum, i::Integer)
end
