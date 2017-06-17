####################################################################################################
# Type Definitions - uses default Constructuors
####################################################################################################
# [Single Objective] SO_FitnessValues   
# [Multi-Objective] MO_FitnessValues
####################################################################################################

type SO_FitnessValues <: FitnessValues
	values::Vector{Real}
end

type MO_FitnessValues <: FitnessValues
	values::Vector{Vector{Real}}
end

####################################################################################################
# FitnessValues methods  (works on all FitnessValues subtypes)
####################################################################################################
#	isoptimal(::FitnessValues, ::OptimumInfo)
#		- determines whether each fitness value is optimal to within some tolerance epsilon
#		- returns the answer as a Bool vector
#	getindex(::FixedLengthChromosomes, i)
#		- allowing indexing directly using [] to select fitness values
#		- i can be a range e.g. 10:20 or 1:3:27, or it can be a vector of indicies
#		- no external type information is retained (i.e. not wrapped in the FitnessValues type)
# 	- duplicate select all FitnessValues(s) i, where i can be a range; returns a FitnessValues object
# 		- calls the appropriate getindex() function (i.e implements the [] operator) 
#		- then wraps the result the appropriate Chromosome subtype C
#	- combine combines two FitnessValues vectors together using a vertical concatenation
#		and constructs a FitnessValues instance of the same type that was passed in (either SO_FitnessValues or MO_FitnessValues)
#	- show / print functions - includes printing all fitness values or selected fitness values
####################################################################################################

function popnsize(fit::FitnessValues)
	size(fit.values,1)
end

function isoptimal(fitness::Vector, opt::OptimumInfoExists)
	(opt.maximize ? opt.fitness - fitness      <= opt.epsilon 
				  : fitness     -  opt.fitness <= opt.epsilon)
end

function isoptimal(fitness::FitnessValues, opt::OptimumInfoExists)
	OT = opt.OptimumType
	OT(isoptimal(fitness.values, opt))
end

function isoptimal(fitness::FitnessValues, opt::NoOptimumInfo)
	UnknownOptimum()
end

function getindex(fit::FitnessValues, i)
	fit.values[i]
end

function duplicate{FV <: FitnessValues}(fit::FV, i)
	FV(getindex(fit,i))
end 

function combine{FV <: FitnessValues}(fit_1::FV, fit_2::FV)
	return FV(vcat(fit_1.values, fit_2.values))
end

import Base.print
function print(fit::FitnessValues, i::Integer)
	print(fit.values[i])
end

function print(fit::FitnessValues)
	for i = 1:size(fit.values, 1)
		print("[", i, "] ")
		print(fit.values[i])
		print("\n")
	end
end

####################################################################################################
# SO_FitnessValues methods
####################################################################################################
#	- currently there is no "non-intersectional" methods  
#		(functions with arguments of more than one type) that are different from MO fitnesses
####################################################################################################



####################################################################################################
# MO_FitnessValues methods
####################################################################################################
#	- currently there is no "non-intersectional" methods  
#		(functions with arguments of more than one type) that are different from SO fitnesses
####################################################################################################


