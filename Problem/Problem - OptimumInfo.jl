####################################################################################################
# Type Definitions - uses default Constructuors
####################################################################################################
# stores information about the optimum if it is known 
# 	- stored in the immutable type Problem
####################################################################################################
abstract OptimumInfoExists <: OptimumInfo

immutable NoOptimumInfo <: OptimumInfo
	OptimumType::Type				# the Type to return from isoptimal - currently can only be KnownOptimum or UnknownOptimum; may change when MO are added

	function NoOptimumInfo()
		new(UnknownOptimum)			# only return type possible if there is no known optimum
	end
end


immutable SO_OptimumInfo <: OptimumInfoExists
	OptimumType::Type				# the Type to return from isoptimal - currently can only be KnownOptimum or UnknownOptimum; may change when MO are added
	maximize::Bool 					# true if finding the maximum; false if finding the minimum - also stored in Problem itself (stored here for convinience)
	fitness::Real 					# the optimal fitness value (may correspond to more than one solution)
	epsilon::Real 					# the maximum distance from the optimal fitness to still be considered a success
end

function SO_OptimumInfo(OptimalType::Type, maximize::Bool, fit::Real)
	zero = convert(typeof(fit), 0)
	SO_OptimumInfo(OptimalType, maximize, fit, zero)
end

function SO_OptimumInfo(maximize::Bool, fit::Real, epsilon::Real)
	SO_OptimumInfo(KnownOptimum, maximize, fit, epsilon)
end

function SO_OptimumInfo(maximize::Bool, fit::Real)
	zero = convert(typeof(fit), 0)
	SO_OptimumInfo(maximize, fit, zero)
end


####################################################################################################
# OptimumInfo methods  (works on all OptimumInfo subtypes unless overridden)
####################################################################################################

function optknown(optinfo::OptimumInfoExists)
	true
end

function optunknown(optinfo::OptimumInfoExists)
	false
end

####################################################################################################
# NoOptimumInfo methods - overides in case of no optimum information (e.g. optimum unkown)
####################################################################################################

function optknown(optinfo::NoOptimumInfo)
	false
end

function optunknown(optinfo::NoOptimumInfo)
	true
end
