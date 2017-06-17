# Information used to setup the population during construction
immutable PopulationInfo
	popn_size::Integer 			# inital population size 
	PopnType::Type 				# currently set to either GeneralPopn or SpatialPopn
end 

function PopulationInfo(popn_size::Integer)
	PopulationInfo(popn_size, GeneralPopn)
end