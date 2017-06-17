################################################################################################
# Uniform (Random) Selection
################################################################################################
#		- not based on fitness, so fitness is not asked for as an arguement - pop.size is used instead
# 		- uses selectlocations(::Symbol, ::Vector, ::Integer) in "Selection - General.jl"
#		- used in GA variants and is the main technique in ES to select for reproduction 
#		- in ES 
#			- fitness based selection is done after reproduction, primarily using truncation selection
#			- count = lambda; 
#			- it is unclear if selection is done with replacement (roulette) or without (sus)
#				 - most of the literature implies selection without replacement -> sus as default
################################################################################################

type UniformSelection <: SelectionType
	sel_type::Symbol
	order::Symbol
end

UniformSelection(;sel_type = :sus, order = :original) = UniformSelection(sel_type, order)

function selectlocations(info::UniformSelection, count::Integer, popn_size::Integer)
		pmf = fill(1.0 / popn_size, popn_size)
		selectloc(info, pmf, count)
end
