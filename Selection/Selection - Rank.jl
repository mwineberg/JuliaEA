import StatsBase.tiedrank

# not yet fully defined (is it a subset of Default or CDFinv?)
type LinearRankSel <: 
	sel_type::Symbol
	order::Symbol
	maximizing::Bool
	slope::Real
end

function LinearRankSelection(slope = 1.0; sel_type = :roulette, order = :original, maximizing = true)
	LinearRankSelection(sel_type, order, maximizing, slope)
end 

################################################################################################
# Rank Selection
################################################################################################
# 	- uses selectlocations(::Symbol, ::Vector, ::Integer) in "Selection - General.jl"
################################################################################################

# DefaultRankSel not yet defined
function selectlocations(info::DefaultRankSel, fit::Vector, count::Integer)
	ranks = (info.maximizing ? tiedrank(fit)
							 : length(fit) - tiedrank(fit) + 1)
	pmf = rank2pmf(info, ranks)
	selectloc(info, pmf, count)
	orderloc!(info, loc)
end

# function not complete
# CDFinvRankSel not yet defined
function selectlocations(info::CDFinvRankSel, fit::Vector, count::Integer)
	ranks = (info.maximizing ? tiedrank(fit)
							 : length(fit) - tiedrank(fit) + 1)
	if(info.sel_type == :sus)
		pmf = rank2pmf(info, ranks)
		selectloc(info, pmf, count)
	else
		sample(sel_type, count)
	orderloc!(info, loc)
end


###################################
# Linear Rank Selection
###################################
#	- slope [0, 1] with 
#		- slope = 0 is the same as uniform selection 
#		- slope = 1 is the maximal (unique) slope 
#			where the probability of selection of the least fit member is 0
###################################

function rank2pmf(info::LinearRankSel, ranks::Vector)
	n = length(ranks)
	denom = n * (n - 1) / 2
	m = info.slope / denom
	b = ((n + 1) * (1 - info.slope) / 2 - 1) / denom
	return m * ranks + b
end

###################################
# Geometric Rank Selection
###################################
