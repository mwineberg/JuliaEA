# population as a matrix
# row is member of the population
# col is the gene

import StatsBase.sample
import StatsBase.weights

################################################################################################
#  Helper function
################################################################################################

# creates a new vector that contains the nomralized values from the current vector
# 	- used to renormalized selection probabilities

function normalize(v::Vector)
	v/sum(v)
end

function normalize!(v::Vector)
	vsum = sum(v)
	for i = 1:length(v)
		v[i] = v[i] / vsum
	end
end

#
# restriction: length(counts) >= length(values)
#	- if length of values is smaller than counts, remaining counts ignored
#	- if length values is greater than counts, an array out of bound error is thrown
# example:
# 	rep([3, 4, 2], ["a", "b", "c"]) returns ["a", "a", "a", "b", "b", "b", "b", "c", "c", "c"]

function rep{T}(counts::Vector, values::AbstractArray{T,1})
	x_length = sum(counts)
	x = Array(T, x_length)
	i = 1
	for (j, value) in enumerate(values)
		for r = 1:(counts[j])
			x[i] = value
			i += 1
		end
	end
	return x
end

function rep(counts::Vector)
	rep(counts, 1:length(counts))
end

################################################################################################
# General Selection Routines: Public Methods
################################################################################################

# Sets the order field in the selection information - currently should be :original, :sorted, :randomized
function setorder!(sel_info::SelectionType, order::Symbol)
	sel_info.order = order
end

function orderloc!(info::SelectionType, loc::Vector)
	orderloc!(info.order, loc)
end

# Called from many other selection routines such as FPS, RS, and uniform select (but not TS)
function selectloc(info::SelectionType, prob::Vector, chr_count::Integer)
	if info.sel_type == :roulette
		roulette(prob, chr_count, order = info.order)
	elseif info.sel_type == :sus
		sus(prob, chr_count, order = info.order)
	elseif info.sel_type == :baker_sus
		sus(prob, chr_count, order = info.order)
	else
		error("selection_type must be :roulette, :sus or :baker_sus")
	end
end

################################################################################################
## Internal Routines
################################################################################################

################################################################################################
# selectorder
################################################################################################
# post-process ordering of the selected indices
# 	- can be used with any selection routine
#	- should be used as a final call
#	- currently used with any selection that calls roullette or sus
#	  as well as in tournament selection
#	- order symbols:
#		[:natural]	  leaves the order as is, i.e. as produced by the original selection algorithm
#		[:randomized] randomly permutes the order of the selected indices (ie uses shuffle)
#		[:sorted]	  sorts the indices in increasing order
################################################################################################

function orderloc!(order::Symbol, selected::Vector)
	if order == :randomized
		return shuffle!(selected)
	elseif  order == :sorted
		return sort!(selected)
	elseif order == :original
		selected
	else
		error("order_type was $order, which is not one of [:randomized, :sorted, :natural].")
	end
end



################################################################################################
# Rooulette Wheel selection
################################################################################################

function roulette(values::Vector, count; order = :original)
	selected = sample(1:length(values), weights(values), count)
	orderloc!(order, selected)
end


################################################################################################
# Stochastic Universal Sampling (SUS)
################################################################################################

# A vectorized implementation originally written to be fast in R
#	- makes explicit and segments the two parts of the sampling: stochastic and non-stochastic
#	- behavior is different then baker_sus (Baker's original SUS)  when count < length(pmf), especially when count | length(pmf)
#		- can be seen when count | length(pmf) and pmf has a discrete uniform distribution (selection "alternates")
#			- not a problem if members are arbitrarily placed, but can be a problem if not (e.g. sorted by fitness)
#		-  with this implementation of SUS, these non-random positional dependencies do not occur
#	- actually the same speed (may be slightly slower, but not significantly) but take more space (1/3 more) then baker_sus
function sus(pmf::Vector, count; order = :original)
	expected_count = count * pmf

	# non-stochastic sampling
	base_count = ifloor(expected_count)
	base_sel = rep(base_count)																		# select the members with expected counts >= 1 the base integer # of times

	# stochastic sampling added ontop of the non-stochastic samples
	rest_count = count - length(base_sel)													# compute how many more are needed
	if (rest_count > 0)																						# if we still need more
		rest_pmf = weights(expected_count - base_count)										# re-compute the prob for each member over the base amount olready selected
		domain   = 1:length(pmf)																					# location domain; i.e locations being sampled (support set of the pmf)
		rest_sel = sample(domain, rest_pmf, rest_count, replace = false)	# randomly select the rest without replacement
		selected = [base_sel, rest_sel]
	else
		selected = base_sel
	end

	orderloc!(order, selected)
end

# The original algorithm from Baker
# 	- assumes pmf is normalized (sums to 1)
function baker_sus(pmf::Vector, count; order = :original)
	selected = Array(Integer, count)						# define output array
	stepsize = 1/count													# the sus unit step-size is an equal fraction of the probability (unit circle)
	i = loc = 0																	# i counts the number of selected locations, final value will be the chr_count
	pmf_sum = 0.0																# loc is the index within the probability vector (the pmf) -> supplies the selected locations
	sel_sum = rand() * stepsize 								# sel_sum start at some random perturbation within the first step-size

	while sel_sum <	1.0
		loc += 1
		pmf_sum += pmf[loc]
		while sel_sum < pmf_sum
			i += 1
			selected[i] = loc
			sel_sum += stepsize
		end
	end

	orderloc!(order, selected)
end


# The original algorithm from Baker written to clarify the intent of the algorithm
#	- same speed (probably slower, but not significantly) but take more space (slightly)
# 	- assumes prob is normalized (sums to 1)
function alt1_sus(pmf::Vector, count; order = :original)
	selected = Array(Integer, count)						# define output array
	cdf = cumsum(pmf)										# create the cdf from the pmf
	stepsize = 1/count										# the sus unit step-size is an equal fraction of the probability (unit circle)
	sum = rand() * stepsize 								# sum start at some random perturbation within the first step-size
	i = loc = 1												# i counts the number of selected locations, final value will be the chr_count
															# loc is the index within the probability vector (the pmf) -> supplies the selected locations
	while sum <	1.0
		while sum < cdf[loc]
			selected[i] = loc
			i   += 1
			sum += stepsize
		end
		loc += 1
	end

	orderloc!(order, selected)
end

function sus(fitness::FitnessValues, count)
	sus(normalize(fitness), count)
end

function sus(popn::Population, count)
	sus(popn.fitness)
end
