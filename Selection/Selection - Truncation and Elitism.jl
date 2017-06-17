################################################################################################
#  Elite / Truncation Selection
################################################################################################
#	- Elitism tends to be used in GA systems while truncation selection is typically used in ES systems
#	- The mechanism is the same. Differences are in how and when they are used:
#			- few members for elitism, frequently only 1 or 2 and the members are passed unchanged into the next populaiton
#			- many more can be selected in an ES system based on the mu parameter set by the user
#			  and used during post-selection after uniform inflation and mutation has been applied
################################################################################################

import StatsBase.tiedrank
import StatsBase.competerank

################################################################################################
# Helper Methods - vector ranking with differnt tie ranking routines
################################################################################################
#  competerank is [1, 2, 2, 4] ranking
#  modcompeterank is [1, 3, 3, 4] ranking
#	 - does not exist in Julia so must be written
#	 - slower than competerank but still fast
#  shufflerank uses [1,2,3,4] ranking but shuffles the ranks of the ties
#	 - does not exist in Julia so must be written
#	 - slower than modcompeterank but still fast
#  tiedrank is [1, 2.5, 2.5, 4] ranking
#	 - built-in
#	 - used to compute modcompeterank
#	 - not directly used for Trunction selection
################################################################################################


# uses [1, 2, 2, 4] ranking and [1, 2.5, 2.5, 4] ranking to create [1, 3, 3, 4] ranking
function modcompeterank(v::Vector, cr = competerank(v))
	tr = tiedrank(v)
	ifloor(2 * tr - cr)
end

function shufflerank(v::Vector, cr = competerank(v))
	ties = Array(Vector,length(cr))

	# initialize tie buckets
	for i = 1:length(ties)
		ties[i] = Array(Integer, 0)
	end

	# peform bucket sort
	for (loc,rnk) in enumerate(cr)
		push!(ties[rnk], loc)
	end

	rand_ties = cr

	for i = 1:length(ties)						# iterate through the ties list
		if length(ties[i]) > 1						# if there are any ties for that particular rank value
			shuffle!(ties[i])							# shuffle the order of the locations of the ties
			for j = 1:length(ties[i])					# iterate through the ties
				rand_ties[ties[i][j]] += j - 1				# find the (shuffled) tie locations in the original list of ranks, re-ranking as you go
			end
		end
	end

	rand_ties
end

################################################################################################
# Tie handling routines for both elitism and trunctation selection
################################################################################################
#	- Uses tie handling for ranking above
#	- if the cut-off falls between non-tied members or between tie sets
#		 the tie breaking method used doies not matter
#		 as all members above the cut-off are selected and all members below are not
#	- if the cut-off falls inside one of the tie set
#		 all ties not in that tie set are unaffected as before
#		 for the tie set with the cut-off, three choices are offered:
#		 1) random	- the ties are randomly assigned ranks within the tie group
#					- the number of selected members is exactly equal to the requested count
#					- default (expecially useful for fixed population size)
#		 2) all: 	- all members of that tied set are selected
#					- the number of selected members will be greater than the requested count
#		 3) none: 	- none of the members of that tied set are selected
#					- the number of selected members will be less than the requested count
#
#	note: 	- elitism and truncation selection utilize a ranking function with a different
#			  tie handling than RS, which uses a [1,2.5,2.5,4] tied ranking method
#			- so saving the ranks from RS would not be helpful and the ranks need to be regenerated
################################################################################################

###################################
# Type Definitions and Constructors
###################################

abstract FinalTies

type AllFinalTies 	 <: FinalTies
	maximizing::Bool
end

type NoFinalTies 	 <: FinalTies
	maximizing::Bool
end

type RandomFinalTies <: FinalTies
	maximizing::Bool
end

AllFinalTies(final_ties::FinalTies) 		= AllFinalTies(final_ties.maximizing)
NoFinalTies(final_ties::FinalTies) 			= NoFinalTies(final_ties.maximizing)
RandomFinalTies(final_ties::FinalTies) 	= RandomFinalTies(final_ties.maximizing)

function FinalTiesChoice(choice::Symbol, maximizing::Bool)
	((choice == :all)    ? AllFinalTies(maximizing)    :
	 (choice == :none)   ? NoFinalTies(maximizing)     :
	 (choice == :random) ? RandomFinalTies(maximizing) :
	 error(""))
end

###################################
# Ranking Methods
###################################

function ranked(info::AllFinalTies, fit::Vector)
	info.maximizing ? length(fit) - modcompeterank(fit) + 1 : competerank(fit)
end

function ranked(info::NoFinalTies, fit::Vector)
	info.maximizing ? length(fit) - competerank(fit) + 1 : modcompeterank(fit)
end

function ranked(info::RandomFinalTies, fit::Vector)
	info.maximizing ? length(fit) - shufflerank(fit) + 1 : shufflerank(fit)
end



################################################################################################
#  Elite / Truncation Selection type definitions and constructors
################################################################################################

###################################
# Type Definitions and Constructors
###################################

type TruncationSelection <: SelectionType
	order::Symbol
	maximizing::Bool
	final_ties::FinalTies
end

typealias Elitism TruncationSelection

function TruncationSelection(;maximizing = true, ties = :random, order = :original)
	TruncationSelection(order, maximizing, FinalTiesChoice(ties, maximizing))
end

#######################################
#  Elite / Truncation Selection Method
#######################################
# The basis idea of the algorithm
#	1) perform sortperm (find the indices that would sort the fitnesses if applied) based on fitness
#	2) select indicies from 1 to popn_size based on the result
#	3) truncate
# ranks are used instead of sortperm to provide more control
#######################################

function selectlocations(info::TruncationSelection, fit::Vector, count::Integer)
	# find the ranks based on fitness (diferent routines would be used based on tied methods and whether we are maximizing or minimizing)
	ranks = ranked(info.final_ties, fit)

	# set the ranks greater than count to 0 (to be eliminated in the next step)
	for i = 1:length(ranks)
		if(ranks[i] > count)
			ranks[i] = 0
		end
	end

	# find() returns the indices of all non-zero elements and since we set any rank greater than count to 0, only the correct indices are returned
	# these locations are then selected from the indices within the population - which is what we want
	popn_size = length(fit)
	loc = [1:popn_size][find(ranks)]
	orderloc!(info, loc)
end
