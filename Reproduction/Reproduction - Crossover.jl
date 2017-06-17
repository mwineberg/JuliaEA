# which location is to be crossed
#	- uniform vs k-pt
#	- implemented as a mask
# How many offspring
#	- one versus two
#		- two offspring only applicable when there is is only two parents
# How many parents
# 	- two versus more than two
# 	- if more than two, sample(1:k) from the k alternative
# Crossover shouldn't have an "inplace!()" version
#	- only make sense for the two parent with two children crossover (when applied to a duplicated subpopulation)
#	- design cost not worth it

abstract Crossover
abstract CrossoverType
abstract CrossoverStructure
abstract CrossoverParents <: CrossoverStructure
abstract CrossoverOffspring <: CrossoverStructure
abstract MultipleOffspring <: CrossoverOffspring

type HomologousCrossover <: Crossover
	xtype::CrossoverType
	parents::CrossoverParents
	offspring::CrossoverOffspring
end

function HomologousCrossover(chr_length::Integer; xover_type = :Uniform, parent_count = 2, offspring_count = 1, args...)
	println(args)
	xtype = (xover_type == :Uniform ? UniformXover(chr_length; args...) :
			 xover_type == :TwoPt	? KptXover(chr_length)				:
			 xover_type == :KPt		? KptXover(chr_length; args...) 	:
			 error("xover_type must be one of {:Uniform, :TwoPt, :KPt}; $xover_type was entered instead"))
	parents = ((parent_count > 2)  ? MultipleParents(parent_count) :
			   (parent_count == 2) ? TwoParents()				   :
			   error("parent_count must be >= 2; $parent_count was entered instead"))
	offspring = ((offspring_count > 1)  ? MultpleOffspring(offspring_count, parent_count, chr_length) :
				 (offspring_count == 1) ? SingleOffspring(chr_length)								   :
				 error("offspring_count must be >= 1; $offspring_count was entered instead"))
	HomologousCrossover(xtype, parents, offspring)
end

type UniformXover <: CrossoverType
	alpha::Real
	length::Integer
	mask::Vector
end

function UniformXover(chr_length::Integer; alpha = 0.25)
	UniformXover(alpha, chr_length, Array(Bool, chr_length))
end

type KptXover <: CrossoverType
	k::Integer
	length::Integer
	mask::Vector
end

function KptXover(chr_length::Integer; k = 2, args...)
	KptXover(k, chr_length, Array(Bool, chr_length))
end

OnePtXover(chr_length::Integer) = KptXover(chr_length, k = 1)

#################################

type SingleOffspring <: CrossoverOffspring
	choice::Matrix
	count::Integer
	length::Integer
end

function SingleOffspring(chr_length::Integer)
	offspring = Array(Integer, 1, chr_length)
	SingleOffspring(offspring, 1, chr_length)
end

type SomeOffspring <: MultipleOffspring
	choice::Matrix
	count::Integer
	length::Integer
	original_order::Vector
	recomb_order::Vector
end

type AllOffspring <: MultipleOffspring
	choice::Matrix
	count::Integer
	length::Integer
	original_order::Vector
	recomb_order::Vector
end

# note: the recomb_order will be shuffled/sub-selected during recomb!()
#	- except when parent.count == 2, where the reverse order ([2,1]) is kept
#	- ocount  means  offspring count
#	- pcount  means  parent count
#	- ocount <= pcount
function MultpleOffspring(ocount::Integer, pcount::Integer, chr_length::Integer)
	if (ocount > pcount) error("Offspring count  must not be greater than parent count. Currently offspring = $ocount > parent = $pcount.") end
	offspring = Array(Integer, ocount, chr_length)
	(ocount == pcount ? AllOffspring(offspring, ocount, chr_length, [1:ocount], [pcount:-1:1])
					  : SomeOffspring(offspring, ocount, chr_length, [1:ocount], [pcount:-1:1]))
end

function getchoices(offspring::CrossoverOffspring)
	offspring.choices
end

#################################
# count::Integer
#		- nnumber of parents used in a single crossover application
# nsets::Integer
#		- number of sets of parent.count in the total number of parents provided for xover i.e. number of crossover applications
#		- each parent.count set produces offspring.count number of offspring
# members::Vector
#		- holds all the parents used for all xover applications
#		- length(members) should equal to nsets * count
#		- actual parent chromosomes are not stored, but rather their location within the current population
# index::Integer
#		- the number of the xovver applications applied so far, used as an index into the set of members used as parents
#		- not set in the constructor, but in setlocations!()

type TwoParents <: CrossoverParents
	count::Integer 						# equals 2 for TwoParents - set in an inner constructor
	nsets::Integer
	index::Integer
	members::Vector
	TwoParents() = new(2)
end

type MultipleParents <: CrossoverParents
	count::Integer 						# how many parent to be crossed; for GA = 2, for ES >= 2 (count != length(members))
	nsets::Integer
	index::Integer
	members::Vector
	MultipleParents(count::Integer) = new(count)
end

function setlocations!(parents::CrossoverParents, new_locations::Vector)
	parents.index = 1
	parents.nsets = div(length(new_locations), parents.count)
	parents.members = new_locations
end

function setlocations!(xover::HomologousCrossover, new_locations::Vector)
	setlocations!(xover.parents, new_locations)
	xover
end

function noparents(parents::CrossoverParents)
	parents.index > parents.nsets
end

function nextparents!(parents::CrossoverParents)
	if noparents(parents) error("There are no parents left - use setlocations!() to add new parents for crossover") end
	pend = parents.index * parents.count
	pbegin = pend - parents.count + 1
	parents.index += 1
	slice(parents.members, pbegin:pend)
end

function parentscount(parents::CrossoverParents)
	parents.count
end

function allparentscount(parents::CrossoverParents)
	size(parents.members, 1)
end

#################################

# function select(sel::Selection, xover::Crossover, fit::Vector)
# end

function createmask!(xover::UniformXover)
	rand!(xover.mask)
	for i = 1:xover.length
		xover.mask[i] = (xover.mask[i] < xover.alpha)
	end
	xover.mask
end

function createmask!(xover::KptXover)
	cross_pts = Array(Integer, xover.k + 1)
	cross_pts[xover.k + 1] = xover.length + 1
	sort!(sample!(2:xover.length, slice(cross_pts, 1:xover.k), replace = false))

	xchoice = false
	xpt = 1

	for i = 1:xover.length
		if i == cross_pts[xpt]
			xchoice = !xchoice
			xpt += 1
		end
		xover.mask[i] = xchoice
	end

	(xover.mask, cross_pts)
end
#################################

function xoverchoice!(offspring::SingleOffspring, parents::TwoParents, cross::Vector)
	for j = 1:offspring.length
		offspring.choice[1,j] = cross[j] ? 2 : 1
	end
	offspring.choice
end

function xoverchoice!(offspring::SingleOffspring, parents::MultipleParents, cross::Vector)
	for j = 1:offspring.length
		offspring.choice[1,j] = cross[j] ? sample(2:parents.count) : 1
	end
	offspring.choice
end

function xoverchoice!(offspring::MultipleOffspring, parents::TwoParents, cross::Vector)
	for j = 1:offspring.length
		offspring.choice[:,j] = (cross[j] ? offspring.recomb_order
										  : offspring.original_order)
	end
	offspring.choice
end

function xoverchoice!(offspring::MultipleOffspring, parents::MultipleParents, cross::Vector)
	for j = 1:offspring.length
		offspring.choice[:,j] = (cross[j] ? recombine!(offspring)
										  : offspring.original_order)
	end
	offspring.choice
end

# recombination at a specific locus in the chromosome
# 	- randomly sample from the parents into the 'k' offspring
# 	- if number of parents equals the number of offspring, perform a random permutation (shuffle)
function recombine!(o::SomeOffspring)
	shuffle!(o.recomb_order)
	mustmove!(o.recomb_order)												# ensure that recombination occurs for every offspring
	slice(o.recomb_order, 1:o.count)
end

function recombine!(o::AllOffspring)
	shuffle!(o.recomb_order)
	mustmove!(o.recomb_order)												# ensure that recombination occurs for every offspring
	o.recomb_order
end

# ensure that the no parent "positional" offspring is keeping the same value as the parent during the recombination
function mustmove!(rorder::Vector)
	pcount = length(rorder)
	top_offset = pcount - 1
	for i = 1:pcount
		if rorder[i] == i
			j = ((sample(1:top_offset) + i - 1) % pcount) + 1		# find a random location among offspring.count not equal to i
			rorder[i] = rorder[j]
			rorder[j] = i
		end
	end
end

#################################

function maxoffspring(x::Crossover)
	x.offspring.count * x.parents.nsets
end

function countcheck(x::Crossover, count::Integer)
	0 < count <= maxoffspring(x)
end

function crossover{T}(x::Crossover, chrs::Matrix{T}, count = maxoffspring(x))
	println("countcheck() = $(countcheck(x,count))")
	if !countcheck(x, count)
		error("count must be between 0 and $(maxoffspring(x)). Instead count = $count.")
	end
	chr_length = chrlength(chrs)
	xover = Array(T, count, chr_length)
	i = 1

	while i <= count
		mask = createmask!(x.xtype)
		parents = nextparents!(x.parents)
		choice = xoverchoice!(x.offspring, x.parents, mask)

		println("mask: \n$mask")
		println("parents: $parents")
		println("choice: \n$choice")

		for j = 1:x.offspring.count
			# choose genes from the parents for offspring j and place them in xover output location i
			for gene = 1:chr_length
				xover[i, gene] = chrs[parents[choice[j,gene]],gene]
			end

			# go to the next xover output location
			i += 1
			if i > count
				break
			end
		end
	end

	xover
end
