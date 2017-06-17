##########################################################################
#	Base Counts
##########################################################################
#	population size    		{copied from population parameters}
#	post_reproduction size = ρ * (population size - elite size)
#	reduction size 	  	   = population size - elite size
#	elite size 		  	   = ξ * population size
#
#	notes:
#		population size = reduction size + elite size
#		post_reproduction size
#			- is the total number of children that are result from reproduction (before reduction)
#			- this is different from the reproduction size (stored in counts.selection),
#			  which is the number of chromosomes that need to be selected for the reproduction process
#			- is called λ in ES (different from the chromosome length λ)
#		"reduction" is reduced from either
#			- post repoduction size 				{(μ,λ) in ES}
#			- post reproduction size + population size 	{(μ+λ) in ES}
#		The base counts do not change from generation to generation so do not need to be updated
##########################################################################

type BaseCounts
	popn::Integer
	post_reproduction::Integer
	reduction::Integer
	elite::Integer
end

function BaseCounts(ea::EA_Parameters)
	reduction = ea.popn_size - ea.ξ
	post_reproduction = ea.ρ * sc.reduction
	BaseCounts(ea.popn_size, post_reproduction, reduction, ea.ξ)
end

##########################################################################
#	Crossover Counts
##########################################################################
#			max number of xovers  = reproduction size / offspring per xover
#			xover applications    = pχ * max number of xovers   or
#			xover parents count   = number of crossover applications * number parents per crossover
#			xover offspring count = crossover applications * offspring per xover     {i.e. offspring count}
#
#	note:
#		- The counts could changed from generation to generation and so must be updated
#			- to avoid duplication of code, the counts are created with empty fields using new() within an inner constructor
#			- the fields are then updated using the same code that updates the counts every generation
#			- update!() checks if nothing has changed between generations; if so the counts are not updated
#		- As max_appl is calculated using only information from base count and EA parameters, which do not change throughout the run,
#		   the field is computed and stored when the type is created and is not included in update!()
##########################################################################

type CrossoverCounts
	max_appl::Integer
	xover_appl::Integer
	all_parents::Integer
	all_offspring::Integer
	function CrossoverCounts(ea::EA_Parameters, bcounts::BaseCounts)
		max_crosses = max_xovers(ea.offpring_per_xover, bcounts.post_reproduction)
		xover_counts = new(max_crosses)
		update!(xover_counts, ea, changed = true)
		xover_counts
	end
end

function max_xovers(offpring_per_xover::Integer, count::Integer)
	div(count, offpring_per_xover)
end

function update!(xc::CrossoverCounts, ea::EA_Parameters; changed = false)
	if changed || update(ea.pχ)
		xc.xover_appl 	 = modificationcount(ea.pχ, xc.max_appl)
		xc.all_parents   = xc.xover_appl * ea.parents_per_xover
		xc.all_offspring = xc.xover_appl * ea.offpring_per_xover
		changed = true
	end
	return changed
end

##########################################################################
#	Reproduction Partitioning
##########################################################################
#	both 		   = crossover+mutation = pμ_chr * offspring count
#	mutation_only  = pμ * (reproduction size - offspring count)
#	crossover_only = offspring count - crossover+mutation
#	neither 	   = directly_copied = reproduction size - crossover+mutation - crossover_only - mutation_only
#	all_crossed    = crossover+mutation + crossover_only
#	all_mutated    = crossover+mutation + mutation_only
#
#	note:
#		- post reproduction size = total = both + crossover_only + mutation_only + directly_copied
#		- The counts could changed from generation to generation and so must be updated
#			- to avoid duplication of code, the counts are created with empty fields using new() within an inner constructor
#			- the fields are then updated using the same code that updates the counts every generation
#			- update!() checks if nothing has changed between generations; if so the counts are not updated
#		- As total is copied from post_reproduction in base count, which does not change throughout the run,
#		   the field is stored when the type is created and is not included in update!()
##########################################################################

type ReproductionPartitions
	total::Integer 					# total number of children that are result from reproduction (before reduction) = counts.base.post_reproduction
	crossed::Integer 				# number of children that are xover offspring (may or may not be mutated)
	not_crossed::Integer 		# number of children that are not xover offspring (may or may not be mutated)
	mutated::Integer 				# number of children that have had mutation applied (may or may not be xover offspring)
	not_mutated::Integer 		# number of children that have not had mutation applied (may or may not be xover offspring)
	both::Integer 					# number of children with both crossover and mutation applied
	xover_only::Integer 		# number of children with only crossover applied (no mutation)
	mutation_only::Integer  # number of children with only mutation applied (are not offpring from crossover)
	neither::Integer 				# number of direct copies from the previous generation - not (necessarily) elite members

	function ReproductionPartitions(ea::EA_Parameters, bcounts::BaseCounts, xcounts::CrossoverCounts)
		partn_counts = new(bcounts.post_reproduction)
		update!(partn_counts, xcounts, ea, changed = true)
		partn_counts
	end
end

function update!(pn::ReproductionPartitions, xcounts::CrossoverCounts, ea::EA_Parameters; changed = false)
	if changed || update(ea.pμ_chr)
		pn.crossed 		 		= xcounts.all_offspring
		pn.not_crossed	 	= pn.total - pn.crossed
		pn.both 		 			= modificationcount(ea.pμ_xchr, pn.crossed)
		pn.mutation_only	= modificationcount(ea.pμ_nxchr, pn.not_crossed)
		pn.xover_only 	 	= pn.crossed - pn.both
		pn.neither 		 		= pn.not_crossed - pn.mutation_only
		pn.mutated 		 		= pn.both + pn.mutation_only
		pn.not_mutated 		= pn.total - pn.mutated
		changed = true
	end
	return changed
end

##########################################################################
#	Selection Counts
##########################################################################
#	population size    		{copied from base count}
#	reduction size 	  		{copied from base count}
#	elite size 		 		{copied from base count}
#	reproduction size 		= counts.xover.all_parents + counts.partition.mutation_only + counts.partition.neither
#
#	note:
#		The reproduction size could changed from generation to generation and so must be updated
#			- to avoid duplication of code, reproduction size count is created with empty fields using new() within an inner constructor
#			- the field is then updated using the same code that updates the reproduction size counts every generation
#			- update!() checks if nothing has changed between generations; if so the reproduction size is not updated
#		As all other fields are copied from base count, which does not change throughout the run,
#		   these fields are stored when the type is created and is not included in update!()
##########################################################################

type SelectionCounts
	popn::Integer
	reduction::Integer
	elite::Integer
	reproduction::Integer
	replication::Integer
	function SelectionCounts(bcounts::BaseCounts, xcounts::CrossoverCounts, pn::ReproductionPartitions)
		scounts = new(bcounts.popn, bcounts.reduction, bcounts.elite)
		update!(scounts, xcounts, pn, changed = true)
		scounts
	end
end

function update!(sc::SelectionCounts, xcounts::CrossoverCounts, pn::ReproductionPartitions; changed = false)
	if changed
		sc.reproduction = xcounts.all_parents + pn.mutation_only
		sc.replication = pn.neither
	end
	return changed
end

##########################################################################
#	EA Counts
##########################################################################
#	Stores all four types of counts
#		base 	   - basic reproduction counts
#		xover 	   - counts used for xover (number of parents used and offspring produced)
#		partitions - partitioning the reproduction counts between xover and mutation and neither
#		selection  - calculating the number of chromosomes to be selected for reproduction
#					     and copying the elite count and reduction count for easy access
##########################################################################

type EA_Counts
	base::BaseCounts
	xover::CrossoverCounts
	partitions::ReproductionPartitions
	selection::SelectionCounts
end

function EA_Counts(ea::EA_Parameters)
	base 	   = BaseCounts(ea)
	xover	   = CrossoverCounts(ea, base)
	partitions = ReproductionPartitions(ea, base, xover)
	selection =  SelectionCounts(base, xover, partitions, changed = true)
	EA_Counts(base, xover, partitions, selection)
end

function update!(pn::ReproductionPartitions, counts::EA_Counts, ea::EA_Parameters; changed = false)
	update!(pn, counts.xover, ea, changed = changed)
end

function update!(sc::SelectionCounts, counts::EA_Counts; changed = false)
	update!(sc, counts.xover, counts.partitions, changed = changed)
end

function update!(counts::EA_Counts, ea::EA_Parameters; force = false)
	counts_changed = update!(counts.xover, ea, changed = force)
	counts_changed = update!(counts.partitions, counts, ea, changed = counts_changed)
	counts_changed = update!(counts.selection, counts, changed = counts_changed)
	return counts_changed
end
