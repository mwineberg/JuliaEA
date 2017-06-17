##########################################################################
# Should turn this into Model View Controller style system
# The Population should evolve by verious step amounts
# The state of the population should be observable so different methods can be called
#		during different states (e.g. selection cannot be done after reproduction but before evaluation)
# Different methods would perform diffeernt sequences of the evolutionary process on the population
# This will also allow statstical queries to more easily be asked
# Different controllers could be developed to run different EA systems
# This will also allow the development of mechansims for interacting populations
##########################################################################

function generational_EA(ea::EA_Parameters, popn_info::PopulationInfo, problem::Problem, stats::EAStats,
						 sel::SelectionType, xover::Xover, mutation::Mutation; verbose = false)
	# variables added: use.goal
	# functions added: optimum.found()
	ea_setup(ea.env)		# includes optimum.found(), update.stats(), and selection.count
	ea_counts = EA_Counts(ea)
	counts = counts.selection

	cp = current_popn = create_population(popn_info, problem.encoding)
	evaluate!(current_popn, problem)

	update!(stats, current_popn)
	display_console_report(0, stats, current_popn, verbose)

	for gen = 1:controls.max_gen
		if goodenough(current_popn) break end

		# perform elitism
		elite_loc  = selectlocations(sel.elite, cp.fitness, counts.elite)		 # select - produce locations
		elite_popn = duplicate(current_popn, elite_loc)							 # create a new popn using copies of selected members for all fields
																				 # no decoding or evaluation needed

		# reproduce with mutation (both xover and random mutations)
		repr_loc = selectlocations(sel.regular, cp.fitness, counts.reproduction) # select	- produce locations
		repr_popn = reproduce(current_popn, repr_loc, ea_counts)			 	 # reproduce - produces chromosome members only
		evaluate!(repr_popn, problem)											 # evaluate  - sets fitness (and optimum_found if applicable)

		# replication with no mutation (neither xover nor random mutations)
		repl_loc = selectlocations(sel.regular, cp.fitness, counts.replication)  # select	- produce locations
		repl_popn = duplicate(current_popn, repr_loc)			    			 # reproduce - produces chromosome members only
																				 # no decoding or evaluation needed

		# combine reproduction and replicated populations together in preparation for reduction selection
		rp = rep_popn = (ea.include_parents ? combine(repr_popn, repl_popn)
											: combine(repr_popn, repl_popn, current_popn))

		# perform reduction selection
		reduced_loc  = selectlocations(sel.reduce, rp.fitness, counts.reduction) # select - produce locations
		reduced_popn = duplicate(rep_popn, reduced_loc)							 # create a new popn using copies of selected members for all fields
																				 # no decoding or evaluation needed

		# change current population to the evolved population and update stats
		current_popn = combine(elite_popn, reduced_popn)

		# update!(stats, current_popn)
		# display_console_report(gen, stats, current_popn, verbose)
		display_console_report(gen, current_popn, verbose)
		update!(counts, ea)
	end

	return (stats, current_popn)
end
