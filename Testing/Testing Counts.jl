#	popn_size 	= 100
#	ρ			= 1.0		{ratio of chromosome count after vs before reproduction} 
#	ξ			= 1			{number of elites}
#	λ 			= 32		{chromosome length; must be > 1 or -1 if variable length}
#	pμ_chr 		= 0.75		{percent/prob that a chromosome will be mutated}
#	pμ_gene		= 0.05		{percent/prob that a given gene within a chromosome will be mutated}
#	pχ 			= 0.75		{percent/prob of the max number of possible crossovers (AKA the probability of crossover)}
# 	parents_per_xover = 2 	{number of parents used per individual crossover application}
#	offpring_per_xover = 1	{number of offpring produced per individual crossover application}

ea0 = EA_Parameters()
	# EA_Parameters(100,1.0,1,32,PercentageParm(0.75),PercentageParm(0.05),PercentageParm(0.75),2,1)
ea1 = EA_Parameters(popn_size = 200, ρ = 2, ξ = 2, λ = 20, pμ_chr = 0.8, pμ_gene = 0.1, pχ = 0.7, parents_per_xover = 3, offpring_per_xover = 2)
	# EA_Parameters(200,2,2,20,PercentageParm(0.8),PercentageParm(0.1),PercentageParm(0.7),3,2)
ea3 = EA_Parameters(popn_size = 200, ρ = 2, ξ = 2, λ = 20, pμ_chr = (0.8, :prob),
					pμ_gene = (0.1, :percent), pχ = (0.7, :prob), parents_per_xover = 3, offpring_per_xover = 2)
	# EA_Parameters(200,2,2,20,ProbabilityParm(0.8),PercentageParm(0.1),ProbabilityParm(0.7),3,2)
ea4 = EA_Parameters(popn_size = 200, ρ = 2, ξ = 2, λ = 20, pμ_chr = (0.8, :prob),
					pμ_gene = (0.1, :percentage), pχ = (0.7, :probability), parents_per_xover = 3, offpring_per_xover = 2)
	# EA_Parameters(200,2,2,20,ProbabilityParm(0.8),PercentageParm(0.1),ProbabilityParm(0.7),3,2)
EA_Parameters(popn_size = 200, ρ = 2, ξ = 2, λ = 20, pμ_chr = (0.8, :prob),
					pμ_gene = (0.1, :pcent), pχ = (0.7, :pr), parents_per_xover = 3, offpring_per_xover = 2)
	# ERROR: :pcent is not recognized in pμ_gene/prob_gene_mutated = (0.1,:pcent). The parameter must be set to either :prob or :percent
EA_Parameters(popn_size = 200, ρ = 2, ξ = 2, λ = 20, pμ_chr = (0.8, :prob),
					pμ_gene = (0.1, :percent), pχ = (0.7, :pr), parents_per_xover = 3, offpring_per_xover = 2)
	# ERROR: :pr is not recognized in pχ/prob_xover = (0.7,:pr). The parameter must be set to either :prob or :percent

# type VariationPartitions
#	total::Integer 			# total = both + xover_only + mutation_only + neither = counts.reproduction.variation = 99
#	crossed::Integer 		# crossed = both + crossover_only = counts.xover.all_offspring = 99 * 0.75 = ifloor(74.25) = 74
#	mutated::Integer 		# mutated = both + mutation_only
#	both::Integer 			# crossed * 0.75 = 74 * 0.75 = ifloor(55.5) = 55
#	xover_only::Integer 	# crossed - both = 74 - 55 = 19
#	mutation_only::Integer  # (total - crossed) * 0.75 = (99 - 74) * 0.75 = ifloor(18.75) = 18
#	neither::Integer 		# total - crossed - mutation_only = 99 - 74 - 18 = 7

ea_counts0 = EA_Counts(ea0)
	# EA_Counts(ReproductionCounts(100,99,99,1),VariationPartitions(99,74,73,55,19,18,7),CrossoverCounts(74,148,74))

# type VariationPartitions
#	total::Integer 			# total = both + xover_only + mutation_only + neither = counts.reproduction.variation = 396
#	crossed::Integer 		# ifloor(div(total, offspring_per_chr) * 0.7) * offspring_per_chr = ifloor(div(396,2) * 0.7) * 2 = 276
#	mutated::Integer 		# mutated = both + mutation_only = 220 + 96 = 316
#	both::Integer 			# crossed * 0.8 = 276 * 0.8 = ifloor(220.8) = 220
#	xover_only::Integer 	# crossed - both = 276 - 220 = 56
#	mutation_only::Integer  # (total - crossed) * 0.8 = (396 - 276) * 0.8 = ifloor(96.0) = 96
#	neither::Integer 		# total - crossed - mutation_only = 396 - 276 - 96 = 24

ea_counts1 = EA_Counts(ea1)
	# EA_Counts(ReproductionCounts(200,396,198,2),VariationPartitions(396,276,316,220,56,96,24),CrossoverCounts(138,414,276))
update!(ea_counts1, ea1)		# false
ea_counts1
	# EA_Counts(ReproductionCounts(200,396,198,2),VariationPartitions(396,276,316,220,56,96,24),CrossoverCounts(138,414,276))

ea_counts3 = EA_Counts(ea3)
update!(ea_counts3, ea3)		# true
ea_counts3