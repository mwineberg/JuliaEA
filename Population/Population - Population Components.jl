####################################################################################################
# Population Component methods  
#	- works on all subtypes of Population Component including:
#		- Chromosome, PhenotypeMembers, FitnessValues, OptimalMembers
####################################################################################################
#	getindex(::PopulationComponent, i)
#		- allowing indexing directly using [] to select phenotype members
#		- i can be a range e.g. 10:20 or 1:3:27, or it can be a vector of indicies
#		- no external type information is retained (i.e. not wrapped in the Phenotype type)
# 	duplicate() select all phenotype(s) i, where i can be a range; returns a Phenotype object
# 		- calls the appropriate getindex() function (i.e implements the [] operator) 
#		- then wraps the result the appropriate Phenotype subtype PT
#	combine() combines two phenotype matrices together using a vertical concatenation 
#	  and constructs a PhenotypeMember instance of the same type that was passed in 
####################################################################################################
#	Note:
#		- ideally all types that inherit from PopulationComponent should have a "members" field
#		- unfortunately, the name "members" only is meaningful for Chromosomes and PhenotypesMembers
#		- for FitnessValues "value" is more meaningful as is "optimial" for OptimalMembers
#		- unfortunately, Julia does not support field name aliases as it does Type aliases
#		- consequently, a more indirect way of accessing the relevant field (the field that holds the values) is used 
#			- the field that holds the values across all members of the population is always the **first** field (and so far the only field)
#			- for any future expansion, the first field of a PopulationComponent must hold the values of the component
####################################################################################################
#	future implementation (once Julia 0.4 is released)
#
#		select(PopulationComponent, i)
#			- would be a "splice" version of duplicate
#			- could increase efficiency, but since not selecting consecutive blocks, 
#				the indexing could require almost as much space as just duplicating, except for Chromosomes
#			- would be used in reduction() part of the main EA routine
####################################################################################################

function getindex(pc::PopulationComponent, i)
	field = names(pc)[1]									# the field name used to store the values for this particular component 
	pc.field[i]
end

function duplicate{PC <: PopulationComponent}(pc::PC, i)
	PC(getindex(pc, i))
end 

function combine{PC <: PopulationComponent}(pc1::C, pc_rest...)
	field = names(c1)[1]									# the field name used to store the values for this particular component
	others = Array(typeof(c1.field), length(pc_rest))		# storage for the field values for the other components (of the same type) being combined
	for i = 1:length(pc_rest)
		others[i] = pc_rest[i].field 						# copy the field values
	end
	PC(vcat(c1.field, others...))							# merge the field values and store in an object of PopulationComponent runtime type
end


