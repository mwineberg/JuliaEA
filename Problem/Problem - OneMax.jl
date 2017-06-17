####################################################################################################
# Type declarations
####################################################################################################

immutable OneMaxSetup <: EvalSetup
end

immutable OneMaxProblem <: Problem
	setup::EvalSetup
	maximize::Bool 					# true if finding the maximum; false if finding the minimum - also stored in OptimumInfo inside Optimum inside Population
	encoding::BinaryEncoding
	optinfo::OptimumInfo
end

####################################################################################################
#  Problem Constructor 
####################################################################################################

function OneMaxProblem(chr_length::Integer, epsilon = 0)
	maximize  = true 
	encoding  = BinaryEncoding(chr_length)
	setup = OneMaxSetup()
	optinfo   = SO_OptimumInfo(maximize, chr_length, epsilon)
	return OneMaxProblem(setup, maximize, encoding, optinfo)
end

####################################################################################################
#  Evaluation of the Chromosome (since no phenotype is used)
####################################################################################################

function evaluate(chr::FixedLengthChromosomes, setup::OneMaxSetup)
	popn_size = popnsize(chr)
	chr_length = chrlength(chr)
	fitness = Array(Integer, popn_size)
	for i = 1:popn_size
		total = 0
		for locus = 1:chr_length
			total = total + chr[i,locus]
		end
		fitness[i] = total
	end

	return SO_FitnessValues(fitness)								# <- need to return appropriate FitnessValues type
end
