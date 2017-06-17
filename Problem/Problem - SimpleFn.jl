####################################################################################################
# Type declarations
####################################################################################################

immutable SimpleFnSetup <: EvalSetup
	fitfn::Function
	args::Tuple
end

immutable SimpleFnProblem <: Problem
	setup :: SimpleFnSetup
	maximize::Bool
	encoding::Encoding
	optinfo::OptimumInfo
end

# Virtual Type SimpleFunctionEncoding: Constructor (combination of RealEncoding, BinaryEncodedVector, GrayEncodedVector)
#	- :BCD for binary coded decimal (BinaryEncodedVector) 
#	- :Gray for Gray coded decimal (GrayEncodedVector)
#	- :Real for simple encoding using a real valued vector (genotype == phenotype i.e. NoPhenotype)

# default shift: the domain is transformed from [0, x] where x is an integer, to [-x/2, x/2 - 1]
# default scale: the domain is transformed from [-x, x] where x is an integer, to [-y.z, y.z] where xy is a single digit from 1..9
# traditionally dim_size = 10 for most binary vector encodings, so [0, 1023] -> [-512, 511] -> [-5.12, 5.11]
function SimpleFunctionEncoding(encode_type::Symbol, dim_count::Integer; 
								max = 5.12, min = -5.12, dim_size = 10, 
								gene_shift = 2^(dim_size - 1), 					# note: x/2 == 2^dim_size/2 == 2^(dim_size - 1) which eliminates the division
								gene_scale = 10^ifloor(log10(gene_shift)))
	if encode_type == :Real
		return RealEncoding(dim_count, max, min)
	else
		return CreateEncodedVector(encode_type, dim_count, dim_size, gene_shift, gene_scale)
	end
end


####################################################################################################
#  Evaluation of the GeneralPhentype
####################################################################################################

function evaluate(chr::FixedLengthChromosomes, fs::SimpleFnSetup)
	fitfn = fs.fitfn
	popn_size = popnsize(chr)
	chr_length = chrlength(chr)
	fitness = Array(FloatingPoint, popn_size)
	for i = 1:popn_size
		soln = reshape(chr.members[i,:], chr_length)
		fitness[i] = fitfn(soln, fs)
	end
	return SO_FitnessValues(fitness)								# <- need to return appropriate FitnessValues type
end


function evaluate(solutions::Vector{Vector}, fs::SimpleFnSetup)
	fitfn = fs.fitfn
	fitness = Array(FloatingPoint, popnsize(solutions))
	for (i, soln) in enumerate(solutions)
		fitness[i] = fitfn(soln, fs)
	end

	return SO_FitnessValues(fitness)								# <- need to return appropriate FitnessValues type
end

function evaluate(chr::GeneralChromosomes, fs::SimpleFnSetup)
	evaluate(chr.members, fs)
end

function evaluate(phtype::GeneralPhenotype, fs::SimpleFnSetup)
	evaluate(phtype.members, fs)
end
