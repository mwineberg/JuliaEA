####################################################################################################
#  The SimpleFn library
####################################################################################################
#  	- The libarary will contain all of the basic optimization test suite problems used in EA
#	- Each function implements
#		- the function itself acting on the phenotype/chromosome vector of real values (decoded/not_decoded)
#		- the FnProblem, including Optimum Information (if known), calling parameters,
#			whether the function maximizes or minimizes, and the encoding information
#	- The encode_type can equal
#		- :BCD for binary coded decimal (BinaryEncodedVector)
#		- :Gray for Gray coded decimal (GrayEncodedVector)
#		- :Real for simple encoding using a real valued vector (genotype == phenotype i.e. NoPhenotype)
#	- Note 1:
#		When encode_TYPE = :Real, max and min are used, o.w. dim_size is used
#	- Note 2:
#		when simple real value encodings are used to create random chromosomes,
#		each dimension can either use the same max and min values (max,min are FloatingPoint)
#		or each dimension has its own maximum and minimum values (max,min are of type Vector of chromosome length)
####################################################################################################
# Implemented Functions
# 	- Rastrigin
####################################################################################################

####################################################################################################
#  Rastrigin Fn and Problem setup
####################################################################################################

function rastrigin(phtype::Vector, fs::SimpleFnSetup)
	n = length(phtype)
	rast = 0
	for dim = 1:n
		rast = rast + phtype[dim]^2 - 10 * cos(2 * pi * phtype[dim])
	end
	return 10 * n + rast
end

function RastigrinProblem(encode_type::Symbol, dim_count::Integer; epsilon = 0.0, args...)
	fnsetup = SimpleFnSetup(rastrigin, ())															# rastrigin has no calling parameters
	maximize = false 																										# rastrigin is a minimization problem
	optinfo = SO_OptimumInfo(maximize, 0.0, epsilon)										# rastrigin has an optimal fitness of 0.0
	encoding = SimpleFunctionEncoding(encode_type, dim_count; args...)
	SimpleFnProblem(fnsetup, maximize, encoding, optinfo)
end
