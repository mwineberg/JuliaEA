##########################################################################
# Types Definition and Constructuors
##########################################################################

immutable  BinaryEncoding <: SimpleEncoding
	chr_length::Integer
	create_chromosome::Function
end

function BinaryEncoding(chr_length::Integer)
	BinaryEncoding(chr_length, random_chromosome)
end

function BinaryEncoding()
	BinaryEncoding(0, random_chromosome)
end

immutable CategoricalEncoding <: SimpleEncoding
	chr_length::Integer
	create_chromosome::Function
	alphabet_size::Integer
end

function CategoricalEncoding(chr_length::Integer, alphabet_size::Integer)
	CategoricalEncoding(chr_length, random_chromosome, alphabet_size)
end

function CategoricalEncoding(alphabet_size::Integer)
	CategoricalEncoding(0, random_chromosome, alphabet_size)
end

immutable  IntegerEncoding <: SimpleEncoding
	chr_length::Integer
	create_chromosome::Function
	range
end

function IntegerEncoding(chr_length::Integer, range)
	IntegerEncoding(chr_length, random_chromosome, range)
end

function IntegerEncoding(range)
	IntegerEncoding(0, random_chromosome, range)
end

immutable RealEncoding <: SimpleEncoding
	chr_length::Integer
	create_chromosome::Function
	max 									# can be an integer (slow because will be converted to a float), float, or vector of chr length
	min 									# can be an integer (slow because will be converted to a float), float, or vector of chr length
	interval 							# can be an integer (slow because will be converted to a float), float, or vector of chr length
end

function RealEncoding(chr_length::Integer, max = 1.0, min = 0.0)
	RealEncoding(chr_length, random_chromosome, max, min, max - min)
end

function RealEncoding(max = 1.0, min = 0.0)
	RealEncoding(0, random_chromosome, max, min, max - min)
end

##########################################################################
# Encoding Methods
##########################################################################

function random_chromosome(encode::BinaryEncoding)
	return rand(0:1, encode.chr_length)
end

function random_chromosome(encode::CategoricalEncoding)
	return rand(0:(encode.alphabet_size - 1), encode.chr_length)
end

# range can be a range e.g. 1:5 or 0:2:10, or it can be an array [0,5,9,3,1]
# future update: weights could be used to bias the sampling by switching to sample(array, weights, ...)
function random_chromosome(encode::IntegerEncoding)
	return rand(encode.range, encode.chr_length)
end

function random_chromosome(encode::RealEncoding)
	return encode.interval .* rand(encode.chr_length) + encode.min
end

##########################################################################
# Decoding Method
##########################################################################
# 	Must return a result of type PhenotypeMember
##########################################################################

function decode(chr::Chromosomes, encoding::SimpleEncoding)
	return NoPhenotype()
end

function decode(popn::Population, encoding::SimpleEncoding)
	return NoPhenotype()
end
