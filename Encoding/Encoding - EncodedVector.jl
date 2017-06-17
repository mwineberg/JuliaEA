##########################################################################
# Types Definition and Constructuors
##########################################################################

immutable BinaryEncodedVector <: EncodedVector
	dim_size::Integer
	dim_count::Integer
	chr_length::Integer
	gene_shift::Real 		# converts natural number (after decoded from binary) to integer/float using linear translation e.g. [0, 1023] becomes [-512, 511]
	gene_scale::Real 		# scales the shifted number; e.g. [-512, 511] becomes [-5.12, 5.11]; scale is often an integer, but can be a float
	create_chromosome::Function
end

function BinaryEncodedVector(dim_count::Integer, dim_size::Integer, gene_shift::Real, gene_scale::Real)
	BinaryEncodedVector(dim_count, dim_size, dim_size * dim_count, gene_shift, gene_scale, random_chromosome)
end

immutable GrayEncodedVector <: EncodedVector
	dim_count::Integer
	dim_size::Integer
	chr_length::Integer
	gene_shift::Real 		# converts natural number (after decoded from binary) to integer e.g. [0, 1023] becomes [-512, 511]
	gene_scale::Real 		# scales integer after shifting to a real; e.g. [-512, 511] becomes [-5.12, 5.11]; scale often an integer, but can be a real
	create_chromosome::Function
end

function GrayEncodedVector(dim_count::Integer, dim_size::Integer, gene_shift::Real, gene_scale::Real)
	GrayEncodedVector(dim_count, dim_size, dim_size * dim_count, gene_shift, gene_scale, random_chromosome)
end

##########################################################################
#  Constructors for EncodedVector;
#		= where appropriate subtype constructor is determined from a parameter
##########################################################################

function CreateEncodedVector(encoding_type::Symbol, dim_size::Integer, dim_count::Integer, gene_shift::Real, gene_scale::Real)
	if(encoding_type == :Gray)
		GrayEncodedVector(dim_size, dim_count, gene_shift, gene_scale)
	elseif (encoding_type == :BCD)
		BinaryEncodedVector(dim_size, dim_count, gene_shift, gene_scale)
	else
		error("Encoding type must be either :Gray or :BCD")
	end
end

##########################################################################
# Encoding Methods
##########################################################################

function random_chromosome(encode::EncodedVector)
	return rand(0:1, encode.chr_length)
end

##########################################################################
# Decoding Display Methods
##########################################################################

function printchr(chr::Union{Matrix, Chromosomes}, e::EncodedVector, i::Integer, dim::Integer)
	gene_start = (dim - 1) * e.dim_size + 1
	gene_end = gene_start + e.dim_size - 1
	print(chr[i, gene_start:gene_end])
end

function printchr(chr::Union{Matrix, Chromosomes}, e::EncodedVector, i::Integer)
	for dim = 1:e.dim_count
		printchr(chr, e, i, dim)
		print("\n")
	end
end

##########################################################################
# Decoding Helper Methods: gray2bin (Gray code to Binary Coded Decimal)
##########################################################################

# conversion of a single gene in a chromosome vector from Gray to BCD
function gray2binary(gray_chrs::Matrix, e::BinaryEncodedVector, i::Integer, dim::Integer)
	bin_chrs = Array(Integer, e.dim_size)
	gk = gray_end = dim * e.dim_size
	bk = bin_end = e.dim_size
	bin_chrs[bin_end] = gray_chrs[i, gray_end]
	while bk > 1
		bk = bk - 1; gk = gk - 1
		bin_chrs[bk] = ifelse(gray_chrs[i, gk] != bin_chrs[bk+1], 1, 0)
	end
	bin_chrs
end

# conversion of a single chromosome in the population from Gray to BCD
function gray2binary(gray_chrs::Matrix, e::EncodedVector, i::Integer)
	bin_chrs = Array(Integer, e.chr_length)
	for dim = 1:e.dim_count
		bin_end     = dim * e.dim_size
		bin_start   = bin_end - e.dim_size + 1
		bin_chrs[bin_start:bin_end] = gray2binary(gray_chrs, e, i, dim)
	end
	return bin_chrs
end

# conversion of the entire population from Gray to BCD
# integrated as one function
#	- for speed
#	- to store as a matrix, not a vector of vectors
function gray2binary(gray_chrs::Matrix, e::EncodedVector)
	popn_size = popnsize(gray_chrs)
	bins_chrs = similar(gray_chrs)
	for dim = 1:e.dim_count
		gene_end   = dim * e.dim_size
		gene_start = gene_end - e.dim_size + 1

		for i = 1:popn_size
			bins_chrs[i, gene_end] = gray_chrs[i, gene_end]

			for k = (gene_end - 1):-1:gene_start
#				print("i = ", i, ", k = ", k, "\n")
				bins_chrs[i, k] = ifelse(gray_chrs[i, k] != bins_chrs[i, k+1], 1, 0)
			end
		end
	end
	return bins_chrs
end

##########################################################################
# Decoding Helper Methods: binary2decimal (BCD to Integer)
##########################################################################

# conversion of a single gene in a chromosome vector from binary to decimal
function binary2decimal(bin_chrs::Matrix, e::BinaryEncodedVector, i::Integer, dim::Integer)
	gene_value = 0
	gene_start = (dim-1) * e.dim_size + 1
	for k = 0:(e.dim_size - 1)
		 gene_value = gene_value + bin_chrs[i, k + gene_start] * 2^k
	end
	return gene_value
end

# conversion of a single chromosome in the population from binary to decimal
function binary2decimal(bin_chrs::Matrix, e::EncodedVector, i::Integer)
	dec_chr = Array(Integer, e.dim_count)
	for dim = 1:e.dim_count
		dec_chr[dim] = binary2decimal(bin_chrs, e, i, dim)
	end
	return dec_chr
end

# conversion of the entire population from binary to decimal
# integrated as one function
#	- for speed
#	- to store as a matrix, not a vector of vectors
function binary2decimal(bin_chrs::Matrix, e::EncodedVector)
	popn_size = popnsize(bin_chrs)
	dec_chr = Array(Integer, popn_size, e.dim_count)
	for dim = 1:e.dim_count
		gene_start = (dim - 1) * e.dim_size + 1

		for i = 1:popn_size
			gene_value = 0

			for k = 0:(e.dim_size - 1)
				 gene_value = gene_value + bin_chrs[i, k + gene_start] * 2^k
			end

			dec_chr[i, dim] = gene_value
		end
	end
	return dec_chr
end

##########################################################################
# Decoding Helper Methods: binary2decimal
##########################################################################

function decimal2float(dec_chrs::Matrix, e::EncodedVector)
	popn_size = size(dec_chrs, 1)
	float_chrs = Array(FloatingPoint, popn_size, e.dim_count)

	for i = 1:popn_size, dim = 1:e.dim_count
		float_chrs[i, dim] = (dec_chrs[i, dim] - e.gene_shift) / e.gene_scale
	end

	return float_chrs
end

##########################################################################
# Decoding Methods
##########################################################################
# 	Must return a result of type PhenotypeMember
##########################################################################

function decode(chr::Chromosomes, encoding::BinaryEncodedVector)
	decimal_chr = binary2decimal(chr.members, encoding)
	float_chr = decimal2float(decimal_chr, encoding)
	return GeneralPhenotype(float_chr)
end

function decode(chr::Chromosomes, encoding::GrayEncodedVector)
	binary_chr = gray2binary(chr.members, encoding)
	decimal_chr = binary2decimal(binary_chr, encoding)
	float_chr = decimal2float(decimal_chr, encoding)
	return GeneralPhenotype(float_chr)
end

function decode(popn::Population, encoding::EncodedVector)
	return decode(popn.chromosomes, encoding)
end
