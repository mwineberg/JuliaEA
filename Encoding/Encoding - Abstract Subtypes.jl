abstract FixedLengthEncoding <: Encoding
abstract VarLengthEncoding <: Encoding
abstract SimpleEncoding <: FixedLengthEncoding
abstract EncodedVector <: FixedLengthEncoding


##########################################################################
# Common Accessor Methods
##########################################################################

function chr_length(se::FixedLengthEncoding)
	se.chr_length
end

function chr_length!(se::FixedLengthEncoding, newLength::Integer)
	se.chr_length = newLength
end

function create_chromosome(se::Encoding)
	se.create_chromosome(se)
end

# used when creating a population in Population()
function chromosomes_constructor(encoding::FixedLengthEncoding)
	FixedLengthChromosomes
end

# used when creating a population in Population()
function chromosomes_constructor(encoding::VarLengthEncoding)
	VarLengthChromosomes
end
