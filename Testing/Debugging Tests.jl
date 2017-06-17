# Testing Encoding - SimpleGenes
be = BinaryEncoding(12)
rc = random_chromosome(be)
ce = CategoricalEncoding(12, 4)
rc = random_chromosome(ce)
ie = IntegerEncoding(12, 0:10:90)
rc = random_chromosome(ie)
re = RealEncoding(12)
rc = random_chromosome(re)

# Testing FixedLengthChromosomes
bfc = FixedLengthChromosomes(24, be)
cfc = FixedLengthChromosomes(24, ce)
ifc = FixedLengthChromosomes(24, ie)
rfc = FixedLengthChromosomes(24, re)
popnsize(cfc)
bf1 = FixedLengthChromosomes(24, be)
bf2 = FixedLengthChromosomes(24, be)
bf3 = combine(bf1,bf2)
bf3[1:24]
bf3[25:48]
duplicate(bf3, 25:32)
duplicate(rfc, 8:18)
ie1 = FixedLengthChromosomes(24, ie)
ie2 = FixedLengthChromosomes(24, ie)
ie3 = combine(ie1, ie2)
popnsize(ie1)
popnsize(ie3)
chrlength(ie3)
chrlength(cfc)

# Testing GeneralChromosomes
bgc = GeneralChromosomes(24, be)
print(bgc)
cgc = GeneralChromosomes(24, ce)
print(cgc)
igc = GeneralChromosomes(24, ie)
print(igc)
rgc = GeneralChromosomes(24, re)
print(rgc)
popnsize(cgc)
chrlength(igc)	  						# => error
chrlength(cgc)  						# => error	
bg1 = GeneralChromosomes(24, be)
bg2 = GeneralChromosomes(24, be)
bg3 = combine(bg1, bg2)
print(bg1)
print(bg2)
print(bg3)
bg4 = combine(bf1, bg2)  				# => error
bg3[1:24]
bg3[25:48]
bg4 = duplicate(bg3, 25:48)
print(bg4)
popnsize(bg4)
ig1 = GeneralChromosomes(24, ie)
ig2 = GeneralChromosomes(24, ie)
ig3 = combine(ig1, ig2)
print(ig1)
print(ig2)
print(ig3)


# Testing create_chromosome
# haven't created a variable length chromosome encoding - test that side later
bf5 = create_chromosomes(24, be)
cf5 = create_chromosomes(24, ce)
if5 = create_chromosomes(24, ie)
rf5 = create_chromosomes(24, re)
pi_be = PopulationInfo(24)
bf6 = create_chromosomes(pi_be, be)
cf6 = create_chromosomes(pi_be, ce)
if6 = create_chromosomes(pi_be, ie)
rf6 = create_chromosomes(pi_be, re)


# Testing PhenotypeMembers
npt = NoPhenotype()
gpt = GeneralPhenotype(rand(18, 9))
popnsize(gpt)
print(gpt)
usesphenotype(npt)
usesphenotype(gpt)
npt[3]
gpt[2:6]
gpt2 = duplicate(gpt, [2:6])
print(gpt2)
gpt == gpt 								# => true
gpt == gpt2 							# => false
gp1 = GeneralPhenotype(rand(4, 6))
gp2 = GeneralPhenotype(rand(4, 6))
gp3 = combine(gp1, gp2)
print(gp1)
print(gp2)
print(gp3)
np1 = duplicate(npt, [2:6])				# => NoPhenotype
np2 = combine(npt, np1)					# => NoPhenotype

# Testing OptimumInfo
SO_OptimumInfo(KnownOptimum, false, 7)		# => SO_OptimumInfo(KnownOptimum,false,7,0)
SO_OptimumInfo(false, 7.3)					# => SO_OptimumInfo(KnownOptimum,false,7.3,0.0)
SO_OptimumInfo(false, 7.3, 0.2)				# => SO_OptimumInfo(KnownOptimum,false,7.3,0.2)
optknown(SO_OptimumInfo(false, 7.3))		# => true
optunknown(SO_OptimumInfo(false, 7.3))		# => false
optunknown(NoOptimumInfo())					# => true
optknown(NoOptimumInfo())					# => false

# Testing Fitness
optinfo = SO_OptimumInfo(false, 0.0, 0.01)
nooptinfo = NoOptimumInfo()
fv = SO_FitnessValues(rand(32))
print(fv)
isopt = isoptimal(fv.values, optinfo)
isopt = isoptimal(fv, optinfo)
isntopt = isoptimal(fv.values, nooptinfo)
isntopt = isoptimal(fv, nooptinfo)
fv[3:5]
fv1 = duplicate(fv, [1,3,6])
print(fv1)
fv2 = combine(fv1, fv1)
print(fv2)

# Testing OptimalMembers
optinfo = SO_OptimumInfo(false, 0.0, 0.05)
fv = SO_FitnessValues(rand(80))
print(fv)
opt = isoptimal(fv, optinfo)
print(opt)
success(opt)
opt[4:8]
opt1 = duplicate(opt, 4:8)			# choose range so no successes
success(opt1)
opt2 = duplicate(opt, 15:20)		# choose range so there are successes
success(opt2)
opt3 = combine(opt1, opt1)
success(opt3)
opt4 = combine(opt1, opt2)
success(opt4)
opt5 = combine(opt2, opt2)
success(opt5)
optknown(opt3)
optunknown(opt3)
noopt = UnknownOptimum()
optknown(noopt)
optunknown(noopt)
success(noopt)
noopt[1:3]
duplicate(noopt, 1:3)
combine(noopt,noopt)

# Testing Problem - OneMax
popn_info = PopulationInfo(24)							# popn_size  = 24
om_prob = OneMaxProblem(12)								# chr_length = 12
bfc = create_chromosomes(popn_info, om_prob.encoding)
print(bfc)
evaluate(bfc, om_prob.eval_args)

# Testing Population
popn_info = PopulationInfo(24)							# popn_size  = 24
om_prob = OneMaxProblem(12)								# chr_length = 12
p1 = create_population(popn_info, om_prob.encoding)
popnsize(p1)
evaluate!(p1, om_prob)
p1
p1[4:8]
p2 = duplicate(p1, 6:12)
p3 = duplicate(p1, 8:20)
p4 = combine(p2, p3)

# Testing Problem - Rastrigin
popn_info = PopulationInfo(24)							# popn_size  = 24
rastprob = RastigrinProblem(:Real, 10)
p1 = create_population(popn_info, rastprob.encoding)
p1 
evaluate!(p1, rastprob)
p1 

# Testing binary2decimal
popn_info = PopulationInfo(24)							# popn_size  = 24
rastprob = RastigrinProblem(:BCD, 10)
p2 = create_population(popn_info, rastprob.encoding)
printchr(p2.chromosomes, rastprob.encoding, 1)
binary2decimal(p2.chromosomes.members, rastprob.encoding, 1, 1)
binary2decimal(p2.chromosomes.members, rastprob.encoding, 1, 2)
binary2decimal(p2.chromosomes.members, rastprob.encoding, 1)
printchr(p2.chromosomes, rastprob.encoding, 2)
binary2decimal(p2.chromosomes.members, rastprob.encoding, 2)
printchr(p2.chromosomes, rastprob.encoding, 10)
binary2decimal(p2.chromosomes.members, rastprob.encoding, 10)
dp2 = binary2decimal(p2.chromosomes.members, rastprob.encoding)
decimal2float(dp2, rastprob.encoding)

# Testing gray2bin
popn_info = PopulationInfo(24)							# popn_size  = 24
rastprob = RastigrinProblem(:BCD, 10)
p2 = create_population(popn_info, rastprob.encoding)
printchr(p2.chromosomes, rastprob.encoding, 1)
gray2bin(p2.chromosomes.members, rastprob.encoding, 1, 1)
gray2bin(p2.chromosomes.members, rastprob.encoding, 1, 2)
gray2bin(p2.chromosomes.members, rastprob.encoding, 1)
gray2bin(p2.chromosomes.members, rastprob.encoding, 10)
gb2 = gray2bin(p2.chromosomes.members, rastprob.encoding)

# Testing evaluation using BCD decoding 
popn_info = PopulationInfo(24)							# popn_size  = 24
rastprob = RastigrinProblem(:BCD, 10)
p2 = create_population(popn_info, rastprob.encoding)
p2 
evaluate!(p2, rastprob)
p2 

# Testing evaluation using Gray decoding 
popn_info = PopulationInfo(24)							# popn_size  = 24
rastprob = RastigrinProblem(:Gray, 10)
p3 = create_population(popn_info, rastprob.encoding)
p3 
evaluate!(p3, rastprob)
p3 