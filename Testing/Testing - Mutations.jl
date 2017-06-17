# mutation testing

mpc1 = MutnPerChr(0.1, 68)							# => CountMutn(6,68)
mpc2 = MutnPerChr(0.1, 68, ptype = :percentage)		# => CountMutn(6,68)
mpc3 = MutnPerChr(0.1, 68, ptype = :probability)	# => BinomialMutn(Binomial(n=68, p=0.1),68)
mpc4 = MutnPerChr(.1)								# => PercentageMutn(.1)
mpc5 = MutnPerChr(0.1, ptype = :percentage)			# => PercentageMutn(.1)
mpc6 = MutnPerChr(0.1, ptype = :probability)		# => ProbabilityMutn(.1)
mpc7 = MutnPerChr(4, 28)							# => CountMutn(4, 28)
mpc8 = MutnPerChr(4, 28, ctype = :count)			# => CountMutn(4, 28)
mpc9 = MutnPerChr(4, 28, ctype = :probability)		# => BinomialMutn(Binomial(n=28, p=0.14285714285714285), 28)
mpc10 = MutnPerChr(4)								# no matching
mpc11 = MutnPerChr(4, ptype = :percentage)			# no matching
mpc12 = MutnPerChr(4, ptype = :probability)			# no matching

count(mpc1)
count(mpc2)
count(mpc3)
count(mpc4, 100)
count(mpc5, 100)
count(mpc6, 100)
count(mpc7)
count(mpc8)
count(mpc9)

randomloc(mpc1)
randomloc(mpc2)
randomloc(mpc3)
randomloc(mpc4, 100)
randomloc(mpc5, 100)
randomloc(mpc6, 100)
randomloc(mpc7)
randomloc(mpc8)
randomloc(mpc9)

import StatsBase.sample
import StatsBase.weights

function popnsize(chromosomes::Array)
	size(chromosomes, 1)
end

function chrlength(chr::Matrix)
	size(chr, 2)
end

function chrlength(chr::Vector)
	popn_size = popnsize(chr)
	chr_length = Array(Integer, popn_size)
	for i = 1:popn_size
		chr_length[i] = length(chr[i])
	end
	chr_length
end

function createpopn(chr_length, range = 0:1)
	pr = Array(Array{Real}, 48)
	for i = 1:length(pr)
		pr[i] = sample(range, rand(Binomial(chr_length, 0.5)))
	end
	pr 
end 

function createrealpopn(chr_length)
	pr = Array(Array{Real}, 48)
	for i = 1:length(pr)
		pr[i] = rand(rand(Binomial(chr_length, 0.5)))
	end
	pr 
end 

function mutdisplay(p::Matrix, mp::Matrix, loc, i)
	println("loc = $(loc[i])")
	println("p   = $(p[i,1:end])")
	println("mp  = $(mp[i,1:end])")
end

function mutdisplay(p::Vector, mp::Vector, loc, i)
	println("loc = $(loc[i])")
	println("p   = $(p[i])")
	println("mp  = $(mp[i])")
end

p1 = reshape(sample(0:1, 24 * 48), 48, 24)
p2 = reshape(sample(1:9, 24 * 48), 48, 24)
p3 = rand(48, 8)
p4 = createpopn(18, 0:1)
p5 = createpopn(18, 1:9)
p6 = createrealpopn(8)
mp1 = copy(p1)
mp2 = copy(p2)
mp3 = copy(p3)
mp4 = deepcopy(p4)
mp5 = deepcopy(p5)
mp6 = deepcopy(p6)

chr_length = 24
mutn_percent = 0.2
bmutn = BinaryMutation(MutnPerChr(mutn_percent, chr_length))
loc, genes, alt = mutate!(bmutn, mp1)
mutdisplay(p1,mp1,loc,1)
bmutn = BinaryMutation(MutnPerChr(mutn_percent, chr_length, ptype = :probability))
loc, genes, alt = mutate!(bmutn, mp1)
mutdisplay(p1,mp1,loc,1)
smutn = SymbolicMutation(MutnPerChr(mutn_percent, chr_length, ptype = :percentage), 1:9)
loc, genes, alt = mutate!(smutn, mp2)
mutdisplay(p2,mp2,loc,1)
chr_length = 8
gmutn = GaussianMutation(MutnPerChr(mutn_percent, chr_length, ptype = :percentage), 1.0)
loc, genes, alt = mutate!(gmutn, mp3)
mutdisplay(p3,mp3,loc,1)
bmutn = BinaryMutation(MutnPerChr(mutn_percent))
loc, genes, alt = mutate!(bmutn, mp4)
mutdisplay(p4,mp4,loc,1)
