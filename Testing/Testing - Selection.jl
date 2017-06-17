# testing normalization
x = [1:4]
y = normalize(x)
x
z = [1.0:4.0]
normalize!(z)
z


# testing SUS
fit = rand(5)
pmf = normalize(fit)
sel = sus(20, pmf, randomize = false)
sel = sus(20, pmf)
sel = sus(20000, pmf);
hsel = hist(sel)
hsel[2]/20000

# testing alt_SUS
rep([5, 3, 7, 2])
fit = rand(5)
pmf = normalize(fit)
sel = alt2_sus(20, pmf, randomize = false)
sel = alt2_sus(20, pmf)
sel = alt2_sus(20000, pmf);
hsel = hist(sel)
hsel[2]/20000

# testing time - both same time, but alt2 takes 33% more space
pmf = normalize(rand(1000));
@time alt2_sus(200000, pmf);
@time sus(200000, pmf);

fit = rand(5)
pmf = normalize(fit)
sel = roulette(pmf, 20)

# testing tournament selector constructors
TournSel(false, 3, 0.2)
TournSel(false, 0.2)
TournSel(false, [0.2, 0.3, 0.4])
TournSel(3, 0.2)
TournSel(0.2)
TournSel([0.2, 0.3, 0.4])
TournSel(false, 3)
TournSel(3)
TournSel(false)
TournSel()	


# testing tournament selection
tsinfo1 = TournSel(10)
tsinfo1s = TournSel(3)
small_fit = rand(6)
selectlocations(tsinfo1s, small_fit, 5, verbose = true)
sl = selectlocations(tsinfo1s, small_fit, 100000);
hist(sl)[2]/100000 										# should be a quadratic distribution
tsinfo2s = TournSel(3, 0.4)
selectlocations(tsinfo2s, small_fit, 5, verbose = true)
setorder!(tsinfo1s, :sorted)
selectlocations(tsinfo1s, small_fit, 24)

tsinfo1 = StrictTournSel(10)
tsinfo1s = StrictTournSel(3)
small_fit = rand(6)
selectlocations(tsinfo1s, small_fit, 5, verbose = true)
sl = selectlocations(tsinfo1s, small_fit, 100000);
hist(sl)[2]/100000 										# should be a quadratic distribution
tsinfo2s = RelaxedTournSel(3, 0.4)
selectlocations(tsinfo2s, small_fit, 5, verbose = true)

# testing time - tournament selection
fit = rand(1000);
@time selectlocations(tsinfo1, fit, 100);
tsinfo2 =  RelaxedTournSel(10, 0.2)
@time selectlocations(tsinfo2, fit, 100);
selectlocations(tsinfo2, small_fit, 5)

# testing UniformSelection
unifsel = UniformSelection()
unifsels = UniformSelection()
setorder!(unifsels, :sorted)
unifselr = UniformSelection()
setorder!(unifselr, :randomized)

selectlocations(unifsel, 12, 24)	# select 12 from a popn_size of 24
selectlocations(unifsels, 12, 24)	# select 12 from a popn_size of 24
selectlocations(unifsel, 24, 24)	# select 12 from a popn_size of 24
selectlocations(unifselr, 24, 24)	# select 12 from a popn_size of 24
selectlocations(unifsel, 25, 24)	# select 12 from a popn_size of 24
selectlocations(unifsels, 25, 24)	# select 12 from a popn_size of 24
selectlocations(unifsel, 48, 24)	# select 12 from a popn_size of 24

unifsel = UniformSelection(:roulette)
unifsels = UniformSelection(:roulette, :sorted)
selectlocations(unifsel, 12, 24)	# select 12 from a popn_size of 24
selectlocations(unifsels, 12, 24)	# select 12 from a popn_size of 24
selectlocations(unifsel, 24, 24)	# select 12 from a popn_size of 24
selectlocations(unifsels, 24, 24)	# select 12 from a popn_size of 24
selectlocations(unifsel, 25, 24)	# select 12 from a popn_size of 24
selectlocations(unifsels, 25, 24)	# select 12 from a popn_size of 24
selectlocations(unifsel, 48, 24)	# select 12 from a popn_size of 24

# testing ranking for truncation selection
fit = rand(1:6,24)
rfit = rank(AllFinalTies(true), fit)
allr = rank(AllFinalTies(false), fit)
rfit = rank(NoFinalTies(true), fit)
none = rank(NoFinalTies(false), fit)
find_final_ties(none, allr, 15)
rfit = rank(RandomFinalTies(true), fit, 15)

@time rfit = rank(AllFinalTies(true), fit);
@time rfit = rank(AllFinalTies(false), fit);
@time rfit = rank(NoFinalTies(true), fit);
@time rfit = rank(NoFinalTies(false), fit);


# testing truncation selection
truncsel = TruncationSelection()
fit = rand(24)
selectlocations(truncsel, fit, 13)
fit = rand(1:6,24)
selectlocations(truncsel, fit, 13)
truncsel = TruncationSelection(:none)

