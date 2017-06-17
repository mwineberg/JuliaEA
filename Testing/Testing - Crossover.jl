# from "Selection - General.jl"
import StatsBase.sample
import StatsBase.sample!
import StatsBase.weights

# from "Chromosomes.jl"
function chrlength(chr::Matrix)
	size(chr, 2)
end

x = [  10 	 20	   30	40	 50;
	   55 	 65	   75	85	 95;
	  100 	200	  300  400  500;
	  555	666	  777  888  999;
	 1000  2000  3000 4000 5000;
	 5432  6543  7654 8765 9876]

xvr21 = HomologousCrossover(size(x,2), xover_type = :Uniform, parent_count = 2, offspring_count = 1)
xvr22 = HomologousCrossover(size(x,2), xover_type = :Uniform, parent_count = 2, offspring_count = 2)
xvr31 = HomologousCrossover(size(x,2), xover_type = :Uniform, parent_count = 3, offspring_count = 1)
xvr32 = HomologousCrossover(size(x,2), xover_type = :Uniform, parent_count = 3, offspring_count = 2)
xvr33 = HomologousCrossover(size(x,2), xover_type = :Uniform, parent_count = 3, offspring_count = 3)

sloc = shuffle([1:6])
setlocations!(xvr21, sloc)
y = crossover(xvr21, x)
setlocations!(xvr22, sloc)
y = crossover(xvr22, x)
setlocations!(xvr31, sloc)
y = crossover(xvr31, x)
setlocations!(xvr32, sloc)
y = crossover(xvr32, x)
setlocations!(xvr33, sloc)
y = crossover(xvr33, x)

xkpt = HomologousCrossover(18, xover_type = :KPt, parent_count = 2, offspring_count = 1, k = 1)
createmask!(xkpt.xtype)
xkpt.xtype