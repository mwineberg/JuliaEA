import StatsBase.mean
import StatsBase.std
import StatsBase.mean_and_std


################################################################################################
## Fitness Proportional Selection - includes fitness scaling/shifts including removing the minimum fitness
## 		  - uses perform.selection in "Selection - General.R"

abstract FitPropSel <: SelectionType
abstract ScaledFitPropSel <: FitPropSel
abstract SigmaScaling <: ScaledFitPropSel

type BasicFitPropSel <: FitPropSel
	sel_type::Symbol
	order::Symbol
end

type BasicScaledFitPropSel <: ScaledFitPropSel
	sel_type::Symbol
	order::Symbol
	shift::Real
end

type RemoveMinFitPropSel <: ScaledFitPropSel
	sel_type::Symbol
	order::Symbol
	shift::Real
end

type TaneseSigmaScaling <: SigmaScaling
	sel_type::Symbol
	order::Symbol
	α::Real 							# α = 1/desired(cv) = desired (mean/std) - Tanese set this to 2
	base_prob::Real
end

type MinFitSigmaScaling <: SigmaScaling
	sel_type::Symbol
	order::Symbol
	α::Real 							# α = 1/desired(cv) = desired (mean/std) - Tanese set this to 2
end


function FPS(;sel_type = :roulette, order = :original, remove_min = false)
	(remove_min ? RemoveMinFitPropSel(sel_type, order, 0.0)
				: BasicFitPropSel(sel_type, order))
end

function FPS(shift::Integer; sel_type = :roulette, order = :original, remove_min = false)
	(remove_min ? RemoveMinFitPropSel(sel_type, order, shift)
				: BasicScaledFitPropSel(sel_type, order, shift))
end

function SigmaScalingFPS(base_prob::Real; factor = 2.0, sel_type = :roulette, order = :original)
	TaneseSigmaScaling(sel_type, order, factor, base_prob)
end

function SigmaScalingFPS(;factor = 2.0, sel_type = :roulette, order = :original)
	MinFitSigmaScaling(sel_type, order, factor)
end

################################################################################################

function selectlocations(info::BasicFitPropSel, fit::Vector, count::Integer; verbose = false)
	normalize!(fit)
	selectloc(info, fit, count)
end

function selectlocations(info::BasicScaledFitPropSel, fit::Vector, count::Integer; verbose = false)
	fit += info.shift
	normalize!(fit)
	selectloc(info, fit, count)
end

function selectlocations(info::RemoveMinFitPropSel, fit::Vector, count::Integer; verbose = false)
	fit += info.shift - min(fit)
	normalize!(fit)
	selectloc(info, fit, count)
end

# Sigma scaling: Tanese version (perhaps slightly modifed)
function selectlocations(info::TaneseSigmaScaling, fit::Vector, count::Integer; verbose = false)
	n = length(fit)
	μ, σ = mean_and_std(fit)
	fit += 	α * σ - μ
	# the denominator that normalizes the fitness in the pmf calculation below is the sum of fitnesses; all fitness terms cancel out
	pmf = fit/(n * info.factor * σ)
	for i = 1:n
		if pmf[i] <= info.base_pr
			pmf[i] = info.base_pr
		end
	end
	normalize!(pmf)
	selectloc(info, pmf, count)
end

# Sigma scaling: an alternative to Tanses solution to the problem of negative fitnesses after shifting
#	-
function selectlocations(info::MinFitSigmaScaling, fit::Vector, count::Integer; verbose = false)
	n = length(fit)
	μ, σ = mean_and_std(fit)
	min_fit = min(fit)
	shift =   α * σ  - μ
	shift =  (min_fit + shift < 0) ? -min_fit : shift
	fit += shift
	normalize!(fit)
	selectloc(info, fit, count)
end
