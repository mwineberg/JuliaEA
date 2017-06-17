Pkg.update()
Pkg.add("Distributions")
using Distributions

ea_path = "/Users/mwineberg/Dropbox/Work/Research/Projects/03 Programs/EAs in Julia/Basic EA/Source"
cd(ea_path)

function load_ea(base_path = "")
	files = ["AbstractTypes.jl",
				#"Aux - VectorComparisons.jl",
			 "Encoding - Abstract Subtypes.jl",
			 "Encoding - Simple Genes.jl",
			 "Encoding - EncodedVector.jl",
			 "Population - PopulationInfo.jl",
			 "Population - Chromosomes.jl",
			 "Population - PhenotypeMembers.jl",
			 "Problem - OptimumInfo.jl",
			 "Population - FitnessValues.jl",
			 "Population - OptimalMembers.jl",
			 "Problem - OneMax.jl",
			 "Problem - SimpleFn.jl",
			 "Problem - SimpleFn library.jl",
			 "Population.jl",
			 "Selection - General.jl",
			 #"Selection - Rank.jl",
			 "Selection - Scaled FPS and Sigma Scaling.jl",
			 "Selection - Tournament.jl",
			 "Selection - Truncation and Elitism.jl",
			 "Selection - Uniform (Random).jl",
			 "Reproduction - Mutation.jl",
			 "Reproduction - Crossover.jl",
			 "EA Parameters",
			 "EA Counts"]

	for file in files
		include("$base_path$file")
	end
end

load_ea()
