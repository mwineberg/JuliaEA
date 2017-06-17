Pkg.update()
Pkg.add("StatsBase")
Pkg.add("Distributions")

import Base: print, getindex      # used throughout "Population - *"
# import Base: <=, <, >=, >, ==     # for "Aux - Vector Comparisons"
import StatsBase.sample           # for "Selection - General"
import StatsBase.weights          # for "Selection - General"
import StatsBase.mean             # for "Selection - Scaled FPS and Sigma Scaling"
import StatsBase.std              # for "Selection - Scaled FPS and Sigma Scaling"
import StatsBase.mean_and_std     # for "Selection - Scaled FPS and Sigma Scaling"
import StatsBase.tiedrank         # for "Selection - Truncation and Elitism"
import StatsBase.competerank      # for "Selection - Truncation and Elitism"
import Distributions.Binomial     # for "Reproduction - Mutation"
