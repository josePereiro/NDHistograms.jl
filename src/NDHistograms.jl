module NDHistograms

import Random: default_rng
import Random: AbstractRNG
using Base.Threads
using Distributions

#! include .
include("0_base.jl")
include("1_threaded.jl")

end
