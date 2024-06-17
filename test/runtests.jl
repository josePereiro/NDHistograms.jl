@time begin
    using NDHistograms
    # using CairoMakie
    using Distributions
    using Base.Threads
    using Random
    using Test
end

# .-- .-. . .- .---. . ...- -- - --. ..- 
@time @testset "NDHistograms.jl" begin
    Random.seed!(124)

    # random distribution
    D = 3
    A = rand(D, D)
    Σ = A' * A
    N = MultivariateNormal(zeros(D), Σ)
    steps = rand([0.05, 0.03, 0.1], D)

    ntasks = 2 * nthreads()
    tasks = map(1:ntasks) do _
        @spawn let
            h = NDHistogram([
                "dim$d" => -100.0:steps[d]:100.0 
                for d in 1:D
            ]...)
            for it in 1:2e5
                x = rand(N)
                count!(h, Tuple(x))
            end
            return h
        end
    end
    h0 = merge!(map(fetch, tasks)...)

    H0 = entropy(N)
    @show H0
    H1 = entropy(h0)
    @show H1
    @test abs(H1 - H0) / abs(H0) < 0.05
end
nothing

