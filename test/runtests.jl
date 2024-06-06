using NDHistograms
using Test

## --.- . -. - .-.- - 
@testset "NDHistograms.jl" begin
        
    ## --------------------------------------------------------
    # _find_nearest
    let
        _v = 1:100
        @test NDHistograms._find_nearest(10, _v) == 10
        _v = [1:100;] .+ 0.1
        @test NDHistograms._find_nearest(10, _v) == 10
    end
    return

    ## --------------------------------------------------------
    # Histogram (threaded)
    let
        Random.seed!(123)

        h0 = rangehistogram(-1e2, 1e2; step = 0.1)

        h_pool = [deepcopy(h0) for th in 1:nthreads()]

        T = 1 # K
        @threads :static for _ in 1:nthreads()
            h = h_pool[threadid()]
            for it in 1:Int(1e5)
                p = rand()
                count!(h, - log(p) * T)
            end
        end
        count!(h0, h_pool...) # reduce

        @show sum(values(h0))
        xs = sort(collect(keys(h0)))
        xs = xs[2:end-1] # avoid edge artifacts
        ws = counts(h0, xs)
        ws = ws ./ maximum(ws)
        
        fs = exp.(-xs/T)
        fs = fs ./ maximum(fs)

        wH = sum(ws .* log.(ws))
        @show wH
        fH = sum(fs .* log.(fs))
        @show fH
        # @show isapprox(wH, fH; rtol = 1e-2)
        @test isapprox(wH, fH; rtol = 1e-2)

        # p = plot()
        # plot!(p, xs, ws; label = "h")
        # plot!(p, xs, fs; label = "f")
        # p
    end

    ## ------------------------------------------------------------
    # Identity histogram
    let
        Random.seed!(123)
        
        # Histogram
        h0 = identity_histogram(Vector{Int})

        # count
        n = 2
        h_pool = [deepcopy(h0) for th in 1:nthreads()]
        @threads :static for _ in 1:nthreads()
            h = h_pool[threadid()]
            for it in 1:Int(1e6)
                comb = rand(1:5, n)
                count!(h, comb)
            end
        end
        count!(h0, h_pool...) # reduce

        @show length(values(h0))
        
        ws = collect(values(h0))
        ws = ws ./ maximum(ws)

        # @show all(isapprox.(ws, 1.0; rtol = 5e-2))
        @test all(isapprox.(ws, 1.0; rtol = 1e-2))

        # plot(ws; label = "", 
        #     ylim = [0, maximum(ws)]
        # )
    end

end
