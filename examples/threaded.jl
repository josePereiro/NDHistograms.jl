@time begin
    using NDHistograms
    using CairoMakie
    using Base.Threads
end

## .-- .-. . .- .---. . ...- -- - --. ..- 
# TODO: add filtering capabilities
@time let
    tasks = map(1:10) do _
        @spawn let
            h = NDHistogram(
                "rand0" => -10.0:0.1:10.0,
                "rand1" => -10.0:0.1:10.0
            )
            for it in 1:1e5
                r0 = randn()
                r1 = r0 * rand()
                count!(h, (r0, r1))
            end
            return h
        end
    end
    global h0 = merge!(fetch.(tasks)...)
end


## .-- .-. . .- .---. . ...- -- - --. ..- 
let
    h1 = marginal(h0, "rand1")
    xs = keys(h1, 1) |> collect
    sidxs = sortperm(xs)
    xs = xs[sidxs]
    ys = values(h1) |> collect
    ys = ys[sidxs]
    scatter(xs, ys)
end