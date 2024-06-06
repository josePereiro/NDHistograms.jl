@time begin
    using NDHistograms
    using CairoMakie
end

## .-- .-. . .- .---. . ...- -- - --. ..- 
@time let
    global h0 = NDHistogram(
        "rand0" => -10.0:0.1:10.0,
        "rand1" => -10.0:0.1:10.0
    )
    for it in 1:1e5
        r0 = randn()
        r1 = r0 * rand()
        count!(h0, (r0, r1))
    end
end

## .-- .-. . .- .---. . ...- -- - --. ..- 
# test rebin
let
    dim = "rand1"
    h1 = marginal(h0, dim)
    h2 = rebin(h1, dim => -10.0:1.0:10.0)
    
    f = Figure()
    ax = Axis(f[1,1]; xlabel = "value", ylabel = "count")
    
    x1s = keys(h1, 1) |> collect
    sidxs = sortperm(x1s)
    x1s = x1s[sidxs]
    y1s = values(h1) |> collect
    y1s = y1s[sidxs]

    x2s = keys(h2, 1) |> collect
    sidxs = sortperm(x2s)
    x2s = x2s[sidxs]
    y2s = values(h2) |> collect
    y2s = y2s[sidxs]

    barplot!(ax, x2s, y2s ./ maximum(y2s); color = :blue)
    scatter!(ax, x1s, y1s ./ maximum(y1s); color = :red)
    lines!(ax, x1s, y1s ./ maximum(y1s); color = :red)

    f
end

## .-- .-. . .- .---. . ...- -- - --. ..- 
## .-- .-. . .- .---. . ...- -- - --. ..- 
## .-- .-. . .- .---. . ...- -- - --. ..- 
# TEST marginal/filter/rand
let
    h1 = filter(h0) do v
        v0 = v[1]
        # filter
        # return true
        return v0 > 1.0 || v0 < -1.0
    end
    h1 = marginal(h1, "rand0")
    
    f = Figure()
    ax = Axis(f[1,1])
    
    xs = keys(h1, 1) |> collect
    sidxs = sortperm(xs)
    xs = xs[sidxs]
    ys = values(h1) |> collect
    ys = ys[sidxs]

    s = [rand(h1)[1] for i in 1:1e5]
    
    # barplot!(ax, xs, ys ./ sum(ys))
    hist!(ax, s; 
        normalization = :probability, 
        color = :values, bins = length(xs) + 13, 
        fillto = 0.0, 
    )
    scatter!(ax, xs, ys ./ sum(ys); color = :red)
    lines!(ax, xs, ys ./ sum(ys); color = :red)
    
    f
end