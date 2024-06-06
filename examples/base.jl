@time begin
    using NDHistograms
    using CairoMakie
end

## .-- .-. . .- .---. . ...- -- - --. ..- 
@time let
    # NDHistogram is just a Dict where each key is a point (a tuple) in
    # a given space and each value is a count of the frequency of such value.

    # define support space (each pair is a dimention name => support)

    global h0 = NDHistogram(
        # support dim called 'rand0'
        "rand0" => -10.0:0.5:10.0, 
        # support dim called 'rand1'
        # Note that this support is not homogeneusly distributed
        "rand1" => [-1 * 10 .^ (-3:0.1:1.0); 10 .^ (-3:0.1:1.0)] |> sort 
    )
    # check 'descretize' to see the supported supports :). 

    for it in 1:1e5
        r0 = randn()
        r1 = r0 * rand()
        count!(h0, (r0, r1)) # count the occurence of (r0, r1)
    end
end

## .-- .-. . .- .---. . ...- -- - --. ..- 
# test rebin
let
    dim = "rand0" # play with this
    # build a new NDHistogram which is the marginal of the given dimentions
    # The dimentions can be a range (ex 1:2)
    h1 = marginal(h0, dim)
    # rebin change the support of a dimention and recompute the histogram
    h2 = rebin(h1, dim => -10.0:1.0:10.0)
    
    # plots
    f = Figure()
    ax = Axis(f[1,1]; xlabel = "value", ylabel = "count")
    
    # The Dict-like interface is use to retrive data
    # 'keys' get all points in the support which has a non zero count
    x1s = keys(h1, 1) |> collect
    sidxs = sortperm(x1s)
    x1s = x1s[sidxs]
    # 'values' get the counts for each point
    y1s = values(h1) |> collect
    y1s = y1s[sidxs]

    x2s = keys(h2, 1) |> collect
    sidxs = sortperm(x2s)
    x2s = x2s[sidxs]
    y2s = values(h2) |> collect
    y2s = y2s[sidxs]

    barplot!(ax, x2s, y2s ./ maximum(y2s); color = :gray)
    scatter!(ax, x1s, y1s ./ maximum(y1s); color = :red)
    lines!(ax, x1s, y1s ./ maximum(y1s); color = :red)

    f
end

## .-- .-. . .- .---. . ...- -- - --. ..- 
# TEST marginal/filter/rand
let
    # Filter build a new NDHistogram but excluding some portion of the
    # original support space
    h1 = filter(h0) do v
        # 'v' is a point in the support space
        v1 = v[1]
        # filter
        return v1 < 1.0 && v1 > -1.0
    end
    
    # marginal 
    dim = "rand1" # play with this
    h1 = marginal(h1, dim)
    
    xs = keys(h1, 1) |> collect
    sidxs = sortperm(xs)
    xs = xs[sidxs]
    ys = values(h1) |> collect
    ys = ys[sidxs]

    # The interface rand(::NDHistogram) is supported
    # It retursn a set of points in the support which must 
    # have a distribution regulated by its weights on the NDHistogram
    s = [rand(h1)[1] for i in 1:1e5]

    # Plots
    f = Figure()
    ax = Axis(f[1,1], title = "non homogeneous support")
    hist!(ax, s; 
        normalization = :probability, 
        color = :values, bins = length(xs) + 2, 
    )
    scatter!(ax, xs, ys ./ sum(ys); color = :red)
    lines!(ax, xs, ys ./ sum(ys); color = :red)
    
    # NOTE: that the support shape is important for this probability. 
    # in the case of 'rand1' the support is not homogeneous...
    # For sampling such case you must rebin first to an homogeneous support
    h2 = rebin(h1, dim => -10.0:0.1:10.0)
    s2 = [rand(h2)[1] for i in 1:1e5]
    
    xs = keys(h2, 1) |> collect
    sidxs = sortperm(xs)
    xs = xs[sidxs]
    ys = values(h2) |> collect
    ys = ys[sidxs]

    ax = Axis(f[1,2], title = "homogeneous support")
    hist!(ax, s2; 
        normalization = :probability, 
        color = :values, bins = length(xs) + 2, 
    )
    scatter!(ax, xs, ys ./ sum(ys); color = :red)
    lines!(ax, xs, ys ./ sum(ys); color = :red)
    
    f
end