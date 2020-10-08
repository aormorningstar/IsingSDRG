
pbc(i::Int64, L::Int64) = mod(i-1, L)+1 # periodic boundaries for a chain with sites in [1, L]

function logspace(x1::Integer, xn::Integer, n::Integer)
    #=
    n values spaced evenly on a log scale in the range [x1, xn].
    =#

    @assert x1 < xn && n > 2 "Invalid inputs."

    lambda = log(xn/x1) / (n - 1)
    delta = log(x1) - lambda

    xi = [round(Integer, exp(lambda*i + delta)) for i in 1:n]

    @assert xi[1] == x1 && xi[end] == xn "logspace failed."

    xi

end

function dict_to_string(d::Dict)
    #=
    Convert a dict of "stringable" objects and keys into a long string.
    =#

    flat_list = collect(Iterators.flatten([[p[1],p[2]] for p in collect(d)]))
    
    join(string.(flat_list), "_")

end

function quantile_probs(n::Integer)
    #=
    Get n probability values in (0.0, 1.0). Endpoints 0 and 1 are excluded.
    =#

    dp = 1.0 / (n + 1.0)

    [i*dp for i in 1:n]

end
