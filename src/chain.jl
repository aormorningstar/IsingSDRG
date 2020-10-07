
mutable struct Chain
  #=
  Represents a disordered Hamiltonian written in terms of Majorana operators. The operators act on
  a one-dimensional lattice of Majorana degrees of freedom with disordered coefficients. Sites get
  eliminated via an RG procedure.
  =#

  L::Int64 # total number of original sites
  N::Int64 # types of terms

  logc::Array{Float64, 2} # log couplings; first dimension is the type of term, second is the site

  l::Vector{Int64} # maps a site to its left-neighboring site
  r::Vector{Int64} # maps a site to its right-neighboring site

  mask::Vector{Bool} # maps site to true if it has not been integrated out yet
  nsites::Int64 # the number of sites remaining = sum(mask)

  function Chain(L::Int64, logc::Array{Float64, 2}, l::Vector{Int64}, r::Vector{Int64},
  mask::Vector{Bool}, nsites::Int64)::Chain
      #=
      default constructor
      =#

      N = 4 # only considering case of XX, XIX, XIIX, XXXX terms
      new(L, N, logc, l, r, mask, nsites)

  end

end

function Chain(L::Int64, dists::Vector{UnivariateDistribution})::Chain
  #=
  Initialize the chain with length L and nonzero couplings sampled from dists.
  =#

  N = 4 # only considering case of XX, XIX, XIIX, XXXX terms

  @assert length(dists) == N "Must provide 4 distributions."

  logc = Array{Float64, 2}(undef, (N, L)) # make space for log couplings

  # init the log couplings randomly
  for i in 1:N

      logc[i, :] .= rand(dists[i], L)

  end

  # set the neighbor relations and start with no sites integrated out
  l = circshift(1:L, 1)
  r = circshift(1:L, -1)
  mask = trues(L)
  nsites = L

  Chain(L, logc, l, r, mask, nsites)

end

function step(c::Chain, i::Int64, di::Int64)::Int64
  #=
  Starting at site i, go di sites to the right (left if delta is negative), skipping over sites
  that have been integrated out.
  =#

  right = delta >= 0 # are we moving to the right?

  # take |di| steps
  for _ in 1:abs(di)

    i = right ? c.r[i] : c.l[i]

  end

  i

end

function findmax(c::Chain, n::Int64=1)::Tuple{Float64, Int64}
  #=
  Find the largest type-n term (t, u, v, f are type 1, 2, 3, 4 terms).
  =#

  @assert 1 < n <= c.N "n must be in [1, N]"

  # step through the chain and find the maximum coupling
  i = findfirst(c.mask)
  imx, mx = i, c.logc[n, i]

  for _ in 1:c.nsites-1

    i = c.r[i]
    logc = c.logc[n, i]

    if logc > mx
      imx, mx = i, logc
    end

  end

  mx, imx

end
