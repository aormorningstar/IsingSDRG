
mutable struct Chain
  #=
  Represents a 1D lattice of Majorana degrees of freedom with disordered terms.
  =#

  L::Int64 # length of the lattice
  logc::Array{Float64, 2} #= log couplings; first index sets type of term;
  second index sets the location (site) =#
  l::Vector{Int64} # maps site to the left-neighboring site
  r::Vector{Int64} # maps site to the right-neighboring site
  mask::Vector{Bool} # maps site to true if it has not been integrated out yet
  nsites::Int64 # the number of sites left = sum(mask)

end

function Chain(L::Int64)
  #=
  Initialize the Chain randomly with length L and nonzero couplings.
  =#

  N = 4 # t, u , v, f terms (n = 1, 2, 3, 4, respectively)
  logc = fill(-Inf, (N, L)) # container to keep log couplings

  # distribution of couplings is exponential, but all n terms are greater than any n-1 term
  dist = Exponential(1)
  cuts = [log(N/(N-(n-1))) for n in N:-1:1] # left boundaries

  # fill the container of log couplings
  jt = ones(Int64, N)
  while sum(jt .<= L) > 0

    c = rand(dist)
    n = findfirst(c .> cuts)

    if jt[n] <= L
      logc[n, jt[n]] = c
      jt[n] += 1
    end

  end

  # set the neighbor relations and start with no sites integrated out
  l = circshift(1:L, 1)
  r = circshift(1:L, -1)
  mask = trues(L)
  nsites = L

  Chain(L, logc, l, r, mask, nsites)

end

#= ----- functions for working with the Chain struct ----- =#

function step(c::Chain, i::Int64, delta::Int64)
  #=
  Starting at site i, go delta sites to the right (goes left if delta is
  negative) not including sites that have been integrated out already.
  =#

  right = delta >= 0 # are we moving to the right?

  # take |di| steps
  for _ in 1:abs(delta)
    i = right ? c.r[i] : c.l[i]
  end

  i

end

function findmax(c::Chain, n::Int64=1)
  #=
  Find the largest type-n term (t, u, v, f are type 1, 2, 3, 4 terms).
  =#

  # step through the chain and find the maximum coupling
  i = findfirst(c.mask)
  imx = i
  mx = c.logc[n, i]

  for _ in 2:c.nsites
    i = c.r[i]
    logc = c.logc[n, i]
    if logc > mx
      imx = i
      mx = logc
    end
  end

  # return both the maximum log coupling and the site
  mx, imx

end
