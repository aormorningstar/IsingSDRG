#
# function analyze(c::Chain)
#   #=
#   An analysis to run on a snapshot of the chain during the RG simulation.
#   Outputs a tuple of:
#     - the fraction of type-n terms that is zero
#     - the fraction of type-1 terms that is smaller than the greatest type-n term
#     - is the strongest term a type-1 term? (1 if yes, 0 if no)
#     - the fraction of type-1 terms that have a stronger type-n (n != 1) term they overlap with
#     - the means of type-n terms
#     - the standard deviations of type-n terms
#   =#
#
#   logc = view(c.logc, :, c.mask)
#   N, nsites = c.N, c.nsites
#
#   @assert (N, nsites) == size(logc), "Chain has inconsistency: probably between mask and nsites."
#
#   maxlogc = vec(maximum(logc, dims=2)) # type-n maximums
#   maxint = maximum(maxlogc[2:end])
#   maxist = argmax(maxlogc) == 1
#   nzero = 0
#   nless = 0
#   intdom = 0
#
#   for i in 1:L
#
#     logc = logcs[:, i]
#
#     nzero += sum(logc[2:end] .== -Inf)
#     nless += logc[1] < maxint
#
#     id = false
#     for n in 2:N
#       js = pbc.(collect(i-n:i+n), L)
#       id = id || maximum(logcs[n, js]) > logc[1]
#     end
#     intdom += id
#
#   end
#
#   nzero/L, nless/L, maxist, intdom/L
#
# end
