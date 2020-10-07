
function update!(c::Chain)
  #=
  Integrate out the two sites participating in the strongest type-1 coupling.
  =#

  N = 4 # t, u , v, f terms (type n = 1, 2, 3, 4, respectively)
  L = c.L

  Omega, imax = findmax(c) # find the largest type-1 coupling
  #= steps away from the "zero" site (imax site) that we'll need =#
  deltas = [-3:3, -2:3, -2:2, -2:2]
  # pack the local couplings that are needed (and their locations) into arrays
  inds = [[step(c, imax, delta) for delta in deltas[n]] for n in 1:N]
  lcs = [c.logc[n, inds[n]] for n in 1:N] # lc stands for log coupling

  zs = [1-deltas[n][1] for n in 1:N] # the indices of the "zero site"

  # to store renormalized couplings
  newinds = [vcat(inds[n][1:zs[n]-1], inds[n][zs[n]+2:end]) for n in 1:N]
  newlcs = [fill(-Inf, length(newinds[n])) for n in 1:N]

  # more convenient names for t, u, v, f couplings
  lt = view(lcs[1], :)
  lu = view(lcs[2], :)
  lv = view(lcs[3], :)
  lf = view(lcs[4], :)
  nlt = view(newlcs[1], :)
  nlu = view(newlcs[2], :)
  nlv = view(newlcs[3], :)
  nlf = view(newlcs[4], :)
  zt, zu, zv, zf = zs

  # t-type couplings (n = 1)
  nlt[zt-3] = max(lt[zt-3], lt[zt-1]+lf[zf-2]-Omega, lv[zv-2]+lv[zv-1]-Omega)
  nlt[zt-2] = max(lt[zt-2], lf[zf-1], lv[zv-2]+lf[zf-2]-Omega, lt[zt-1]+lv[zv-1]
    -Omega, lu[zu-1]+lu[zu]-Omega)
  nlt[zt-1] = max(lf[zf], lt[zt-1]+lt[zt+1]-Omega, lu[zu]+lu[zu+1]-Omega)
  nlt[zt] = max(lt[zt+2], lf[zf+1], lv[zv+2]+lf[zf+2]-Omega, lt[zt+1]+lv[zv+1]
    -Omega, lu[zu+1]+lu[zu+2]-Omega)
  nlt[zt+1] = max(lt[zt+3], lt[zt+1]+lf[zf+2]-Omega, lv[zv+1]+lv[zv+2]-Omega)

  # u-type couplings (n = 2)
  nlu[zu-2] = max(lu[zu-2], lu[zu-1]+lf[zf-2]-Omega, lv[zv-2]+lu[zu]-Omega)
  nlu[zu-1] = max(lu[zu-1]+lt[zt+1]-Omega, lv[zv-1]+lu[zu+1]-Omega)
  nlu[zu] = max(lt[zt-1]+lu[zu+2]-Omega, lu[zu]+lv[zv+1]-Omega)
  nlu[zu+1] = max(lu[zu+3], lu[zu+2]+lf[zf+2]-Omega, lu[zu+1]+lv[zv+2]-Omega)

  # v-type couplings (n = 3)
  nlv[zv-2] = lv[zv-2]+lt[zt+1]-Omega
  nlv[zv-1] = max(lu[zu-1]+lu[zu+2]-Omega, lv[zv-1]+lv[zv+1]-Omega)
  nlv[zv] = lt[zt-1]+lv[zv+2]-Omega

  # f-type couplings (n = 4)
  nlf[zf-2] = lf[zf-2]+lt[zt+1]-Omega
  nlf[zf-1] = -Inf #= max(lu[zu-1]+lu[zu]+lf[zf+1]-2*Omega, lf[zf-1]+lu[zu+1]+lu[zu+2]
    -2*Omega, lt[zt-1]+lv[zv-1]+lf[zf+1]-2*Omega, lf[zf-1]+lt[zt+1]+lv[zv+1]
    -2*Omega, lu[zu-1]+lf[zf]+lu[zu+2]-2*Omega, lv[zv-1]+lf[zf]+lv[zv+1]
    -2*Omega, lu[zu-1]+lf[zf]+lv[zv+1]-2*Omega, lv[zv-1]+lf[zf]+lu[zu+2]
    -2*Omega) =#
  nlf[zf] = lt[zt-1]+lf[zf+2]-Omega

  # delete two sites from the chain
  i1, i2, i3, i4 = [step(c, imax, dlta) for dlta in -1:2]
  c.l[i2], c.r[i2], c.l[i3], c.r[i3] = fill(-1, 4) #= sites that get integrated
  out are marked with -1 =#
  c.mask[[i2, i3]] .= false
  c.nsites -= 2
  # i1 and i4 are now neighbors; i2 and i3 have been integrated out
  c.r[i1] = i4
  c.l[i4] = i1

  # insert renormalized log couplings
  for n in 1:N
    c.logc[n, newinds[n]] .= newlcs[n]
  end

  nothing

end

function update!(c::Chain, nsteps::Int64)
  #=
  Update the chain multiple times.
  =#

  for step in 1:nsteps
    update!(c)
  end

  nothing

end
