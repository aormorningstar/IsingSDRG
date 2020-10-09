
function update!(c::Chain)::Nothing
  #=
  Integrate out the two sites participating in the strongest XX coupling. Refer to XX, XIX, XIIX,
  XXXX types by t, u, v, f, respectively.
  =#

  N, L = c.N, c.L

  lt0, imax = findmax(c) # find the largest XX (t type) coupling

  #= steps away from the "zero" site (imax site) that we'll need for each type =#
  dis = [-3:3, -2:3, -2:2, -2:2]

  # pack the local couplings that are needed (and their locations) into arrays
  inds = [[step(c, imax, di) for di in dis[n]] for n in 1:N]
  lcs = [c.logc[n, inds[n]] for n in 1:N] # lc stands for log coupling

  zs = [1-dis[n][1] for n in 1:N] # the indices of the "zero site" for each type

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

  # t type couplings, XX, n = 1
  nlt[zt-3] = max(lt[zt-3], lt[zt-1]+lf[zf-2]-lt0, lv[zv-2]+lv[zv-1]-lt0)
  nlt[zt-2] = max(lt[zt-2], lf[zf-1], lv[zv-2]+lf[zf-2]-lt0, lt[zt-1]+lv[zv-1]
    -lt0, lu[zu-1]+lu[zu]-lt0)
  nlt[zt-1] = max(lf[zf], lt[zt-1]+lt[zt+1]-lt0, lu[zu]+lu[zu+1]-lt0)
  nlt[zt] = max(lt[zt+2], lf[zf+1], lv[zv+2]+lf[zf+2]-lt0, lt[zt+1]+lv[zv+1]
    -lt0, lu[zu+1]+lu[zu+2]-lt0)
  nlt[zt+1] = max(lt[zt+3], lt[zt+1]+lf[zf+2]-lt0, lv[zv+1]+lv[zv+2]-lt0)

  # u type couplings, XIX, n = 2
  nlu[zu-2] = max(lu[zu-2], lu[zu-1]+lf[zf-2]-lt0, lv[zv-2]+lu[zu]-lt0)
  nlu[zu-1] = max(lu[zu-1]+lt[zt+1]-lt0, lv[zv-1]+lu[zu+1]-lt0)
  nlu[zu] = max(lt[zt-1]+lu[zu+2]-lt0, lu[zu]+lv[zv+1]-lt0)
  nlu[zu+1] = max(lu[zu+3], lu[zu+2]+lf[zf+2]-lt0, lu[zu+1]+lv[zv+2]-lt0)

  # v type couplings, XIIX, n = 3
  nlv[zv-2] = lv[zv-2]+lt[zt+1]-lt0
  nlv[zv-1] = max(lu[zu-1]+lu[zu+2]-lt0, lv[zv-1]+lv[zv+1]-lt0)
  nlv[zv] = lt[zt-1]+lv[zv+2]-lt0

  # f type couplings, XXXX, n = 4
  nlf[zf-2] = lf[zf-2]+lt[zt+1]-lt0
  nlf[zf-1] = max(lu[zu-1]+lu[zu]+lf[zf+1]-2*lt0, lf[zf-1]+lu[zu+1]+lu[zu+2]
    -2*lt0, lt[zt-1]+lv[zv-1]+lf[zf+1]-2*lt0, lf[zf-1]+lt[zt+1]+lv[zv+1]
    -2*lt0, lu[zu-1]+lf[zf]+lu[zu+2]-2*lt0, lv[zv-1]+lf[zf]+lv[zv+1]
    -2*lt0, lu[zu-1]+lf[zf]+lv[zv+1]-2*lt0, lv[zv-1]+lf[zf]+lu[zu+2]
    -2*lt0) # -Inf
  nlf[zf] = lt[zt-1]+lf[zf+2]-lt0

  # delete two sites from the chain
  i1, i2, i3, i4 = [step(c, imax, di) for di in -1:2]
  c.l[i2], c.r[i2], c.l[i3], c.r[i3] = fill(-1, 4) # mark deleted sites with -1
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

function update!(c::Chain, nsteps::Int64)::Nothing
  #=
  Update the chain multiple times.
  =#

  for step in 1:nsteps
    update!(c)
  end

  nothing

end
