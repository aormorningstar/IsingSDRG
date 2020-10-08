using HDF5
using Dates
using ArgParse
using Distributions
using StatsBase
using Statistics

include("../src/IsingSDRG.jl")
include("./utils.jl")

function run(; Li::Int64, Lf::Int64, nt::Int64, np::Int64, alpha::Float64, kwargs...)
    #=
    The "main" function. Run the simulations. Output the data.
    =#

    @assert alpha > 0 "alpha > 0 must be true"
    @assert Li > Lf "initial number of sites must be larger than final"

    # set up integer times spaced approx evenly on a log scale
    T = div(Li - Lf, 2)
    ts = logspace(1, T, nt)
    dts = diff(ts)

    # initialize chain
    d1 = Exponential(alpha) # t ~ t^{-(1+alpha)}, so logt ~ exp(-alpha logt)
    a, b = -1.1, -0.9 # default params
    d2 = Uniform(a, b) # a default dist for all other couplings
    c = Chain(Li, [d1, d2, d2, d2])

    # track some points of the empirical cumulative probablity distribution
    p = quantile_probs(np)
    q_data = Array{Float64, 3}(undef, (np, c.N, nt)) # container for the data

    # RG updates
    for (i, dt) in enumerate(dts)

        update!(c, dt)

        # calculate quantiles for all types of couplings
        for j in 1:c.N
            quantile!(view(q_data, :, j, i), c.logc[j, c.mask], p) # sorts the copied data in place
        end

    end

    ts, p, q_data

end

# parse command line arguments
aps = ArgParseSettings()

@add_arg_table aps begin

    "--Li" # length of system initially
        arg_type = Int64
        default = 1000000

    "--Lf" # number of sites left after the RG is run
        arg_type = Int64
        default = 1000

    "--nt" # number of snapshots (in RG time) to take
        arg_type = Int64
        default = 100

    "--np" # number of quantiles to record (evaluations of the empirical CDF)
        arg_type = Int64
        default = 100

    "--alpha" # controls power law for the initial distribution of the type-1 couplings
        arg_type = Float64
        default = 0.1

    "--testing" # store data in local testing directory
        action = :store_true

end

const ad = parse_args(aps, as_symbols=true)
const testing = ad[:testing]

# run the simulation
t_dat, p_dat, q_dat = run(; ad...)
data = Dict("time" => t_dat, "cdf" => p_dat, "logcoupling" => q_dat)

# set up data file name
local_data_dir = "../data/testing/"
della_data_dir = "/scratch/gpfs/aorm/ising_sdrg/"
data_dir = testing ? local_data_dir : della_data_dir
timestamp = Dates.format(now(), "yyyy_mm_dd_HH_MM_SS")

data_file_name = join([timestamp, string(rand(Int64))[2:7], dict_to_string(ad)], "_")
data_file_path = join([data_dir, data_file_name, ".hdf5"])

# write to file
data_file = h5open(data_file_path, "w")
for ky in keys(ad)
    strky = string(ky)
    if strky != "testing"
        attrs(data_file)[strky] = ad[ky]
    end
end

for ky in keys(data)
    write(data_file, ky, data[ky])
end

close(data_file)