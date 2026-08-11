## Power Systems Flexibility from District Heating Networks
# Implementation of the Combined Heat and Power Dispatch of
#   L. Mitridati and J. A. Taylor, "Power systems flexibility from district
#   heating networks", PSCC 2018, doi:10.23919/PSCC.2018.8442617
#
# Three models are built on the same data:
#   - Section II    the original bilinear MINLP
#   - Section III-B the convex relaxation (McCormick + convex quadratic)
#   - Section IV-A  the conventional economic dispatch, as a benchmark
#
# Run with:  julia --project=. chpd.jl

## Step 0: Activate environment - ensure consistency across computers
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

## Step 1: input data
using CSV
using DataFrames
using YAML

data = YAML.load_file(joinpath(@__DIR__, "data", "network.yaml"))
ts = CSV.read(joinpath(@__DIR__, "data", "timeseries.csv"), DataFrame)
println("Data loaded: $(length(data["dhnNodes"]))-node DHN, $(length(data["buses"]))-bus power system")

## Step 2: create model & pass data to model
using JuMP

# One time step to start with, no pipeline time delays. See NOTES.md.
NSTEPS = 1

include(joinpath(@__DIR__, "solvers.jl"))
include(joinpath(@__DIR__, "init_model.jl"))
include(joinpath(@__DIR__, "build_chpd_minlp.jl"))
include(joinpath(@__DIR__, "build_chpd_misocp.jl"))
include(joinpath(@__DIR__, "build_ced.jl"))
include(joinpath(@__DIR__, "results.jl"))

function setup(optimizer; nsteps::Int=NSTEPS)
    m = Model(optimizer)
    define_sets!(m, data, ts; nsteps=nsteps)
    process_time_series_data!(m, data, ts)
    process_parameters!(m, data)
    return m
end

## Step 3: build the models
m_minlp = setup(nonconvex_solver())
build_chpd_minlp!(m_minlp)

m_misocp = setup(convex_solver())
build_chpd_misocp!(m_misocp)

m_ced = setup(lp_solver())
build_ced!(m_ced)

## Step 4: solve
println("\n================ solving ================")
optimize!(m_ced);    report_status(m_ced, "Section IV-A (CED)")
optimize!(m_minlp);  report_status(m_minlp, "Section II (MINLP)")
optimize!(m_misocp); report_status(m_misocp, "Section III-B (relaxed)")

## Step 5: validation
println("\n================ validation of the one-step case ================")
check_lower_bound(m_misocp, m_minlp)
report_envelope_widths(m_misocp)
check_eq36(m_misocp)
check_mccormick(m_misocp)

# The relaxed solution is a bound, not a dispatch. The Section II solution is
# the one that has to make physical sense, so both get checked.
println("\n>>> physical checks on the Section II solution (must all close)")
check_physics(m_minlp)
println("\n>>> the same checks on the relaxed solution (these are allowed to drift)")
check_physics(m_misocp)

println("\n--- CHPD versus CED ---")
zc = objective_value(m_ced)
zh = objective_value(m_minlp)
println("  CED  : $(round(zc, digits=2)) \$")
println("  CHPD : $(round(zh, digits=2)) \$")
println("  difference: $(round((zh - zc) / zc * 100, digits=2)) %")
println("  (at a single time step the DHN cannot store anything, so the CHPD",
        " is expected to cost slightly MORE than the CED: it pays for the",
        " water pumps and the pipe heat losses that the CED ignores.)")
