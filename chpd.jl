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
# Run with:  julia --project=. chpd.jl        (one time step)
#            julia --project=. chpd.jl 3      (the first three time steps)

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
# Pass a number on the command line to run more steps, and "delays" to switch
# on the pipeline storage of Eqs. (6)-(13).
NSTEPS = (isempty(ARGS) || !occursin(r"^\d+$", ARGS[1])) ? 1 : parse(Int, ARGS[1])
DELAYS = "delays" in ARGS

# "from=N" starts the horizon at hour N instead of hour 1. Useful because the
# interesting hours of this case study are the evening ones, where the wind
# drops away while the heat load peaks, and a short window has to be placed
# over them to show anything.
FROM = let a = findfirst(s -> startswith(s, "from="), ARGS)
    a === nothing ? 1 : parse(Int, split(ARGS[a], '=')[2])
end
println("Horizon: hours $FROM..$(FROM + NSTEPS - 1), ",
        DELAYS ? "pipeline time delays ON" : "no time delays")

include(joinpath(@__DIR__, "solvers.jl"))
include(joinpath(@__DIR__, "init_model.jl"))
include(joinpath(@__DIR__, "delays.jl"))
include(joinpath(@__DIR__, "build_chpd_minlp.jl"))
include(joinpath(@__DIR__, "build_chpd_misocp.jl"))
include(joinpath(@__DIR__, "build_ced.jl"))
include(joinpath(@__DIR__, "results.jl"))

tsw = ts[FROM:min(FROM + NSTEPS - 1, nrow(ts)), :]

function setup(optimizer; nsteps::Int=NSTEPS)
    m = Model(optimizer)
    define_sets!(m, data, tsw; nsteps=nsteps)
    process_time_series_data!(m, data, tsw)
    process_parameters!(m, data)
    return m
end

## Step 3: build the models
m_minlp = setup(nonconvex_solver())
build_chpd_minlp!(m_minlp; delays=DELAYS)

m_misocp = setup(convex_solver())
build_chpd_misocp!(m_misocp; delays=DELAYS)

# The delay model is hard enough that a solver can spend its whole budget
# looking for a first feasible point. The no-delay version of the same model is
# cheap and its solution is always feasible with every delay set to zero, so it
# is solved first and handed over as a starting point. See warm_start_from!.
if DELAYS
    println("\nwarm start: solving the no-delay model first")
    m_warm = setup(nonconvex_solver())
    build_chpd_minlp!(m_warm; delays=false)
    safe_optimize!(m_warm, "  no-delay warm start")
    if solved(m_warm)
        warm_start_from!(m_minlp, m_warm)
        warm_start_from!(m_misocp, m_warm)
    else
        println("  warm start unavailable, carrying on without one")
    end
end

m_ced = setup(lp_solver())
build_ced!(m_ced)

## Step 4: solve
println("\n================ solving ================")
safe_optimize!(m_ced, "Section IV-A (CED)")
safe_optimize!(m_minlp, "Section II (MINLP)")
safe_optimize!(m_misocp, "Section III-B (relaxed)")

## Step 5: validation
println("\n================ validation over $NSTEPS time step$(NSTEPS == 1 ? "" : "s") ================")
check_lower_bound(m_misocp, m_minlp)

if solved(m_misocp)
    report_envelope_widths(m_misocp)
    check_eq36(m_misocp)
    check_mccormick(m_misocp)
end

# The relaxed solution is a bound, not a dispatch. The Section II solution is
# the one that has to make physical sense, so both get checked.
if solved(m_minlp)
    println("\n>>> physical checks on the Section II solution (must all close)")
    check_physics(m_minlp)
end
if solved(m_misocp)
    println("\n>>> the same checks on the relaxed solution (these are allowed to drift)")
    check_physics(m_misocp)
end

println("\n--- CHPD versus CED ---")
zc = objective_value(m_ced)
println("  CED  : $(round(zc, digits=2)) \$")
for (label, mm) in [("CHPD, Section II   ", m_minlp), ("CHPD, Section III-B", m_misocp)]
    solved(mm) || continue
    zh = objective_value(mm)
    println("  $label: $(round(zh, digits=2)) \$",
            "   difference vs CED: $(round((zh - zc) / zc * 100, digits=2)) %")
end
if !DELAYS
    println("  (with the delays off the DHN cannot store anything, so the CHPD is",
            " expected to cost MORE than the CED: it pays for the water pumps,",
            " which the CED does not model. The saving in the paper is a storage",
            " effect and needs Eqs. (6)-(13).)")
end
