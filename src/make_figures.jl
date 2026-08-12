## Figures for the report.
#
# Two figures, both about the single-time-step result that the assignment asks
# for. Run with:
#
#   julia --project=. src/make_figures.jl          (hour 20, the default)
#   julia --project=. src/make_figures.jl 1        (any other hour)
#
# Output lands in results/ as PNG.

const ROOT = normpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(ROOT)

using CSV, DataFrames, YAML, JuMP, Plots

include(joinpath(@__DIR__, "solvers.jl"))
include(joinpath(@__DIR__, "init_model.jl"))
include(joinpath(@__DIR__, "delays.jl"))
include(joinpath(@__DIR__, "build_chpd_minlp.jl"))
include(joinpath(@__DIR__, "build_chpd_misocp.jl"))

HOUR = isempty(ARGS) ? 20 : parse(Int, ARGS[1])

data = YAML.load_file(joinpath(ROOT, "data", "network.yaml"))
ts = CSV.read(joinpath(ROOT, "data", "timeseries.csv"), DataFrame)
tsw = ts[HOUR:HOUR, :]

outdir = joinpath(ROOT, "results")
mkpath(outdir)
gr()

function setup(optimizer)
    m = Model(optimizer)
    define_sets!(m, data, tsw; nsteps=1)
    process_time_series_data!(m, data, tsw)
    process_parameters!(m, data)
    return m
end

println("Building the Section II reference at hour $HOUR ...")
m2 = setup(nonconvex_solver())
build_chpd_minlp!(m2)
optimize!(m2)
z2 = objective_value(m2)
println("  Section II: $(round(z2, digits=3)) \$")

## ---------------------------------------------------------------- figure 1
# Where the heat and the electricity actually come from at this hour. This is
# the dispatch the single-step validation is checking.

Qv = value.(m2.ext[:variables][:Q])
Pv = value.(m2.ext[:variables][:P])
IHS = m2.ext[:sets][:IHS]
IE = m2.ext[:sets][:IE]

heat_names = collect(IHS)
heat_vals = [Qv[j, 1] for j in heat_names]
elec_names = collect(IE)
elec_vals = [Pv[g, 1] for g in elec_names]

# drop units that produce nothing, they only add clutter
hk = heat_vals .> 1e-6
ek = elec_vals .> 1e-6

p1 = bar(heat_names[hk], heat_vals[hk],
         title="Heat production", ylabel="MW", legend=false,
         color=:steelblue, bar_width=0.5)
p2 = bar(elec_names[ek], elec_vals[ek],
         title="Electricity production", ylabel="MW", legend=false,
         color=:indianred, bar_width=0.5)

fig1 = plot(p1, p2, layout=(1, 2), size=(760, 320),
            plot_title="Section II dispatch, hour $HOUR",
            titlefontsize=10, plot_titlefontsize=11, guidefontsize=9,
            tickfontsize=8, left_margin=5Plots.mm, bottom_margin=5Plots.mm)
savefig(fig1, joinpath(outdir, "fig1_dispatch_hour$HOUR.png"))
println("  wrote fig1_dispatch_hour$HOUR.png")

for (n, v) in zip(heat_names, heat_vals)
    println("    heat  $n = $(round(v, digits=2)) MW")
end
for (n, v) in zip(elec_names, elec_vals)
    println("    elec  $n = $(round(v, digits=2)) MW")
end

## ---------------------------------------------------------------- figure 2
# The relaxation gap as a function of how wide a box the McCormick envelopes
# are built on. `widen` interpolates the heat-exchanger and station flow bounds
# from the load-implied ones (0) back to the raw equipment limits (1). The
# Section II model is unchanged throughout, so every point is a gap against the
# same reference.

println("Sweeping the envelope box width ...")
widens = collect(0.0:0.1:1.0)
gaps = Float64[]
boxw = Float64[]

for w in widens
    mr = setup(convex_solver())
    build_chpd_misocp!(mr; widen=w)
    optimize!(mr)
    zr = objective_value(mr)

    mfHES = mr.ext[:variables][:mfHES]
    i1 = first(mr.ext[:sets][:IHES])
    width = upper_bound(mfHES[i1, 1]) - lower_bound(mfHES[i1, 1])

    push!(gaps, (z2 - zr) / z2 * 100)
    push!(boxw, width)
    println("  widen $(round(w, digits=1)): box $(round(width, digits=1)) kg/s",
            " -> relaxed $(round(zr, digits=1)) \$, gap $(round(gaps[end], digits=2)) %")
end

fig2 = plot(boxw, gaps,
            marker=:circle, markersize=5, linewidth=2, color=:darkorange,
            legend=false,
            xlabel="width of the heat-exchanger flow box [kg/s]",
            ylabel="relaxation gap [%]",
            title="Relaxation gap against McCormick box width, hour $HOUR",
            size=(720, 380), titlefontsize=11, guidefontsize=9, tickfontsize=8,
            left_margin=5Plots.mm, bottom_margin=5Plots.mm)
annotate!(fig2, boxw[1], gaps[1],
          text("  load-implied bounds\n  (what the model uses)", 8, :left, :bottom))
savefig(fig2, joinpath(outdir, "fig2_gap_vs_width_hour$HOUR.png"))
println("  wrote fig2_gap_vs_width_hour$HOUR.png")
