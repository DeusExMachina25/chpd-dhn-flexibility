using Test
using JuMP, YAML, CSV, DataFrames

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "solvers.jl"))
include(joinpath(ROOT, "init_model.jl"))
include(joinpath(ROOT, "build_chpd_minlp.jl"))
include(joinpath(ROOT, "build_chpd_misocp.jl"))

silent(; nonconvex=false) = nonconvex ? nonconvex_solver() : convex_solver()

@testset "McCormick envelope" begin
    # With x pinned to one of its own bounds the envelope is exact, so the
    # feasible w collapses to the single point x*y.
    for xfix in (1.0, 3.0)
        m = Model(silent())
        @variable(m, 1 <= x <= 3)
        @variable(m, 2 <= y <= 5)
        @variable(m, w)
        mccormick!(m, w, x, y, 1.0, 3.0, 2.0, 5.0)
        @constraint(m, x == xfix)
        @constraint(m, y == 4.0)

        @objective(m, Min, w); optimize!(m); wmin = objective_value(m)
        @objective(m, Max, w); optimize!(m); wmax = objective_value(m)
        @test wmin ≈ xfix * 4.0 atol = 1e-6
        @test wmax ≈ xfix * 4.0 atol = 1e-6
    end

    # In the interior the envelope is a genuine relaxation: it must contain
    # the true product but is allowed to be wider.
    m = Model(silent())
    @variable(m, 1 <= x <= 3)
    @variable(m, 2 <= y <= 5)
    @variable(m, w)
    mccormick!(m, w, x, y, 1.0, 3.0, 2.0, 5.0)
    @constraint(m, x == 2.0)
    @constraint(m, y == 3.5)
    @objective(m, Min, w); optimize!(m); wmin = objective_value(m)
    @objective(m, Max, w); optimize!(m); wmax = objective_value(m)
    @test wmin <= 2.0 * 3.5 + 1e-9
    @test wmax >= 2.0 * 3.5 - 1e-9
end

@testset "Two-node network against the analytic optimum" begin
    # Geometry: node 1 (CHP) --pipe--> node 2 (heat load), one time step.
    #
    # With the delays switched off the time delay is zero, and the loss bracket
    # of Eq. (5) is multiplied by that delay, so there is no thermal loss:
    # TS2 = TS1 and TR1 = TR2. Then
    #     LH = c*mf*(TS2 - TR2) = c*mf*(TS1 - TR2)
    #     Q  = c*mf*(TS1 - TR1) = c*mf*(TS1 - TR2)
    # so Q = LH exactly - whatever the CHP puts in arrives at the load - and
    # the only freedom left is how much water to move. Cost is increasing in
    # the pump load, hence in mf, so the optimum takes the widest temperature
    # difference the bounds allow: TS1 at its ceiling, TR2 at its floor.
    #
    # (An earlier version of this test carried a loss term and expected
    # Q > LH. That was wrong; see NOTES.md section 4.)
    data = YAML.load_file(joinpath(@__DIR__, "testcase.yaml"))
    ts = CSV.read(joinpath(@__DIR__, "timeseries.csv"), DataFrame)

    c = 4.182e-3
    LH = 100.0
    TS1 = 373.0                # upper bound of the supply temperature
    TR1 = 313.0                # lower bound of the return temperature
    TR2 = TR1                  # no loss, so the two ends agree
    mf_exp = LH / (c * (TS1 - TR2))
    Q_exp = LH

    # Only the original model has to reproduce this. The relaxation is checked
    # separately, below - a relaxation is not supposed to land on the exact
    # dispatch, only to bound it.
    m = Model(silent(nonconvex=true))
    define_sets!(m, data, ts; nsteps=1)
    process_time_series_data!(m, data, ts)
    process_parameters!(m, data)
    build_chpd_minlp!(m)
    optimize!(m)

    @test termination_status(m) in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)

    TS = value.(m.ext[:variables][:TS]); TR = value.(m.ext[:variables][:TR])
    mfS = value.(m.ext[:variables][:mfS]); Q = value.(m.ext[:variables][:Q])

    @test TS[1, 1] ≈ TS1 atol = 1e-3
    @test TR[1, 1] ≈ TR1 atol = 1e-3
    @test TR[2, 1] ≈ TR2 atol = 1e-3
    @test mfS[1, 1] ≈ mf_exp rtol = 1e-4
    @test Q["CHP1", 1] ≈ Q_exp rtol = 1e-4

    # With no delay there are no losses, so production must equal the load.
    @test Q["CHP1", 1] ≈ LH rtol = 1e-9

    # And the heat equation must hold on the recovered numbers.
    @test Q["CHP1", 1] ≈ c * mfS[1, 1] * (TS[1, 1] - TR[1, 1]) rtol = 1e-6

    println("  Section II: mf = $(round(mfS[1,1], digits=3)) kg/s (expected $(round(mf_exp, digits=3))),",
            " Q = $(round(Q["CHP1",1], digits=3)) MW (expected $(round(Q_exp, digits=3)))")
end

@testset "Relaxation bounds the original from below" begin
    data = YAML.load_file(joinpath(@__DIR__, "testcase.yaml"))
    ts = CSV.read(joinpath(@__DIR__, "timeseries.csv"), DataFrame)

    function solve_with(builder, opt)
        m = Model(opt)
        define_sets!(m, data, ts; nsteps=1)
        process_time_series_data!(m, data, ts)
        process_parameters!(m, data)
        builder(m)
        optimize!(m)
        return objective_value(m)
    end

    z_orig = solve_with(build_chpd_minlp!, silent(nonconvex=true))
    z_rel = solve_with(build_chpd_misocp!, silent())

    # The whole point of Section III-B: a bound, cheap to compute, never above
    # the thing it relaxes.
    @test z_rel <= z_orig + 1e-6
    println("  relaxed $(round(z_rel, digits=2)) \$ <= original $(round(z_orig, digits=2)) \$",
            "   (gap $(round((z_orig - z_rel) / z_orig * 100, digits=2)) %)")
end
