using Test
using JuMP, Gurobi, YAML, CSV, DataFrames

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "init_model.jl"))
include(joinpath(ROOT, "build_chpd_minlp.jl"))
include(joinpath(ROOT, "build_chpd_misocp.jl"))

silent(; nonconvex=false) = nonconvex ?
    optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => 2, "OutputFlag" => 0) :
    optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0)

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
    # Geometry: node 1 (CHP) --pipe--> node 2 (heat load). One time step, so
    # the outlet temperature is just the inlet times (1 - gamma).
    #
    # Mass balance forces a single flow mf through everything, and Eqs. (19),
    # (20) give TS2 = TS1*(1-g) and TR1 = TR2*(1-g). Then
    #     LH = c*mf*(TS2 - TR2) = c*mf*(TS1*(1-g) - TR2)
    #     Q  = c*mf*(TS1 - TR1) = c*mf*(TS1 - TR2*(1-g))
    # so, eliminating mf,
    #     Q = LH * (TS1 - TR2*(1-g)) / (TS1*(1-g) - TR2).
    # dQ/dTS1 < 0 and dQ/dTR2 > 0, so the cheapest dispatch pushes the supply
    # temperature to its upper bound and the return temperature to its lower
    # bound. Everything else follows.
    data = YAML.load_file(joinpath(@__DIR__, "testcase.yaml"))
    ts = CSV.read(joinpath(@__DIR__, "timeseries.csv"), DataFrame)

    g = 2 * 1.0 * 3600.0 / (963.0 * 4182.0 * 0.30)   # gamma of the pipe
    c = 4.182e-3
    LH = 100.0
    TS1 = 373.0        # upper bound of the supply temperature
    TR2 = 313.0        # lower bound of the return temperature
    mf_exp = LH / (c * (TS1 * (1 - g) - TR2))
    Q_exp = LH * (TS1 - TR2 * (1 - g)) / (TS1 * (1 - g) - TR2)

    for (name, builder, opt) in [("Section II", build_chpd_minlp!, silent(nonconvex=true)),
                                 ("Section III-B", build_chpd_misocp!, silent())]
        m = Model(opt)
        define_sets!(m, data, ts; nsteps=1)
        process_time_series_data!(m, data, ts)
        process_parameters!(m, data)
        builder(m)
        optimize!(m)

        @test termination_status(m) in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)

        TS = value.(m.ext[:variables][:TS]); TR = value.(m.ext[:variables][:TR])
        mfS = value.(m.ext[:variables][:mfS]); Q = value.(m.ext[:variables][:Q])

        @test TS[1, 1] ≈ TS1 atol = 1e-3
        @test TR[2, 1] ≈ TR2 atol = 1e-3
        @test mfS[1, 1] ≈ mf_exp rtol = 1e-4
        @test Q["CHP1", 1] ≈ Q_exp rtol = 1e-4

        # The loss identity Q - LH = c*mf*gamma*(TS1 + TR2), derived by hand
        # from the same three equations.
        @test Q["CHP1", 1] - LH ≈ c * mfS[1, 1] * g * (TS[1, 1] + TR[2, 1]) rtol = 1e-6

        println("  $name: mf = $(round(mfS[1,1], digits=3)) kg/s, Q = $(round(Q["CHP1",1], digits=3)) MW")
    end
end

@testset "Relaxation is a valid lower bound" begin
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
    @test z_rel <= z_orig + 1e-6
end
