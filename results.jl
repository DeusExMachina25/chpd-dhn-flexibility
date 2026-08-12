## Validation of the solved models
# The point of these checks is to answer three questions:
#   1. is the Section III-B model a valid relaxation of Section II?
#   2. is that relaxation tight, i.e. do the dropped equalities hold anyway?
#   3. does the solution make physical sense?
# At a single time step all of this is small enough to be checked by hand,
# which is exactly why we start there.

# A model counts as solved only if there is an actual primal point to read.
# has_values alone is not enough: a solver that hits its time limit with a dual
# bound but no incumbent still reports values, and they come back as NaN.
solved(m::Model) = termination_status(m) != MOI.OPTIMIZE_NOT_CALLED &&
                   primal_status(m) == MOI.FEASIBLE_POINT &&
                   isfinite(objective_value(m))

# Solving can fail outright rather than just fail to converge - the usual cause
# here is the size-limited Gurobi licence refusing a model with a few thousand
# delay binaries. That is worth reporting and carrying on with, not crashing
# the whole run, because the other models still have something to say.
function safe_optimize!(m::Model, label::String)
    try
        optimize!(m)
    catch err
        msg = sprint(showerror, err)
        if occursin("size-limited", msg)
            # The licence, not the model, is the problem. Hand it to SCIP,
            # which has no size cap, and say so rather than doing it silently.
            println(rpad(label, 22), " too large for the size-limited Gurobi",
                    " licence, retrying with Juniper")
            try
                set_optimizer(m, uncapped_solver())
                optimize!(m)
            catch err2
                println(rpad(label, 22), " FAILED TO SOLVE (Juniper)")
                println("    ", first(split(sprint(showerror, err2), '\n')))
                return :failed
            end
            return report_status(m, label * " [Juniper]")
        end
        println(rpad(label, 22), " FAILED TO SOLVE")
        println("    ", first(split(msg, '\n')))
        return :failed
    end
    return report_status(m, label)
end

function report_status(m::Model, label::String)
    st = termination_status(m)
    msg = rpad(label, 22) * " status: " * string(st)
    if !solved(m) && st == MOI.TIME_LIMIT
        msg *= "   no feasible point found; dual bound only: " *
               "$(round(objective_bound(m), digits=3)) \$"
    elseif solved(m)
        msg *= "   objective: $(round(objective_value(m), digits=3)) \$"
        # When the time limit bites the incumbent is not proven optimal, so
        # say how far the solver still had to go.
        if st != MOI.OPTIMAL && st != MOI.LOCALLY_SOLVED
            msg *= "   bound: $(round(objective_bound(m), digits=3)) \$" *
                   "   (MIP gap $(round(relative_gap(m) * 100, digits=2)) %)"
        end
    end
    println(msg)
    return st
end

# 1. The relaxation has to be a lower bound of the original minimisation.
#
# If either model stopped on the time limit the comparison still works, but it
# has to be made between the right quantities: the relaxation's valid lower
# bound is its dual bound, and the original's valid upper bound is its
# incumbent. Comparing two unproven incumbents would say nothing.
function check_lower_bound(m_relaxed::Model, m_original::Model)
    if !solved(m_relaxed) || !solved(m_original)
        println("\n--- Relaxation quality ---")
        println("  skipped: one of the models has no solution")
        return NaN
    end
    exact = termination_status(m_relaxed) in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED) &&
            termination_status(m_original) in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)

    zr = exact ? objective_value(m_relaxed) : objective_bound(m_relaxed)
    zo = objective_value(m_original)
    gap = (zo - zr) / abs(zo) * 100

    println("\n--- Relaxation quality ---")
    println("  Section III-B (relaxed) : $(round(zr, digits=3)) \$",
            exact ? "" : "   (dual bound, time limit hit)")
    println("  Section II   (original) : $(round(zo, digits=3)) \$",
            exact ? "" : "   (incumbent, time limit hit)")
    println("  gap                     : $(round(gap, digits=4)) %")
    if zr > zo + 1e-6
        println("  !! the relaxation is ABOVE the original - envelopes or Eq. (36) are wrong")
    else
        println("  ok: valid lower bound")
    end
    return gap
end

# 2a. Is Eq. (36) tight, i.e. does the relaxed pressure loss sit on the
# original quadratic equality (2)?
function check_eq36(m::Model)
    T = m.ext[:sets][:T]; IP = m.ext[:sets][:IP]
    pfrom = m.ext[:sets][:pfrom]; pto = m.ext[:sets][:pto]
    phi = m.ext[:parameters][:phi]
    mfS = value.(m.ext[:variables][:mfS]); mfR = value.(m.ext[:variables][:mfR])
    prS = value.(m.ext[:variables][:prS]); prR = value.(m.ext[:variables][:prR])

    println("\n--- Exactness of Eq. (36), pressure loss ---")
    worst = 0.0
    for pp in IP, t in T
        sS = (prS[pfrom[pp], t] - prS[pto[pp], t]) - phi[pp] * mfS[pp, t]^2
        sR = (prR[pto[pp], t] - prR[pfrom[pp], t]) - phi[pp] * mfR[pp, t]^2
        worst = max(worst, abs(sS), abs(sR))
        println("  pipe $pp, t=$t : supply slack $(round(sS, digits=2)) Pa, return slack $(round(sR, digits=2)) Pa")
    end
    println("  worst slack: $(round(worst, digits=2)) Pa  (0 means the relaxation is exact here)")
    return worst
end

# 2b. Do the bilinear equalities that McCormick replaced still hold?
function check_mccormick(m::Model)
    T = m.ext[:sets][:T]
    IHS = m.ext[:sets][:IHS]; IHES = m.ext[:sets][:IHES]
    INmixS = m.ext[:sets][:INmixS]; INmixR = m.ext[:sets][:INmixR]
    SPminus = m.ext[:sets][:SPminus]; SPplus = m.ext[:sets][:SPplus]
    p = m.ext[:parameters]; c = p[:c]; rho = p[:rho]
    LH = m.ext[:timeseries][:LH]

    TS = value.(m.ext[:variables][:TS]); TR = value.(m.ext[:variables][:TR])
    TSout = value.(m.ext[:variables][:TSout]); TRout = value.(m.ext[:variables][:TRout])
    mfS = value.(m.ext[:variables][:mfS]); mfR = value.(m.ext[:variables][:mfR])
    mfHS = value.(m.ext[:variables][:mfHS]); mfHES = value.(m.ext[:variables][:mfHES])
    prS = value.(m.ext[:variables][:prS]); prR = value.(m.ext[:variables][:prR])
    Q = value.(m.ext[:variables][:Q]); Lpump = value.(m.ext[:variables][:Lpump])

    println("\n--- Exactness of the McCormick envelopes ---")
    res = Dict{String,Float64}()

    res["Eq. (21) HES heat [MW]"] = maximum(abs(LH[i, t] - c * mfHES[i, t] *
        (TS[p[:HESnode][i], t] - TR[p[:HESnode][i], t])) for i in IHES, t in T)

    res["Eq. (24) HS heat [MW]"] = maximum(abs(Q[j, t] - c * mfHS[j, t] *
        (TS[p[:HSnode][j], t] - TR[p[:HSnode][j], t])) for j in IHS, t in T)

    res["Eq. (27) pump [MW]"] = maximum(abs(Lpump[j, t] - p[:pressure_to_MW] / (rho * p[:etaPump][j]) *
        mfHS[j, t] * (prS[p[:HSnode][j], t] - prR[p[:HSnode][j], t])) for j in IHS, t in T)

    res["Eq. (16) mixing [kg.K/s]"] = isempty(INmixS) ? 0.0 :
        maximum(abs(TS[n, t] * sum(mfS[pp, t] for pp in SPminus[n]) -
                    sum(mfS[pp, t] * TSout[pp, t] for pp in SPminus[n])) for n in INmixS, t in T)

    res["Eq. (17) mixing [kg.K/s]"] = isempty(INmixR) ? 0.0 :
        maximum(abs(TR[n, t] * sum(mfR[pp, t] for pp in SPplus[n]) -
                    sum(mfR[pp, t] * TRout[pp, t] for pp in SPplus[n])) for n in INmixR, t in T)

    for (k, v) in sort(collect(res))
        println("  ", rpad(k, 28), " max residual: ", round(v, sigdigits=4))
    end
    return res
end

# How loose can each envelope possibly be? For w = x*y over a box, the worst
# case distance between the McCormick hull and the true surface is
# (xhi-xlo)*(yhi-ylo)/4, reached in the middle of the box. This says up front
# how much error the bounds allow, before any solving happens.
function report_envelope_widths(m::Model)
    T = m.ext[:sets][:T]
    IHS = m.ext[:sets][:IHS]; IHES = m.ext[:sets][:IHES]
    INmixS = m.ext[:sets][:INmixS]; SPminus = m.ext[:sets][:SPminus]
    p = m.ext[:parameters]; c = p[:c]
    mfHS = m.ext[:variables][:mfHS]; mfHES = m.ext[:variables][:mfHES]
    dTHS = m.ext[:variables][:dTHS]; dTHES = m.ext[:variables][:dTHES]

    println("\n--- How loose the envelopes are allowed to be ---")
    for i in IHES, t in T
        dx = upper_bound(mfHES[i, t]) - lower_bound(mfHES[i, t])
        dy = upper_bound(dTHES[i, t]) - lower_bound(dTHES[i, t])
        println("  Eq. (21) $i t=$t: mf range $(round(dx)) kg/s, dT range $(round(dy)) K",
                " -> up to $(round(c * dx * dy / 4, digits=1)) MW of slack")
    end
    for j in IHS, t in T
        dx = upper_bound(mfHS[j, t]) - lower_bound(mfHS[j, t])
        dy = upper_bound(dTHS[j, t]) - lower_bound(dTHS[j, t])
        println("  Eq. (24) $j t=$t: mf range $(round(dx)) kg/s, dT range $(round(dy)) K",
                " -> up to $(round(c * dx * dy / 4, digits=1)) MW of slack")
    end
    for n in INmixS, pp in SPminus[n]
        dx = p[:TSmax][n] - p[:TSmin][n]
        dy = p[:mfSmax][pp] - p[:mfSmin][pp]
        println("  Eq. (16) node $n pipe $pp: T range $(round(dx)) K, mf range $(round(dy)) kg/s",
                " -> up to $(round(dx * dy / 4)) kg.K/s of slack")
    end
    return nothing
end

# 3. Physical sanity of the dispatch.
function check_physics(m::Model)
    T = m.ext[:sets][:T]
    IN = m.ext[:sets][:IN]; IP = m.ext[:sets][:IP]; IB = m.ext[:sets][:IB]
    IHS = m.ext[:sets][:IHS]; IHES = m.ext[:sets][:IHES]
    ICHP = m.ext[:sets][:ICHP]; IHP = m.ext[:sets][:IHP]
    SPplus = m.ext[:sets][:SPplus]; SPminus = m.ext[:sets][:SPminus]
    SHS = m.ext[:sets][:SHS]; SHES = m.ext[:sets][:SHES]
    SE = m.ext[:sets][:SE]; SHP = m.ext[:sets][:SHP]
    SHSbus = m.ext[:sets][:SHSbus]; SB = m.ext[:sets][:SB]
    p = m.ext[:parameters]; c = p[:c]
    LE = m.ext[:timeseries][:LE]; LH = m.ext[:timeseries][:LH]
    COP = m.ext[:timeseries][:COP]

    mfS = value.(m.ext[:variables][:mfS]); mfR = value.(m.ext[:variables][:mfR])
    mfHS = value.(m.ext[:variables][:mfHS]); mfHES = value.(m.ext[:variables][:mfHES])
    TSin = value.(m.ext[:variables][:TSin]); TSout = value.(m.ext[:variables][:TSout])
    TRin = value.(m.ext[:variables][:TRin]); TRout = value.(m.ext[:variables][:TRout])
    Q = value.(m.ext[:variables][:Q]); P = value.(m.ext[:variables][:P])
    LHP = value.(m.ext[:variables][:LHP]); Lpump = value.(m.ext[:variables][:Lpump])
    theta = value.(m.ext[:variables][:theta])

    println("\n--- Physical checks ---")

    # Nodal mass balance, Eq. (14)
    mS = maximum(abs(sum(mfS[pp, t] for pp in SPminus[n]; init=0.0) + sum(mfHS[j, t] for j in SHS[n]; init=0.0)
                     - sum(mfS[pp, t] for pp in SPplus[n]; init=0.0) - sum(mfHES[i, t] for i in SHES[n]; init=0.0))
                 for n in IN, t in T)
    mR = maximum(abs(sum(mfR[pp, t] for pp in SPplus[n]; init=0.0) + sum(mfHES[i, t] for i in SHES[n]; init=0.0)
                     - sum(mfR[pp, t] for pp in SPminus[n]; init=0.0) - sum(mfHS[j, t] for j in SHS[n]; init=0.0))
                 for n in IN, t in T)
    println("  mass balance residual   supply $(round(mS, sigdigits=3)) kg/s, return $(round(mR, sigdigits=3)) kg/s")

    # Heat balance: production = load + thermal losses in the pipes
    for t in T
        prod = sum(Q[j, t] for j in IHS)
        load = sum(LH[i, t] for i in IHES)
        loss = sum(c * mfS[pp, t] * (TSin[pp, t] - TSout[pp, t]) for pp in IP) +
               sum(c * mfR[pp, t] * (TRin[pp, t] - TRout[pp, t]) for pp in IP)
        println("  t=$t heat: produced $(round(prod, digits=3)) MW = load $(round(load, digits=3))",
                " + losses $(round(loss, digits=3))  -> residual $(round(prod - load - loss, sigdigits=3))")
    end

    # CHP feasible region and heat pump COP
    for j in ICHP, t in T
        fuel = p[:fuelE][j] * P[j, t] + p[:fuelH][j] * Q[j, t]
        println("  t=$t $j: P=$(round(P[j,t], digits=2)) MW, Q=$(round(Q[j,t], digits=2)) MW,",
                " P-r*Q=$(round(P[j,t] - p[:r][j]*Q[j,t], digits=2)) (>=0),",
                " fuel=$(round(fuel, digits=1))/$(p[:Fmax][j]) MW")
    end
    for j in IHP, t in T
        println("  t=$t $j: Q=$(round(Q[j,t], digits=2)) MW, LHP=$(round(LHP[j,t], digits=2)) MW,",
                " Q-COP*LHP=$(round(Q[j,t] - COP[j,t]*LHP[j,t], sigdigits=3))")
    end

    # Power balance and line loading
    pb = maximum(abs(LE[b, t] + sum(LHP[j, t] for j in SHP[b]; init=0.0) + sum(Lpump[j, t] for j in SHSbus[b]; init=0.0)
                     - sum(P[g, t] for g in SE[b]; init=0.0)
                     - sum(p[:B][(b, mm)] * (theta[mm, t] - theta[b, t]) for mm in SB[b]; init=0.0))
                 for b in IB, t in T)
    over = maximum(abs(p[:B][(b, mm)] * (theta[mm, t] - theta[b, t])) / p[:fmax][(b, mm)]
                   for b in IB for mm in SB[b] for t in T)
    println("  power balance residual  $(round(pb, sigdigits=3)) MW")
    println("  worst line loading      $(round(over * 100, digits=1)) % of the limit")

    # Pump consumption, the part of the cost the CED never sees
    println("  total water pump load   $(round(sum(Lpump[j, t] for j in IHS, t in T), digits=3)) MW")
    return nothing
end
