## Validation of the solved models
# The point of these checks is to answer three questions:
#   1. is the Section III-B model a valid relaxation of Section II?
#   2. is that relaxation tight, i.e. do the dropped equalities hold anyway?
#   3. does the solution make physical sense?
# At a single time step all of this is small enough to be checked by hand,
# which is exactly why we start there.

function report_status(m::Model, label::String)
    st = termination_status(m)
    println(rpad(label, 22), " status: ", st,
            st == MOI.OPTIMAL || st == MOI.LOCALLY_SOLVED ?
            "   objective: $(round(objective_value(m), digits=3)) \$" : "")
    return st
end

# 1. The relaxation has to be a lower bound of the original minimisation.
function check_lower_bound(m_relaxed::Model, m_original::Model)
    zr = objective_value(m_relaxed)
    zo = objective_value(m_original)
    gap = (zo - zr) / abs(zo) * 100
    println("\n--- Relaxation quality ---")
    println("  Section III-B (relaxed) : $(round(zr, digits=3)) \$")
    println("  Section II   (original) : $(round(zo, digits=3)) \$")
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

    res["Eq. (27) pump [MW]"] = maximum(abs(Lpump[j, t] - p[:Pa_to_MW] / (rho * p[:etaPump][j]) *
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
    mS = maximum(abs(sum(mfS[pp, t] for pp in SPminus[n]) + sum(mfHS[j, t] for j in SHS[n])
                     - sum(mfS[pp, t] for pp in SPplus[n]) - sum(mfHES[i, t] for i in SHES[n]))
                 for n in IN, t in T)
    mR = maximum(abs(sum(mfR[pp, t] for pp in SPplus[n]) + sum(mfHES[i, t] for i in SHES[n])
                     - sum(mfR[pp, t] for pp in SPminus[n]) - sum(mfHS[j, t] for j in SHS[n]))
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
    pb = maximum(abs(LE[b, t] + sum(LHP[j, t] for j in SHP[b]) + sum(Lpump[j, t] for j in SHSbus[b])
                     - sum(P[g, t] for g in SE[b])
                     - sum(p[:B][(b, mm)] * (theta[mm, t] - theta[b, t]) for mm in SB[b]))
                 for b in IB, t in T)
    over = maximum(abs(p[:B][(b, mm)] * (theta[mm, t] - theta[b, t])) / p[:fmax][(b, mm)]
                   for b in IB, mm in SB[b], t in T)
    println("  power balance residual  $(round(pb, sigdigits=3)) MW")
    println("  worst line loading      $(round(over * 100, digits=1)) % of the limit")

    # Pump consumption, the part of the cost the CED never sees
    println("  total water pump load   $(round(sum(Lpump[j, t] for j in IHS, t in T), digits=3)) MW")
    return nothing
end
