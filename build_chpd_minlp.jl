## Section II: the Combined Heat and Power Dispatch, as written in the paper
# This is the original mixed integer nonlinear model. The quadratic equality
# (2) and the bilinear equalities (16), (17), (21), (24) and (27) are kept
# exactly as they appear in the paper, which makes the model non-convex. It is
# built so that build_chpd_misocp! can relax it afterwards.
#
# `delays` switches the pipeline energy storage on. With delays = false the
# time delay is zero everywhere, Eqs. (5)-(13) collapse into a single loss
# factor, and no binary variables are created. That is the one-step case.

function build_chpd_minlp!(m::Model; delays::Bool=false)
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()

    # Extract sets
    T = m.ext[:sets][:T]
    IN = m.ext[:sets][:IN]
    IP = m.ext[:sets][:IP]
    IHS = m.ext[:sets][:IHS]
    IHES = m.ext[:sets][:IHES]
    ICHP = m.ext[:sets][:ICHP]
    IHP = m.ext[:sets][:IHP]
    IE = m.ext[:sets][:IE]
    IB = m.ext[:sets][:IB]
    pfrom = m.ext[:sets][:pfrom]
    pto = m.ext[:sets][:pto]
    SPplus = m.ext[:sets][:SPplus]
    SPminus = m.ext[:sets][:SPminus]
    SHS = m.ext[:sets][:SHS]
    SHES = m.ext[:sets][:SHES]
    SE = m.ext[:sets][:SE]
    SHP = m.ext[:sets][:SHP]
    SHSbus = m.ext[:sets][:SHSbus]
    SB = m.ext[:sets][:SB]
    INminus = m.ext[:sets][:INminus]
    INplus = m.ext[:sets][:INplus]
    INmixS = m.ext[:sets][:INmixS]
    INmixR = m.ext[:sets][:INmixR]

    # Extract time series
    LE = m.ext[:timeseries][:LE]
    LH = m.ext[:timeseries][:LH]
    COP = m.ext[:timeseries][:COP]
    Wavail = m.ext[:timeseries][:Wavail]

    # Extract parameters
    p = m.ext[:parameters]
    c = p[:c]
    rho = p[:rho]
    gamma = p[:gamma]
    phi = p[:phi]
    HSnode = p[:HSnode]
    HESnode = p[:HESnode]

    if delays
        error("Time delays (Eqs. 5-13) are not implemented yet - phase 7.")
    end

    # Global temperature bounds, used for the pipe temperature variables and
    # later for the McCormick envelopes.
    TSlo = minimum(p[:TSmin][n] for n in IN); TShi = maximum(p[:TSmax][n] for n in IN)
    TRlo = minimum(p[:TRmin][n] for n in IN); TRhi = maximum(p[:TRmax][n] for n in IN)
    m.ext[:parameters][:TSlo] = TSlo; m.ext[:parameters][:TShi] = TShi
    m.ext[:parameters][:TRlo] = TRlo; m.ext[:parameters][:TRhi] = TRhi

    ##### Create variables
    # Mass flow rates in the pipes, bounded by Eq. (1)
    mfS = m.ext[:variables][:mfS] = @variable(m, [pp=IP, t=T],
        lower_bound=p[:mfSmin][pp], upper_bound=p[:mfSmax][pp], base_name="mfS")
    mfR = m.ext[:variables][:mfR] = @variable(m, [pp=IP, t=T],
        lower_bound=p[:mfRmin][pp], upper_bound=p[:mfRmax][pp], base_name="mfR")

    # Nodal pressures, bounded by Eq. (15)
    prS = m.ext[:variables][:prS] = @variable(m, [n=IN, t=T],
        lower_bound=p[:prSmin][n], upper_bound=p[:prSmax][n], base_name="prS")
    prR = m.ext[:variables][:prR] = @variable(m, [n=IN, t=T],
        lower_bound=p[:prRmin][n], upper_bound=p[:prRmax][n], base_name="prR")

    # Nodal temperatures, bounded by Eq. (18)
    TS = m.ext[:variables][:TS] = @variable(m, [n=IN, t=T],
        lower_bound=p[:TSmin][n], upper_bound=p[:TSmax][n], base_name="TS")
    TR = m.ext[:variables][:TR] = @variable(m, [n=IN, t=T],
        lower_bound=p[:TRmin][n], upper_bound=p[:TRmax][n], base_name="TR")

    # Pipe inlet and outlet temperatures
    TSin = m.ext[:variables][:TSin] = @variable(m, [pp=IP, t=T],
        lower_bound=TSlo, upper_bound=TShi, base_name="TSin")
    TSout = m.ext[:variables][:TSout] = @variable(m, [pp=IP, t=T],
        lower_bound=TSlo * (1 - maximum(values(gamma))), upper_bound=TShi, base_name="TSout")
    TRin = m.ext[:variables][:TRin] = @variable(m, [pp=IP, t=T],
        lower_bound=TRlo, upper_bound=TRhi, base_name="TRin")
    TRout = m.ext[:variables][:TRout] = @variable(m, [pp=IP, t=T],
        lower_bound=TRlo * (1 - maximum(values(gamma))), upper_bound=TRhi, base_name="TRout")

    # Mass flow rates at the stations, bounded by Eqs. (26) and (22).
    #
    # The bounds in the data file are the physical limits of the equipment, but
    # (21) and (24) imply much narrower ones: a station moving heat X across a
    # temperature difference that (18) confines to [dTlo, dThi] can only be
    # running a flow between X/(c*dThi) and X/(c*dTlo). Since the heat load is
    # known, that pins the heat exchanger flow to a narrow band. These bounds
    # are valid for the original model too - they are implied by its own
    # constraints - but they matter most for the relaxation, whose envelopes
    # are only as tight as the box they are built on.
    dTlo(n) = p[:TSmin][n] - p[:TRmax][n]
    dThi(n) = p[:TSmax][n] - p[:TRmin][n]

    mfHESlo = Dict((i, t) => max(p[:mfHESmin][i], LH[i, t] / (c * dThi(HESnode[i]))) for i in IHES, t in T)
    mfHEShi = Dict((i, t) => min(p[:mfHESmax][i], LH[i, t] / (c * dTlo(HESnode[i]))) for i in IHES, t in T)
    mfHShi = Dict(j => min(p[:mfHSmax][j], p[:Qmax][j] / (c * dTlo(HSnode[j]))) for j in IHS)

    mfHS = m.ext[:variables][:mfHS] = @variable(m, [j=IHS, t=T],
        lower_bound=p[:mfHSmin][j], upper_bound=mfHShi[j], base_name="mfHS")
    mfHES = m.ext[:variables][:mfHES] = @variable(m, [i=IHES, t=T],
        lower_bound=mfHESlo[(i, t)], upper_bound=mfHEShi[(i, t)], base_name="mfHES")

    # Heat production of the heat stations, bounded by Eq. (25)
    Q = m.ext[:variables][:Q] = @variable(m, [j=IHS, t=T],
        lower_bound=0.0, upper_bound=p[:Qmax][j], base_name="Q")

    # Electricity production, bounded by Eq. (33)
    P = m.ext[:variables][:P] = @variable(m, [g=IE, t=T],
        lower_bound=0.0, upper_bound=p[:Pmax][g], base_name="P")

    # Electricity consumption of the heat pumps and of the water pumps
    LHP = m.ext[:variables][:LHP] = @variable(m, [j=IHP, t=T],
        lower_bound=0.0, upper_bound=p[:LHPmax][j], base_name="LHP")
    Lpump = m.ext[:variables][:Lpump] = @variable(m, [j=IHS, t=T],
        lower_bound=0.0, base_name="Lpump")

    # Voltage angles
    theta = m.ext[:variables][:theta] = @variable(m, [b=IB, t=T], base_name="theta")

    ##### Objective, Eq. (35)
    # Marginal cost of the electricity-only generators, plus the fuel cost of
    # the CHPs from Eq. (34). Heat pumps carry no direct cost: they pay for
    # their heat through the electricity they consume.
    m.ext[:objective] = @objective(m, Min,
        sum(p[:alphaE][g] * P[g, t] for g in IE, t in T)
        + sum(p[:alphaF][j] * (p[:fuelE][j] * P[j, t] + p[:fuelH][j] * Q[j, t]) for j in ICHP, t in T)
    )

    ##### Constraints
    # Eq. (2) - pressure loss due to friction. Quadratic equality, non-convex.
    m.ext[:constraints][:eq2S] = @constraint(m, [pp=IP, t=T],
        prS[pfrom[pp], t] - prS[pto[pp], t] == phi[pp] * mfS[pp, t]^2)
    m.ext[:constraints][:eq2R] = @constraint(m, [pp=IP, t=T],
        prR[pto[pp], t] - prR[pfrom[pp], t] == phi[pp] * mfR[pp, t]^2)

    # Eqs. (3)-(4) - inlet temperature of a pipe is the temperature of the node
    # it leaves from (supply), respectively arrives at (return).
    m.ext[:constraints][:eq3] = @constraint(m, [pp=IP, t=T], TSin[pp, t] == TS[pfrom[pp], t])
    m.ext[:constraints][:eq4] = @constraint(m, [pp=IP, t=T], TRin[pp, t] == TR[pto[pp], t])

    # Eq. (5) with zero time delay - the water crosses the pipe within one time
    # step and only loses heat. See NOTES.md for what this drops.
    m.ext[:constraints][:eq5S] = @constraint(m, [pp=IP, t=T],
        TSout[pp, t] == TSin[pp, t] * (1 - gamma[pp]))
    m.ext[:constraints][:eq5R] = @constraint(m, [pp=IP, t=T],
        TRout[pp, t] == TRin[pp, t] * (1 - gamma[pp]))

    # Eq. (14) - nodal mass balance. Supply: what arrives through the pipes and
    # what the heat stations inject equals what leaves through the pipes and
    # what the heat exchanger stations draw. Return is the mirror image.
    m.ext[:constraints][:eq14S] = @constraint(m, [n=IN, t=T],
        sum(mfS[pp, t] for pp in SPminus[n]) + sum(mfHS[j, t] for j in SHS[n])
        == sum(mfS[pp, t] for pp in SPplus[n]) + sum(mfHES[i, t] for i in SHES[n]))
    m.ext[:constraints][:eq14R] = @constraint(m, [n=IN, t=T],
        sum(mfR[pp, t] for pp in SPplus[n]) + sum(mfHES[i, t] for i in SHES[n])
        == sum(mfR[pp, t] for pp in SPminus[n]) + sum(mfHS[j, t] for j in SHS[n]))

    # Eqs. (16)-(17) - the nodal temperature is the mass-flow weighted mix of
    # the outlet temperatures of the arriving pipes. Bilinear equality.
    # Only needed where two or more pipes meet.
    m.ext[:constraints][:eq16] = @constraint(m, [n=INmixS, t=T],
        TS[n, t] * sum(mfS[pp, t] for pp in SPminus[n])
        == sum(mfS[pp, t] * TSout[pp, t] for pp in SPminus[n]))
    m.ext[:constraints][:eq17] = @constraint(m, [n=INmixR, t=T],
        TR[n, t] * sum(mfR[pp, t] for pp in SPplus[n])
        == sum(mfR[pp, t] * TRout[pp, t] for pp in SPplus[n]))

    # Eqs. (19)-(20) - at a node with a single arriving pipe the mix above is
    # exact and linear, so we use it instead of the bilinear form.
    m.ext[:constraints][:eq19] = @constraint(m, [n=INminus, t=T],
        TS[n, t] == TSout[only(SPminus[n]), t])
    m.ext[:constraints][:eq20] = @constraint(m, [n=INplus, t=T],
        TR[n, t] == TRout[only(SPplus[n]), t])

    # Eq. (21) - heat exchanger stations are the heat loads. Bilinear equality.
    m.ext[:constraints][:eq21] = @constraint(m, [i=IHES, t=T],
        LH[i, t] == c * mfHES[i, t] * (TS[HESnode[i], t] - TR[HESnode[i], t]))

    # Eq. (23) - the pressure difference at a HES has to exceed a threshold so
    # that water actually flows through it.
    m.ext[:constraints][:eq23] = @constraint(m, [i=IHES, t=T],
        prS[HESnode[i], t] - prR[HESnode[i], t] >= p[:dprHESmin][i])

    # Eq. (24) - heat output of a heat station. Bilinear equality.
    m.ext[:constraints][:eq24] = @constraint(m, [j=IHS, t=T],
        Q[j, t] == c * mfHS[j, t] * (TS[HSnode[j], t] - TR[HSnode[j], t]))

    # Eq. (27) - the water pump at a heat station works against the pressure
    # difference between the supply and the return network. Bilinear equality.
    m.ext[:constraints][:eq27] = @constraint(m, [j=IHS, t=T],
        Lpump[j, t] == p[:Pa_to_MW] / (rho * p[:etaPump][j]) *
                       mfHS[j, t] * (prS[HSnode[j], t] - prR[HSnode[j], t]))

    # Eqs. (28)-(29) - feasible operating region of an extraction CHP.
    m.ext[:constraints][:eq28] = @constraint(m, [j=ICHP, t=T], P[j, t] >= p[:r][j] * Q[j, t])
    m.ext[:constraints][:eq29] = @constraint(m, [j=ICHP, t=T],
        p[:fuelE][j] * P[j, t] + p[:fuelH][j] * Q[j, t] <= p[:Fmax][j])

    # Eq. (30) - heat pumps.
    m.ext[:constraints][:eq30] = @constraint(m, [j=IHP, t=T], Q[j, t] == COP[j, t] * LHP[j, t])

    # Eq. (31) - nodal power balance of the linearised power flow.
    m.ext[:constraints][:eq31] = @constraint(m, [b=IB, t=T],
        LE[b, t] + sum(LHP[j, t] for j in SHP[b]) + sum(Lpump[j, t] for j in SHSbus[b])
        == sum(P[g, t] for g in SE[b]) + sum(p[:B][(b, mm)] * (theta[mm, t] - theta[b, t]) for mm in SB[b]))

    # Eq. (32) - line flow limits.
    m.ext[:constraints][:eq32] = @constraint(m, [b=IB, mm=SB[b], t=T],
        -p[:fmax][(b, mm)] <= p[:B][(b, mm)] * (theta[mm, t] - theta[b, t]) <= p[:fmax][(b, mm)])

    # Not in the paper, but a linearised power flow needs an angle reference.
    m.ext[:constraints][:refangle] = @constraint(m, [t=T], theta[p[:refbus], t] == 0)

    # Wind is capped by what is available in that hour, not by its capacity.
    m.ext[:constraints][:wind] = @constraint(m, [t=T], P["W1", t] <= Wavail[t])

    return m
end
