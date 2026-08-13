## Section IV-A: the Conventional Economic Dispatch used as a benchmark
# The district heating network is not modelled at all. Everything in Eqs.
# (1)-(27) disappears and is replaced by a single heat balance (38); the power
# balance (31) loses its water pump term and becomes (39). The CHP and heat
# pump characteristics and the power system constraints are untouched.
#
# With no network, heat produced in an hour must be consumed in that hour, so
# this model has no storage. That is the flexibility the CHPD is meant to add.

function build_ced!(m::Model)
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()

    # Extract sets
    T = m.ext[:sets][:T]
    IHS = m.ext[:sets][:IHS]
    IHES = m.ext[:sets][:IHES]
    ICHP = m.ext[:sets][:ICHP]
    IHP = m.ext[:sets][:IHP]
    IE = m.ext[:sets][:IE]
    IB = m.ext[:sets][:IB]
    SE = m.ext[:sets][:SE]
    SHP = m.ext[:sets][:SHP]
    SB = m.ext[:sets][:SB]

    # Extract time series
    LE = m.ext[:timeseries][:LE]
    LH = m.ext[:timeseries][:LH]
    COP = m.ext[:timeseries][:COP]
    Wavail = m.ext[:timeseries][:Wavail]

    # Extract parameters
    p = m.ext[:parameters]

    ##### Create variables
    Q = m.ext[:variables][:Q] = @variable(m, [j=IHS, t=T],
        lower_bound=0.0, upper_bound=p[:Qmax][j], base_name="Q")
    P = m.ext[:variables][:P] = @variable(m, [g=IE, t=T],
        lower_bound=0.0, upper_bound=p[:Pmax][g], base_name="P")
    LHP = m.ext[:variables][:LHP] = @variable(m, [j=IHP, t=T],
        lower_bound=0.0, upper_bound=p[:LHPmax][j], base_name="LHP")
    theta = m.ext[:variables][:theta] = @variable(m, [b=IB, t=T], base_name="theta")

    ##### Objective, Eq. (35)
    m.ext[:objective] = @objective(m, Min,
        sum(p[:alphaE][g] * P[g, t] for g in IE, t in T)
        + sum(p[:alphaF][j] * (p[:fuelE][j] * P[j, t] + p[:fuelH][j] * Q[j, t]) for j in ICHP, t in T)
    )

    ##### Constraints
    # Eq. (38) - heat balance, replacing the whole heating network.
    m.ext[:constraints][:eq38] = @constraint(m, [t=T],
        sum(LH[i, t] for i in IHES) == sum(Q[j, t] for j in IHS))

    # Eqs. (28)-(29) - extraction CHP, unchanged.
    m.ext[:constraints][:eq28] = @constraint(m, [j=ICHP, t=T], P[j, t] >= p[:r][j] * Q[j, t])
    m.ext[:constraints][:eq29] = @constraint(m, [j=ICHP, t=T],
        p[:fuelE][j] * P[j, t] + p[:fuelH][j] * Q[j, t] <= p[:Fmax][j])

    # Eq. (30) - heat pumps, unchanged.
    m.ext[:constraints][:eq30] = @constraint(m, [j=IHP, t=T], Q[j, t] == COP[j, t] * LHP[j, t])

    # Eq. (39) - power balance without the water pumps.
    m.ext[:constraints][:eq39] = @constraint(m, [b=IB, t=T],
        LE[b, t] + sum(LHP[j, t] for j in SHP[b])
        == sum(P[g, t] for g in SE[b]) + sum(p[:B][(b, mm)] * (theta[mm, t] - theta[b, t]) for mm in SB[b]))

    # Eqs. (32)-(33) - power system, unchanged.
    m.ext[:constraints][:eq32] = @constraint(m, [b=IB, mm=SB[b], t=T],
        -p[:fmax][(b, mm)] <= p[:B][(b, mm)] * (theta[mm, t] - theta[b, t]) <= p[:fmax][(b, mm)])

    m.ext[:constraints][:refangle] = @constraint(m, [t=T], theta[p[:refbus], t] == 0)
    m.ext[:constraints][:wind] = @constraint(m, [t=T], P["W1", t] <= Wavail[t])

    return m
end
