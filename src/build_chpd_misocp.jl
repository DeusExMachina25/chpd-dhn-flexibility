## Section III-B: convex relaxation of the CHPD
# Two things happen here, both on top of the Section II model:
#
#   III-A  the bilinear equalities (16), (17), (21), (24) and (27) are replaced
#          by McCormick envelopes, which is where the auxiliary variables come
#          from;
#   III-B  the quadratic equality (2) is replaced by the convex quadratic
#          inequality (36).
#
# With the time delays switched off there are no binaries left, so the result
# is a convex QCP. With delays on it is the MISOCP of the paper.

# The four McCormick inequalities for w = x*y over [xlo,xhi] x [ylo,yhi].
function mccormick!(m::Model, w, x, y, xlo, xhi, ylo, yhi)
    @constraint(m, w >= xlo * y + x * ylo - xlo * ylo)
    @constraint(m, w >= xhi * y + x * yhi - xhi * yhi)
    @constraint(m, w <= xhi * y + x * ylo - xhi * ylo)
    @constraint(m, w <= xlo * y + x * yhi - xlo * yhi)
    return w
end

# Remove a constraint block that was added by build_chpd_minlp!.
function drop_constraint!(m::Model, key::Symbol)
    haskey(m.ext[:constraints], key) || return m
    for con in m.ext[:constraints][key]
        delete(m, con)
    end
    delete!(m.ext[:constraints], key)
    return m
end

function build_chpd_misocp!(m::Model; delays::Bool=false, widen::Float64=0.0,
                            mfHScap::Union{Nothing,Float64}=nothing)
    # Build the original model first, then relax it. `widen` is the analysis
    # knob documented in build_chpd_minlp!; it defaults to no change.
    build_chpd_minlp!(m; delays=delays, widen=widen, mfHScap=mfHScap)

    # Extract sets
    T = m.ext[:sets][:T]
    IN = m.ext[:sets][:IN]
    IP = m.ext[:sets][:IP]
    IHS = m.ext[:sets][:IHS]
    IHES = m.ext[:sets][:IHES]
    pfrom = m.ext[:sets][:pfrom]
    pto = m.ext[:sets][:pto]
    SPplus = m.ext[:sets][:SPplus]
    SPminus = m.ext[:sets][:SPminus]
    INmixS = m.ext[:sets][:INmixS]
    INmixR = m.ext[:sets][:INmixR]

    # Extract time series
    LH = m.ext[:timeseries][:LH]

    # Extract parameters
    p = m.ext[:parameters]
    c = p[:c]
    rho = p[:rho]
    phi = p[:phi]
    gamma = p[:gamma]
    HSnode = p[:HSnode]
    HESnode = p[:HESnode]

    # Extract variables
    mfS = m.ext[:variables][:mfS]
    mfR = m.ext[:variables][:mfR]
    prS = m.ext[:variables][:prS]
    prR = m.ext[:variables][:prR]
    TS = m.ext[:variables][:TS]
    TR = m.ext[:variables][:TR]
    TSout = m.ext[:variables][:TSout]
    TRout = m.ext[:variables][:TRout]
    mfHS = m.ext[:variables][:mfHS]
    mfHES = m.ext[:variables][:mfHES]
    Q = m.ext[:variables][:Q]
    Lpump = m.ext[:variables][:Lpump]

    # Bounds of the pipe outlet temperatures, needed by the envelopes. These
    # are set by build_chpd_minlp!, which knows whether delays are on.
    TSoutlo = p[:TSoutlo]; TSouthi = p[:TShi]
    TRoutlo = p[:TRoutlo]; TRouthi = p[:TRhi]

    ## Eq. (2) -> Eq. (36): convex quadratic relaxation of the pressure loss.
    # Written with the square on the small side so that the constraint is
    # convex; flipping this would silently give back a non-convex problem.
    drop_constraint!(m, :eq2S)
    drop_constraint!(m, :eq2R)
    m.ext[:constraints][:eq36S] = @constraint(m, [pp=IP, t=T],
        phi[pp] * mfS[pp, t]^2 <= prS[pfrom[pp], t] - prS[pto[pp], t])
    m.ext[:constraints][:eq36R] = @constraint(m, [pp=IP, t=T],
        phi[pp] * mfR[pp, t]^2 <= prR[pto[pp], t] - prR[pfrom[pp], t])

    ## Eq. (21) -> McCormick. The product is mass flow times temperature
    # difference. We envelope the difference as a single variable rather than
    # the two products separately, which gives a tighter relaxation.
    drop_constraint!(m, :eq21)
    dTHES = m.ext[:variables][:dTHES] = @variable(m, [i=IHES, t=T],
        lower_bound=p[:TSmin][HESnode[i]] - p[:TRmax][HESnode[i]],
        upper_bound=p[:TSmax][HESnode[i]] - p[:TRmin][HESnode[i]], base_name="dTHES")
    @constraint(m, [i=IHES, t=T], dTHES[i, t] == TS[HESnode[i], t] - TR[HESnode[i], t])

    # The envelope is built on whatever bounds the variables actually carry,
    # which are the tightened ones computed in build_chpd_minlp!, not the raw
    # equipment limits from the data file.
    wHES = m.ext[:variables][:wHES] = @variable(m, [i=IHES, t=T], base_name="wHES")
    for i in IHES, t in T
        mccormick!(m, wHES[i, t], mfHES[i, t], dTHES[i, t],
                   lower_bound(mfHES[i, t]), upper_bound(mfHES[i, t]),
                   lower_bound(dTHES[i, t]), upper_bound(dTHES[i, t]))
    end
    m.ext[:constraints][:eq21mc] = @constraint(m, [i=IHES, t=T], LH[i, t] == c * wHES[i, t])

    ## Eq. (24) -> McCormick, same structure at the heat stations.
    drop_constraint!(m, :eq24)
    dTHS = m.ext[:variables][:dTHS] = @variable(m, [j=IHS, t=T],
        lower_bound=p[:TSmin][HSnode[j]] - p[:TRmax][HSnode[j]],
        upper_bound=p[:TSmax][HSnode[j]] - p[:TRmin][HSnode[j]], base_name="dTHS")
    @constraint(m, [j=IHS, t=T], dTHS[j, t] == TS[HSnode[j], t] - TR[HSnode[j], t])

    wHS = m.ext[:variables][:wHS] = @variable(m, [j=IHS, t=T], base_name="wHS")
    for j in IHS, t in T
        mccormick!(m, wHS[j, t], mfHS[j, t], dTHS[j, t],
                   lower_bound(mfHS[j, t]), upper_bound(mfHS[j, t]),
                   lower_bound(dTHS[j, t]), upper_bound(dTHS[j, t]))
    end
    m.ext[:constraints][:eq24mc] = @constraint(m, [j=IHS, t=T], Q[j, t] == c * wHS[j, t])

    ## Eq. (27) -> McCormick. Product of mass flow and pressure difference.
    drop_constraint!(m, :eq27)
    dprHS = m.ext[:variables][:dprHS] = @variable(m, [j=IHS, t=T],
        lower_bound=p[:prSmin][HSnode[j]] - p[:prRmax][HSnode[j]],
        upper_bound=p[:prSmax][HSnode[j]] - p[:prRmin][HSnode[j]], base_name="dprHS")
    @constraint(m, [j=IHS, t=T], dprHS[j, t] == prS[HSnode[j], t] - prR[HSnode[j], t])

    wPump = m.ext[:variables][:wPump] = @variable(m, [j=IHS, t=T], base_name="wPump")
    for j in IHS, t in T
        mccormick!(m, wPump[j, t], mfHS[j, t], dprHS[j, t],
                   lower_bound(mfHS[j, t]), upper_bound(mfHS[j, t]),
                   lower_bound(dprHS[j, t]), upper_bound(dprHS[j, t]))
    end
    m.ext[:constraints][:eq27mc] = @constraint(m, [j=IHS, t=T],
        Lpump[j, t] == p[:pressure_to_MW] / (rho * p[:etaPump][j]) * wPump[j, t])

    ## Eqs. (16)-(17) -> McCormick. Two families of products per mixing node:
    # the nodal temperature times each arriving flow, and each arriving flow
    # times its own outlet temperature.
    drop_constraint!(m, :eq16)
    drop_constraint!(m, :eq17)

    mixS = [(n, pp) for n in INmixS for pp in SPminus[n]]
    mixR = [(n, pp) for n in INmixR for pp in SPplus[n]]

    if !isempty(mixS)
        wTSmf = m.ext[:variables][:wTSmf] = @variable(m, [k=mixS, t=T], base_name="wTSmf")
        wmfTSout = m.ext[:variables][:wmfTSout] = @variable(m, [k=mixS, t=T], base_name="wmfTSout")
        for (n, pp) in mixS, t in T
            mccormick!(m, wTSmf[(n, pp), t], TS[n, t], mfS[pp, t],
                       p[:TSmin][n], p[:TSmax][n], p[:mfSmin][pp], p[:mfSmax][pp])
            mccormick!(m, wmfTSout[(n, pp), t], mfS[pp, t], TSout[pp, t],
                       p[:mfSmin][pp], p[:mfSmax][pp], TSoutlo, TSouthi)
        end
        m.ext[:constraints][:eq16mc] = @constraint(m, [n=INmixS, t=T],
            sum(wTSmf[(n, pp), t] for pp in SPminus[n]) == sum(wmfTSout[(n, pp), t] for pp in SPminus[n]))
    end

    if !isempty(mixR)
        wTRmf = m.ext[:variables][:wTRmf] = @variable(m, [k=mixR, t=T], base_name="wTRmf")
        wmfTRout = m.ext[:variables][:wmfTRout] = @variable(m, [k=mixR, t=T], base_name="wmfTRout")
        for (n, pp) in mixR, t in T
            mccormick!(m, wTRmf[(n, pp), t], TR[n, t], mfR[pp, t],
                       p[:TRmin][n], p[:TRmax][n], p[:mfRmin][pp], p[:mfRmax][pp])
            mccormick!(m, wmfTRout[(n, pp), t], mfR[pp, t], TRout[pp, t],
                       p[:mfRmin][pp], p[:mfRmax][pp], TRoutlo, TRouthi)
        end
        m.ext[:constraints][:eq17mc] = @constraint(m, [n=INmixR, t=T],
            sum(wTRmf[(n, pp), t] for pp in SPplus[n]) == sum(wmfTRout[(n, pp), t] for pp in SPplus[n]))
    end

    return m
end
