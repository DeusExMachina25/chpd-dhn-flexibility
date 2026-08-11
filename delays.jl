## Eqs. (6)-(13): pipeline time delays, and the energy storage they carry
#
# This is where the flexibility in the paper actually comes from. Water pushed
# into a pipe now arrives at the far end tau steps later, so heat produced in a
# cheap hour can be consumed in an expensive one. tau itself is a decision, and
# it depends on the mass flow, which is what makes the model non-convex before
# the big-M reformulation.
#
# The chain is:
#   (6)      tau is the smallest number of steps whose accumulated mass fills
#            the pipe, rho*pi*R^2*L
#   (7),(8)  that "smallest sigma such that" is turned into binaries u and a
#            linear expression for tau
#   (9)-(12) binaries v select the single sigma equal to tau, and Ttilde picks
#            up the correspondingly delayed and cooled inlet temperature
#   (13)     the outlet temperature is the sum over sigma of Ttilde
#
# add_delays! is called once for the supply network and once for the return
# network; they are structurally identical, only the variables differ.

# History before the start of the horizon. The pipe is not empty at t = 1, so
# Eq. (7) needs mass flows and Eq. (10) needs inlet temperatures for the steps
# before it. The authors do the same: they hold the flow at its minimum and the
# temperature at the cold end of the band.
function delay_history(m::Model, network::Symbol)
    p = m.ext[:parameters]
    IP = m.ext[:sets][:IP]
    mfhist = Dict{Tuple{Int,Int},Float64}()
    Thist = Dict{Tuple{Int,Int},Float64}()
    for pp in IP, k in (1 - p[:taumax][pp]):0
        mfhist[(pp, k)] = network === :supply ? p[:mfSmin][pp] : p[:mfRmin][pp]
        Thist[(pp, k)] = network === :supply ? p[:TSlo] : p[:TRhi]
    end
    return mfhist, Thist
end

function add_delays!(m::Model, network::Symbol)
    T = m.ext[:sets][:T]
    IP = m.ext[:sets][:IP]
    p = m.ext[:parameters]
    dt = p[:dt]
    gamma = p[:gamma]
    taumax = p[:taumax]

    sup = network === :supply
    mf = sup ? m.ext[:variables][:mfS] : m.ext[:variables][:mfR]
    Tin = sup ? m.ext[:variables][:TSin] : m.ext[:variables][:TRin]
    Tout = sup ? m.ext[:variables][:TSout] : m.ext[:variables][:TRout]
    tag = sup ? "S" : "R"

    mfhist, Thist = delay_history(m, network)

    # Volume of each pipe expressed as a mass, the right hand side of Eq. (6)
    V = Dict(pp => p[:rho] * pi * p[:R][pp]^2 * p[:L][pp] for pp in IP)

    # Accumulated mass over the last sigma+1 steps, mixing history and
    # variables. This is the left hand side of Eqs. (6) and (7).
    accum(pp, t, sigma) = dt * sum(t - k >= 1 ? mf[pp, t - k] : mfhist[(pp, t - k)]
                                   for k in 0:sigma)

    # Index set: one entry per pipe, time step and candidate delay.
    idx = [(pp, t, s) for pp in IP for t in T for s in 0:taumax[pp]]

    ## Eq. (7): u = 1 exactly when sigma steps of flow have filled the pipe.
    # Big-M from the widest the residual can get, computed per pipe rather than
    # picked out of the air - a loose M here is the classic way to get a model
    # that solves fast and answers the wrong question.
    Mmass = Dict(pp => max(V[pp], (taumax[pp] + 1) * p[sup ? :mfSmax : :mfRmax][pp] * dt - V[pp]) * 1.01
                 for pp in IP)

    u = m.ext[:variables][Symbol("u$tag")] =
        @variable(m, [k=idx], binary=true, base_name="u$tag")

    m.ext[:constraints][Symbol("eq7$(tag)a")] = @constraint(m, [k=idx],
        accum(k[1], k[2], k[3]) - V[k[1]] <= Mmass[k[1]] * u[k])
    m.ext[:constraints][Symbol("eq7$(tag)b")] = @constraint(m, [k=idx],
        accum(k[1], k[2], k[3]) - V[k[1]] >= Mmass[k[1]] * (u[k] - 1))

    # Accumulated mass only grows with sigma, so u must be monotone. Not in the
    # paper, but implied by (7), and it cuts the search space a lot.
    m.ext[:constraints][Symbol("eq7$(tag)mono")] = @constraint(m,
        [pp=IP, t=T, s=1:taumax[pp]], u[(pp, t, s)] >= u[(pp, t, s - 1)])

    ## Eq. (8): tau is the number of sigma values that did NOT yet fill the
    # pipe. With u monotone 0,0,...,0,1,1,...,1 this picks the switch point.
    tau = m.ext[:variables][Symbol("tau$tag")] =
        @variable(m, [pp=IP, t=T], lower_bound=0, upper_bound=taumax[pp], base_name="tau$tag")

    m.ext[:constraints][Symbol("eq8$tag")] = @constraint(m, [pp=IP, t=T],
        tau[pp, t] == taumax[pp] - sum(u[(pp, t, s)] for s in 0:taumax[pp]) + 1)

    ## Eqs. (9)-(12): v selects the single sigma that equals tau, and Ttilde
    # carries the delayed, cooled inlet temperature for that sigma alone.
    Thi = sup ? p[:TShi] : p[:TRhi]
    MT = Thi * 1.01

    v = m.ext[:variables][Symbol("v$tag")] =
        @variable(m, [k=idx], binary=true, base_name="v$tag")
    Ttil = m.ext[:variables][Symbol("Ttil$tag")] =
        @variable(m, [k=idx], lower_bound=0.0, upper_bound=Thi, base_name="Ttil$tag")

    # Eq. (9): Ttilde is zero unless its sigma is the selected one.
    m.ext[:constraints][Symbol("eq9$tag")] = @constraint(m, [k=idx], Ttil[k] <= MT * v[k])

    # Eq. (10): when selected, Ttilde equals the inlet temperature sigma steps
    # ago times the loss over that travel time, the bracket of Eq. (5).
    delayed(pp, t, s) = t - s >= 1 ? Tin[pp, t - s] : Thist[(pp, t - s)]

    m.ext[:constraints][Symbol("eq10$(tag)a")] = @constraint(m, [k=idx],
        Ttil[k] - delayed(k[1], k[2], k[3]) * (1 - k[3] * gamma[k[1]]) <= MT * (1 - v[k]))
    m.ext[:constraints][Symbol("eq10$(tag)b")] = @constraint(m, [k=idx],
        Ttil[k] - delayed(k[1], k[2], k[3]) * (1 - k[3] * gamma[k[1]]) >= MT * (v[k] - 1))

    # Eq. (11): v can only be on for the sigma that equals tau.
    m.ext[:constraints][Symbol("eq11$(tag)a")] = @constraint(m, [k=idx],
        k[3] - tau[k[1], k[2]] <= taumax[k[1]] * (1 - v[k]))
    m.ext[:constraints][Symbol("eq11$(tag)b")] = @constraint(m, [k=idx],
        k[3] - tau[k[1], k[2]] >= taumax[k[1]] * (v[k] - 1))

    # Eq. (12): exactly one sigma is selected.
    m.ext[:constraints][Symbol("eq12$tag")] = @constraint(m, [pp=IP, t=T],
        sum(v[(pp, t, s)] for s in 0:taumax[pp]) == 1)

    ## Eq. (13): the outlet temperature is the selected Ttilde.
    m.ext[:constraints][Symbol("eq13$tag")] = @constraint(m, [pp=IP, t=T],
        Tout[pp, t] == sum(Ttil[(pp, t, s)] for s in 0:taumax[pp]))

    return m
end
