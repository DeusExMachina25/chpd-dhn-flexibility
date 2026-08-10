## Sets, time series and parameters for the Combined Heat and Power Dispatch
# Reference: L. Mitridati and J. A. Taylor, "Power systems flexibility from
# district heating networks", PSCC 2018. Equation numbers below refer to that
# paper.

# Step 2a: create sets
function define_sets!(m::Model, data::Dict, ts::DataFrame; nsteps::Int=24)
    m.ext[:sets] = Dict()

    # Time steps
    T = m.ext[:sets][:T] = 1:nsteps

    # District heating network
    IN = m.ext[:sets][:IN] = [n["id"] for n in data["dhnNodes"]]
    IP = m.ext[:sets][:IP] = [p["id"] for p in data["pipes"]]

    # Heat stations, split by technology. A CHP is both a heat station and an
    # electricity generator, which is exactly where the coupling comes from.
    IHS = m.ext[:sets][:IHS] = [j["id"] for j in data["heatStations"]]
    ICHP = m.ext[:sets][:ICHP] = [j["id"] for j in data["heatStations"] if j["type"] == "chp"]
    IHP = m.ext[:sets][:IHP] = [j["id"] for j in data["heatStations"] if j["type"] == "hp"]

    # Heat exchanger stations = the heat loads
    IHES = m.ext[:sets][:IHES] = [i["id"] for i in data["heatExchangerStations"]]

    # Electricity system
    IE = m.ext[:sets][:IE] = [g["id"] for g in data["electricityGenerators"]]
    IB = m.ext[:sets][:IB] = data["buses"]

    # Pipe end points, in the supply direction
    pfrom = m.ext[:sets][:pfrom] = Dict(p["id"] => p["from"] for p in data["pipes"])
    pto = m.ext[:sets][:pto] = Dict(p["id"] => p["to"] for p in data["pipes"])

    # SP+_n / SP-_n : pipes starting / ending at node n (supply direction)
    m.ext[:sets][:SPplus] = Dict(n => [p for p in IP if pfrom[p] == n] for n in IN)
    m.ext[:sets][:SPminus] = Dict(n => [p for p in IP if pto[p] == n] for n in IN)

    # Stations connected to each DHN node
    m.ext[:sets][:SHS] = Dict(n => [j["id"] for j in data["heatStations"] if j["dhnNode"] == n] for n in IN)
    m.ext[:sets][:SHES] = Dict(n => [i["id"] for i in data["heatExchangerStations"] if i["dhnNode"] == n] for n in IN)

    # IN- / IN+ : nodes with a single arriving / departing pipe. Equations (19)
    # and (20) hold exactly there, so we can skip the bilinear mixing (16)/(17).
    m.ext[:sets][:INminus] = [n for n in IN if length(m.ext[:sets][:SPminus][n]) == 1]
    m.ext[:sets][:INplus] = [n for n in IN if length(m.ext[:sets][:SPplus][n]) == 1]

    # Nodes where two or more pipes actually mix. Only these need the bilinear
    # constraints (16)/(17); a node with no arriving pipe has a free
    # temperature, set by the heat station or the heat exchanger sitting on it.
    m.ext[:sets][:INmixS] = [n for n in IN if length(m.ext[:sets][:SPminus][n]) >= 2]
    m.ext[:sets][:INmixR] = [n for n in IN if length(m.ext[:sets][:SPplus][n]) >= 2]

    # Units connected to each electricity bus
    m.ext[:sets][:SE] = Dict(b => [g["id"] for g in data["electricityGenerators"] if g["bus"] == b] for b in IB)
    m.ext[:sets][:SHP] = Dict(b => [j["id"] for j in data["heatStations"] if j["type"] == "hp" && j["bus"] == b] for b in IB)
    m.ext[:sets][:SHSbus] = Dict(b => [j["id"] for j in data["heatStations"] if j["bus"] == b] for b in IB)

    # Lines, as an ordered list of (from, to) pairs, and the buses connected to
    # each bus (S^B_n in the paper)
    L = m.ext[:sets][:L] = [(l["from"], l["to"]) for l in data["lines"]]
    m.ext[:sets][:SB] = Dict(b => vcat([n for (f, t) in L for n in (t,) if f == b],
                                       [n for (f, t) in L for n in (f,) if t == b]) for b in IB)

    return m
end

# Step 2b: add time series
function process_time_series_data!(m::Model, data::Dict, ts::DataFrame)
    T = m.ext[:sets][:T]
    IB = m.ext[:sets][:IB]
    IHES = m.ext[:sets][:IHES]
    IHP = m.ext[:sets][:IHP]

    m.ext[:timeseries] = Dict()

    # Electricity load, split over the load buses using the shares in the yaml
    share = Dict(b => 0.0 for b in IB)
    for l in data["electricityLoads"]
        share[l["bus"]] += l["share"]
    end
    m.ext[:timeseries][:LE] = Dict((b, t) => share[b] * ts.ElectricityLoad[t] for b in IB, t in T)

    # Available wind power [MW]
    Pw = only([g["Pmax"] for g in data["electricityGenerators"] if g["id"] == "W1"])
    m.ext[:timeseries][:Wavail] = Dict(t => Pw * ts.WindAF[t] for t in T)

    # Heat load at each heat exchanger station [MW]. There is a single HES in
    # this case study, so it carries the whole profile.
    m.ext[:timeseries][:LH] = Dict((i, t) => ts.HeatLoad[t] / length(IHES) for i in IHES, t in T)

    # Coefficient of performance of the heat pumps, Eq. (30)
    m.ext[:timeseries][:COP] = Dict((j, t) => ts.COP[t] for j in IHP, t in T)

    return m
end

# Step 2c: process input parameters
function process_parameters!(m::Model, data::Dict)
    IN = m.ext[:sets][:IN]
    IP = m.ext[:sets][:IP]
    IHS = m.ext[:sets][:IHS]
    IHES = m.ext[:sets][:IHES]
    ICHP = m.ext[:sets][:ICHP]
    IHP = m.ext[:sets][:IHP]
    IE = m.ext[:sets][:IE]

    m.ext[:parameters] = Dict()
    p = m.ext[:parameters]

    # Physical constants
    cst = data["constants"]
    p[:rho] = cst["rho"]
    p[:c] = cst["c"]                # [MJ/(kg.K)], so Q[MW] = c * mf * dT
    p[:dt] = cst["dt"]              # [s]
    p[:Pa_to_MW] = cst["Pa_to_MW"]

    # DHN nodes: temperature and pressure bounds, Eqs. (15) and (18)
    nodes = Dict(n["id"] => n for n in data["dhnNodes"])
    for (key, field) in [(:TSmin, "TSmin"), (:TSmax, "TSmax"), (:TRmin, "TRmin"), (:TRmax, "TRmax"),
                         (:prSmin, "prSmin"), (:prSmax, "prSmax"), (:prRmin, "prRmin"), (:prRmax, "prRmax")]
        p[key] = Dict(n => nodes[n][field] for n in IN)
    end

    # Pipes
    pipes = Dict(pp["id"] => pp for pp in data["pipes"])
    for (key, field) in [(:R, "R"), (:L, "L"), (:mu, "mu"), (:phi, "phi"),
                         (:mfSmin, "mfSmin"), (:mfSmax, "mfSmax"), (:mfRmin, "mfRmin"), (:mfRmax, "mfRmax")]
        p[key] = Dict(pp => pipes[pp][field] for pp in IP)
    end

    # Per-time-step thermal loss factor of a pipe, i.e. the bracket in Eq. (5)
    # evaluated for a travel time of sigma steps. gamma[p] is the loss for one
    # step of travel; the factor for sigma steps is (1 - sigma * gamma[p]).
    p[:gamma] = Dict(pp => 2 * p[:mu][pp] * p[:dt] / (p[:rho] * 4182.0 * p[:R][pp]) for pp in IP)

    # Maximum time delay of a pipe, from the pipe volume and the minimum flow.
    # This is the tau-bar in Eq. (6); it is only used when delays are enabled.
    p[:taumax] = Dict(pp => Int(ceil(p[:rho] * pi * p[:R][pp]^2 * p[:L][pp] /
                                     (p[:mfSmin][pp] * p[:dt]))) for pp in IP)

    # Heat stations, Eqs. (24)-(27)
    hs = Dict(j["id"] => j for j in data["heatStations"])
    p[:Qmax] = Dict(j => hs[j]["Qmax"] for j in IHS)
    p[:mfHSmin] = Dict(j => hs[j]["mfmin"] for j in IHS)
    p[:mfHSmax] = Dict(j => hs[j]["mfmax"] for j in IHS)
    p[:etaPump] = Dict(j => hs[j]["etaPump"] for j in IHS)
    p[:HSnode] = Dict(j => hs[j]["dhnNode"] for j in IHS)
    p[:HSbus] = Dict(j => hs[j]["bus"] for j in IHS)

    # Heat exchanger stations, Eqs. (21)-(23)
    hes = Dict(i["id"] => i for i in data["heatExchangerStations"])
    p[:mfHESmin] = Dict(i => hes[i]["mfmin"] for i in IHES)
    p[:mfHESmax] = Dict(i => hes[i]["mfmax"] for i in IHES)
    p[:dprHESmin] = Dict(i => hes[i]["dprmin"] for i in IHES)
    p[:HESnode] = Dict(i => hes[i]["dhnNode"] for i in IHES)

    # Electricity generators, Eq. (33) and the objective (35)
    gen = Dict(g["id"] => g for g in data["electricityGenerators"])
    p[:Pmax] = Dict(g => gen[g]["Pmax"] for g in IE)
    p[:alphaE] = Dict(g => gen[g]["cost"] for g in IE)
    p[:Ebus] = Dict(g => gen[g]["bus"] for g in IE)

    # Extraction CHP, Eqs. (28), (29), (34)
    chp = Dict(j["id"] => j for j in data["chpUnits"])
    p[:r] = Dict(j => chp[j]["r"] for j in ICHP)
    p[:fuelE] = Dict(j => chp[j]["fuelE"] for j in ICHP)
    p[:fuelH] = Dict(j => chp[j]["fuelH"] for j in ICHP)
    p[:Fmax] = Dict(j => chp[j]["Fmax"] for j in ICHP)
    p[:alphaF] = Dict(j => chp[j]["fuelCost"] for j in ICHP)

    # Heat pumps, Eq. (30)
    hp = Dict(j["id"] => j for j in data["heatPumps"])
    p[:LHPmax] = Dict(j => hp[j]["LHPmax"] for j in IHP)

    # Lines, Eq. (32)
    p[:B] = Dict((l["from"], l["to"]) => l["B"] for l in data["lines"])
    p[:fmax] = Dict((l["from"], l["to"]) => l["fmax"] for l in data["lines"])
    # susceptance is symmetric, add the reverse direction for convenience
    for l in data["lines"]
        p[:B][(l["to"], l["from"])] = l["B"]
        p[:fmax][(l["to"], l["from"])] = l["fmax"]
    end
    p[:refbus] = data["referenceBus"]

    return m
end
