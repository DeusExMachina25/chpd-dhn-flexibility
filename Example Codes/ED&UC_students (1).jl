## Practical Session II: ED & UC problems
# Author: Kenneth Bruninx, Sebastian Gonzato
# Last update: November 2, 2020

## Step 0: Activate environment - ensure consistency accross computers
using Pkg
Pkg.activate(@__DIR__) # @__DIR__ = directory this script is in
Pkg.instantiate() # If a Manifest.toml file exist in the current project, download all the packages declared in that manifest. Else, resolve a set of feasible packages from the Project.toml files and install them.

##  Step 1: input data
using CSV
using DataFrames
using YAML

data = YAML.load_file(joinpath(@__DIR__, "data.yaml"))
ts = CSV.read(joinpath(@__DIR__, "timeseries.csv"), DataFrame)
println("Done")

## Step 2: create model & pass data to model
using JuMP
using Gurobi
m = Model(optimizer_with_attributes(Gurobi.Optimizer))

# Step 2a: create sets
function define_sets!(m::Model, data::Dict, ts::DataFrame)
    m.ext[:sets] = Dict()

    # Time steps
    J = m.ext[:sets][:J] = 1:24

    # Dispatchable generators per type
    IDtype = m.ext[:sets][:IDtype] = [id for id in keys(data["dispatchableGenerators"])]
    # Dispatchable generators per unit
    ID = Array{Union{Nothing,String}}(nothing,0)
    for idtype in IDtype, i in 1:data["dispatchableGenerators"][idtype]["numberOfUnits"]
       ID = m.ext[:sets][:ID] = push!(ID,string(idtype,"_$(i)"))
    end

    # Variable generators
    IV = m.ext[:sets][:IV] = [iv for iv in keys(data["variableGenerators"])]

    # All variable and dispatchable generators, per unit
    I = m.ext[:sets][:I] = union(m.ext[:sets][:IV],m.ext[:sets][:ID])

    # return model
    return m
end

# Step 2b: add time series
function process_time_series_data!(m::Model, ts::DataFrame)
    # extract the relevant sets
    IV = m.ext[:sets][:IV]
    J = m.ext[:sets][:J]

    # create dictionary to store time series
    m.ext[:timeseries] = Dict()
    m.ext[:timeseries][:AF] = Dict()

    # add time series to dictionary
    m.ext[:timeseries][:D] = ts.Load[1:24]
    m.ext[:timeseries][:AF][IV[1]] = ts.LFW[1:24]
    m.ext[:timeseries][:AF][IV[2]] = ts.LFS[1:24]

    # return model
    return m
end

# step 2c: process input parameters
function process_parameters!(m::Model, data::Dict)
    # extract sets
    ID = m.ext[:sets][:ID]
    IV = m.ext[:sets][:IV]

    # Create parameter dictonary
    m.ext[:parameters] = Dict()

    # general parameters
    m.ext[:parameters][:VOLL] = data["valueOfLostLoad"] # VOLL

    # parameters of dispatchable generators per unit
    d = data["dispatchableGenerators"]
    m.ext[:parameters][:MUT] = Dict(i => d[SubString(i,1:length(i)-2)]["minUpTime"] for i in ID)

    # parameters of variable generators
    d = data["variableGenerators"]
    m.ext[:parameters][:GmaxV] = Dict(i => d[i]["installedCapacity"] for i in IV)

    # parameters storage

    # return model
    return m
end

# call functions
define_sets!(m, data, ts)
process_time_series_data!(m, ts)
process_parameters!(m, data)

## Step 3: construct your model
function build_basic_ED_model!(m::Model)
    # Create m.ext entries "variables", "expressions" and "constraints"
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()

    # Extract sets
    I = m.ext[:sets][:I]
    J = m.ext[:sets][:J]

    # Extract time series data

    # Extract parameters

    # create variables (g is already created for you, don't change this)
    g = m.ext[:variables][:g] = @variable(m, [i=I,j=J], lower_bound=0, base_name="generation")

    # Create affine expressions (= linear combinations of variables)

    # Objective
    obj = m.ext[:objective] = @objective(m, Min,
    #     ENTER OBJECTIVE HERE and remove 0 below
    0
    )

    # Constraints
    # Formulate constraints (market clearing constraint "con2a" is given - don't change this name, but add your energy not served variable and demand parameter)
    con2a = m.ext[:constraints][:con2a] = @constraint(m, [j=J],
        sum(g[i,j] for i in I) == 0
    )
end

function build_basic_ED_storage_model!(m::Model)
    # Build the basic model
    build_basic_ED_model!(m)

    # Extract sets

    # Extract time series

    # Extract parameters

    # Extract variables
    g = m.ext[:variables][:g]

    # create variables (state-of-charge "e", charging "c" and discharging "d" are given, don't change these names)
    c =  m.ext[:variables][:c] = @variable(m, [j=J], lower_bound=0, base_name="charging")
    d =  m.ext[:variables][:d] = @variable(m, [j=J], lower_bound=0, base_name="discharging")
    e =  m.ext[:variables][:e] = @variable(m, [j=J], lower_bound=0, base_name="SOC")

    # Remove market clearing constraint using "delete!"

    # Add new market clearing constraint

    # Add storage model with cyclic boundary conditions
    # m.ext[:constraints][:con4a] = @constraint(m, [j=J[2:end-1]],
    #     ENTER CONSTRAINT HERE
    # )
    # m.ext[:constraints][:con4b] = @constraint(m, [j=J[1]],
    #     ENTER CONSTRAINT HERE
    # )
    # m.ext[:constraints][:con4c] = @constraint(m, [j=J[end]],
    #     ENTER CONSTRAINT HERE
    # )

    return m
end

function build_basic_UC_model!(m::Model)
    # Build the basic model
    build_basic_ED_model!(m)

    # Extract sets

    # Extract parameters

    # Extract variables

    # Create variables (commitment status "z" and start-up "v" are given, don't change this)
    z = m.ext[:variables][:z] = @variable(m, [i=ID,j=J], binary=true, base_name="commitment")
    v = m.ext[:variables][:v] = @variable(m, [i=ID,j=J], binary=true, base_name="start_up")

    # Remove constraints & expressions that need to be adapted

    # Create affine expressions (= linear combinations of variables)

    # Formulate objective

    # Formulate constraints

    return m
end

function build_basic_UC_technical_constraints_model!(m::Model)
    # Build the basic model
    build_basic_UC_model!(m)

    # Extract sets

    # Extract parameters
    MUT = m.ext[:parameters][:MUT]

    # Extract variables
    z = m.ext[:variables][:z]
    v = m.ext[:variables][:v]

    # Additional technical constraints (Minimum up-time constraint "con3d" is given, don't change this name)
    con3d = m.ext[:constraints][:con3d] = Dict()
    for i in ID
        for j in J[MUT[i]:end]
            con3d[i,j] = @constraint(m,
                z[i,j] >= sum(v[i,j-jj] for jj in 0:MUT[i]-1)
            )
        end
    end

    return m
end

# Build model: uncomment the relevant function-call, depending on which model you're running
build_basic_ED_model!(m)
# build_basic_ED_storage_model!(m)
# build_basic_UC_model!(m)
# build_basic_UC_technical_constraints_model!(m)

## Step 4: solve
optimize!(m)
# check termination status
print(
    """

    Termination status: $(termination_status(m))

    """
)
# check objective value
@show value(m.ext[:objective])
@show value(m.ext[:objective])/sum(m.ext[:timeseries][:D])

## Step 5: Visualization
using Plots
using StatsPlots

# sets
J = m.ext[:sets][:J]
I = m.ext[:sets][:I]
ID = m.ext[:sets][:ID]

# parameters
D = m.ext[:timeseries][:D]

# variables/expressions
Ivec = [i for  i in I]
g = value.(m.ext[:variables][:g])
gvec = [g[i,j] for  i in I, j in J]
if haskey(m.ext[:variables],:z)
    z = value.(m.ext[:variables][:z])
    uc = [z[i,j] for i in ID, j in J]
    uc = [zeros(2,24); uc] # plotting purposes -- added zero rows for wind and solor to make legend match
else
    λ = dual.(m.ext[:constraints][:con2a])
    λvec = [λ[j] for j in J]
end

# Check if storage is present in the model
if haskey(m.ext[:variables],:e)
    d = value.(m.ext[:variables][:d])
    c = value.(m.ext[:variables][:c])
    dcvec = [d[j]-c[j] for j in J]
    gvec = [transpose(dcvec); gvec]
    Ivec =["Storage"; Ivec]
end

if haskey(m.ext[:variables],:z)
# UC schedule
p1 = groupedbar(transpose(uc[:,:]),bar_position = :stack,xlabel="Timesteps [-]",ylabel="Commitment [# units]",label="",legend=:outertopright);
else
# electricity price price
p1 = plot(J,λvec, xlabel="Timesteps [-]",ylabel="λ [EUR/MWh]",label="",legend=:outertopright);
end
# dispatch
p2 = groupedbar(transpose(gvec[:,:]),bar_position = :stack,label=permutedims(Ivec),legend=:outertopright);
plot!(J,D, xlabel="Timesteps [-]",ylabel="Generation [MWh]", linewidth = 2,label="");
plot(p1, p2, layout = (1,2))
plot!(size=(900,400))
if haskey(m.ext[:variables],:z) # unit commitment
    if haskey(m.ext[:constraints],:con3d) # technical constraints in the model
    savefig(joinpath(@__DIR__,"BasicUCTechCon"))
    else
    savefig(joinpath(@__DIR__,"BasicUC"))
    end
else
    if haskey(m.ext[:variables],:e) # storage in the model
    savefig(joinpath(@__DIR__, "BasicEDStorage.png"))
    else
    savefig(joinpath(@__DIR__,"BasicED"))
    end
end
