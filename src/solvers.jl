## Solver selection
#
# Without a Gurobi licence we fall back to Ipopt and HiGHS. Section III-B is
# convex, so Ipopt's local optimum is the global one; Section II is not, so its
# objective is then only an upper bound, not the optimum.
#
# Juniper is the escape hatch for when the size-limited licence (2000 variables
# and constraints) refuses the delay model, which happens from about 8 hours.
# SCIP is the natural choice but its Windows binary segfaults inside optimize!
# here, even on a two-variable MILP. Juniper is pure Julia and branches over the
# integers with Ipopt at each node. Same caveat as above: it is global on the
# convex Section III-B and only a heuristic on Section II.
using JuMP, Gurobi, Ipopt, HiGHS, Juniper

const HAS_GUROBI = try
    Gurobi.Env(output_flag=0)
    true
catch
    @warn "No Gurobi licence found, falling back to Ipopt and HiGHS. See NOTES.md."
    false
end

# Section II over 24 hours with the delay binaries does not close in any
# reasonable time. When the limit bites we report the incumbent and the bound
# rather than pretending the number is optimal.
const TIMELIMIT = 600.0

# Juniper solves an NLP at every node, so it needs more room than Gurobi.
const TIMELIMIT_UNCAPPED = 1800.0

nonconvex_solver() = HAS_GUROBI ?
    optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => 2, "OutputFlag" => 0,
                              "TimeLimit" => TIMELIMIT) :
    optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)

convex_solver() = HAS_GUROBI ?
    optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0,
                              "TimeLimit" => TIMELIMIT) :
    optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)

lp_solver() = HAS_GUROBI ?
    optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0) :
    optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)

# Used when Gurobi refuses a model for being too large.
uncapped_solver() = optimizer_with_attributes(Juniper.Optimizer,
    "nl_solver" => optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0),
    "mip_solver" => optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false),
    "time_limit" => TIMELIMIT_UNCAPPED,
    "log_levels" => [])
