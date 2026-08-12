## Solver selection
# Gurobi is the natural choice here: it handles the non-convex Section II
# directly (spatial branch and bound, NonConvex=2) and the convex quadratic
# Section III-B, and it is the only one of the three that will cope with the
# MISOCP once the time delay binaries are switched on.
#
# It needs a licence though. Without one we fall back to Ipopt and HiGHS:
#   - Section III-B is convex, so Ipopt's local optimum is the global one and
#     nothing is lost;
#   - Section II is not, so its objective is a local solution and should be
#     read as an upper bound on the optimum, not as the optimum;
#   - the CED is a plain LP, HiGHS solves it exactly.

# The escape hatch, for when the size-limited Gurobi licence (2000 variables
# and constraints) refuses the delay model of Eqs. (6)-(13), which it does from
# about 8 hours onwards.
#
# What is needed is a solver that handles integers together with nonlinearity.
# Ipopt does nonlinear but not integers, HiGHS does integers but not quadratic
# constraints, Clarabel does cones but not integers. SCIP was tried first and
# is the natural choice on paper, but its Windows binary segfaults inside
# optimize! on this machine - even on a two-variable MILP - so it is unusable
# here. Juniper is the working alternative: pure Julia, no native binary, and
# it does branch and bound over the integers with Ipopt on each node.
#
# The distinction that matters when reading the results:
#   - Section III-B is CONVEX once the integers are fixed, so Juniper's
#     branch and bound returns a genuine global optimum of the MISOCP.
#   - Section II is NOT convex, so Juniper is a heuristic there and its
#     objective is a valid upper bound, not a certified optimum. That makes
#     the Section III-B lower bound more important, not less.
using JuMP, Gurobi, Ipopt, HiGHS, Juniper

const HAS_GUROBI = try
    Gurobi.Env(output_flag=0)
    true
catch
    @warn "No Gurobi licence found, falling back to Ipopt and HiGHS. See NOTES.md."
    false
end

# Time limit for the two hard models. Section II over 24 hours with the delay
# binaries is a non-convex MINLP and does not close in any reasonable time -
# which is the whole reason the paper proposes a relaxation. When the limit
# bites we report the incumbent and the bound rather than pretending the number
# is optimal.
const TIMELIMIT = 600.0

# Juniper gets longer: it is branch and bound written in Julia solving an NLP
# at every node, so it is far slower per node than Gurobi and needs the room.
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
