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

using JuMP, Gurobi, Ipopt, HiGHS

const HAS_GUROBI = try
    Gurobi.Env(output_flag=0)
    true
catch
    @warn "No Gurobi licence found, falling back to Ipopt and HiGHS. See NOTES.md."
    false
end

nonconvex_solver() = HAS_GUROBI ?
    optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => 2, "OutputFlag" => 0) :
    optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)

convex_solver() = HAS_GUROBI ?
    optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0) :
    optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)

lp_solver() = HAS_GUROBI ?
    optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0) :
    optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)
