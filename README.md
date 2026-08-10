# Power Systems Flexibility from District Heating Networks

Julia/JuMP implementation of the Combined Heat and Power Dispatch (CHPD) of

> L. Mitridati and J. A. Taylor, *Power systems flexibility from district
> heating networks*, 2018 Power Systems Computation Conference (PSCC), Dublin,
> Ireland, pp. 1–7. doi:[10.23919/PSCC.2018.8442617](https://doi.org/10.23919/PSCC.2018.8442617)

Three models are built on the same data set:

| Model | Paper section | What it is |
|---|---|---|
| `build_chpd_minlp!` | II | the original formulation: bilinear equalities and a quadratic pressure loss, so non-convex |
| `build_chpd_misocp!` | III-B | the convex relaxation: McCormick envelopes (III-A) on the bilinear terms plus the convex quadratic pressure loss (36) |
| `build_ced!` | IV-A | conventional economic dispatch, the benchmark with no heating network at all |

The Section III-B model is the deliverable; Section II exists so that the
relaxation can be checked against the problem it relaxes.

## How to run it

1. Install [Julia](https://julialang.org/downloads/) 1.11 or later and make
   sure a Gurobi installation with a valid licence is on the machine
   (`GUROBI_HOME` set). Gurobi is used because Section II is non-convex and
   Section III-B is a convex quadratically constrained problem.

2. From this directory, instantiate the environment. This downloads and
   precompiles everything listed in `Project.toml`; the first run takes a few
   minutes.

   ```bash
   julia --project=. -e "using Pkg; Pkg.instantiate()"
   ```

3. Run the whole study end to end:

   ```bash
   julia --project=. chpd.jl
   ```

   This builds the three models, solves them, and prints the validation
   report: objective values, the relaxation gap, how exact the relaxation
   turned out to be, and a set of physical checks on the solution.

4. Run the tests:

   ```bash
   julia --project=. test/runtests.jl
   ```

   The tests check the McCormick helper against a product it must reproduce
   exactly, and solve a two-node network whose optimal dispatch is derived by
   hand in the test file.

## Directories

```
.
├── chpd.jl               entry point: loads data, builds and solves the three models, validates
├── init_model.jl         define_sets!, process_time_series_data!, process_parameters!
├── build_chpd_minlp.jl   Section II, the original non-convex CHPD
├── build_chpd_misocp.jl  Section III-B, McCormick envelopes + convex quadratic relaxation
├── build_ced.jl          Section IV-A, the conventional economic dispatch benchmark
├── results.jl            validation checks and reporting
├── data/                 case study input: network.yaml (topology, units) and timeseries.csv (loads, wind, COP)
├── results/              generated tables and figures
├── test/                 unit tests and a small analytically solvable network
├── report/               the assignment report
├── NOTES.md              equation-to-code map, assumptions, and problems encountered
└── Project.toml          the Julia environment this code runs in
```

## Current scope

The model currently runs a **single time step with the pipeline time delays
switched off**, which is the first validation stage asked for in the
assignment. With one period the delay variables of Eq. (6) can only be zero, so
Eqs. (6)–(13) and all their binary variables disappear and Eq. (5) reduces to a
pure heat loss factor.

An important consequence for reading the results: with one time step the
network cannot store anything, so the CHPD is expected to cost slightly *more*
than the CED, because it pays for water pump electricity and pipe heat losses
that the CED does not model. The 3.51 % saving reported in the paper comes from
shifting heat production across a 24 hour horizon.

The case study data in `data/` is reconstructed from the paper's Figures 2 and
3 and is **not** the authors' own data; every choice is justified in
`NOTES.md`. Running against the authors' online appendix data, enabling the
time delays and extending to 24 hours is the next step.
