# Power Systems Flexibility from District Heating Networks

Julia/JuMP implementation of the Combined Heat and Power Dispatch (CHPD) of

> L. Mitridati and J. A. Taylor, *Power systems flexibility from district
> heating networks*, 2018 Power Systems Computation Conference (PSCC), Dublin,
> Ireland, pp. 1–7. doi:[10.23919/PSCC.2018.8442617](https://doi.org/10.23919/PSCC.2018.8442617)

A district heating network stores thermal energy simply because water takes time
to travel down a pipe. That storage lets combined heat and power units follow the
electricity system instead of the heat load. This repository builds the paper's
dispatch model, which co-optimises a DC optimal power flow with the hydraulic and
thermal state of the heating network, and the convex relaxation that makes it
tractable.

The assignment deliverable is the **Section II** model and its **Section III-B**
relaxation, validated at a single time step. Section IV-A and the pipeline time
delays are additional work, marked as such below.

| Builder | Paper section | What it is | Scope |
|---|---|---|---|
| `build_chpd_minlp!` | II | the original formulation: bilinear equalities and a quadratic pressure loss, so non-convex | required |
| `build_chpd_misocp!` | III-B | the convex relaxation: McCormick envelopes (III-A) on the bilinear terms plus the convex quadratic pressure loss (36) | required |
| `build_ced!` | IV-A | conventional economic dispatch, the benchmark with no heating network at all | extension |
| `add_delays!` | II, Eqs. (6)–(13) | pipeline time delays, i.e. the storage itself | extension |

Section II exists so that the relaxation can be checked against the problem it
relaxes.

## Prerequisites

Install [Julia](https://julialang.org/downloads/) 1.11 or later.

**No commercial licence is required.** `src/solvers.jl` picks what is available:

| situation | solver used |
|---|---|
| licensed Gurobi present | Gurobi (`NonConvex=2` for Section II) |
| no Gurobi licence | Ipopt for Sections II and III-B, HiGHS for the CED |
| model too large for a size-limited Gurobi licence | Juniper (Ipopt + HiGHS), automatically |

Gurobi is preferred because Section II is non-convex and it does spatial branch
and bound. The fallbacks are not equivalent, and the difference matters when
reading the output: **Section III-B is convex**, so any of them returns its
global optimum; **Section II is not**, so under Ipopt or Juniper its objective is
a valid upper bound rather than a certified optimum. Where both stacks could run,
they agreed (3390.427 $ to six significant figures on the earlier reconstructed
case).

## How to run it

From the repository root, instantiate the environment. This installs the exact
versions pinned in `Manifest.toml`; the first run takes a few minutes to
precompile, later runs start in seconds.

```bash
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

Then run the study:

```bash
julia --project=. src/chpd.jl
```

This builds the models, solves them, and prints the validation report: objective
values, the relaxation gap, how exact the relaxation turned out to be, and a set
of physical checks on the solution.

Everything is driven by command-line arguments — no file needs editing to change
what is solved. A bare integer sets the number of time steps, `from=N` starts the
horizon at hour `N`, and `delays` enables Eqs. (6)–(13).

| Command | What it solves |
|---|---|
| `julia --project=. src/chpd.jl` | one time step at hour 1, no delays — the default |
| `julia --project=. src/chpd.jl 1 from=20` | one time step at hour 20, the wind-scarce hour reported in the report |
| `julia --project=. src/chpd.jl 24` | the full 24-hour horizon, no delays |
| `julia --project=. src/chpd.jl 6 from=18 delays` | six hours from hour 18, with pipeline storage on |

The last pair is the storage result: the same six hours spanning the evening wind
drought, with and without the delays.

Run the tests with:

```bash
julia --project=. test/runtests.jl
```

The tests check the McCormick helper against a product it must reproduce exactly,
and solve a two-node network whose optimal dispatch is derived by hand in the
test file.

### Expected runtime and output

| Run | Runtime | Notes |
|---|---|---|
| single time step | seconds | after precompilation |
| 24 h, no delays | under a minute | |
| 6 h with delays | a few minutes | branch and bound over the delay binaries |
| tests | under a minute | 14 tests, all passing |

As a reference for checking a reproduction, `src/chpd.jl 1 from=20` prints:

```
Section IV-A (CED)      status: OPTIMAL   objective: 4377.174 $
Section II (MINLP)      status: OPTIMAL   objective: 5223.621 $
Section III-B (relaxed) status: OPTIMAL   objective: 5208.639 $
  gap                     : 0.2868 %
  ok: valid lower bound
```

together with mass, heat and power balance residuals at or near machine
precision (1e-14 to 1e-11), CHP1 at its 250 MW fuel limit, and the heat pump
exactly on `Q = COP·L`.

The delay model of Eqs. (6)–(13) exceeds the 2000 variable/constraint cap of the
free size-limited Gurobi licence from roughly 8 hours onward; the run detects
this and switches to Juniper rather than failing.

## Repository layout

```
.
├── src/          all model code; src/chpd.jl is the single entry point
├── data/         case study input: network.yaml (topology, units) and
│                 timeseries.csv (loads, wind availability, COP)
├── test/         unit tests and a small analytically solvable network
├── results/      generated output (git-ignored except .gitkeep)
├── report/       the assignment report
├── appendix/     the authors' online appendix and reference code, the data source
├── NOTES.md      equation-to-code map, assumptions, and problems encountered
├── Project.toml  the environment, with version bounds
└── Manifest.toml the exact dependency tree this was produced with
```

Inside `src/`:

| File | Contents |
|---|---|
| `chpd.jl` | entry point: loads data, builds and solves the models, validates |
| `init_model.jl` | `define_sets!`, `process_time_series_data!`, `process_parameters!` |
| `build_chpd_minlp.jl` | Section II, the original non-convex CHPD |
| `build_chpd_misocp.jl` | Section III-B, McCormick envelopes + convex quadratic relaxation |
| `build_ced.jl` | Section IV-A, the conventional economic dispatch benchmark |
| `delays.jl` | Eqs. (6)–(13), the pipeline time delays and their warm start |
| `results.jl` | validation checks and reporting |
| `solvers.jl` | solver selection and the fallbacks described above |

## Data and results

`data/` holds the **authors' own case study**, taken from their online appendix
(doi:[10.5281/zenodo.1195508](https://doi.org/10.5281/zenodo.1195508)) and
reference implementation, both in `appendix/`. That code is the only complete and
self-consistent source: the appendix PDF's tables are ambiguous once the PDF is
de-laid-out, and its load and wind profiles appear there only as figures. Section
6 of `NOTES.md` records the provenance of every parameter and the places where
their appendix PDF and their code disagree.

What has been established:

- **The Section III-B relaxation is a valid lower bound** on Section II at every
  horizon tested, and Eq. (36) comes out tight — the relaxed pressure loss sits on
  the original quadratic surface, because pump cost pushes it there. The remaining
  gap is entirely McCormick's, and tracks the width of the bounds the envelopes
  are built on.
- **The physical checks close at machine precision** on the Section II solution:
  mass, heat and DC power balance, the CHP operating region, and the heat pump
  conversion.
- **The CED benchmark reproduces the paper exactly**: 19276.98 $ against the
  19277 $ reported in Section IV-B. That single number depends on the data, the
  units, the topology and the CED formulation all being right at once. *(extension)*
- **Pipeline storage does what the paper says.** Over hours 18–23, which span the
  evening collapse in wind availability at peak heat load, switching on Eqs.
  (6)–(13) saves 305.34 $ (2.64 %) against the identical model with them off. The
  paper reports 3.51 % over a full day. *(extension)*

With one time step the network cannot store anything, so the CHPD is expected to
cost *more* than the CED there — it pays for water pumping, which the CED does not
model at all. That is not a discrepancy with the paper.

The paper's full-horizon 18600 $ CHPD figure is not reproduced. It needs 24 hours
with the delays on, which exceeds the available size-limited solver licence;
reproducing it was outside the scope of this assignment. See `NOTES.md` §7.
