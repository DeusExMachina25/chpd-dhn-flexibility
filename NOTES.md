# Implementation notes

Working notes for the assignment on Mitridati & Taylor, *Power Systems Flexibility
from District Heating Networks*, PSCC 2018. This file records where every
equation lives in the code, every assumption made because the paper does not
give the number, and every problem hit along the way. The report is written
from it.

---

## 1. Where each equation lives

| Eq. | Meaning | File |
|---|---|---|
| (1) | mass flow bounds in the pipes | `build_chpd_minlp.jl`, variable bounds on `mfS`, `mfR` |
| (2) | pressure loss, quadratic equality | `build_chpd_minlp.jl` `:eq2S`, `:eq2R` |
| (3), (4) | pipe inlet temperatures | `:eq3`, `:eq4` |
| (5) | pipe outlet temperature, Taylor heat loss | `:eq5S`, `:eq5R` (zero-delay case only, see §3.1) |
| (6)–(13) | time delays, big-M linearisation | `delays.jl`, `add_delays!` |
| (14) | nodal mass balance | `:eq14S`, `:eq14R` |
| (15) | nodal pressure bounds | variable bounds on `prS`, `prR` |
| (16), (17) | nodal temperature mixing, bilinear | `:eq16`, `:eq17` |
| (18) | nodal temperature bounds | variable bounds on `TS`, `TR` |
| (19), (20) | single-pipe node shortcut | `:eq19`, `:eq20` |
| (21) | HES heat load, bilinear | `:eq21` |
| (22) | HES mass flow bounds | variable bounds on `mfHES` |
| (23) | HES pressure difference threshold | `:eq23` |
| (24) | HS heat output, bilinear | `:eq24` |
| (25) | HS heat bounds | variable bounds on `Q` |
| (26) | HS mass flow bounds | variable bounds on `mfHS` |
| (27) | water pump power, bilinear | `:eq27` |
| (28), (29) | extraction CHP region | `:eq28`, `:eq29` |
| (30) | heat pump COP | `:eq30` |
| (31) | DC power balance | `:eq31` |
| (32) | line flow limits | `:eq32` |
| (33) | generator bounds | variable bounds on `P` |
| (34), (35) | cost and objective | `m.ext[:objective]` |
| (36) | convex quadratic pressure loss | `build_chpd_misocp.jl` `:eq36S`, `:eq36R` |
| McCormick (III-A) | envelopes of (16), (17), (21), (24), (27) | `build_chpd_misocp.jl` `:eq16mc` … `:eq27mc` |
| (37) | outer approximation, MILP | not implemented, optional extra |
| (38), (39) | CED benchmark | `build_ced.jl` `:eq38`, `:eq39` |

---

## 2. Case study: SUPERSEDED, kept as a record

> **This section describes the reconstruction used in phases 1–6, which is no
> longer what the code runs.** Section 6 has the authors' real data. It is kept
> because the report's "challenges" section needs the before/after, and because
> the reconstruction is what the model was first proven on. Do not quote the
> numbers below as results.

Reconstructed from Fig. 2 and Section IV-A. **These are not the paper's
numbers**, the paper puts its data in the online appendix
(doi:10.5281/zenodo.1195508), which we deliberately did not use until the model
itself was proven. Results in this phase therefore validate the *model*, not the
paper's results.

### Topology

DHN, 3 nodes, supply direction:

```
 node 1  (CHP1, heat station)  --- pipe 1 (11 km) --->  node 3  (HES1, the heat load)
 node 2  (HP1,  heat station)  --- pipe 2 ( 8 km) --->  node 3
```

The return network runs the other way along the same node pairs.

This layout was chosen deliberately: node 3 is the only node where two supply
pipes meet, so it is the only place where the bilinear mixing equation (16) is
actually needed. Nodes 1 and 2 have a single arriving *return* pipe, so (20)
applies exactly there and no envelope is needed. Nodes 1 and 2 have no arriving
supply pipe and node 3 has no arriving return pipe, so those three temperatures
are free decision variables set by the stations, which is physically right: a
heat station chooses its supply temperature.

Power system: 6 buses, G1 at bus 1 (reference), CHP1 at bus 2, W1 at bus 4,
HP1 at bus 5, loads at buses 3, 4 and 6.

### Reconstructed parameters and why

| Parameter | Value | Reasoning |
|---|---|---|
| ρ | 963 kg/m³ | water at ~90 °C |
| c | 4182 J/(kg·K) | water |
| Δt | 3600 s | the paper uses 24 hourly steps |
| supply temperature | 353–373 K | 80–100 °C, normal DH supply band. Started at 70–100 °C and narrowed: see §4 on envelope tightness |
| return temperature | 313–323 K | 40–50 °C, normal DH return band, narrowed for the same reason |
| pipe radius | 0.36 m / 0.30 m | sized for ~1.5 m/s at the nominal flow, ≈ DN700/DN600 |
| pipe length | 11 km / 8 km | gives a travel time of roughly 2 h at nominal flow, which is what makes pipeline storage worth modelling at all |
| μ | 1.0 W/(m²·K) | insulated DH pipe; gives ≈ 0.5 % temperature loss per hour of travel |
| φ | 1.4 / 2.0 Pa/(kg/s)² | ≈ 5 bar drop across the long pipe at nominal flow |
| pressure bounds | supply 2–16 bar, return 1–8 bar | typical DH transmission levels |
| Δpr min at HES | 1 bar | enough to drive flow through the exchanger |
| heat load | 84–124 MW | shape from Fig. 3b: morning and evening peaks, midday dip |
| electricity load | 124–224 MW | shape from Fig. 3a: daytime peak |
| wind | 100 MW capacity, AF 0.33–0.92 | anti-correlated with the electricity load, so surplus occurs at night, that is the situation the paper's storage argument depends on |
| G1 | 150 MW, 30 $/MWh | marginal thermal unit |
| CHP1 | ηE 0.35, ηH 0.45, r 0.15, F̄ 500 MW, 12 $/MWh_fuel | typical extraction unit |
| HP1 | 25 MW electric, COP 2.6–3.3 | COP varies with ambient temperature, lowest at night |

Sanity of the scale: at nominal conditions the heat load of ~100 MW needs
`mf = Q/(c·ΔT) = 100/(4.182e-3 · 40) ≈ 600 kg/s`, which at 0.36 m radius is
about 1.5 m/s, a realistic DH transmission velocity. Pump consumption comes
out below 1 MW against a 200 MW power system, which is also the right order of
magnitude.

---

## 3. Assumptions and interpretations

These are the places where the paper is ambiguous or silent, and what was
decided. Worth raising with the TA.

### 3.1 One step means zero time delay

The assignment says to validate one step first. With a single time period the
delay variables τ of (6) can only be zero, so (6)–(13) and their binaries
collapse.

**Corrected in phase 7:** at τ = 0 the loss bracket of (5) is exactly 1, so
there is no thermal loss at all, and not the `T_out = T_in · (1 − 2μΔt/(ρ c R))`
that was written here first. See §6. Zero delay also has to be *earned*: it
requires `mf·Δt ≥ ρπR²L`, which is now imposed.

Consequence, and it matters for reading the results: **with one step the DHN
has no storage**. The whole benefit the paper reports comes from shifting heat
production in time. So at one step the CHPD should cost slightly *more* than
the CED, not less: it pays for pump electricity and pipe heat losses that the
CED simply does not model. The paper's 3.51 % saving is a 24-hour effect and
must not be expected here.

### 3.2 Eq. (29) and (34): efficiencies or their inverses

The nomenclature calls ηH, ηE "heat/electricity fuel efficiency", but they
appear as `ηE·P + ηH·Q ≤ F̄` with F̄ a fuel consumption, and again inside the
cost (34). Read literally as efficiencies this is dimensionally wrong. The
coefficients are therefore implemented as *fuel per unit of output*, i.e. 1/η,
and named `fuelE`, `fuelH` in the yaml and the code. The formulas are then
exactly the paper's and the units work out.

### 3.3 Eq. (14) in the return network

The extracted text of (14) uses the same pipe sets `SP+`/`SP−` for the supply
and the return balance. Since a return pipe runs opposite to its supply pipe,
that reads as the mirror image of the supply balance, and it is implemented
that way: at node n, return flow arrives from the pipes that *depart* in the
supply direction. The heat exchanger injects into the return network and the
heat station withdraws from it. This is the only reading that conserves mass.

### 3.4 Angle reference

A linearised power flow needs `θ = 0` at a reference bus. The paper does not
state one; bus 1 is used. Without it the model is unbounded in θ (though the
objective is not affected).

### 3.5 Wind

The paper does not give a wind cost. Wind is modelled at zero marginal cost
with an hourly availability cap, so curtailment is free and shows up as
`P_W1 < available`.

### 3.6 Envelope bounds

McCormick tightness is entirely determined by the variable bounds it is given.
The envelope for (21) and (24) is built on the *temperature difference* as a
single variable rather than on `TS` and `TR` separately, one product instead
of two, and a much narrower interval. Where equations (19)/(20) apply, the
bilinear constraint is dropped entirely rather than relaxed, because those are
exact.

---

## 4. Problems hit and how they were handled

*(filled in as they come up; this section feeds the "challenges" part of the
report)*

- **Gurobi licence.** Gurobi 11.0 was installed but initially unlicensed
  (`Gurobi Error 10009`), so the first working version of the code ran on Ipopt
  (Sections II and III-B) and HiGHS (the CED). The licence was activated part
  way through and Gurobi is now used for everything. `solvers.jl` still
  detects which of the two is available and picks accordingly, so the project
  runs either way.

  That detour turned into a useful cross-check rather than wasted effort. Ipopt
  is an interior point NLP solver and finds a *local* solution of the
  non-convex Section II; Gurobi with `NonConvex=2` does spatial branch and
  bound and returns a certified *global* one. Both land on **3390.427 $**,
  agreeing to six significant figures. Two unrelated algorithms reaching the
  same point is decent evidence that the Section II formulation has no stray
  local optima on this case, and that the objective reported below is the real
  optimum rather than whatever Ipopt happened to walk into.
- **Julia precompilation** on this machine is slow (several minutes for the
  JuMP/Gurobi/Plots stack). Not a modelling problem, but worth knowing before
  blaming the model for a hang.

- **The McCormick relaxation is loose, and bound width is the whole story.**
  The first run gave a relaxed objective 50 % below the original. Nothing was
  wrong with the formulation: the mass flow bounds in the data file were the
  physical limits of the equipment (`mfHES` anywhere in 40–1600 kg/s) while the
  actual flows sit near 600 kg/s. For `w = x*y` on a box, the McCormick hull
  can sit as far as `(xhi-xlo)*(yhi-ylo)/4` from the true surface, which for
  those bounds is 83 MW of slack on a 98 MW load, so the relaxation could pretend
  to serve the load without moving any water.

  Two rounds of tightening, both rigorous:

  1. Realistic operating bands in `network.yaml` (supply 80–100 °C, return
     40–50 °C, flows within the pipes' real range). Gap 50.0 % → 25.6 %.
  2. Bounds implied by the model's own constraints, computed in
     `build_chpd_minlp!`: a station moving heat `X` across a temperature
     difference that (18) confines to `[dTlo, dThi]` must run a flow between
     `X/(c*dThi)` and `X/(c*dTlo)`. Since the heat load is a known parameter,
     this pins the heat exchanger flow to 390–781 kg/s instead of 250–1400.
     Gap 25.6 % → 23.2 %. The envelopes are built on whatever bounds the
     variables actually carry, so this propagates automatically.

  `report_envelope_widths` in `results.jl` prints the worst-case slack each
  envelope is allowed, which is what made the diagnosis quick.

  23 % is still loose, and no amount of box tightening fixes that: plain
  McCormick is a single linear hull over a wide box. Closing it further needs
  piecewise McCormick (partition the mass flow range, one binary per segment),
  which is the standard remedy and would make the one-step problem a genuine
  MISOCP rather than the convex QCP it currently is. Noted for the report as
  the obvious extension.

---

## 5. One-step results: SUPERSEDED, kept as a record

> **Run on the reconstructed data of §2 and with the Eq. (5) bug described in
> §6.** The relaxation-gap story and the test narrative below still stand; the
> objective values and every "pipe losses" figure do not. Section 6 has the
> current numbers.

Reconstructed data, `t = 1`, no time delays. `julia --project=. chpd.jl`.

| model | objective | solver | status |
|---|---|---|---|
| Section IV-A, CED | 3204.32 $ | Gurobi (LP) | optimal |
| Section II, original | 3390.43 $ | Gurobi `NonConvex=2` | optimal (global) |
| Section III-B, relaxed | 2605.11 $ | Gurobi (convex QCP) | optimal |

**The relaxation is a valid lower bound**, 2605.11 ≤ 3390.43, gap 23.2 %.

**Eq. (36) is exact.** Slack on the relaxed pressure loss inequality is 0.17 Pa
at worst, i.e. it holds with equality at the optimum. This is expected rather
than lucky: pump cost (27) is increasing in the pressure difference, so the
objective pushes the constraint down onto the original quadratic surface. The
convex quadratic part of III-B costs nothing in accuracy here: all of the gap
comes from the McCormick side.

**Physical checks on the Section II solution** all close at machine precision:
mass balance 0 kg/s, heat balance 6.3e-11 MW (104.655 MW produced = 98 MW load
+ 6.655 MW pipe losses), DC power balance 3.3e-11 MW, CHP inside its (28)/(29)
region, heat pump exactly on `Q = COP*LHP`, worst line at 36.8 % of its limit.

**CHPD vs CED: 3390.43 $ against 3204.32 $, i.e. the CHPD costs 5.81 % *more*.**
This is the expected result at a single time step and not a contradiction of
the paper. The DHN's value is storage, and one period has nowhere to store
anything; what remains is the cost the CED ignores: water pumping (0.146 MW)
and 6.655 MW of pipe heat losses. The paper's 3.51 % saving is a 24-hour
effect that only appears once the pipeline can be charged and discharged, so it
cannot be reproduced or expected here. Reproducing it is Phase 7.

### Tests

`julia --project=. test/runtests.jl` runs 14 tests, all passing:

- `mccormick!` reproduces a known bilinear product exactly at the corners of
  the box and brackets it in the interior.
- A two-node network solved by hand: supply temperature at its 373 K ceiling,
  mass flow 427.759 kg/s, heat output 107.333 MW, matched to 1e-4 relative.
  This one earned its place. The first version of the hand derivation assumed
  the return temperature `TR2` would sit on its own 313 K lower bound and the
  test failed at 427.759 vs 413.866. The model was right: the water arrives at
  the heat station already cooled, `TR1 = TR2*(1-g)`, so it is the bound on
  `TR1` at the *far* end of the pipe that binds first, giving
  `TR2 = TRmin/(1-g) = 314.876 K`. The derivation was corrected, not the code.
- The relaxation bounds the original from below on the test network too.

---

## 6. Phase 7: the authors' own data

`appendix/` holds the online appendix (Zenodo doi:10.5281/zenodo.1195508, MD5
`a7388358c72b622a3397f99a620f8403`, verified) and their reference Python
implementation. `data/network.yaml` and `data/timeseries.csv` are now built
from those, replacing the reconstruction of section 2.

### Provenance

Everything comes from `step 1.a initialization MILP MINI 3.py`, which is the
only complete and self-consistent source: the appendix PDF gives the tables but
its column order is ambiguous once the PDF is de-laid-out, and the load and
wind profiles appear there only as figures.

Topology, which differs from the reconstruction: the DHN nodes and the
electricity buses are the **same** six nodes, and the heating network is a
**line**, not a Y:

```
 node 1 (HP1, W1) --pipe 1--> node 2 (CHP1, G2) --pipe 2--> node 3 (HES1, load)
```

G1 sits at bus 6. Lines: 1-2, 2-3, 3-6, 6-5, 5-4, 4-1, 3-5.

A consequence worth stating: **every node has at most one arriving pipe**, so
Eqs. (19)/(20) apply everywhere and the bilinear mixing equations (16)/(17) are
never needed. The McCormick envelopes for them are still implemented and still
correct, they simply have nothing to relax on this network. Only (21), (24) and
(27) are actually relaxed.

Key parameters: pipes R = 0.8 m, L = 500 m, mu = 20 W/(m2K), mf 50-300 kg/s;
supply 90-120 degC, return 30-60 degC; rho = 1000 kg/m3; the heat coefficient
is their `1e-6 * 1.1704 * 3600 = 4.21344e-3`, i.e. water at 4213 J/(kg.K).
CHP1: r = 0.6, rho_E = 2.4, rho_H = 0.25, Fmax = 250 MW, 12.5 $/MWh_fuel,
Qmax = 100 MW. HP1: COP 2.5, Qmax 150 MW. G1 180 MW at 11 $/MWh, G2 100 MW at
33 $/MWh, W1 500 MW at 0.0001 $/MWh. Loads split 20/40/40 over buses 3, 4, 5.

### Where the appendix PDF and their code disagree

The PDF is the published record, the code is what actually produces numbers.
The code was followed, and these were noted rather than silently reconciled:

| | appendix PDF | their code |
|---|---|---|
| pressure bounds | 0-100 kPa | 50-5000 kPa |
| minimum HES flow | 50 kg/s | 150 kg/s |
| water density | 988 kg/m3 | 1000 kg/m3 |
| units present | G1, W1, CHP1, HP1 | also G2 and a heat-only HO1 |

HO1 is omitted here because their own data sets its maximum heat output to
zero, so it can never produce anything.

Their objective also carries two numerical devices that are not in the paper:
elastic mass balance slacks penalised at 1000, and a 1e-5 regularisation on the
time delays. Neither is reproduced here; our mass balances are hard equalities
and close at machine precision without help.

### Two real bugs the authors' data exposed

**1. Eq. (5) was charging thermal loss at zero delay.** The loss bracket in (5)
is `(1 - 2*mu*tau*dt/(rho*c*R))` and it is multiplied by the delay `tau`. At
`tau = 0` the bracket is exactly 1, so a pipe crossed within one time step
loses nothing. The code instead applied a full step of loss. On the
reconstruction this merely inflated the losses; on the authors' data it made
the model **infeasible**, and the infeasibility is what exposed it: their pipes
lose 4.27 % of the *absolute* temperature per step, so two in series drop a
393.15 K supply to 360.3 K, below the 363.15 K floor of Eq. (18). Fixed to
`TSout = TSin`. All the "pipe losses" reported in section 5 were an artefact of
this and are withdrawn.

Losing exactly nothing across a real pipe is of course not physical. It is an
artefact of discretising the delay into whole time steps, and it is what the
paper's Eq. (5) says. Worth a paragraph in the report's critical reflection,
together with the fact that the loss is proportional to the *absolute*
temperature, which implicitly places the ambient at 0 K.

**2. The zero-delay assumption needed a consistency bound.** Asserting
`tau = 0` is only legitimate if Eq. (6) would actually return 0, i.e. if a
whole pipe volume is pushed through within one step: `mf*dt >= rho*pi*R^2*L`.
For these pipes that is `mf >= 279.25 kg/s` against a bound of 50-300. Without
it the no-delay model quietly allows flows that could not possibly clear the
pipe in one step. The bound is now imposed whenever delays are off.

Two smaller things found in the same pass: `gamma` was computed with a
hard-coded `4182.0` instead of the `c` in the data file (harmless then, wrong
the moment the data changed), and `Pa_to_MW` was renamed `pressure_to_MW`
because the authors' pressures are in kPa, not Pa.

### Results on the authors' data, delays still off

| model | 24 h objective |
|---|---|
| Section IV-A, CED | **19276.98 $** |
| Section II, original | 21041.73 $ |
| Section III-B, relaxed | 20749.78 $ (gap 1.39 %) |

**The CED reproduces the paper's 19277 $ exactly.** That is the strongest
evidence so far that the data, the units, the network and the CED formulation
are all right, since it is a single number that depends on nearly all of them.

The CHPD at 21041.73 $ does not yet match the paper's 18600 $, and should not:
with the delays off there is no pipeline storage, so the model has no
flexibility to exploit and simply pays for water pumping on top of the CED. The
24-hour problem is currently 24 independent single-period problems. Closing
that gap is what implementing Eqs. (6)-(13) is for.

The relaxation gap is now only 1.39 %, against 23 % on the reconstruction. The
reason is the flow bound above: it pins the pipe flows into 279-300 kg/s, and
McCormick tightness is driven by exactly that width. Eq. (36) is tight to
within 3.3 Pa.

---

## 7. Time delays and pipeline storage (Eqs. 6-13)

Implemented in `delays.jl`, one `add_delays!` call per network. The chain:

- **(6)-(7)** `u[p,t,s] = 1` iff the mass pushed through the pipe over the last
  `s+1` steps reaches a pipe volume `rho*pi*R^2*L`. Big-M is computed per pipe
  from the widest the residual can get, not picked out of the air.
- **(8)** with `u` monotone in `s` (added explicitly - implied by (7), and it
  cuts the search space a lot), `tau = taumax + 1 - sum(u)` picks the switch
  point.
- **(9)-(12)** `v[p,t,s] = 1` for the single `s` equal to `tau`, and
  `Ttilde[p,t,s]` then carries `Tin[p,t-s] * (1 - s*gamma)`.
- **(13)** the outlet temperature is the sum of `Ttilde` over `s`.

Steps before the start of the horizon need history. As the authors do, the flow
is held at its minimum and the temperature at the cold end of the band.

`taumax` is 6 hours for these pipes: `rho*pi*R^2*L / (mf_min*dt)` with
`mf_min = 50 kg/s`.

### The licence ceiling

The available Gurobi licence is the **size-limited** one (2000 variables and
constraints). The delay model needs roughly 224 constraints per time step -
`u`, `v` and `Ttilde` over 7 candidate delays, 2 pipes, 2 networks - so it fits
up to **6 hours** and fails at 8 with `Gurobi Error 10010`. `safe_optimize!`
now reports this and carries on instead of crashing.

This is a tooling limit, not a modelling one. It is exactly the situation the
paper describes in a different form: the delayed model does not scale, which is
why they relax it.

### Storage does what the paper says it does

Same 6-hour window, hours 18-23, chosen because it straddles the evening event
in this data: the wind availability factor collapses from 0.68 to 0.06 and 0.04
at hours 20 and 21 while the heat load peaks at 111 MW at hour 20.

| hours 18-23 | Section II CHPD |
|---|---|
| delays off, no storage | 11583.52 $ |
| delays on, storage available | **11278.18 $** |

**305.34 $, or 2.64 %, and it comes from nothing but Eqs. (6)-(13).** The two
runs are the same data, the same network and the same units; the only
difference is whether the model is allowed to push heat into the pipes early
and take it out later. The paper reports 3.51 % over 24 hours, so a 2.64 %
saving over a 6-hour window containing the single sharpest wind drought of the
day is the right order of magnitude.

Run them yourself with:

```
julia --project=. chpd.jl 6 from=18 delays
julia --project=. chpd.jl 6 from=18
```

### What is still open

The 24-hour CHPD cannot be solved with the current licence, so the paper's
18600 $ is **not** reproduced. Only the CED half of the comparison is
(19276.98 $ against 19277 $, §6).

One honest observation while it is open: over hours 18-23 our CHPD costs *more*
than the CED (14.3 % with storage), where the paper has it costing less. That
is not obviously a bug. The CED has no heating network at all, so it lets the
cheap heat pump at node 1 serve the whole load; the CHPD has to push that heat
down two pipes capped at 300 kg/s, which at the 90 K maximum temperature
difference carries at most `c*300*90 = 113.8 MW` - against a 111 MW heat load
at hour 20. The network is genuinely close to binding, which forces the more
expensive CHP on and is a real effect the CED cannot see. Whether the 24-hour
picture reverses, as the paper reports, needs the full run.

### Working around the licence: SCIP, then Juniper

The size-limited Gurobi licence caps the delay model at 6 hours, so the 24-hour
run needs a solver without a size cap. What is required is one that handles
**integers together with nonlinearity**, which rules out most of the free
options individually:

| solver | integers | nonlinear constraints | usable here |
|---|---|---|---|
| Ipopt | no | yes | no |
| HiGHS | yes | no (linear/quadratic objective only) | no |
| Clarabel | no | conic | no |
| SCIP | yes | yes, incl. non-convex | crashes, see below |
| Juniper | yes | yes, via Ipopt per node | **yes** |

SCIP was tried first and is the right answer on paper: it does spatial branch
and bound, so it covers the non-convex Section II as well as the MISOCP of
Section III-B. On this machine its Windows binary **segfaults inside
`optimize!`**, and not on our model but, on a two-variable MILP. It loads and
constructs an optimizer fine, then dies during solve. It was removed from
`Project.toml` rather than left as a broken dependency.

Juniper is the working alternative. It is pure Julia, so there is no native
binary to go wrong, and it branches on the integers with Ipopt solving the
continuous problem at each node. `safe_optimize!` catches Gurobi's
`Error 10010` and retries with it automatically, announcing the switch.

**This changes what the two objectives mean, and the report has to say so:**

- **Section III-B is convex** once the integers are fixed, so branch and bound
  with a convex NLP at each node returns a genuine global optimum of the
  MISOCP. Nothing is lost.
- **Section II is not convex**, so Juniper is a heuristic there. Its objective
  is a valid *upper* bound on the true optimum, not the optimum itself.

The second point makes the Section III-B lower bound *more* useful, not less:
with Gurobi unavailable at this size, the relaxation is the only thing
producing a certified bound on the optimal cost at all. That is a fair
practical illustration of why the paper bothers with a relaxation, arrived at
for a more mundane reason than theirs.
