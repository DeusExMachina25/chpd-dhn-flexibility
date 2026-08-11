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
| (5) | pipe outlet temperature, Taylor heat loss | `:eq5S`, `:eq5R` (zero delay, see §3.1) |
| (6)–(13) | time delays, big-M linearisation | not yet implemented, phase 7 |
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

## 2. Case study

Reconstructed from Fig. 2 and Section IV-A. **These are not the paper's
numbers** — the paper puts its data in the online appendix
(doi:10.5281/zenodo.1195508), which we deliberately do not use until the model
itself is proven. Results in this phase therefore validate the *model*, not the
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
are free decision variables set by the stations — which is physically right, a
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
| wind | 100 MW capacity, AF 0.33–0.92 | anti-correlated with the electricity load, so surplus occurs at night — that is the situation the paper's storage argument depends on |
| G1 | 150 MW, 30 $/MWh | marginal thermal unit |
| CHP1 | ηE 0.35, ηH 0.45, r 0.15, F̄ 500 MW, 12 $/MWh_fuel | typical extraction unit |
| HP1 | 25 MW electric, COP 2.6–3.3 | COP varies with ambient temperature, lowest at night |

Sanity of the scale: at nominal conditions the heat load of ~100 MW needs
`mf = Q/(c·ΔT) = 100/(4.182e-3 · 40) ≈ 600 kg/s`, which at 0.36 m radius is
about 1.5 m/s — a realistic DH transmission velocity. Pump consumption comes
out below 1 MW against a 200 MW power system, which is also the right order of
magnitude.

---

## 3. Assumptions and interpretations

These are the places where the paper is ambiguous or silent, and what was
decided. Worth raising with the TA.

### 3.1 One step means zero time delay

The assignment says to validate one step first. With a single time period the
delay variables τ of (6) can only be zero, so (6)–(13) and their binaries
collapse and (5) reduces to a pure loss factor
`T_out = T_in · (1 − 2μΔt/(ρ c R))`.

Consequence, and it matters for reading the results: **with one step the DHN
has no storage**. The whole benefit the paper reports comes from shifting heat
production in time. So at one step the CHPD should cost slightly *more* than
the CED, not less — it pays for pump electricity and pipe heat losses that the
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
single variable rather than on `TS` and `TR` separately — one product instead
of two, and a much narrower interval. Where equations (19)/(20) apply, the
bilinear constraint is dropped entirely rather than relaxed, because those are
exact.

---

## 4. Problems hit and how they were handled

*(filled in as they come up — this section feeds the "challenges" part of the
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
  those bounds is 83 MW of slack on a 98 MW load — the relaxation could pretend
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

  23 % is still loose, and no amount of box tightening fixes that — plain
  McCormick is a single linear hull over a wide box. Closing it further needs
  piecewise McCormick (partition the mass flow range, one binary per segment),
  which is the standard remedy and would make the one-step problem a genuine
  MISOCP rather than the convex QCP it currently is. Noted for the report as
  the obvious extension.

---

## 5. One-step results

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
convex quadratic part of III-B costs nothing in accuracy here — all of the gap
comes from the McCormick side.

**Physical checks on the Section II solution** all close at machine precision:
mass balance 0 kg/s, heat balance 6.3e-11 MW (104.655 MW produced = 98 MW load
+ 6.655 MW pipe losses), DC power balance 3.3e-11 MW, CHP inside its (28)/(29)
region, heat pump exactly on `Q = COP*LHP`, worst line at 36.8 % of its limit.

**CHPD vs CED: 3390.43 $ against 3204.32 $, i.e. the CHPD costs 5.81 % *more*.**
This is the expected result at a single time step and not a contradiction of
the paper. The DHN's value is storage, and one period has nowhere to store
anything; what remains is the cost the CED ignores — water pumping (0.146 MW)
and 6.655 MW of pipe heat losses. The paper's 3.51 % saving is a 24-hour
effect that only appears once the pipeline can be charged and discharged, so it
cannot be reproduced or expected here. Reproducing it is Phase 7.

### Tests

`julia --project=. test/runtests.jl` — 14 tests, all passing:

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
