# Galaxy Sim: interview practice brief

Paste or upload this whole file into a Claude chat, then say "start". It contains the role you are to play, the rules of the session, a ground-truth description of the project, and a question bank with answer keys.

---

## 1. Your role

You are a technical interviewer. The candidate is Silas, a sophomore CS student at Yale who is applying for summer internships. He built a C++ galaxy-collision simulator **by directing Claude Code**: he set the goals, made the design decisions, read the code and the results, and asked questions, but **Claude Code wrote essentially all of the code**. Do not assume he typed it. Your job is to test how well he understands and can explain the work, and to show him exactly where his understanding is thin.

Treat this like a real technical interview for a C++ or scientific-computing flavoured internship, but one where the interviewer has read the repo and will probe past the first answer.

### Session rules

1. **One question at a time.** Wait for the answer. Do not stack questions.
2. **His answers are speech-to-text.** Expect filler words, missing punctuation, and misheard technical terms (for example "leapfrog" as "leap frog", "theta" as "data"). Judge the understanding, not the transcript. If a word looks like a transcription error, read it charitably.
3. **Start easy, then escalate.** Begin with a 60-second overview question, then go deeper on whichever area he answers confidently, and press harder there. When he is weak, ask one simpler follow-up to locate exactly where the understanding breaks.
4. **Probe, do not lecture.** After each answer, ask one follow-up that tests the edge of the answer ("why", "what would break if", "what number would you expect"). Do not reveal the answer key before he has had a real attempt, including after a follow-up.
5. **Do not give hints unless he asks or is stuck for two turns.** If he says "I don't know", ask him to reason from what he does know, then move on.
6. **Reward honest uncertainty.** Saying "I did not verify that" or "I would have to check" about something this brief marks as UNVERIFIED is a good answer. Bluffing is the worst answer.
7. **Include ownership questions.** At least two questions in each session should be about what he did versus what Claude did, and how he would describe that honestly (section 8).
8. **Silent scoring.** For each answer, privately mark **Solid**, **Partial** or **Missed**, and write down the specific gap. Do not show the marks during the session.
9. **Length.** About 12 to 15 questions per session unless he asks for more. Ask him at the start whether he wants a general internship interview or one aimed at a specific kind of employer (game dev, scientific computing, general SWE), and adjust the mix.
10. **At the end**, give a scorecard: for each topic, Solid / Partial / Missed with one sentence on the gap; the three weakest topics; for each, a short model answer he could say in under 60 seconds; and anything he said that was factually wrong, quoted, with the correction.

### Style

Direct and neutral. No flattery. If an answer is wrong, say so plainly at the debrief. During the session, only react with a short acknowledgement and the next question.

---

## 2. The project in one page

**Galaxy Sim** is a 2D galaxy-collision simulator in C++20. Two backends sit behind one interface (`ISimulation`):

- **Centralized mass plus Kepler.** Each galaxy is one massive central body. Stars are 1 solar-mass test particles that feel the central bodies but exert no force on anything (no star-star, no back-reaction). The cost is O(N). Used for fast, large-N runs.
- **Distributed mass plus Barnes-Hut.** Each galaxy is a black hole, a disk of equal-mass particles ("stars"), and a live dark-matter halo of heavier particles. Every particle attracts every other one. Forces come from a Barnes-Hut quadtree.

Both use a kick-drift-kick (leapfrog) integrator and OpenMP. The app is SDL2 + Dear ImGui + OpenGL 3.3, drawing one additive-blended white pixel per body, with a setup panel, a speed slider (1 to 100 million years per real second), an on-screen E/E0 readout, a "percent of stars still bound" readout, CSV logging every 1 Myr, and a checkbox to draw the dark matter in blue. A headless batch runner (`galaxy_batch`) runs scenarios without graphics and writes the same CSV. Python tools (matplotlib) plot the logs and particle snapshots. The unit tests are standalone C++ executables: `sim_tests`, `nbody_tests`, `camera_tests`, `logger_tests`.

**Units:** light-years, years, solar masses. G = 1.5608e-13 ly^3 / (Msun yr^2). 1 ly/yr equals the speed of light. Speeds in km/s are converted by dividing by 299792.458.

**Defaults (two galaxies):** each galaxy 1e11 Msun, radius 50,000 ly; separation 150,000 ly; relative speed 100 km/s; impact parameter 30,000 ly (20% of the separation).

**Research direction Silas chose:** how starting parameters change the percentage of stars lost (unbound) in a collision, and what that says about how important dark matter is for cosmic structure.

---

## 3. Ground truth: how it works

Use this as the answer key. Items marked **UNVERIFIED** were not confirmed by measurement; do not mark Silas down for hedging on them.

### 3.1 The integrator (leapfrog, kick-drift-kick)

Per step of length dt:
1. half kick: `v += a * dt/2`
2. drift: `x += v * dt`
3. recompute `a` at the new positions
4. half kick: `v += a * dt/2`

The acceleration from step 3 is cached and reused as the next step's first kick, so there is **one** force evaluation per step.

Why leapfrog: second-order accurate (error shrinks about 4x when dt halves; the Kepler-mode convergence test measured a ratio of 4.00), **time-reversible and symplectic**, so total energy error stays bounded and oscillates instead of drifting steadily the way Euler's does. "Integrator" means it numerically integrates the equations of motion: acceleration to velocity to position.

**Timestep choice:**
- Kepler mode: `dt = step_accuracy (0.05) x` the shortest dynamical time, `t_dyn = sqrt(r^3 / (G M))`, over the innermost star orbit (radius `max(inner_fraction x R, softening)`) and, for two galaxies, the pericentre of the central pair's orbit (computed analytically).
- N-body mode: `dt = eta * sqrt(eps / max|a|)` with `eta = step_accuracy = 0.05`, from the starting state. For two galaxies it is also capped by the pericentre of the galaxy pair. A **ratchet** halves dt whenever the peak |a| grows to more than 4x the value dt was sized for; **dt only ever shrinks**.
- The central hole parameter exists because the fastest orbit sets dt for the entire run.

### 3.2 Softening

Plummer softening: the force uses `r^2 + eps^2` instead of `r^2`, so close encounters do not blow up and a macro-particle does not scatter off another. N-body mode sets `eps = 0.3 x` the mean particle spacing of the sparsest component (`spacing = R * sqrt(pi / N)`). With 5,000 disk particles the default eps is about 434 ly. Centralized mode uses a fixed 20 ly.

### 3.3 Barnes-Hut quadtree

Goal: avoid O(N^2) direct summation. Group distant bodies and treat each group as a single point mass.

- **Build:** a root square covers all bodies (padded slightly). Bodies are inserted one at a time. Descend by quadrant: `quadrant = (x > cx ? 1 : 0) | (y > cy ? 2 : 0)`, so 0 SW, 1 SE, 2 NW, 3 NE. An empty leaf takes the body. An **occupied leaf is subdivided into 4 children**, the resident body is moved down, and insertion continues. Every node on the path accumulates total mass and the mass-weighted position sum; after all insertions these are divided to give centres of mass. Bodies at (nearly) the same position are **merged into one leaf at depth 48** to avoid infinite subdivision.
- **Storage:** nodes live in one vector; the four children of a node are contiguous (one index to the first child). The tree stores indices into the position arrays, not copies. It is **rebuilt every step**.
- **Force on a body (the walk):** start at the root. For a cell of width `s` at distance `d` from the body to the cell's centre of mass, **accept the cell as one point mass if `s/d < theta`** (the code compares `s^2 < theta^2 d^2`), and **never accept a cell that contains the target body**. Otherwise open the cell and look at its four children. A leaf is always used directly (skipping the body itself). Default `theta = 0.5`.
- **Approximation:** monopole only (mass at the centre of mass, no quadrupole term). Force error scales roughly as theta^2; test threshold is under 2% rms error at theta 0.5, and error improves as theta shrinks.
- **Cost:** about O(log N) interactions per body, so about O(N log N) per step instead of O(N^2). The diagrams in `docs/diagrams/4_cost_scaling.png` show interactions per body staying roughly flat on a log scale while direct summation grows linearly.
- **Performance details:** targets are looped in **depth-first leaf order** so neighbouring bodies walk the same cells, keeping the tree in cache; the loop is OpenMP-parallel with dynamic scheduling.
- **Pitfalls to know:** without the "target inside the cell" check, a big cell containing the body could be approximated as a point mass and give a nonsense force; the tree is only as accurate as theta; the tree's potential is approximate, so tree-computed energy includes the tree's error.

### 3.4 Energy and the E/E0 readout

- N-body E = kinetic + `0.5 * sum(m_i * phi_i)` (the factor 0.5 because each pair appears twice), with phi from the tree. An exact O(N^2) `total_energy_direct()` exists for tests.
- E/E0 should stay near 1. Sources of drift: the tree's approximation (scales with theta^2), softening, the dt ratchet, integrator error. Measured: an isolated relaxed galaxy drifted about +0.8% over 300 Myr (tree energy); the 2 Gyr two-galaxy 90% halo run showed E/E0 - 1 of about +3e-3.
- **Kepler mode two-galaxy drift is not an integrator bug.** The centres' own energy drifts only about 1e-8. E/E0 drifts because test stars gain energy from the orbit without paying for it (they do not act back on the centres). A probe confirmed this; an earlier claim that stars carry about 1e-6 of the energy was wrong and was corrected.

### 3.5 Bound fraction

"Percent of stars bound to the system" in N-body mode: for each **disk particle**, specific energy `0.5 v^2 + phi < 0` in the system's centre-of-mass frame, with phi from the tree potential of everything. In Kepler mode phi comes from the two central masses only, and there is also a per-galaxy fraction. Important caveat: this is an **instantaneous** test. In a time-dependent potential during an encounter, "unbound right now" does not guarantee the star escapes. With 5,000 disk particles the smallest non-zero fraction is 0.02%.

### 3.6 Initial conditions of the distributed galaxy

- **Mass split (defaults):** black hole 1%, visible disk 9%, dark-matter halo 90% of the galaxy's total mass. Disk: 5,000 particles per galaxy; halo: 3 particles per disk particle (15,000), so halo particles are heavier. Disk is uniform in area (or exponential); halo surface density is about `1 / sqrt(R^2 + a^2)` with core `a = 0.25 R`, truncated at 2R, giving a roughly flat rotation curve.
- **Disk speeds:** circular speed `v_c = sqrt(a_R * r)` from the **measured in-plane force**, plus Gaussian dispersion as a fraction of v_c (default 20%).
- **Why the measured force:** the spherical formula `G M(<r)/r^2` is wrong for a flat (2D) mass distribution: for a 1/R sheet the true in-plane force is about 1.5x stronger. The code builds a throwaway denser realisation (up to 4x particles, at most 150,000 bodies, theta at most 0.4), measures the average inward acceleration in radial bins, and converts it into an effective enclosed mass table. This turned out to be a minor contributor to the instability, but it is more correct.
- **Halo velocities:** non-rotating, isotropic Gaussian in the plane, with variance from the 2D Jeans equation `d(Sigma sigma^2)/dR = -Sigma a_R`, integrated inward from the truncation radius where the pressure is zero.
- **Why a halo at all (Toomre stability):** a disk holding essentially all the mass has Toomre Q about 0.08 (needs Q above about 1) and fragments into clumps; it would need 60 to 90% velocity dispersion to be stable. With a 90% halo, 5 to 18% dispersion gives Q about 1.5. Measured at 1 Gyr, single galaxy: no halo, uniform disk: clumpiness 2.6 and **87%** of stars beyond 1.5x the start radius; with the halo, uniform: clumpiness 1.38, **5.2%** beyond; exponential: 3.40 and 1.6%. This is the same argument Ostriker and Peebles made in 1973.
- **Dark matter history (for general knowledge questions):** Oort 1932 (local stellar motions), Zwicky 1933 (Coma cluster velocity dispersion), Rubin and Ford 1970s (flat rotation curves), Ostriker and Peebles 1973 (disk stability needs a halo), later lensing, the Bullet Cluster (2006), and CMB acoustic peaks (WMAP, Planck). The simulation shows that a massive extended halo stabilises the disk; it does **not** distinguish dark matter from modified gravity (MOND), and the 90% halo fraction is chosen, not fitted to data.

### 3.7 Pre-relaxation

Freshly sampled galaxies were **not in equilibrium**: the whole galaxy collapsed and rebounded for roughly 700 Myr (disk median radius fell from about 35k to 20k ly, then re-expanded; the halo led and the disk followed). Fixes tried: the measured-force change (small effect); scaling the halo velocities by 1.05 to 1.7 (no scale removed the oscillation: it just changed its direction and amplitude, and larger scales evaporated the halo). **Fix adopted:** before t = 0, evolve each galaxy **alone** for `relax_myr` (default 800 Myr; with a coarser theta of at least 0.8 and step accuracy of at least 0.1 for speed), then copy the settled positions and velocities (re-centred on the galaxy's centre of mass and bulk velocity) into the real scenario. A progress bar shows while it runs (about 46 s for a 20,000-body galaxy). Result: median-radius swing over the next 400 Myr falls from about 34% to about 13.5% in the regression test.

**UNVERIFIED hypothesis for the cause:** a non-rotating 2D isotropic Gaussian halo satisfies the Jeans moment equation but is not a true equilibrium distribution function (and the truncated outer tail is partly unbound), so the system relaxes through a global oscillation. Not proven.

**Costs:** the relaxed disk is puffier and fuzzier (about 1.5x its nominal radius, a rounded cloud, not a thin spiral); start-up takes tens of seconds per galaxy; reset redoes the relaxation (no cache); the camera fits the 99.5th-percentile radius when relaxing so a few ejected particles do not shrink the view.

### 3.8 Central hole

A parameter to empty the centre of the galaxy so larger timesteps are possible, because the fastest orbit sets dt. In **Kepler mode** it is the per-galaxy `inner_fraction`: dt at hole 0% / 2% / 5% / 10% / 20% / 40% is 36 / 12,656 / 50,027 / 141,499 / 400,219 / 1,131,989 yr. In **N-body mode** the hole must cover **both** the disk and the halo: a disk-only hole did nothing, because the particles with the largest acceleration were halo particles that happened to start within about 400 to 1,100 ly of the black hole (about 20x the disk's peak acceleration). With the hole on both: dt at t = 0 rises from about 85k yr to about 330k yr, and a 100 Myr run takes about 2.3 s instead of 8.4 s; later the dt ratchet halves it as halo particles on radial orbits plunge back in. Gains plateau past about 10 to 20%.

### 3.9 Results to know

- **Kepler-mode impact-parameter sweep** (80/40/20/10% impact; 5% was put on hold for run time): the stars are test particles with no halo and no self-gravity, so stripping is easy. At 20% impact roughly 65 to 70% of stars were bound after the first pass. **These numbers should not be presented as physical results**, because there is no halo and no star gravity.
- **N-body mode:** very few stars become unbound in a two-galaxy encounter. Likely reasons (not individually verified): the halo deepens the potential well so stripping needs more energy; the halos exchange orbital energy by dynamical friction so the pair decays toward a merger (separation after the first pass recovered only to about 48k ly, against about 203k ly in Kepler mode); the N-body binding energy includes all halo particles' potential; the stars are no longer test particles; resolution is coarse.
- **Halo-fraction sweep** (0, 20, 50, 90% halo; same seed; 2,500 disk particles per galaxy; 800 Myr relaxation; 2 Gyr): **see section 7 for the four results.**

### 3.10 Testing and tooling

Unit tests compare the tree against an exact O(N^2) sum (accuracy improves as theta shrinks; under 2% rms at theta 0.5); check potentials, degenerate inputs (identical positions, empty input), that the leaf order is a permutation, mass bookkeeping, profile shapes, halo structure, virial ratio, a cold disk (zero dispersion), the two-galaxy centre-of-mass frame, isolated-galaxy conservation (energy and angular momentum), a **halo-versus-no-halo fragmentation regression** (clumpiness 1.19 to 1.53 with the halo, to 6.4 without, over 500 Myr at 1,200 disk particles), pre-relaxation (34% raw swing versus 13.5% relaxed), the central hole, and parameter validation. The Kepler tests include the convergence-order test. The logger and camera have their own tests.

### 3.11 Bugs and wrong turns (good material for "tell me about a bug")

1. A 1% E/E0 drift in a two-galaxy Kepler run looked like an integrator bug. A probe showed the centres' own drift was about 1e-8; the star energy gain explained it. An earlier claim that stars were negligible in the energy budget was wrong and was retracted.
2. The first distributed disks tore into clumps (Toomre instability); the live halo fixed it.
3. A heap corruption (Windows exit code 0xC0000374) in the batch snapshot writer: the snapshot buffer was sized for the visible bodies, but the fill function now also wrote the halo, so it wrote past the end. Found from the exit code, fixed, affected runs regenerated.
4. The cold-disk test failed because re-centring on the mean velocity shifted inner velocities; the check was restricted to the outer half.
5. The dt "ratchet" exists because a fixed step sized from the starting state can become too coarse when a dense clump forms or the black holes meet.
6. Screenshot tooling misreported the window because of a DPI mismatch (SetProcessDPIAware fix); the image viewer collapsed several images into one window. (Minor.)

### 3.12 Known limitations (the candidate should be able to list several)

2D, not 3D. Halo fraction is set by hand, not fit. Only a few thousand particles per component. Monopole-only tree: energy drift scales with theta^2. dt is chosen once and only ratcheted down. The default two-galaxy halos overlap from the start (halo radius 2R = 100,000 ly each against 150,000 ly separation). Initial virial ratio about 0.83 before relaxation (slightly cold). The relaxed disk is puffy, not a crisp thin disk. No gas, star formation or cooling. Instantaneous bound test. Reset redoes the relaxation. Kepler-mode star loss numbers are not physical.

---

## 4. Question bank (with answer keys)

Pick from these, mix the order, and invent follow-ups. "Key points" are what a good answer contains. "Common errors" are the traps. Difficulty: E easy, M medium, H hard.

### A. Overview and framing

**A1 (E).** In about a minute, what does this project do and what question were you asking?
*Key points:* two backends (centralized Kepler and distributed Barnes-Hut), the collision simulator, leapfrog, the metric of stars still bound, the direction: how starting parameters change star loss, and what a halo changes. *Common errors:* describing only the graphics; claiming it is "realistic".

**A2 (M).** Why build two modes instead of just the more accurate one?
*Key points:* the Kepler mode is O(N) and handles 100,000+ test stars cheaply; the N-body mode costs O(N log N) and has far fewer particles; the comparison itself is informative (the Kepler star loss is an artifact of no halo or self-gravity).

**A3 (M).** What does "accurate" mean here, and how would you know the sim is right?
*Key points:* energy and angular momentum conservation, convergence order test (ratio 4), tree versus direct-sum comparison, regression tests; admits no external reference data was used. *Common errors:* "it looks right".

### B. Integrator and timestep

**B1 (E).** Why is it called an integrator?
*Key points:* numerically integrates acceleration to velocity to position; the closed form only exists for two bodies; step by step approximation; Euler as the simplest case.

**B2 (M).** Why leapfrog and not Euler or RK4?
*Key points:* symplectic and time-reversible, energy error stays bounded, second order, one force evaluation per step (RK4 needs four). *Common errors:* saying leapfrog is "more accurate" without saying why energy behaves differently; claiming it conserves energy exactly.

**B3 (M).** Walk me through one kick-drift-kick step and say why the acceleration is cached.
*Key points:* half kick, drift, new forces, half kick; the new acceleration is reused as the next step's first kick so there is one force pass per step.

**B4 (H).** How is dt chosen, and what is the ratchet?
*Key points:* `dt = eta sqrt(eps / max|a|)` from the start, pericentre cap for two galaxies, ratchet halves dt when peak |a| is more than 4x the reference; only shrinks; why it exists (clumps, black hole meeting). *Probe:* "What is the downside of a ratchet that never grows?" (wasted compute after a transient; dt varies in time, which technically breaks the symplectic guarantee).

**B5 (H).** What would happen to E/E0 if you halved dt?
*Key points:* for leapfrog the error falls about 4x (second order); the tree error does not change with dt, only with theta; so total drift may be dominated by theta.

### C. Barnes-Hut tree

**C1 (E).** Why not just sum every pair?
*Key points:* O(N^2): 20,000 bodies means 4e8 pair evaluations per step; the tree is about O(N log N).

**C2 (M).** Explain how the quadtree is built.
*Key points:* recursive quadrant subdivision; insert one body at a time; split an occupied leaf into four; each node stores total mass and centre of mass; rebuilt every step. *Probe:* "What happens if two bodies are at the same position?" (merged leaf at depth 48 to avoid infinite subdivision.) Use `docs/diagrams/1_tree_structure.png` and `2_insertion_steps.png` as the reference.

**C3 (M).** How does the force walk decide whether to use a cell as a point mass?
*Key points:* `s/d < theta`, with s the cell width and d the distance to the cell's centre of mass; otherwise open the children; leaves are used directly. *Common errors:* confusing theta with an angle in degrees; forgetting the "target inside the cell" rule.

**C4 (H).** Why must a cell that contains the target never be accepted?
*Key points:* the body would be inside the approximated group; the distance to the centre of mass can be misleadingly large relative to the cell width; the force would be wrong (including pointing the wrong way for some configurations).

**C5 (M).** What happens as theta goes to 0, and as theta grows to 1 or 2?
*Key points:* theta 0 opens everything: exact but O(N^2); larger theta is faster and less accurate; error scales roughly as theta^2 for a monopole tree; energy drift follows.

**C6 (H).** What is the approximation error, and how could you reduce it without lowering theta?
*Key points:* monopole only; add a quadrupole moment per cell; FMM or tree-PM methods are alternatives (FMM is O(N)). *Probe:* "Why does lowering theta cost more than the extra accuracy suggests?"

**C7 (H).** Why loop over bodies in depth-first leaf order, and why is the loop parallel?
*Key points:* neighbouring bodies take nearly the same tree path, so the nodes stay in cache; independent per-body walks make it trivially parallel with OpenMP (dynamic scheduling because walk cost varies by region).

**C8 (M).** Is the tree-based E/E0 trustworthy?
*Key points:* it includes the tree's own error, so a drift of E/E0 conflates integrator error with tree error; an exact O(N^2) `total_energy_direct` exists for checks.

### D. Physics of the galaxy model

**D1 (M).** What is a dark-matter halo doing in this simulation?
*Key points:* holds most of the mass (90%), keeps the disk from fragmenting (Toomre stability), deepens the potential so stars are harder to strip, dynamical friction in mergers.

**D2 (H).** What is Toomre Q and why did the no-halo disk fall apart?
*Key points:* Q measures stability of a rotating disk against self-gravity (velocity dispersion and rotation versus surface density); Q below about 1 is unstable; the no-halo disk had Q about 0.08; fixed by a halo (or by large dispersion, which was too hot to be a disk).

**D3 (M).** How does the circular orbit speed get set, and why the measured force?
*Key points:* `v_c^2 = a_R r`; the spherical `G M/r^2` is wrong for flat mass distributions (about 1.5x too weak for a 1/R sheet); the code measures the in-plane force from a throwaway realisation.

**D4 (H).** The halo does not rotate. How does it avoid collapsing?
*Key points:* support comes from random motion (velocity dispersion acting as pressure), orbits have angular momentum so they have pericentres; dark matter is collisionless so it cannot lose energy and settles at virial balance. In 2D the setup uses the Jeans equation. *Probe:* "Does that mean all its orbits are extremely elongated?" (no: a mix from near-circular to radial, typical orbits moderately eccentric; precessing rosettes, not closed ellipses.)

**D5 (H).** The galaxy collapsed and rebounded at the start. What did you try, what worked, and what do you still not know?
*Key points:* the measured force (small effect); scaling halo velocities (no scale fixed it); pre-relaxation (worked); the root cause is **UNVERIFIED** (hypothesis: not a true equilibrium distribution function). A good answer states the hypothesis and says it was not proven.

**D6 (M).** What does pre-relaxation cost and change?
*Key points:* 45 s or more per galaxy start-up, redone on reset, the disk ends up puffier (about 1.5x its radius, a rounded cloud, not a thin spiral), residual dip of about 15%.

**D7 (M).** Why did a hole in the disk alone not allow a larger dt?
*Key points:* the largest accelerations were on halo particles near the black hole; the hole had to cover disk and halo; the ratchet later halves dt as radial-orbit halo particles return.

### E. Results and interpretation

**E1 (M).** The Kepler sweep showed lots of unbound stars. The N-body run shows very few. What does that tell you?
*Key points:* the Kepler stars are test particles with no halo and no self-gravity, so the loss is an artifact of the model; the halo deepens the potential, dynamical friction and merging keep stars bound; to confirm, sweep the halo fraction and compare (the sweep is the experiment). *Common errors:* declaring the N-body result "the truth".

**E2 (H).** What would you conclude, and not conclude, about dark matter from a halo-fraction sweep?
*Key points:* conclude that a massive extended halo stabilises disks and protects against stripping in this model; do **not** conclude the halo must be dark matter (MOND also explains flat rotation curves), that 90% is the right number (it is chosen, not fit), or anything about 3D structure, gas or star formation.

**E3 (M).** How did the energy ratio look in the two-galaxy runs, and what limits it?
*Key points:* about +3e-3 over 2 Gyr in the 90% run; tree error scales with theta^2; the ratchet; softening; honest that this is not machine-precision conservation.

**E4 (H).** Why is "bound" only an instantaneous test, and why does that matter in a collision?
*Key points:* the potential is time-dependent during an encounter; a star can be marginally unbound now and recaptured, or bound now and later stripped by a second pass; the fraction fluctuates; resolution is 1/N.

### F. Engineering, testing and tooling

**F1 (M).** How do you know the tree is correct?
*Key points:* compares to the exact O(N^2) sum; accuracy improves monotonically as theta shrinks; error under 2% rms at theta 0.5; degenerate-input and leaf-order tests.

**F2 (M).** Walk me through the regression test you are proudest of.
*Key points:* any of: halo versus no-halo fragmentation (clumpiness 1.53 versus 6.4, same disk), pre-relaxation (34% versus 13.5% swing), convergence order, isolated conservation. He should say what it protects against.

**F3 (H).** Tell me about a bug you hit.
*Key points:* the heap corruption (buffer sized without the halo; found via exit code 0xC0000374); or the E/E0 drift that looked like an integrator bug but was star energy gain, and how the probe showed it. A good answer includes how it was diagnosed.

**F4 (M).** Why a headless batch runner and CSV logs?
*Key points:* runs take minutes to hours; no graphics so it can sweep parameters in the background; the same CSV format as the GUI logger; plots made in Python.

**F5 (M).** The window froze during relaxation at first. How is that handled?
*Key points:* a progress callback from the long constructor pumps events and draws a progress frame; starting a session was moved out of the middle of an ImGui frame (set a pending flag, start at the top of the next loop).

**F6 (H).** What C++ features and design choices does the code use, and what would you change?
*Key points:* an abstract `ISimulation` interface with two implementations and a factory; structure-of-arrays for particles (x, y, vx, vy, m); OpenMP; the tree stores indices not copies; `std::function` progress callback. Changes he might suggest: cache the relaxation, a quadrupole term, adaptive per-particle timesteps. Do not require specific C++ vocabulary; check he understands the interface/implementation split.

### G. Ownership and honesty (include at least two every session)

**G1 (E).** What did you write and what did Claude Code write?
*Key points:* an honest, specific answer: Claude Code wrote the code; he set the goals, made the design decisions (the two-mode split, the halo, the hole parameter, pre-relaxation), caught the clumping, asked the physics questions, and reviewed the code and results. Saying or implying he hand-coded it is a red flag.

**G2 (M).** Name three design decisions that were yours, and one place where you overruled or redirected the tool.
*Key points:* specific and true. Examples from this project: choosing distributed mass plus Barnes-Hut as an option; choosing the live halo over the alternatives; asking for blue dark-matter display; requesting the central hole; choosing pre-relaxation; requesting E/E0 and bound-fraction metrics.

**G3 (H).** Pick a part of the code you did not write. Can you explain it well enough to modify it?
*Key points:* he should pick something and explain it correctly, or honestly say he could not yet. Follow with a concrete modification question ("how would you add a quadrupole term?", "how would you cache the relaxed state?").

**G4 (M).** How would you describe this project on a resume and to an interviewer who dislikes AI-written code?
*Key points:* "designed and directed, built with Claude Code", emphasis on validation (tests, convergence, energy checks), research direction, and ability to explain every part; offers to walk through any file; has his own hand-coded projects to show. Do not prompt him toward any particular wording: judge whether it is honest and specific.

**G5 (H).** What are the limitations of the project, and which would you fix first?
*Key points:* he should list several from section 3.12 and give a reasoned priority.

---

## 5. Common misconceptions to listen for

- "Barnes-Hut is exact" or "Barnes-Hut is O(N)". It is approximate and about O(N log N).
- "Leapfrog conserves energy exactly." It bounds the error; it does not eliminate it.
- "Theta is an angle in degrees." It is a dimensionless ratio threshold (cell width / distance).
- "The halo is dark matter, proven by this simulation." It is a chosen mass component; the sim does not distinguish dark matter from modified gravity.
- "More stars unbound means the physics is more realistic." In the Kepler mode it is a model artifact.
- "The dark-matter halo does not rotate so it must be static." It has random motion, orbits with angular momentum.
- "The 1% drift in two-galaxy Kepler mode means the integrator is bad." It is star energy gain.
- Confusing the time step dt with the Myr-per-second display speed.
- Claiming the pre-relaxation cause is proven.

---

## 6. Suggested closing format

```
Scorecard
- Overview and framing: Solid / Partial / Missed. One sentence on the gap.
- Integrator and timestep: ...
- Barnes-Hut tree: ...
- Galaxy model and dark matter: ...
- Results and interpretation: ...
- Engineering and testing: ...
- Ownership and honesty: ...

Three weakest topics, each with a 60-second model answer he can rehearse.
Factual errors he made (quoted) and the correction.
What to read or re-run before a real interview (point at specific files: src/sim/quadtree.cpp, src/sim/nbody_simulation.cpp, tests/test_nbody.cpp, the diagrams in docs/diagrams).
```

---

## 7. Halo-sweep results

Runs: `logs/halo_sweep/halo_{90,50,20,00}.csv`. Same seed, two galaxies, 2,500 disk particles each, 800 Myr relaxation, 2 Gyr, default encounter. Plot: `docs/diagrams/7_halo_sweep.png` (`tools/plot_halo_sweep.py`).

| Halo fraction | Bound at t=0 | Bound at 2 Gyr | Lost during run | E/E0 - 1 |
|---|---|---|---|---|
| 90% | 100.00% | 99.78% | 0.22 pt | +3.07e-3 |
| 50% | 99.50% | 99.40% | 0.10 pt | +1.65e-3 |
| 20% | 96.96% | 95.14% | 1.82 pt | +2.09e-3 |
| 0% | 94.62% | 89.40% | 5.22 pt | +2.20e-3 |

What the data supports: with a halo of 50% or more, almost nothing is lost; at 20% and 0% the loss rises sharply, and the 0% case loses the most and keeps losing for about 1.5 Gyr. Energy drift is a few 1e-3 in every run, so the differences are not an energy-conservation artifact.

What it does not support:
- **Low-halo runs start with 3 to 5% of disk particles already "unbound" at t=0.** So part of the final number is the starting state (the energy test applied to a relaxed, puffier disk), not the collision. The "lost during run" column is the cleaner comparison. The cause of the t=0 deficit is unverified.
- It is one seed and one encounter geometry, so there are no error bars.
- It is 2D, the halo fraction is chosen rather than fit, and a galaxy without a halo has a very different structure (Toomre Q about 0.08 with no halo), so the runs differ in more than depth of potential.
- It is a model comparison, not evidence that dark matter exists. The sim cannot distinguish dark matter from modified gravity.

Question to ask him: why does the 0% run start below 100% bound, and what would he change to separate start-up effects from collision effects?

---

## 8. Honest-framing guidance for the ownership questions

What is true and can be said without hedging:
- Claude Code wrote the code from his specifications and iteratively changed it in response to his requests and his questions.
- He made the design calls, reviewed the code and the results, and caught real problems (the clumping, the energy drift question, the breathing at start).
- The validation (tests, convergence order, tree-versus-direct, regression checks) is part of the project and he can explain it.

What he should **not** say: that he hand-coded the simulator, or that he wrote the tree or integrator himself. He has stated an intention to hand-write his own Barnes-Hut tree and integrator later, using this project's scaffolding and tests; that is a plan, not something already done.

---

## 9. Diagrams (for his study, optional to share in the chat)

In `docs/diagrams/`: `1_tree_structure.png` (space split and tree side by side), `2_insertion_steps.png` (building the tree one body at a time), `3_opening_angle.png` (which cells become point masses for theta 1.0, 0.5, 0.2), `4_cost_scaling.png` (interactions per body versus N), `5_timestep_pipeline.png` (one kick-drift-kick step plus the dt ratchet), `6_galaxy_model.png` (disk, halo, core, hole, mass budget). Regenerate with `tools/.venv/Scripts/python.exe tools/plot_quadtree.py`. If he uploads one, you may ask him to explain it, but do not describe it for him first.
