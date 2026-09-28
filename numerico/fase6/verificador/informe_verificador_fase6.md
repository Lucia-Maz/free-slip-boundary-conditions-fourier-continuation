# Independent verification of Phase 6: the order-N deformable surface in SPECTER

**Role:** adversarial verifier. **Date:** 2026-09-28.

**Code checked:** `src/boundary/fsorder.f90` as of 15:00, md5 `4dcba4d7…`. This version already has the new stopping rule (it ≥ 6, 1e2·fs_tol, best of the last 4 against the best before them).

**Also read:** the record `numerico/fase6/REGISTRO_fase6.md` (including P13), `DECISIONES.md` entries [D-51]–[P-04], and `verificacion/intent_fase6.txt`.

**What I did not do:** I did not run SPECTER, `make` or any gate. Every computation below is my own light Python, in
`/tmp/claude-1000/-home-lucia-GW-AI-proyecto-final-piv/cfa53625-6433-4683-95b6-41a034e9a63e/scratchpad/verificador/`.

**How the runs were done:**
- One run at a time, each under `ulimit -v 2500000`.
- Peak memory 155 MB (`c3b`); everything else stayed under 45 MB.
- The repository was not modified. The author's replica was imported with bytecode writing off and its outputs redirected to my scratch. The repo's `salidas/` files keep their 14:44 timestamps.
- sympy and mpmath are **not installed**. The "symbolic" checks use exact rational power-series arithmetic (`fractions.Fraction`) instead. Those comparisons are exact, not floating point.

## Summary of verdicts

| claim | verdict |
|---|---|
| C1 model = Taylor polynomials | **verified** (exact). Holds for the header formulas and for a transliteration of `fs_functionals`, under amplitude counting |
| C2, N = 2: spurious mode λ = ν/η² for η < 0 | **verified**. The local analysis survives the full 2D problem: the root is exactly λ = ν(1/η0² − k²) |
| C2, N = 3: neutral modes | **verified for uniform η**; **refuted with slope terms**. A frozen slope s0 gives Re λ ≈ ν(2s0²/η0² − k²) > 0 |
| C3 replica fidelity (a, b, c) and reproduction | **verified**: exact reproduction, and my independent step agrees |
| C3, "N = 2 unstable only for η/dz ≤ −0.70" | **refuted as a threshold**. The onset depends on ν·dt/dz²; at the gate's parameters it is −0.6665 |
| C4 depth-shift degeneracy O(η0³) | **verified**, inviscid (exact) and viscous (numerical) |
| C5 preconditioner signs and fixed points | **verified** |
| C6 real-coefficient Anderson | **verified** |
| C7 Nyquist convention | **verified** |
| C8 guard arithmetic | **verified**; but the guard margin is insufficient (bug B1) |
| C9 roundoff-floor estimate | **refuted**: a 1/nz is missing. The origin of the floor stays unexplained to within about 10× |

---

## C1: the model

**Method** (`c1_modelo.py`, `c1b_codigo.py`). I derived the exact conditions myself, for a graph surface with upward normal n ∝ (−s, 1) and tangents t_a = e_a + s_a e_z:
- **Kinematic:** η_t = w − s·u_h (material surface).
- **Tangential:** t_a·S·n = 0, i.e. S_az − s_b S_ab + s_a S_zz − s_a s_b S_bz = 0.
- **Normal:** p = gη − γκ + ν n·S·n/|n|², with κ = [(1+s_y²)η_xx − 2 s_x s_y η_xy + (1+s_x²)η_yy]/(1+|s|²)^{3/2}. This is Young–Laplace, with p the kinematic pressure minus the rest hydrostatics. It agrees with the Wang–Tice–Kim form quoted in the task.

Every field is evaluated at z = h + η as its full Taylor series. The local data are random rationals, and everything is truncated power series in ε with Fraction coefficients. I compared, coefficient by coefficient, against:
- (i) the header formulas, implemented from the text;
- (ii) a line-by-line transliteration of the pointwise algebra of `fs_functionals`: `trz`, `gph(:,:,1..6)`, B_N, and the κ flux for N ≥ 3.

**Results:**
- **Amplitude counting** (η, its derivatives and all fields ∝ ε), 3D and 2D, N = 1..4: the first mismatch is at degree N+1 in every one of K, T_x, T_y and Φ. Degrees 0..N agree **exactly** (zero difference in rationals).
- The same holds for the transliterated Fortran algebra (N = 1, 2, 3).
- **Negative control:** flipping the sign of the s_b S_ab term makes T_x and T_y fail at degree 2. So the slope signs are right.

**Wording caveat.** Under "deformation-only" counting (η ~ ε, velocity O(1)), the formulas are Taylor polynomials of degree **N−1** in η for every velocity-dependent term: the first mismatch is at degree N. The curvature is kept to degree N (odd N).

This is exactly equivalent to "order N in the amplitude", which is what the header says. But the record (§1 "desarrolladas… hasta orden N" in η; [D-51] "desarrollo en la deformación con el campo de velocidad completo") reads as if the ηᴺ term were kept. It is not: the truncation error is O(ηᴺ ∂ᴺv) relative to the leading term.

**Verdict: verified.**

## C2: N = 2 ill-posed for η < 0; N = 3 neutral

**Method** (`c2_continuo.py`, `c2b_pendiente.py`).

**(i) Scalar model.** The condition is κ·P_{N−1}(ηκ) = 0, with P the truncated exponential.
- N = 2: κ = −1/η, admissible (Re κ > 0) only for η < 0, with λ = ν/η² > 0. There is no admissible mode for η > 0.
- N = 3: ηκ = −1 ± i and λη²/ν = ∓2i (neutral).
- By-product: N = 4 is unstable again. P₃ has a real root at x = −1.596, giving λη²/ν = +2.547. This holds for every even N, which is relevant to [P-04].

**(ii) Full linearized 2D problem.**
- Setup: a half-space below z = h; base at rest with frozen elevation η0 and frozen slope s0; perturbation e^{ikx+λt}.
- Modes: a potential mode A e^{k(z−h)} with p = −λφ, and a rotational mode B e^{m(z−h)} with m² = k² + λ/ν.
- Rows: the order-N kinematic, tangential (with its slope terms) and normal (κ_N and B_N) conditions, all linearized.
- Roots found by Newton; the number of roots with Re λ > 0 counted by the argument principle on a rectangle [1e-3, R] × [−R, R], with R = 50ν/η0².

**N = 2 results:**
- For uniform η0 the root is **exactly** λ = ν(1/η0² − k²), with relative deviation 0 for k = √2 and 10, and η0 = −0.01 and −0.04.
- Reason: Tr₁ = (1 + η0∂z) annihilates e^{−(z−h)/η0}. So that rotational mode is invisible to **all three** transferred rows, not only the tangential one. The ∂x w term, continuity and the normal/kinematic coupling cannot remove it.
- A frozen slope s0 = 0.05–0.1 shifts the root by 0.5–11% only (e.g. 99.99 − 0.56i at η0 = −0.01).
- Argument principle: 1 root in the right half-plane for η0 < 0, 0 for η0 > 0, with and without slope.

**Nuance for N = 2.** For fixed η0 the growth rate is bounded (by ν/η0²), so this is not Hadamard ill-posedness. It is a spurious instability whose rate diverges as η → 0⁻. That is why it appears as soon as the grid resolves |η|, and why N = 2 with η < 0 cannot converge under z-refinement. The record's "mal planteado" is fair in that loose sense.

**N = 3 results:**
- Uniform η0: exact roots λ = ∓2iν/η0² − νk², i.e. damped by νk². No right-half-plane roots.
- **With a frozen slope, the pair becomes unstable:**

| k | η0 | s0 | Re λ (root) | estimate ν(2s0²/η0² − k²) |
|---|---|---|---|---|
| √2 | −0.01 | 0.05 | +0.477 | 0.480 |
| √2 | −0.01 | 0.10 | +1.975 | 1.980 |
| √2 | −0.02 | 0.05 | +0.102 | 0.105 |
| √2 | −0.04 | 0.05 | +0.008 | 0.011 |
| 3 | −0.01 | 0.05 | +0.40 | 0.41 |

- Argument principle: 2 roots in the right half-plane for η0 < 0 with s0 = 0.05.
- Mechanism: the −s_a s_b Tr₀[S_zb] term shifts the zero of P₂(mη0) by about s0².

So N = 3's neutrality does **not** survive the slope terms of the full tangential condition. The growth is weak (rate ~ νs0²/η0², versus N = 2's ν/η0²).

**Limits of this result:**
- The unstable horizontal wavenumbers satisfy k < √2|s0/η0|. That is the scale on which η itself varies, so frozen coefficients are only marginally justified there. Whether real waves excite this is **unclear**.
- The author's 1D replica cannot see it: it has no slope and no horizontal wavenumber.
- It would equally affect the [P-04] idea of taking only the tangential transfer to order 3.

**Verdicts:** N = 2 part **verified**; "N = 3 neutral survives the full 2D problem" **refuted** in the local (frozen-slope) analysis, with practical impact unclear.

## C3: replica fidelity and the stability map

**Method** (`c3_replica.py`, `c3b_main_y_umbral.py`).

**(a) Index conventions.** I wrote a literal 1-based transliteration of the Fortran z-pipeline:
- `load_dirichlet_tables`: A·Qᵀ after `Q = TRANSPOSE(Q)`;
- `load_neumann_tables`: Q(d,:)·Qnᵀ with Q not transposed, and neu(d)·dx/dxp;
- the `boun .eq. 6` branch of `neumann_reconstruct`, weights on f(nz−Cz−d+k);
- the continuation loop of `fftp1d_real_to_complex_z`: dir(ii,jj)·f(n−C−d+jj) + dir(C−ii+1,jj)·f(d−jj+1);
- the kz grid of `specter.fpp`;
- `fs_ph` with the odd-m Nyquist drop;
- the `fs_gen_setup` probe.

I also confirmed the plan sign: the forward z plan has sign −1, so the +i phases in `fs_ph` are the matching inverse.

**Agreement with the replica:**
- Neumann weights identical (difference 0).
- D2 and the trace weights agree to 1e-13 and 2e-12 relative (roundoff).
- The neu weights are (−0.12, 0.64, −1.44, 1.92) with neu(d) = 0.48 dz.

**(b) RK2 structure.** `hd_rkstep1.f90` copies u⁰ into C1–C3. `specter.fpp` loops `DO o = ord,1,-1` with ORD = 2. `hd_rkstep2.f90` does u^(o) = C + dt/o·RHS(current u), then `v_imposebc_and_project`. The replica's `for o in (2,1): u = bc(u0 + DT/o*NU*D2@u)` is the same thing.

**(c) Trace response to a unit datum.** β_m·dz^(m−1) = 0.4800, 1.0052, 1.3967, 1.1515 at nz = 128 and 256, identical from the literal transliteration and from the replica. At nz = 64, β₃dz² = 1.1517 (the record says "equal"; the difference is trivial). The coordinator's quote of the Fortran print matches, but I did not see the print myself.

Consistency check: at η0 = 0.04, σ = β₂η + β₃η²/2 = 5.698 + 9.584 = 15.28. That equals the measured gain −15.28 (quoted).

**Reproduction:**
- The author's `main()`, run with outputs redirected, gives ρ **bit-identical** to the repo JSON (max |Δρ| = 0).
- My own RK2 step, built from the transliterated operators, reproduces the replica's ρ to 8 digits at 16 test points.
- q = 8 interior extrapolation: ρ = 2.770 at η = dz and 39.27 at 4dz.
- Stable sets on the author's grid:
  - q8, N = 2: [−0.4, 0.5], so "unstable for |η| > 0.5dz" ✓.
  - q8, N = 3: [−1.75, 0.3].
  - N = 3 spectral: stable everywhere.

**Break-tests**, with D = ν·dt/dz² the diffusion number (the scaled problem depends only on D and r = η/dz; ν = 0.1, dt = 1e-4 reproduces ν = 0.01, dt = 1e-3 exactly):
- **N = 3 spectral: stable for r ∈ [−8, 9]** at every stable D tried: nz = 64, 128, 256; dt from 1e-4 to 1e-3; ν = 1e-3, 1e-2, 1e-1. **Verified.**
- The nz = 256, dt = 1e-3 and nz = 128, dt = 2e-3 cases are meaningless: the base diffusion step is itself unstable there (ρ(η=0) = 9.33 and 1.033).
- **The N = 2 onset is not at −0.70, and it moves with D.** Bisection with my step:

| D | nz, dt | unstable for r ≤ | 1+σ at onset |
|---|---|---|---|
| 0.010 | 128, 1e-4 | −0.7114 | 0.006 |
| 0.026 | 128, 2.5e-4 | −0.7045 | 0.016 |
| 0.052 | 128, 5e-4 | −0.6925 | 0.033 |
| **0.104** | **128, 1e-3 (gate)** | **−0.6665** | **0.069** |
| 0.130 | 128, 1.25e-3 | −0.6521 | 0.089 |
| 0.132 | 256, 2.5e-4 (record's nz = 256 runs) | −0.6509 | 0.091 |
| 0.177 | 128, 1.7e-3 | −0.6219 | 0.131 |

  On the author's 0.1-step grid this shows up as "≤ −0.70". [H-24]'s "en la malla el cruce está donde 1 + σ = 0, en η = −0,716 Δz" and §1 / [D-51]'s "N = 2 estable con η > −0,72 Δz" are true only as D → 0.
- **Perturbed Neumann weight:** scaling neu(d) by 0.9 or 1.1 moves −1/β₂ to −0.796 or −0.651, and the onset tracks it (−0.76, −0.62). The replica's logic is consistent.
- **Minor:** the record says the replica measures the continuum rate "exactamente" when the layer is resolved. In rate, ln ρ versus ν dt/η² differs by 2.4%, 4.9% and 14% at r = −4, −3, −2.

**Verdicts:** fidelity and reproduction **verified**; the N = 2 threshold statement **refuted** (see B1).

## C4: depth-shift degeneracy

**Inviscid, exact** (`c4_profundidad.py`, Fraction series, random rational T = tanh kh). The transferred relation is R_N = (T + y_N)/(1 + T y_N), with y_N = odd/even truncations of sinh/cosh. Transferred minus exact tanh(kh + x):
- N = 1: −(1−T²)x;
- N = 2: **+(1−T²)x³/3**;
- N = 3: **−(1−T²)x³/6**;
- N = 4: +(1−T²)x⁵/30.

So N = 2 and N = 3 are both O(η0³), and N = 3 has half the error with the opposite sign.

The numerator of S₂/C₂ − tanh(kh+x) is x cosh x − sinh x = x³/3 + x⁵/30 + …. **[H-25] writes "−(ε³/2)/(…)"**: the coefficient is wrong (1/3, not 1/2), and so is the sign for "transferred − exact". The conclusion is unaffected (B2).

**Viscous, numerical.** Finite depth h = 1, no-slip bottom. A 5×5 determinant with modes e^{k(z−h)}, e^{−kz}, e^{m(z−h)}, e^{−mz}. Order N means each mode's top entries are multiplied by P_{N−1}(μη0); the exact problem at depth h + η0 means multiplying by e^{μη0}. Parameters: ν = 0.01, k = √2, g = 1, γ = 0.2. The exact root at depth 1 is λ = −0.05325 − 1.30999i.

|λ_N − λ_exact(h+η0)| for η0 = 0.0025, 0.005, 0.01, 0.02, 0.04:

| N | errors | log₂ slopes |
|---|---|---|
| 1 | 6.5e-4 … 9.9e-3 | 1.00 … 0.96 |
| 2 | 2.66e-9, 2.12e-8, 1.68e-7, 1.31e-6, 1.00e-5 | 2.99, 2.98, 2.97, 2.93 |
| 3 | about half of N = 2 | 2.99 … 2.93 |

- Viscosity changes the magnitude by 15–20% relative to the inviscid |δω|, but not the order.
- Even at ν = 0.1, where the boundary layers overlap (e^{−Re m·h} = 0.10), the slope is 2.98.
- Why it holds: the top rotational mode's factor is a column scaling, and the ±k factors agree through O(x²).

**Verdict: verified**, and it holds with viscosity at these parameters.

## C5: preconditioner

**Code read of the signs.** The datum is ikx(w* − 2W) − Dx, and the traces respond as β_m. So:
- ∂F_Dx/∂Dx_in = −(β₂η + β₃η²/2) = −σ;
- ∂F_Dx/∂W = −2σ ikx.

The code computes `ph(5) = (F − 2σ·ph(1) + σ·ph(2))/(1+σ)`, with ph(1) = to_phys(ikx(W_new − W)), ph(2) = to_phys(Dx_in), ph(3) = Dy_in and ph(4) = ∂yΔW; likewise for y. The signs are consistent.

**Fixed-point equivalence** (my own argument). At a fixed point of G, ΔW = 0. Write e = F − Dx ∈ V. Then G = P[Dx + e/(1+σ)] = Dx + P[e/(1+σ)]. So G's fixed points require P[e/(1+σ)] = 0, i.e. Σ_grid |e|²/(1+σ) = ⟨e, P[e/(1+σ)]⟩ = 0. This uses that P (`fs_to_spec`) is the orthogonal projection in the padded-grid inner product. When 1+σ > 0 at every padded-grid point, this forces e = 0, so the fixed points coincide.

The projection is what makes positivity on the **padded grid** the right condition, and that is exactly where the guard evaluates min(1+σ).

**Numerics** (`c5_c6_c7.py`, literal `fs_to_phys` / `fs_to_spec`):
- Round trip exact (8e-16).
- L = P[·/(1+σ)] on V has real positive spectrum: minimum eigenvalue 0.28–0.83 for N = 3 up to 4dz and for N = 2 within ±0.6dz.
- With a sign-changing σ, L has negative eigenvalues and becomes nearly singular for some scaling (σ_min/σ_max = 5.7e-6), so G can have fixed points F does not.
- On a model pass with the code's datum convention, Picard on G converges. With the W-term sign flipped it diverges or stagnates.

**Also checked:** the padded grid 2nx de-aliases the cubic products of N = 3 exactly (M = 4K + 4 > 4K with K = nx/2 − 1).

**Verdict: verified.**

## C6: Anderson with real coefficients

`fs_anderson` implements:
- Gm = Re(dRᴴdR) and bv = Re(dRᴴr), both MPI-SUMmed;
- a 1e-13·trace Tikhonov term;
- a real solve (`fs_rsolve`, partial pivoting; a singular pivot falls back to plain Picard);
- x_new = F − Σγ_m dF_m (type II).

History: slot 1 is the newest. dR₁ = r − ahr₁ and dR_m = ahr_{m−1} − ahr_m, with mk = min(nhist, mhist). fs_nhist is reset every substage, so only this substage's slots are used. The pushed F is the preconditioned G, consistently.

**Numerics:**
- A literal transliteration reproduces a textbook type-II Anderson (real least squares in ℝ^{2n}) to 9.4e-14 over 18 passes, on an ℝ-linear, non-ℂ-linear map x ↦ Ax + B x̄ + c.
- On a contractive such map, after 18 passes: real γ reaches 8e-5, complex γ 5e-4.
- No index or sign bugs found.

**Verdict: verified.**

## C7: Nyquist convention

- z(top) = (nz − Cz − 1)·dz = Lz exactly: node 102 (0-based) at nz = 128, z/dz = 102.
- For a real column, the symmetric real interpolant (Nyquist term (F_{N/2}/N)·cos(πz/dz)) has at the node exactly the derivatives given by `fs_ph` with the Nyquist term dropped for odd m: m = 1 gives 145.0392 in both, m = 3 gives −7.480065e6 in both. Even m keeps the Nyquist term.
- Keeping kz = −nz/2 for odd m adds −49i (m = 1) and 5.1e6·i (m = 3).
- For complex columns of a real 3D field, keeping it breaks the Hermitian pairing of the (0, ±ky) traces: 6.7e6 versus 5e-9 when dropped.

**Verdict: verified.**

## C8: guard arithmetic

- With β₂dz = 1.3967 and β₃dz² = 1.1515:
  - 1 + σ = 0 at η = −0.71597 dz;
  - the 0.05 guard fires at η ≤ −0.68017 dz;
  - at N = 3, min(1 + σ) = 1 − β₂²/(2β₃) = 0.15294, at η = −1.2129 dz.
- σ uses the physical η: `fs_to_phys` divides by nx·ny, and the round trip is exact.

**Verdict: arithmetic verified.** The guard's **adequacy** fails; see B1.

## C9: the roundoff floor at N = 3

**The record's arithmetic.** ε·max|û|·kz³ = 2.2e-16 × 4.6e-5 × 3.29e7 = **3.4e-13**, written as "~1e-12". It **omits the 1/nz** in `fs_ph`: the trace is Σ_k û_k(ikz)³e^{ikz·Lz}/nz.

**With the 1/nz**, white per-coefficient noise ε·max|û| gives ε·max|û|·(Σ_k|fs_ph(k,3)|²)^{1/2} = **1.1e-14** rms (coherent upper bound ε·max|û|·Σ|fs_ph| = 8.1e-14). Times η²/2 = 8e-4, that is 9e-18 (at most 6.5e-17) in Dx. That is 12–90× below the reported 8e-16. So the stated estimate is not consistent to an order of magnitude, and the "exactamente el piso" agreement is fortuitous.

**Measured noise** (`c9_redondeo.py`):
- Pipeline: the literal Fortran z-pipeline (`neumann_reconstruct` → continuation loop → FFT → trace sum), on a wave-like column calibrated to max|û| = 4.6e-5.
- Measurement: rms of second differences of the m = 3 trace under ±h perturbations of the datum, h = 1e-9 to 1e-3 |g|.
- Result: **2.5e-14 to 9e-14**. For comparison, the m = 2 trace gives 2e-17.
- The dominant source is the **FC continuation's cancellation**: max row sum of |dir| = 9.5e3, max |dir| = 3.5e3. With only the continuation done in extended precision, the m = 3 noise falls to 1–4e-15.
- Implied Dx noise: about 2–7e-17 rms, still about 10–40× below 8e-16. Max-versus-rms over samples could close part of this.
- The quoted Σ_k|vx_k fs_ph(k,3)| = 0.751 only implies summation roundoff ε·0.751 ≈ 1.7e-16 in the trace. That is negligible and does not rescue the estimate.

**Verdict: refuted as stated.**
- The attribution to the m = 3 trace is plausible only to within about one order of magnitude.
- The real amplifier in my replica is the conditioning of the continuation, not "white coefficient roundoff".
- The floor itself is not explained by the given arithmetic.
- P13's new estimate, ε·kz³·(η²/2)·(nz/k)·C ≈ 2e-9, has a free constant C, so it cannot be checked independently.

---

## Adversarial code read (`fs_general_imposebc`, `fs_general_pass`, `fs_anderson`)

**Checked and correct:**
- Module state between passes: `fs_s1..3` are saved after the bottom BC, with vx and vy in the (z,ky,kx) domain and vz in (kz,ky,kx), and restored at every pass. `pr` is OUT from `sol_project_fs` on every pass. The bottom noslip uses the previous substage's `pr` before saving, consistently.
- The accepted state is really the output of the last pass, including at fs_maxit: fs_trn and pr come from that pass, and so does `fs_ptr_guess = fx(4)/fac`. `fs_wtop` takes the actual trace, not W_new.
- `fs_maxit ≥ 2` is enforced in `v_setup`.
- MPI: `fs_to_phys` ALLREDUCE; esc and res MAX; Anderson Gram SUM; the redundant solve is identical on all ranks; σ and the guard are identical on all ranks without a reduction.
- Array sections passed to explicit-shape dummies are contiguous.
- `dder` for m = 0..2 suffices (the pressure transfer uses m ≤ N−1 ≤ 2).
- The fs_eta / fs_etap roles (tangential at η^(o), explicit E and pressure transfer at η^(o+1)) are consistent with the header.
- The esc(4) mean exclusion is correct.
- The new stopping rule's indices are valid (it ≥ 6: rhist(it−3:it) against rhist(1:it−4)).

**Found:**
- **B1 (moderate), guard margin.** The 0.05 guard does not cover the discrete instability of the N = 2 implicit closure.
  - At the gate's parameters (nz = 128, ν = 0.01, dt = 1e-3, D = 0.104), the replica step is unstable for η ∈ (−0.716, −0.6665] dz, with ρ = 1.018, 1.145, 1.44 and 1.90 per step at −0.667, −0.670, −0.675 and −0.680 dz.
  - The guard fires only at η ≤ −0.6802 dz. The window (−0.6802, −0.6665] dz is therefore **unguarded** and blows up within tens of steps.
  - It widens with D: (−0.680, −0.651] at nz = 256, dt = 2.5e-4; (−0.680, −0.622] at D = 0.177.
  - The onset in 1+σ grows from 0.006 to 0.13 as D goes from 0.01 to 0.18. A safe guard needs 1+σ ≳ 0.15, or a D-dependent threshold.
  - None of the record's SPECTER tests falls in the window (−0.51dz was stable; −1.02dz was aborted), so there is no contradiction with its data.
  - Irrelevant for the cell (η/dz ~ 1e-3).
- **B2 (documentation), [H-25].** The residual is +(ε³/3 + ε⁵/30 + …)/(C₂·cosh(kh+ε)) for transferred − exact, not "−(ε³/2)/(…)".
- **B3 (documentation), header of fsorder.f90.**
  - Lines 36–39 state G(x) = (F(x) + σx)/(1+σ). This omits the −2σ∂x(W_new − W) term that the code (≈ lines 707–724) and the record use.
  - Lines 43–47 say the run is "stopped" at −0.72 dz, while the guard is at −0.68 dz and the real onset is at −0.667 dz at the gate's D.
  - The abort message ("ill-posed where eta < -0.7 dz") conflates the continuum (any η < 0) with the grid onset.
- **B4 (minor), silent non-convergence.** Hitting fs_maxit without convergence produces no warning. Only column 12 of the diagnostic shows it, and only for the substage just before each cstep.

**Note, not a bug:** the stagnation exit keeps the last pass rather than the best of the last four. It is bounded by res ≤ 1e2·fs_tol, which is by design.

```claims
[{"status": "verified", "text": "C1 The header formulas and a transliteration of fs_functionals' pointwise algebra are exactly the degree-N Taylor polynomials of the exact graph conditions under amplitude counting (exact Fraction series, 3D and 2D, N=1..4; first mismatch at degree N+1; slope-sign negative control fails at degree 2). Under deformation-only counting the velocity terms are degree N-1 in eta, so the record's 'hasta orden N en eta' is loose wording."},
 {"status": "verified", "text": "C2 (N=2) Scalar model: kappa=-1/eta, lambda=nu/eta^2>0 only for eta<0. Full linearized 2D problem: root exactly lambda=nu(1/eta0^2-k^2), because Tr_1 annihilates e^{-(z-h)/eta0} in all three rows; argument principle gives 1 unstable root for eta0<0 and 0 for eta0>0; frozen slope 0.05-0.1 shifts it by <=11%. Not Hadamard ill-posedness (the rate is bounded for fixed eta0) but a spurious instability whose rate diverges as eta->0-."},
 {"status": "refuted", "text": "C2 (N=3) Neutrality holds in the scalar model (lambda=-/+2i nu/eta^2) and for uniform eta0 in 2D (lambda=-/+2i nu/eta0^2 - nu k^2). It does not survive the slope terms: with a frozen slope s0, Re(lambda) ~ nu(2 s0^2/eta0^2 - k^2) > 0 (e.g. +0.48 at eta0=-0.01, s0=0.05, k=sqrt2; 2 roots in the right half-plane). Growth is weak, and whether real waves excite it is unclear, since frozen coefficients are marginal at the unstable k."},
 {"status": "verified", "text": "C3 (fidelity) Literal 1-based transliteration of load_*_tables, neumann_reconstruct(boun=6), the continuation loop, the kz grid, fs_ph and the probe matches the replica (weights identical, D2 and traces to 1e-12); RK2 structure matches hd_rkstep1/2 with ORD=2; beta_m dz^(m-1)=(0.4800,1.0052,1.3967,1.1515) at nz=128 and 256 (1.1517 at nz=64); the author's main() reproduces the repo JSON bit-for-bit; my own step matches rho to 8 digits; q=8 gives rho 2.770 and 39.27; N=3 spectral stable for r in [-8,9] at every stable nu dt/dz^2."},
 {"status": "refuted", "text": "C3 (threshold) 'N=2 unstable only for eta/dz <= -0.70' is an artifact of the 0.1-step grid. The onset depends on D=nu dt/dz^2: -0.7114 (D=0.01), -0.6665 (D=0.104, the gate's parameters), -0.6509 (nz=256, dt=2.5e-4), -0.6219 (D=0.177). [H-24]'s 'crossing at 1+sigma=0, -0.716 dz' and '[D-51] stable for eta > -0.72 dz' hold only as D->0."},
 {"status": "verified", "text": "C4 Exact series: transferred minus exact is +(1-T^2)x^3/3 for N=2 and -(1-T^2)x^3/6 for N=3 (N=4: x^5/30), so both are O(eta0^3). Viscous finite-depth dispersion relation (nu=0.01, k=sqrt2, g=1, gam=0.2, h=1, no-slip bottom): |lambda_N - lambda_exact| slopes 2.99-2.93 for N=2 and N=3, N=3 about half of N=2; slope still 2.98 at nu=0.1."},
 {"status": "verified", "text": "C5 Signs are consistent: the datum's -Dx gives gain -sigma, compensated by +sigma Dx; the datum's -2ikxW gives -2 sigma dx(W_new-W); ph(1..6) are used correctly. Fixed points coincide when 1+sigma>0 on the padded grid (P is the orthogonal grid projection, so P[e/(1+sigma)]=0 implies e=0); numerically L=P[./(1+sigma)] is positive definite (min eigenvalue 0.28-0.83) and becomes indefinite or near-singular when 1+sigma changes sign."},
 {"status": "verified", "text": "C6 fs_anderson solves Re(dR^H dR) gamma = Re(dR^H r) with 1e-13 trace regularization and a real solve, using the type-II update; the history differences and indexing are correct; a literal transliteration equals textbook real-gamma Anderson to 9.4e-14; complex gamma converges worse (5e-4 vs 8e-5) on an R-linear map."},
 {"status": "verified", "text": "C7 z(top)=Lz exactly (node nz-Cz-1). The symmetric real interpolant's odd derivatives at that node equal fs_ph with the Nyquist term dropped; keeping kz=-nz/2 adds spurious imaginary parts (5.1e6 i at m=3) and breaks the Hermitian pairing of the (0,+-ky) traces."},
 {"status": "verified", "text": "C8 Arithmetic: 1+sigma=0 at -0.71597 dz, the 0.05 guard fires at eta <= -0.68017 dz, and at N=3 min(1+sigma)=0.15294 at -1.2129 dz; sigma uses the physical eta (fs_to_phys divides by nx*ny); the quoted gain 15.28 equals sigma at eta0=0.04."},
 {"status": "refuted", "text": "C9 The estimate eps*max|u|*kz^3 = 3.4e-13 (written ~1e-12) omits the 1/nz in fs_ph. With it, white coefficient noise gives 1.1e-14 (at most 8e-14) in the m=3 trace and 9e-18 to 6.5e-17 in Dx, 12-90x below the 8e-16 floor. Measured replica-pipeline noise is 2.5-9e-14, dominated by the FC continuation's cancellation (row sums |dir|~9.5e3), giving Dx ~2-7e-17. So the floor's origin is not established."},
 {"status": "refuted", "text": "B1 (moderate) The guard 1+sigma<=0.05 leaves an unguarded unstable window at N=2. At the gate's parameters (D=0.104) the step is unstable for eta in (-0.716,-0.6665] dz (rho up to 1.90 per step at -0.680 dz), but the guard fires only at <= -0.6802 dz. The window widens with dt (to -0.651 at nz=256, dt=2.5e-4, and -0.622 at D=0.177). It needs 1+sigma >~ 0.15 or a D-dependent threshold."},
 {"status": "refuted", "text": "B2 (doc) [H-25] states the residual as -(eps^3/2)/(...). The exact numerator is eps*cosh(eps) - sinh(eps) = eps^3/3 + eps^5/30, so transferred minus exact = +(eps^3/3)/(...). The O(eta0^3) conclusion stands."},
 {"status": "refuted", "text": "B3 (doc) The fsorder.f90 header states G=(F+sigma x)/(1+sigma), omitting the -2 sigma dx(W_new-W) term the code uses, and says the run is stopped at -0.72 dz, whereas the guard is at -0.68 dz and the real onset at the gate's dt is -0.667 dz."},
 {"status": "unclear", "text": "B4 (minor) Hitting fs_maxit without convergence is silent: no warning, and diagnostic column 12 only reports the substage just before each cstep. It is a robustness gap, not demonstrated to have caused any wrong result."}]
```
