# Phase 6 verification, round 2

**Role:** adversarial verifier. **Date:** 2026-09-28, 16:40–17:00.

**Versions checked:**

| file | time | md5 (prefix) |
|---|---|---|
| `src/boundary/fsorder.f90` | 16:34:35 | `fe8e4706` |
| `boundary_mod.fpp` | 16:34:35 | `7cfbcd52` |
| `numerico/fase6/REGISTRO_fase6.md` | 16:54:41 | `77f79509` |
| `informe/superficie_orden_N.tex` | 16:54:41 | `17f34894` |
| `informe/superficie_orden_N.pdf` | 16:55:45 | — |
| `DECISIONES.md` | 16:36:59 | `e8e5d9a0` |

The documents changed while I was checking, so I re-grepped every item flagged below at 16:58–16:59.

**Constraints followed:** no SPECTER, `make` or gates; no repository edits; one light Python process at a time under `ulimit -v 2500000` (peak 40 MB).

**New scripts** (in this folder): `r2_guardia.py` (+ `.out`), `r2_n3_limite.py`, `r2_ord.py`, `r2_orden_trazas.py`, `r2_orden_poli.py`. They reuse my round-1 step (`c3_replica.MiPaso`), vectorized over columns and cross-checked against it to 1e-10 at 10 points.

## Summary

| item | verdict |
|---|---|
| 1a guard covers the N = 2 window up to the explicit limit | **verified** for N = 2 (RK orders 2–4). **New gap at N = 3** within ~1.5% of the explicit limit (N1) |
| 1b quantities `fs_nu`, `fs_dt`, `z(2)−z(1)` | **verified** |
| 1c no firing in the physical regime or at N = 3 | **verified** |
| 2 header | **verified** (one inherited nit) |
| 3 fsmaxit counter | **verified**; small residual gap (last output interval) |
| 4 documentation | **not clean**: two refuted statements remain, one stale guard sentence, gate result claimed before the run finished, minor issues |
| N2 (new) "∂z³u trace error is O(Δz²)" | **refuted beyond nz ≈ 128**: the error grows as Δz⁻³, and it is a deterministic defect of the FC tables |

---

## 1a. Guard coverage (N = 2)

**Largest stable D = ν·dt/dz² at η = 0** (RK2 with this FC-Gram Laplacian):

| nz | D_max |
|---|---|
| 64 | 0.2103 |
| 128 | 0.2048 |
| 256 | 0.2034 |
| limit | 2/π² = 0.2026 |

**Onset versus guard at nz = 128.** The unstable set in r = η/dz was scanned in steps of 0.005 on [−0.715, 0), and at 0.05–0.5 on (0, 12]. Its upper edge was then refined by bisection. The guard fires for 1 + σ ≤ sthr = 0.05 + D.

| D/D_max | D | unstable r | 1+σ at onset | sthr | margin |
|---|---|---|---|---|---|
| 0.05 | 0.0102 | [−0.715, −0.7115] | 0.0062 | 0.0602 | 0.054 |
| 0.10 | 0.0205 | [−0.715, −0.7070] | 0.0126 | 0.0705 | 0.058 |
| 0.25 | 0.0512 | [−0.715, −0.6929] | 0.0323 | 0.1012 | 0.069 |
| 0.50 | 0.1024 | [−0.715, −0.6674] | 0.0679 | 0.1524 | 0.085 |
| 0.70 | 0.1434 | [−0.715, −0.6443] | 0.1001 | 0.1934 | 0.093 |
| 0.85 | 0.1741 | [−0.715, −0.6240] | 0.1285 | 0.2241 | 0.096 |
| 0.95 | 0.1946 | [−0.715, −0.6072] | 0.1519 | 0.2446 | 0.093 |
| 0.99 | 0.2028 | [−0.715, −0.5973] | 0.1657 | 0.2527 | 0.087 |
| 0.999 | 0.2046 | [−0.715, −0.5805] | 0.1892 | 0.2546 | 0.065 |

- The unstable set is always contiguous, and nothing is unstable for η > 0 up to r = 12.
- At nz = 64 and 256 (0.5, 0.95 and 0.99 of D_max) the margins are 0.074–0.093.
- The coordinator's own onsets (0.013, 0.069, 0.091, 0.131) agree with these.
- **ORD = 3 and 4** (SPECTER's scheme u^(o) = u⁰ + dt/o·RHS with o = ORD..1): D_max = 0.2573 and 0.2852.
  - The onset in 1+σ depends only on D/D_max: it is 0.0323, 0.0679, 0.1091, 0.1394 and 0.1579 at 0.25, 0.5, 0.75, 0.9 and 0.97, identical for ORD 2, 3 and 4.
  - So sthr = 0.05 + D is conservative there too (margins 0.08–0.17).
- **Verdict: N = 2 is covered at every stable D, with a margin of at least 0.054 in 1+σ.** At the gate's D = 0.104 the guard fires at η ≤ −0.606 dz, against an onset of −0.6665 dz.

**N1 (new): N = 3 has an unguarded unstable band just below the explicit limit.** The N = 3 threshold is 0.05, which can never fire, because 1+σ ≥ 0.153 for every η. Yet the step becomes unstable near η ≈ −1.2 dz, which is where 1+σ is smallest:

| nz | first unstable at D/D_max | D | ρ at 0.995 D_max |
|---|---|---|---|
| 64 | 0.968 | 0.2035 | 1.058 per step |
| 128 | 0.985 | 0.2018 | 1.020 per step |
| 256 | 0.991 | 0.2016 | 1.008 per step |

- At 0.999 D_max (nz = 128) the band is r ∈ [−2.5, −0.75].
- Nothing is unstable at N = 3 for D ≤ 0.97 D_max (ORD 2–4) or for η > 0.
- Practical relevance is low. The band needs ν dt/dz² ≳ 0.20, i.e. dt ≥ 1.94e-3 at the gate's ν and nz, whereas the runs use 1e-3 or 2.5e-4.
- But REGISTRO P3's "N = 3 es estable en todos los casos" holds only up to about 0.97 D_max. A cheap fix is to also require ν dt/dz² ≤ 0.19 (or ≤ 0.95·D_max) whenever fsorder = 3.

## 1b. Quantities used by the guard

- `fs_nu = nu` and `fs_dt = dt` are set in `v_fs_init` (vboundary.f90:846–847).
- It is called at specter.fpp:981 with the global `nu` (the coefficient of `laplak` in `hd_rkstep2`) and the full, fixed `dt`. `dt = dt/mult` at :367 happens earlier, and dt is not changed inside the RK loop.
- So D is built from the full step, as in the replica, and is dimensionless.
- `z(2)−z(1)` is the physical spacing Lz/(nz−Cz−1). `grid`'s `z` is allocated and filled before `v_fs_init`.
- `USE grid, ONLY: z` shadows nothing:
  - the host module `boundary` only uses `fprecision`;
  - no used module exports another `z` (`kes` has `kz`; `var` has `pi` and `im`);
  - `planfc%z` is a derived-type component, not an entity;
  - there is no local `z`.

**Verdict: verified.** One caveat: the formula is calibrated on the ORD = 2 replica, but it is conservative for ORD 3 and 4 (1a).

## 1c. Physical regime and N = 3

- **N = 3:** sthr = 0.05 < min_η(1+σ) = 1 − β₂²/(2β₃) = 0.153, so it never fires.
- **N = 2, physical regime:** with η/dz ~ 1e-3, 1+σ ≥ 0.9986. That is far above sthr ≤ 0.05 + D_max ≈ 0.255; firing at η = 0 would need D ≥ 0.95, far beyond explicit stability.
- **N = 2 with η > 0, and fsgen = 1:** 1+σ ≥ 1 and sg = 0 respectively, so it never fires.

**Verdict: verified.**

## 2. Header of fsorder.f90 (lines 1–55)

- G = (F − 2σ∂x(W_new − W) + σDx)/(1+σ) now matches the code (lines 734–735).
- "Unstable where 1 + σ falls to ~0.7 ν dt/dz²": the measured onset/D ratio is 0.61–0.74 for D ≤ 0.18, rising to 0.80 at 0.98 D_max. That is fair.
- "η = −0.667 dz at 0.104" ✓; "−0.716 only as dt → 0" ✓; "stopped where 1 + σ ≤ 0.05 + ν dt/dz²" ✓ (N = 2).
- "N = 3 neutral for uniform η (weakly growing with slope), 1 + σ > 0.15" ✓.
- **Nit (inherited):** "extrapolating from interior nodes … unstable for |eta| > ~0.5 dz" is true at N = 2. At N = 3 with q = 8 it is stable down to −1.75 dz.
- N1 is new and naturally not in the header.

**Verdict: verified.**

## 3. The fsmaxit counter (fs_nmaxit / fs_resmaxit)

- Both are module variables initialized to 0 (boundary_mod.fpp:48 and :51).
- After normal completion of `DO it = 1,fs_maxit`, the DO variable is fs_maxit + 1 (Fortran: it keeps its last incremented value). After `EXIT` it is ≤ fs_maxit. So `it .gt. fs_maxit` counts exactly the non-converged substages; both EXIT paths are excluded, including convergence on the last pass. fs_maxit ≥ 2 is enforced.
- All the exit decisions use the MPI_ALLREDUCEd `res`, so the counters are identical on every rank. `fs_order_diagnostic` prints only on rank 0 and resets outside that IF, i.e. on every rank. All ranks call it, since it contains an MPI_REDUCE.
- It is called only under `IF (fs_general)` inside `planbc%bczend .eq. 2` (vboundary.f90), which is the same path as `fs_general_imposebc`. So the counter is meaningful.

**Residual gap (minor):** GLOBALOUTPUT runs at the start of a step when `timec == cstep` (specter.fpp:1139–1142), and there is no call after the last step. Non-converged substages in the final (< cstep steps) interval are therefore never reported. The same holds when `bench ≠ 0`.

**Verdict: verified**, with that gap.

## 4. Documentation

**Correct now:**
- +ε³/3 with numerator ε cosh ε − sinh ε ([H-25], `.tex` §3(iii)).
- The onset table: P3 and [H-24] agree with my scan to 1e-3 in r.
- "orden N en la amplitud … la velocidad hasta η^(N−1)" (REGISTRO §1, [D-51], `.tex` §1).
- C9 is now "origen cuantitativo no establecido", with the 1/nz credited (P4, [H-26], `.tex`).
- The N = 3 slope result is attributed to me and marked **no reproducido aparte** in REGISTRO §1 and P14, [H-24] and [V-22].
- The B1–B4 descriptions in P14 and [V-22] are accurate.
- `.tex` numbers and filled placeholders (orders 1.00/1.00 and 2.33/3.02, errors 6.9e-7 to 2.8e-5) match `salidas/tablas/fase6_escalera_no_lineal.json` and its log.
- Fase 2 6/6 (15:09 log) and Fase 5 8/8 (16:22 log): both were measured before the 16:34 edit. That edit only touches the general path, which these gates do not use.
- The NUL-byte note is supported. There is exactly one NUL right after each `prterun` abort banner, at the same place in both finished logs (bytes 3747/5237 and 3748/5239). The red log, which has no aborts, has none. A delayed-allocation hole would be a run of NULs, not isolated bytes followed by valid text.
- The PDF (16:55) is current with the `.tex`.

**Still wrong or unsupported:**

1. **The refuted threshold is still asserted.**
   - REGISTRO P2 table: "traza espectral implícita (final) | ≤ −0,70".
   - DECISIONES [H-23]: "estable para todo η > −0,7 Δz a N = 2".
   - At the replica's own D = 0.104, the band (−0.70, −0.6665] dz is unstable (ρ = 1.02–1.9 per step).
2. **Stale guard sentence.** REGISTRO §8: "A N = 2, la corrida aborta sola si η < −0,68 Δz". With the new guard it is 1+σ ≤ 0.05 + ν dt/Δz², i.e. η ≤ −0.61 Δz at the gate's D, as §3.6 correctly says.
3. **Gate result for the corrected code claimed before the run finished.**
   - REGISTRO §6 says "Tres corridas completas, todas con 4/10 … la definitiva … (`logs/fase6_2026-09-28.log`)", and [V-22] says "El código corregido se volvió a pasar por la puerta".
   - At 16:59 that log had only reached H5, with no RESULTADO line, and the gate process was still running. At 17:02 it had reached H6, still unfinished.
   - By my reading, none of H1–H10 can change outcome with the edits: H7, H8 and H10 abort earlier but abort anyway, and nothing else uses the guard. So 4/10 is the expected result, but it had not yet been measured.
   - [V-21] still cites `fase6_2026-09-28.log` for the 9 min 21 s run, which is now `fase6_segunda_4de10_2026-09-28.log`.
   - "Los dos logs tienen … NUL" now refers to three logs.
4. **The O(Δz²) claim is refuted (N2, below).** It appears in REGISTRO §1 item 4 and §5 ("N = 3 cae ×3,4, que es el O(Δz²) de ∂z³"), [H-26] "Precisión", and `.tex` §3(iv).
5. **"El piso escala como η²"** (P13, [H-26], `.tex` §3(iv)). This is expected from the mechanism (the m = 3 trace enters Dx times η²/2), but no measurement at a second η is cited. It should carry the § mark.
6. **Pass counts disagree.**
   - REGISTRO §1: "N = 3 … pasadas ≤ 7".
   - `.tex` box: "4–12 pasadas".
   - P13 says 8–16 in the depth staircase, and the nonlinear-staircase log shows ≤ 13 at ε = 0.04.
7. **`.tex` §3(ii) and §4:**
   - The N = 3 slope result is attributed ("según el verificador") but not marked as not independently reproduced, unlike the REGISTRO and DECISIONES.
   - "una lectura adversarial del código sin errores" is unqualified, although that read produced B1 and B4. [V-22]'s qualified wording ("sin errores: estado entre pasadas, …") is fine.
8. **Minor:**
   - "cinco cosas" or "cinco correcciones" next to six bullets (P14, [V-22] title).
   - "η > −0,67 Δz" (REGISTRO §1 and §2, [D-51], `.tex` box): the onset is −0.6665, so (−0.67, −0.6665] is unstable. Also the code now aborts from −0.606 at that D, so what runs is η > −0.61 Δz.
   - [H-25] gives the nz = 128 mesh floor as 2.2e-6, while REGISTRO §5 gives 1.2e-6. These are pre-existing and I did not check them.
   - Stale line references in REGISTRO §3: `:562` → 566, `:140` → 144, `:108` → 112, `:127` → 131, `:715-726` → about 718–736, `:684-706` → about 690–716. The `.tex` listings are correct.
   - P3 still says the replica "mide exactamente la tasa del continuo"; it agrees to 2.4–14% in rate.
   - The P3 onset table and the "≈ 0,7 D" relation have no script in the repo. `estabilidad_cierre.py` (14:44) runs at fixed dt, yet the guard comment cites it. The D = 0.021 row is not in my round-1 outputs either.
   - [D-51]'s title and first paragraph still say "orden N en η", and "no es un desarrollo en la amplitud del flujo" sits next to "orden N en la amplitud". The two can be reconciled, but read together they are confusing.

**Verdict:** not clean. Items 1–3 must be fixed; 4–8 are minor.

## N2 (new): the spectral ∂z³u trace does not converge beyond nz ≈ 128

**Method** (`r2_orden_trazas.py`, `r2_orden_poli.py`). I used the round-1 literal transliteration of the Fortran z-pipeline: values at the physical nodes → the continuation loop → forward FFT → Σ û·`fs_ph`. The same computation was repeated fully in extended precision (continuation, DFT by matrix, trace sum), for nz = 64…1024.

**The m = 3 trace does not converge.** The first row uses the gate's H2 profiles.

| profile | error | nz = 64 | 128 | 256 | 512 | 1024 |
|---|---|---|---|---|---|---|
| sin(1.2z) | relative | 2.0e-3 | **2.9e-4** | 4.0e-4 | 3.2e-3 | 2.4e-2 |
| 1 + 0.3z² (exact 0) | absolute | 1.3e-6 | 2.4e-5 | 2.6e-4 | 2.6e-3 | 2.2e-2 |
| z² (exact 0) | absolute | 1.2e-6 | 2.0e-5 | 2.4e-4 | 2.2e-3 | 1.7e-2 |
| exp(0.6z) | relative | 4.5e-4 | 1.5e-4 | 1.0e-3 | 9.1e-3 | — |
| cosh(1.1z) | relative | 1.4e-3 | 2.2e-4 | 2.4e-4 | 2.0e-3 | — |

- The sin(1.2z) value at nz = 128 is the same 2.9e-4 the gate measured in H2.
- **The extended-precision results are the same.** So this is not roundoff. It is a deterministic defect of the continuation with the stored tables (C = 25, d = 5), of about 2e-11·|f|, amplified by Δz⁻³.
- An oscillatory profile, 0.3cos(3z+1), does converge (orders 1.93, 2.00, …).
- The m = 2 trace also stalls near 1e-7 at nz ≥ 256 for the smooth profiles. That is negligible for N = 2, since it is multiplied by η.

**Consequences:**
- The truncation part is O(Δz²) only up to nz ≈ 128, the gate's resolution; beyond that the defect dominates and the error grows.
- N = 3's m = 3 error cannot be reduced by z-refinement with these tables.
- The staircase's ×3.4 improvement from nz = 128 to 256 is not explained by the m = 3 trace, whose error does not fall over that range for these profiles.

This strengthens the "∂z³u de borde más precisa" item in §10. The impact enters multiplied by η²/2, so it is irrelevant for the cell.

```claims
[{"status": "verified", "text": "1a For N=2 the new guard sthr = 0.05 + nu dt/dz^2 covers the whole unstable window at every D up to 0.999 D_max (D_max at eta=0: 0.2103/0.2048/0.2034 for nz=64/128/256, -> 2/pi^2). Margin in 1+sigma is 0.054-0.096 (0.085 at the gate's D, where the guard fires at -0.606 dz vs onset -0.6665 dz). The unstable set is contiguous, nothing is unstable for eta>0, and the onset depends only on D/D_max, so the guard is conservative for ORD=3 and 4 as well."},
 {"status": "verified", "text": "1b fs_nu=nu (the Laplacian coefficient) and fs_dt=dt (the full step, fixed before v_fs_init) are set in v_fs_init; z(2)-z(1) is the physical dz = Lz/(nz-Cz-1); D is dimensionless and built as in the replica; USE grid, ONLY: z shadows nothing (no other used module exports z, planfc%z is a component, no local z). The formula is ORD=2-calibrated but conservative for ORD=3 and 4."},
 {"status": "verified", "text": "1c At N=3, sthr=0.05 < min_eta(1+sigma)=0.153, so it never fires. At N=2 with eta/dz~1e-3, 1+sigma>=0.9986 >> sthr <= 0.255; firing at eta=0 would need D>=0.95. Nor does it fire for eta>0 or fsgen=1 (sigma=0)."},
 {"status": "verified", "text": "2 The fsorder.f90 header now has G=(F-2 sigma dx(W_new-W)+sigma Dx)/(1+sigma), consistent with lines 734-735; onset '~0.7 nu dt/dz^2' (measured ratio 0.61-0.74 for D<=0.18), -0.667 dz at D=0.104, -0.716 only as dt->0, and the stop at 1+sigma <= 0.05+D are all correct. Inherited nit: 'interior extrapolation unstable for |eta|>~0.5 dz' holds for N=2 only (N=3 with q=8 is stable down to -1.75 dz)."},
 {"status": "verified", "text": "3 fs_nmaxit/fs_resmaxit start at 0; after normal completion of DO it=1,fs_maxit, it = fs_maxit+1, so 'it .gt. fs_maxit' counts exactly the non-converged substages and no EXIT path is counted; decisions use the allreduced res, so the counters agree on all ranks; the reset runs on every rank; the diagnostic is called only on the general path. Residual gap: substages after the last cstep output (or with bench/=0) are never reported."},
 {"status": "refuted", "text": "4 The documentation is not clean. Correct: +eps^3/3, the onset table, 'orden N en la amplitud', C9 as 'origen cuantitativo no establecido', and the N=3 slope result attributed and marked unreproduced in REGISTRO and DECISIONES. Still wrong: the refuted threshold persists (REGISTRO P2 table '<= -0,70' and [H-23] 'estable para todo eta > -0,7 dz'); REGISTRO section 8 still says the run aborts at eta < -0,68 dz; section 6 and [V-22] report the corrected-code gate run as complete with 4/10, but at 16:59 its log had only reached H5 (4/10 is expected by reasoning, not yet measured), and [V-21] still cites the renamed log."},
 {"status": "refuted", "text": "N1 (new) At N=3 the guard (0.05) can never fire, yet the replica step is unstable near eta ~ -1.2 dz once D >= 0.968/0.985/0.991 D_max for nz=64/128/256 (D ~ 0.202), with rho up to 1.058/1.020/1.008 per step at 0.995 D_max. It is stable for D <= 0.97 D_max (ORD 2-4). Low practical relevance (needs dt >= 1.94e-3 at the gate's nu and nz), but 'N = 3 es estable en todos los casos' only holds below that; fix by capping nu dt/dz^2 <= ~0.19 at N=3."},
 {"status": "refuted", "text": "N2 (new) The claim that the FC-Gram error of the m=3 spectral trace is O(dz^2) ([H-26], REGISTRO sections 1 and 5, .tex 3(iv)) fails beyond nz~128. With the gate's own H2 profiles the error is: sin(1.2z) 2.9e-4 (nz=128, matching H2) -> 4.0e-4 -> 3.2e-3 -> 2.4e-2 (nz=1024); 1+0.3z^2 and z^2 (exact 0) grow about x10 per doubling from nz=64. Results are identical in extended precision, so this is a deterministic continuation defect (~2e-11|f|/dz^3), not roundoff. N=3 cannot be improved by z-refinement with these tables, and the staircase's x3.4 is not explained by the trace."},
 {"status": "unclear", "text": "N3 (minor docs) 'El piso escala como eta^2' (P13, [H-26], .tex) cites no second-eta measurement and should be marked as inferred. Pass counts disagree (REGISTRO section 1 '<=7', .tex '4-12', P13 '8-16', staircase log '<=13'). The .tex does not mark the N=3 slope result as unreproduced and calls the code read 'sin errores' although it produced B1/B4. Also: 'cinco' next to six bullets; '-0,67' vs onset -0.6665 and guard -0.606; [H-25] 2.2e-6 vs section 5 1.2e-6 mesh floor; stale line references in REGISTRO section 3; 'mide exactamente'; no repo script for the D-dependent onset table."}]
```
