#!/usr/bin/env python3
"""Chequeos reproducibles y valores numéricos para ``notes.tex``.

Se ejecuta en un solo proceso y mantiene matrices menores que 500 x 500.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("capa_local", HERE / "capa.py")
capa = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = capa
SPEC.loader.exec_module(capa)


def derived_values() -> dict[str, object]:
    """Única fuente de los decimales impresos en la nota."""

    free = capa.tasas_decaimiento("noslip-libre", 1.0, 1.0, 4)
    wall = capa.tasas_decaimiento("noslip-noslip", 1.0, 1.0, 4)
    h6 = 6.0e-3
    return {
        "tasas_libre_nu_sobre_h2": free.tolist(),
        "tasas_noslip_nu_sobre_h2": wall.tolist(),
        "alfa_libre_factor": capa.alfa(1.0, 1.0, "noslip-libre"),
        "alfa_noslip_factor": capa.alfa(1.0, 1.0, "noslip-noslip"),
        "cociente_alfas": capa.alfa(1.0, 1.0, "noslip-noslip")
        / capa.alfa(1.0, 1.0, "noslip-libre"),
        "superficie_sobre_promedio": capa.cociente_superficie_promedio(),
        "superficie_exceso_porcentaje": 100.0
        * (capa.cociente_superficie_promedio() - 1.0),
        "beta_promedio": float(np.pi**2 / 8.0),
        "beta_superficie": float(np.pi / 4.0),
        "gap_libre_factor": float(free[1] - free[0]),
        "gap_noslip_factor": float(wall[1] - wall[0]),
        "h6_alfa_sobre_nu_menos_m2": capa.alfa(
            1.0, h6, "noslip-libre"
        ),
        "h6_ejemplo_nu_1e_6_alfa_smenos1": capa.alfa(
            1.0e-6, h6, "noslip-libre"
        ),
        "h6_ejemplo_nu_1e_6_tau_s": 1.0
        / capa.alfa(1.0e-6, h6, "noslip-libre"),
    }


def _fd_rates(bordes: str, nz: int, n: int = 4) -> np.ndarray:
    """Referencia independiente: autovalores de una matriz FD densa."""

    dz = 1.0 / nz
    m = nz - 1 if bordes == "noslip-noslip" else nz
    A = np.diag(np.full(m, 2.0))
    if m > 1:
        A += np.diag(np.full(m - 1, -1.0), 1)
        A += np.diag(np.full(m - 1, -1.0), -1)
    if bordes == "noslip-libre":
        A[-1, -2] = -2.0
    vals = np.sort(np.linalg.eigvals(A).real / dz**2)
    return vals[:n]


def test_spectra_against_fd() -> str:
    details = []
    for bordes in ("noslip-noslip", "noslip-libre"):
        exact = capa.tasas_decaimiento(bordes, 1.0, 1.0, 4)
        errs = []
        for nz in (120, 240):
            fd = _fd_rates(bordes, nz)
            errs.append(float(np.max(np.abs(fd - exact) / exact)))
        order = np.log2(errs[0] / errs[1])
        assert 1.95 < order < 2.05, (bordes, errs, order)
        details.append(f"{bordes}: {errs[0]:.3e}->{errs[1]:.3e}, p={order:.3f}")
    return "; ".join(details)


def test_mode_integrals_and_reduction() -> str:
    z = np.linspace(0.0, 1.0, 200_001)
    f = 0.5 * np.pi * np.sin(0.5 * np.pi * z)
    mean = np.trapezoid(f, z)
    beta = np.trapezoid(f * f, z)
    surface_ratio = f[-1] / mean
    assert abs(mean - 1.0) < 1e-11
    assert abs(beta - np.pi**2 / 8.0) < 1e-11
    assert abs(surface_ratio - capa.cociente_superficie_promedio()) < 1e-11
    bottom_gradient = np.pi**2 / 4.0
    assert abs(bottom_gradient - capa.alfa(1.0, 1.0, "noslip-libre")) < 1e-14
    return f"<f>={mean:.12f}, <f^2>={beta:.12f}, f(1)/<f>={surface_ratio:.12f}"


def _measured_solver_error(bordes: str, nz: int) -> float:
    t_final = 0.3
    _, u = capa.solver_1d(1.0, 1.0, nz, t_final, bordes)
    measured = -np.log(np.max(np.abs(u))) / t_final
    return abs(measured - capa.alfa(1.0, 1.0, bordes))


def test_solver_second_order() -> str:
    details = []
    for bordes in ("noslip-noslip", "noslip-libre"):
        e48 = _measured_solver_error(bordes, 48)
        e96 = _measured_solver_error(bordes, 96)
        order = np.log2(e48 / e96)
        assert e48 > 0.0 and e96 > 0.0
        assert 1.95 < order < 2.05, (bordes, e48, e96, order)
        details.append(f"{bordes}: {e48:.3e}->{e96:.3e}, p={order:.3f}")
    return "; ".join(details)


def test_laplace_dn() -> str:
    cases = [
        (0.0, 1.0, 0.7, -0.3),
        (1.0e-10, 2.0, -0.4, 2.1),
        (1.0, 1.0, 1.0, 0.0),
        (3.7, 2.0, -0.4, 2.1),
        (50.0, 1.0, 1.0, 1.0),
        (350.0, 2.0, 0.5, -1.0),
        (2.5, 0.7, 1.0 + 2.0j, -0.3 + 0.8j),
    ]
    residual = 0.0
    for k, Lz, g0, g1 in cases:
        sol = capa.coef_laplace_dn(k, Lz, g0, g1)
        residual = max(residual, abs(sol.phi(0.0) - g0), abs(sol.dphi(Lz) - g1))
        if k == 0.0:
            assert sol.C1 == g1 and sol.C2 == g0
        z = np.linspace(0.0, Lz, 19)
        assert np.all(np.isfinite(sol.phi(z)))
        assert np.all(np.isfinite(sol.dphi(z)))
        if k > 0.0 and sol.C1 is not None and sol.C2 is not None:
            r = np.exp(-k * Lz)
            assert abs(r * sol.C1 + sol.C2 - g0) < 2e-5
            assert abs(k * (sol.C1 - r * sol.C2) - g1) < 2e-5
    assert residual < 2e-13, residual
    return f"7 casos (incluye kL=700 y k=1e-10), residuo máximo={residual:.3e}"


def test_failures() -> str:
    bad_calls = (
        lambda: capa.tasas_decaimiento("libre-noslip", 1.0, 1.0, 1),
        lambda: capa.tasas_decaimiento("noslip-libre", 0.0, 1.0, 1),
        lambda: capa.tasas_decaimiento("noslip-libre", 1.0, 1.0, 0),
        lambda: capa.solver_1d(1.0, 1.0, 1, 0.1, "noslip-libre"),
        lambda: capa.solver_1d(1.0, 1.0, 20, -0.1, "noslip-libre"),
        lambda: capa.coef_laplace_dn(-1.0, 1.0, 0.0, 0.0),
        lambda: capa.coef_laplace_dn(1.0, 0.0, 0.0, 0.0),
        lambda: capa.coef_laplace_dn(1.0, 1.0, [0.0], 0.0),
    )
    for call in bad_calls:
        try:
            call()
        except (TypeError, ValueError):
            continue
        raise AssertionError("una entrada inválida no produjo excepción")
    sol = capa.coef_laplace_dn(1.0, 1.0, 0.0, 0.0)
    try:
        sol.phi(1.1)
    except ValueError:
        pass
    else:
        raise AssertionError("se aceptó z fuera del intervalo")
    return "9 rutas inválidas rechazadas"


TESTS = (
    ("espectros vs FD independiente", test_spectra_against_fd),
    ("integrales modal y reducción", test_mode_integrals_and_reduction),
    ("solver FD de orden dos", test_solver_second_order),
    ("Laplace DN y estabilidad", test_laplace_dn),
    ("validación de entradas", test_failures),
)


def main() -> int:
    if len(sys.argv) == 2 and sys.argv[1] == "--values":
        print(json.dumps(derived_values(), indent=2, sort_keys=True))
        return 0
    if len(sys.argv) != 1:
        print("uso: python out/checks.py [--values]", file=sys.stderr)
        return 2

    print("CHEQUEOS PROPIOS — CAPA DELGADA")
    passed = 0
    for name, test in TESTS:
        try:
            detail = test()
        except Exception as exc:
            print(f"  [FALLA] {name}: {type(exc).__name__}: {exc}")
        else:
            print(f"  [PASA]  {name}: {detail}")
            passed += 1
    print(f"RESULTADO: {passed}/{len(TESTS)}")
    return 0 if passed == len(TESTS) else 1


if __name__ == "__main__":
    raise SystemExit(main())
