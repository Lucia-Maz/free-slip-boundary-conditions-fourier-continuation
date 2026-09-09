#!/usr/bin/env python
"""
Puerta de aceptación de la Fase 1 — estructura vertical de una capa delgada.

Este archivo lo escribe la humana y el equipo de agentes NO puede editarlo: queda
fijado por hash con --acceptance-guard. Mientras alguno de estos chequeos falle, el
equipo no puede declararse terminado.

PRINCIPIO DE DISEÑO: toda verdad de referencia se CALCULA acá adentro, a partir de
una discretización en diferencias finitas construida en este mismo archivo. No hay
ninguna fórmula analítica escrita a mano contra la cual comparar. Eso importa: si la
fórmula que motivó el proyecto (alfa = nu*pi^2/(4h^2)) estuviera equivocada, este test
lo detectaría en vez de fijar el error. Los chequeos C6 se verifican por residuo:
el test no necesita conocer la solución para comprobar que se cumplen los bordes.

Qué tiene que existir cuando esto pase: un módulo `out/capa.py` que exponga

    tasas_decaimiento(bordes, nu, h, n)      -> ndarray de n tasas, ascendente
    cociente_superficie_promedio()           -> float
    alfa(nu, h, bordes)                      -> float, la tasa del modo más lento
    solver_1d(nu, h, nz, t_final, bordes)    -> (z, u) con u evolucionado hasta t_final
                                                desde el modo más lento de amplitud 1
    coef_laplace_dn(k, Lz, g0, g1)           -> objeto con .phi(z) y .dphi(z)

con `bordes` una de las dos cadenas exactas: "noslip-noslip" o "noslip-libre".

Se corre con:
    /home/lucia/miniforge3/envs/piv-dt/bin/python verificacion/test_aceptacion_fase1.py
"""

import os
import pathlib
import importlib.util
import sys

import numpy as np

TOL_TASAS = 1e-4        # las tasas de FD tienen error O(dz^2); con nz=400 sobra
TOL_COCIENTE = 1e-3
TOL_RESIDUO = 1e-10


# ---------------------------------------------------------------------------
# Verdad de referencia, construida acá: el operador nu * d2/dz2 en [0, h]
# ---------------------------------------------------------------------------

def operador_difusion(bordes, nu, h, nz):
    """Matriz densa de -nu*d2/dz2 con los bordes pedidos, en una grilla uniforme.

    Fondo (z=0) siempre no-deslizante: u = 0, nodo excluido de las incógnitas.
    Tope (z=h):
      - "noslip-noslip": u = 0, nodo excluido.
      - "noslip-libre" : du/dz = 0, nodo incluido, cerrado con nodo fantasma
                          u[nz+1] = u[nz-1], que es la derivada centrada nula.
    Devuelve (A, z) con z las coordenadas de las incógnitas.
    """
    dz = h / nz
    if bordes == "noslip-noslip":
        m = nz - 1
        z = np.arange(1, nz) * dz
    elif bordes == "noslip-libre":
        m = nz
        z = np.arange(1, nz + 1) * dz
    else:
        raise ValueError(f"bordes desconocido: {bordes!r}")

    A = np.zeros((m, m))
    for i in range(m):
        A[i, i] = 2.0
        if i > 0:
            A[i, i - 1] = -1.0
        if i < m - 1:
            A[i, i + 1] = -1.0
    if bordes == "noslip-libre":
        # última fila: (u[nz-1] - 2 u[nz] + u[nz+1])/dz^2 con u[nz+1] = u[nz-1]
        A[m - 1, m - 2] = -2.0
    return nu * A / dz**2, z


_CACHE_AUTOVALORES = {}


def _autovalores(bordes, nu, h, nz):
    clave = (bordes, nu, h, nz)
    if clave not in _CACHE_AUTOVALORES:
        A, _ = operador_difusion(bordes, nu, h, nz)
        _CACHE_AUTOVALORES[clave] = np.sort(np.linalg.eigvals(A).real)
    return _CACHE_AUTOVALORES[clave]


def tasas_referencia(bordes, nu, h, n, nz=400):
    """Tasas de referencia con extrapolación de Richardson.

    El operador en diferencias finitas converge a orden 2 (medido: 2.000 para los
    dos pares de bordes), así que (4*lam(2nz) - lam(nz))/3 cancela el término en
    dz^2 y deja O(dz^4). Sin esto, la referencia con nz=400 tiene error relativo
    (n*pi)^2*dz^2/12, que en el cuarto modo vale 8.2e-5 y deja la tolerancia de
    1e-4 con un factor 1.2 de margen: el chequeo pasaría siempre, pero por poco y
    sin margen si alguien pide más modos. Con Richardson el error baja a ~7e-10.
    """
    a = _autovalores(bordes, nu, h, nz)[:n]
    b = _autovalores(bordes, nu, h, 2 * nz)[:n]
    return (4 * b - a) / 3


def modo_lento_referencia(bordes, nu, h, nz=400):
    A, z = operador_difusion(bordes, nu, h, nz)
    lam, vec = np.linalg.eig(A)
    i = int(np.argmin(lam.real))
    u = vec[:, i].real
    if u[np.argmax(np.abs(u))] < 0:
        u = -u
    return z, u / np.max(np.abs(u))


# ---------------------------------------------------------------------------
# Carga del módulo que el equipo tiene que escribir
# ---------------------------------------------------------------------------

def cargar_capa():
    cand = []
    env = os.environ.get("CAPA")
    if env:
        cand.append(pathlib.Path(env))
    cand.append(pathlib.Path.cwd() / "out" / "capa.py")
    raiz = pathlib.Path(__file__).resolve().parents[1]
    cand.append(raiz / "out" / "capa.py")
    cand.extend(sorted(raiz.glob("jobs/*/out/capa.py")))
    for p in cand:
        if p.is_file():
            spec = importlib.util.spec_from_file_location("capa", p)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            return mod, p
    raise FileNotFoundError(
        "no se encontró out/capa.py; se buscó en: " + ", ".join(str(p) for p in cand)
    )


# ---------------------------------------------------------------------------
# Los chequeos
# ---------------------------------------------------------------------------

def c1_tasas_noslip(capa):
    """Las cuatro tasas más lentas con no-deslizamiento en las dos paredes."""
    nu, h, n = 1.0, 1.0, 4
    ref = tasas_referencia("noslip-noslip", nu, h, n)
    got = np.asarray(capa.tasas_decaimiento("noslip-noslip", nu, h, n), dtype=float)
    if got.shape != ref.shape:
        return False, f"forma {got.shape}, se esperaba {ref.shape}"
    err = np.max(np.abs(got - ref) / ref)
    return err < TOL_TASAS, f"error relativo máximo {err:.2e} (tol {TOL_TASAS:.0e})"


def c2_tasas_libre(capa):
    """Las cuatro tasas más lentas con fondo no-deslizante y tope libre."""
    nu, h, n = 1.0, 1.0, 4
    ref = tasas_referencia("noslip-libre", nu, h, n)
    got = np.asarray(capa.tasas_decaimiento("noslip-libre", nu, h, n), dtype=float)
    if got.shape != ref.shape:
        return False, f"forma {got.shape}, se esperaba {ref.shape}"
    err = np.max(np.abs(got - ref) / ref)
    return err < TOL_TASAS, f"error relativo máximo {err:.2e} (tol {TOL_TASAS:.0e})"


def c3_factor_cuatro(capa):
    """El cociente entre las dos fricciones, que es el resultado que da título al trabajo.

    No se compara contra un 4 escrito a mano: se exige que el cociente que sale del
    módulo coincida con el que sale de la referencia numérica de este archivo.
    """
    nu, h = 1.0, 1.0
    ref = tasas_referencia("noslip-noslip", nu, h, 1)[0] / \
          tasas_referencia("noslip-libre", nu, h, 1)[0]
    got = capa.alfa(nu, h, "noslip-noslip") / capa.alfa(nu, h, "noslip-libre")
    err = abs(got - ref) / ref
    return err < 1e-3, f"cociente {got:.6f}, referencia numérica {ref:.6f}"


def c4_cociente_superficie(capa):
    """u en la superficie sobre u promediada en la vertical, para el modo más lento libre.

    Es el sesgo geométrico del PIV de superficie: las partículas flotan arriba.
    """
    z, u = modo_lento_referencia("noslip-libre", 1.0, 1.0)
    z_full = np.concatenate(([0.0], z))
    u_full = np.concatenate(([0.0], u))
    promedio = np.trapezoid(u_full, z_full) / (z_full[-1] - z_full[0])
    ref = u_full[-1] / promedio
    got = float(capa.cociente_superficie_promedio())
    err = abs(got - ref) / ref
    return err < TOL_COCIENTE, f"cociente {got:.6f}, referencia numérica {ref:.6f}"


def c5_solver_converge(capa):
    """El solver 1D tiene que ser un solver, y tiene que converger a orden 2.

    Se mide la tasa de decaimiento que produce el solver a dos resoluciones y se
    exige que el error contra la tasa analítica caiga como dz^2. Una función que
    devolviera la respuesta analítica sin integrar nada daría error nulo y cociente
    indefinido: eso también se rechaza.
    """
    nu, h, t_final = 1.0, 1.0, 0.3
    exacta = capa.alfa(nu, h, "noslip-libre")
    errores = []
    for nz in (48, 96):
        z, u = capa.solver_1d(nu, h, nz, t_final, "noslip-libre")
        z = np.asarray(z, dtype=float)
        u = np.asarray(u, dtype=float)
        if u.shape != z.shape:
            return False, f"nz={nz}: z y u tienen formas distintas"
        if not np.all(np.isfinite(u)):
            return False, f"nz={nz}: el campo tiene valores no finitos"
        amp = np.max(np.abs(u))
        if not (0.0 < amp < 1.0):
            return False, (f"nz={nz}: amplitud final {amp:.4g}; se esperaba una "
                           "decaída desde amplitud 1")
        medida = -np.log(amp) / t_final
        errores.append(abs(medida - exacta))
    if errores[1] <= 0 or errores[0] <= 0:
        return False, "error exactamente nulo: el solver no está integrando"
    orden = np.log2(errores[0] / errores[1])
    ok = 1.7 < orden < 2.3
    return ok, f"orden medido {orden:.2f} (errores {errores[0]:.2e} -> {errores[1]:.2e})"


def c6_laplace_dirichlet_neumann(capa):
    """La rama que le falta a laplace_z de SPECTER: Dirichlet abajo, Neumann arriba.

    Se verifica por RESIDUO: phi(0) = g0 y dphi(Lz) = g1. El test no necesita saber
    la fórmula. Incluye k = 0 (el modo constante en x,y) y k*Lz grande, donde una
    implementación ingenua con exponenciales crecientes desborda.
    """
    casos = [
        (0.0, 1.0, 0.7, -0.3),
        (1.0, 1.0, 1.0, 0.0),
        (3.7, 2.0, -0.4, 2.1),
        (50.0, 1.0, 1.0, 1.0),
        (350.0, 2.0, 0.5, -1.0),   # k*Lz = 700: exp(700) desborda en doble precisión
    ]
    for k, Lz, g0, g1 in casos:
        sol = capa.coef_laplace_dn(k, Lz, g0, g1)
        v0, d1 = float(sol.phi(0.0)), float(sol.dphi(Lz))
        if not (np.isfinite(v0) and np.isfinite(d1)):
            return False, f"k={k}, Lz={Lz}: valores no finitos (phi(0)={v0}, dphi(Lz)={d1})"
        e0 = abs(v0 - g0) / max(1.0, abs(g0))
        e1 = abs(d1 - g1) / max(1.0, abs(g1))
        if max(e0, e1) > TOL_RESIDUO:
            return False, (f"k={k}, Lz={Lz}: residuo phi(0)={e0:.2e}, "
                           f"dphi(Lz)={e1:.2e} (tol {TOL_RESIDUO:.0e})")
    return True, f"{len(casos)} casos, residuos por debajo de {TOL_RESIDUO:.0e}"


CHEQUEOS = [
    ("C1  tasas, no-slip / no-slip", c1_tasas_noslip),
    ("C2  tasas, no-slip / libre", c2_tasas_libre),
    ("C3  cociente entre las dos fricciones", c3_factor_cuatro),
    ("C4  superficie sobre promedio vertical", c4_cociente_superficie),
    ("C5  el solver 1D converge a orden 2", c5_solver_converge),
    ("C6  Laplace Dirichlet-Neumann, por residuo", c6_laplace_dirichlet_neumann),
]


def main():
    try:
        capa, ruta = cargar_capa()
    except FileNotFoundError as e:
        print("PUERTA DE ACEPTACIÓN — FASE 1")
        print(f"\n  no se puede evaluar: {e}")
        print("\nRESULTADO: 0/%d — falta el módulo." % len(CHEQUEOS))
        return 1

    print("PUERTA DE ACEPTACIÓN — FASE 1")
    print(f"módulo: {ruta}\n")
    pasan = 0
    for nombre, fn in CHEQUEOS:
        try:
            ok, detalle = fn(capa)
        except Exception as e:
            ok, detalle = False, f"{type(e).__name__}: {e}"
        print(f"  [{'PASA' if ok else 'FALLA'}]  {nombre}\n            {detalle}")
        pasan += bool(ok)
    print(f"\nRESULTADO: {pasan}/{len(CHEQUEOS)}")
    return 0 if pasan == len(CHEQUEOS) else 1


if __name__ == "__main__":
    sys.exit(main())
