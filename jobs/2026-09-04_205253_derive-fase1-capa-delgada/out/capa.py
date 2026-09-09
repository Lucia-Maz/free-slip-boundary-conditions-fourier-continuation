"""Estructura vertical de una capa viscosa y extensión Laplace DN.

Todas las derivadas en ``coef_laplace_dn`` están orientadas en la dirección
positiva de z, igual que la variable ``b`` construida por ``laplace_z`` en
SPECTER.  El módulo no modifica ni depende del código Fortran de SPECTER.
"""

from __future__ import annotations

import operator

import numpy as np
from scipy.sparse import diags, lil_matrix
from scipy.sparse.linalg import expm_multiply


_BORDES = ("noslip-noslip", "noslip-libre")


def _bordes_validos(bordes: str) -> str:
    if bordes not in _BORDES:
        raise ValueError(
            f"bordes desconocido: {bordes!r}; opciones válidas: {_BORDES}"
        )
    return bordes


def _positivo_finito(nombre: str, valor: float) -> float:
    try:
        valor = float(valor)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{nombre} debe ser un escalar real") from exc
    if not np.isfinite(valor) or valor <= 0.0:
        raise ValueError(f"{nombre} debe ser positivo y finito")
    return valor


def _entero_minimo(nombre: str, valor: int, minimo: int) -> int:
    if isinstance(valor, (bool, np.bool_)):
        raise TypeError(f"{nombre} debe ser entero, no booleano")
    try:
        entero = operator.index(valor)
    except TypeError as exc:
        raise TypeError(f"{nombre} debe ser entero") from exc
    if entero < minimo:
        raise ValueError(f"{nombre} debe ser >= {minimo}")
    return entero


def tasas_decaimiento(bordes: str, nu: float, h: float, n: int) -> np.ndarray:
    """Devuelve las ``n`` menores tasas positivas, en orden ascendente.

    Para no-slip/no-slip los números de onda son ``m*pi/h``, ``m=1,...``;
    para no-slip/libre son ``(m+1/2)*pi/h``, ``m=0,...``.
    """

    bordes = _bordes_validos(bordes)
    nu = _positivo_finito("nu", nu)
    h = _positivo_finito("h", h)
    n = _entero_minimo("n", n, 1)

    if bordes == "noslip-noslip":
        indices = np.arange(1, n + 1, dtype=float)
    else:
        indices = np.arange(n, dtype=float) + 0.5
    return nu * (np.pi * indices / h) ** 2


def cociente_superficie_promedio() -> float:
    """Razón u(h)/<u> del modo fundamental no-slip/libre."""

    return float(np.pi / 2.0)


def alfa(nu: float, h: float, bordes: str) -> float:
    """Coeficiente de fricción lineal: tasa del modo vertical más lento."""

    return float(tasas_decaimiento(bordes, nu, h, 1)[0])


def _operador_fd(bordes: str, nu: float, h: float, nz: int):
    """Construye -nu*d2/dz2 en los nodos incógnita (orden espacial 2)."""

    dz = h / nz
    if bordes == "noslip-noslip":
        m = nz - 1
    else:
        m = nz

    A = diags(
        diagonals=(
            -np.ones(max(m - 1, 0)),
            2.0 * np.ones(m),
            -np.ones(max(m - 1, 0)),
        ),
        offsets=(-1, 0, 1),
        shape=(m, m),
        format="lil",
        dtype=float,
    )
    if bordes == "noslip-libre":
        # u_{nz+1}=u_{nz-1}: derivada centrada nula en el tope.
        A[m - 1, m - 2] = -2.0
    return (nu / dz**2) * A.tocsr()


def solver_1d(
    nu: float, h: float, nz: int, t_final: float, bordes: str
) -> tuple[np.ndarray, np.ndarray]:
    """Integra la difusión 1D desde el modo lento de amplitud unitaria.

    ``nz`` es el número de intervalos uniformes. Se discretiza ``-nu*d2/dz2``
    con orden dos y se aplica la exponencial del operador semidiscreto; por
    tanto, el error que queda es espacial, no una fórmula analítica insertada.
    Se devuelven ambos extremos del intervalo.
    """

    bordes = _bordes_validos(bordes)
    nu = _positivo_finito("nu", nu)
    h = _positivo_finito("h", h)
    nz = _entero_minimo("nz", nz, 2)
    try:
        t_final = float(t_final)
    except (TypeError, ValueError) as exc:
        raise TypeError("t_final debe ser un escalar real") from exc
    if not np.isfinite(t_final) or t_final < 0.0:
        raise ValueError("t_final debe ser no negativo y finito")

    z = np.linspace(0.0, h, nz + 1)
    if bordes == "noslip-noslip":
        z_int = z[1:-1]
        u0 = np.sin(np.pi * z_int / h)
    else:
        z_int = z[1:]
        u0 = np.sin(0.5 * np.pi * z_int / h)
    u0 /= np.max(np.abs(u0))

    A = _operador_fd(bordes, nu, h, nz)
    # traceA evita la estimación estocástica de traza de expm_multiply.
    u_int = expm_multiply(-t_final * A, u0, traceA=-t_final * A.diagonal().sum())

    u = np.empty_like(z)
    u[0] = 0.0
    if bordes == "noslip-noslip":
        u[1:-1] = u_int
        u[-1] = 0.0
    else:
        u[1:] = u_int
    return z, u


def _escalar_finito(nombre: str, valor):
    arr = np.asarray(valor)
    if arr.ndim != 0:
        raise TypeError(f"{nombre} debe ser escalar")
    valor = arr.item()
    if not np.isfinite(valor):
        raise ValueError(f"{nombre} debe ser finito")
    return valor


class SolucionLaplaceDN:
    """Un modo de Laplace con Dirichlet abajo y derivada +z arriba."""

    __slots__ = ("k", "Lz", "g0", "g1", "C1", "C2")

    def __init__(self, k, Lz, g0, g1, C1, C2):
        # Clase explícita: la puerta carga el módulo sin registrarlo en
        # sys.modules, cosa que dataclasses no admite con anotaciones diferidas.
        self.k = k
        self.Lz = Lz
        self.g0 = g0
        self.g1 = g1
        self.C1 = C1
        self.C2 = C2

    @staticmethod
    def _salida(valor, era_escalar: bool):
        return np.asarray(valor).item() if era_escalar else valor

    def _validar_z(self, z):
        arr = np.asarray(z, dtype=float)
        if not np.all(np.isfinite(arr)):
            raise ValueError("z debe contener valores finitos")
        # Las exponenciales están normalizadas para el dominio físico.
        tol = 16.0 * np.finfo(float).eps * max(1.0, self.Lz)
        if np.any(arr < -tol) or np.any(arr > self.Lz + tol):
            raise ValueError("z debe pertenecer a [0, Lz]")
        return np.clip(arr, 0.0, self.Lz), arr.ndim == 0

    def phi(self, z):
        """Evalúa phi sin exponenciales crecientes en 0 <= z <= Lz."""

        z, era_escalar = self._validar_z(z)
        if self.k == 0.0:
            ans = self.g0 + self.g1 * z
        elif self.k * self.Lz < 0.5:
            # Forma regular cuando C1 y C2 se vuelven grandes por k -> 0.
            kz = self.k * z
            K = self.k * self.Lz
            ans = self.g0 * (np.cosh(kz) - np.tanh(K) * np.sinh(kz))
            ans += self.g1 * z * np.sinc(1j * kz / np.pi).real / np.cosh(K)
        else:
            ans = self.C1 * np.exp(self.k * (z - self.Lz))
            ans += self.C2 * np.exp(-self.k * z)
        return self._salida(ans, era_escalar)

    def dphi(self, z):
        """Evalúa dphi/dz en la dirección positiva de z."""

        z, era_escalar = self._validar_z(z)
        if self.k == 0.0:
            ans = self.g1 + np.zeros_like(z)
        elif self.k * self.Lz < 0.5:
            kz = self.k * z
            K = self.k * self.Lz
            ans = self.k * self.g0 * (
                np.sinh(kz) - np.tanh(K) * np.cosh(kz)
            )
            ans += self.g1 * np.cosh(kz) / np.cosh(K)
        else:
            ans = self.k * (
                self.C1 * np.exp(self.k * (z - self.Lz))
                - self.C2 * np.exp(-self.k * z)
            )
        return self._salida(ans, era_escalar)


def coef_laplace_dn(k: float, Lz: float, g0, g1) -> SolucionLaplaceDN:
    """Coeficientes estables para phi(0)=g0, phi'(Lz)=g1.

    Para k>0 se usa

        phi = C1*exp(k*(z-Lz)) + C2*exp(-k*z),

    con un denominador ``1 + exp(-2*k*Lz)``. Para k=0 la solución es
    lineal: ``C1=g1`` es la pendiente y ``C2=g0`` el intercepto, siguiendo
    la rama de modo cero de SPECTER.
    """

    try:
        k = float(k)
    except (TypeError, ValueError) as exc:
        raise TypeError("k debe ser un escalar real") from exc
    if not np.isfinite(k) or k < 0.0:
        raise ValueError("k debe ser no negativo y finito")
    Lz = _positivo_finito("Lz", Lz)
    g0 = _escalar_finito("g0", g0)
    g1 = _escalar_finito("g1", g1)

    if k == 0.0:
        # En la rama lineal de SPECTER, coef(1) es la pendiente y coef(2)
        # el valor inferior.
        return SolucionLaplaceDN(k, Lz, g0, g1, g1, g0)

    r = np.exp(-k * Lz)
    den = 1.0 + r * r
    C1 = (g1 / k + r * g0) / den
    C2 = (g0 - r * g1 / k) / den
    return SolucionLaplaceDN(k, Lz, g0, g1, C1, C2)


__all__ = [
    "tasas_decaimiento",
    "cociente_superficie_promedio",
    "alfa",
    "solver_1d",
    "coef_laplace_dn",
]
