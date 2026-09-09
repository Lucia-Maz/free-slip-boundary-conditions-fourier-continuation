#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Parámetros de la celda y de la campaña, en un único lugar.

Regla de diseño, y es una decisión de Lucía: el espesor de la capa, la viscosidad y la
geometría del forzado **cambian entre experiencias**. No son propiedades del problema,
son condiciones de una corrida. Por eso viven acá y en ningún otro lado: ningún script,
derivación ni tabla del informe debe escribir 6e-3 en el medio de una cuenta.

Cada campo lleva su procedencia. Un valor sin procedencia es un defecto, no un valor.

Uso:
    from celda import CAMPANA_02_06_25 as celda
    alfa = celda.alfa()
"""

from dataclasses import dataclass, field, asdict
import math


@dataclass(frozen=True)
class Fluido:
    """Propiedades del fluido de trabajo."""

    nu: float                     # viscosidad cinemática [m^2/s]
    descripcion: str
    procedencia_nu: str

    def __post_init__(self):
        if self.nu <= 0:
            raise ValueError("nu tiene que ser positiva")


@dataclass(frozen=True)
class Forzado:
    """Geometría del arreglo de imanes que produce la fuerza de Lorentz.

    El modelo es: imanes cilíndricos de diámetro `diametro` con sus ejes verticales,
    dispuestos en una red cuadrada de paso `paso` con la polaridad alternada como un
    tablero de ajedrez, y una corriente horizontal aproximadamente uniforme que
    atraviesa la celda. La fuerza por unidad de volumen es J x B, y con B casi vertical
    bajo cada imán el patrón de fuerza reproduce el patrón de polaridades.
    """

    paso: float                   # separación entre imanes vecinos [m]
    diametro: float               # diámetro de cada imán [m]
    patron: str                   # "tablero" (polaridad alternada en red cuadrada)
    procedencia: str

    def k_fundamental(self) -> float:
        """Número de onda del armónico fundamental del patrón de polaridades [1/m].

        Para el signo (-1)^(i+j) sobre una red cuadrada de paso a, el primer armónico
        no nulo está en k = (pi/a)(±1, ±1): es DIAGONAL, y su longitud de onda es
        a*sqrt(2), no a. Este valor es una propiedad de la red, no del imán; la forma
        del imán sólo pesa las amplitudes. `numerico/fase4/p01_validez.py` lo verifica
        numéricamente contra la transformada del patrón real, con imanes de diámetro
        finito, en vez de darlo por sentado.
        """
        if self.patron != "tablero":
            raise NotImplementedError("sólo está modelado el patrón de tablero")
        return math.pi * math.sqrt(2.0) / self.paso


@dataclass(frozen=True)
class Camara:
    """Lo que fija la resolución espacial y temporal de la medición."""

    fps: float
    escala_px_por_m: float
    lado_px: int
    cuadros: int
    procedencia: str

    def duracion(self) -> float:
        return self.cuadros / self.fps

    def campo_de_vision(self) -> float:
        return self.lado_px / self.escala_px_por_m


@dataclass(frozen=True)
class Celda:
    """Una campaña de medición: la capa, el fluido, el forzado y la cámara."""

    nombre: str
    h: float                      # espesor de la capa [m]
    procedencia_h: str
    fluido: Fluido
    forzado: Forzado
    camara: Camara
    notas: tuple = field(default_factory=tuple)

    # ---- cantidades derivadas -------------------------------------------------

    def alfa(self, bordes: str = "noslip-libre") -> float:
        """Coeficiente de fricción de fondo de la clausura de un modo [1/s].

        Fase 1, verificada numéricamente: con fondo no-deslizante y tope libre de
        tensiones alfa = pi^2 nu / (4 h^2); con no-deslizamiento en las dos caras,
        pi^2 nu / h^2. El cociente es exactamente 4.
        """
        base = math.pi ** 2 * self.fluido.nu / self.h ** 2
        if bordes == "noslip-libre":
            return base / 4.0
        if bordes == "noslip-noslip":
            return base
        raise ValueError("bordes debe ser 'noslip-libre' o 'noslip-noslip'")

    def tau_alfa(self, bordes: str = "noslip-libre") -> float:
        """Tiempo de fricción 1/alfa [s]."""
        return 1.0 / self.alfa(bordes)

    def tau_difusion(self) -> float:
        """Tiempo de difusión viscosa a través del espesor, h^2/nu [s]."""
        return self.h ** 2 / self.fluido.nu

    def t_olvido_modal(self) -> float:
        """1/(lambda_1 - lambda_0) con tope libre [s].

        Cuánto hay que esperar para que el segundo modo vertical pierda memoria y el
        perfil sea el fundamental. lambda_m = nu ((m+1/2) pi / h)^2, así que la brecha
        es 2 pi^2 nu / h^2.
        """
        return self.h ** 2 / (2.0 * math.pi ** 2 * self.fluido.nu)

    def sesgo_piv_superficie(self) -> float:
        """u(h)/<u> para el modo fundamental con tope libre: pi/2.

        Las partículas trazadoras flotan, así que el PIV mide la velocidad de la
        superficie y no el promedio vertical. Para pasar de lo medido al promedio hay
        que dividir por este factor.
        """
        return math.pi / 2.0

    def a_dict(self) -> dict:
        return asdict(self)


# --------------------------------------------------------------------------------
# La campaña del 2 de junio de 2025
# --------------------------------------------------------------------------------

AGUA = Fluido(
    nu=1.0e-6,
    descripcion="agua a ~20 C",
    procedencia_nu=(
        "SUPUESTO NO VERIFICADO, decidido explícitamente por Lucía el 2026-09-09: se "
        "usa la viscosidad del agua y no la del electrolito real (KNO3 al 16 % m/m). "
        "Entra a primer orden en alfa, que va como nu/h^2, así que el error que "
        "introduce es del mismo orden que su desviación relativa respecto del agua. "
        "No se midió ni se buscó en tablas."
    ),
)

CAMPANA_02_06_25 = Celda(
    nombre="02-06-25",
    h=6.0e-3,
    procedencia_h=(
        "Lucía, 2026-09-09, zanjando la discrepancia entre los 5 mm de README.md y los "
        "~6 mm del manuscrito de sinh-Poisson: para esta campaña es 6 mm. Es un valor "
        "que cambia entre experiencias, por eso está acá y no dentro de ninguna cuenta."
    ),
    fluido=AGUA,
    forzado=Forzado(
        paso=5.0e-2,
        diametro=1.0e-2,
        patron="tablero",
        procedencia=(
            "Lucía, 2026-09-09: 'los imanes tienen un diámetro de 1 cm, están en una "
            "grid tipo checkerboard separados por 5 cm'. Corrige el valor heredado de "
            "l = 2 cm que se venía usando en [P-01] como escala de inyección."
        ),
    ),
    camara=Camara(
        fps=60.0,
        escala_px_por_m=3800.0,
        lado_px=1024,
        cuadros=3072,
        procedencia=(
            "Cabeceras .cih de la campaña y la grilla de calibración med_S0001; "
            "resumido en README.md y verificado en codigo/01_inventario.py."
        ),
    ),
    notas=(
        "Electrolito: KNO3 al 16 % m/m. Partículas trazadoras de 100 um, que flotan.",
        "Corrientes de forzado: 0,9-1,0 A (16 V), 1,9-2,0 A y 2,1-2,2 A (32 V).",
    ),
)


# --------------------------------------------------------------------------------
# Velocidades medidas, con su procedencia y su estado de verificación
# --------------------------------------------------------------------------------

VELOCIDADES_MEDIDAS = {
    "forzado_2.15A": {
        "u_rms_superficie_m_s": 1.3214710285395503e-3,
        "ic68_m_s": (1.3060363837019684e-3, 1.3456906789877758e-3),
        "que_es": (
            "Ordenada al origen del barrido en Dt sobre med_S0008, es decir la "
            "velocidad cuadrática media del campo sobre los vectores válidos, "
            "extrapolada a ruido nulo."
        ),
        "procedencia": "salidas/tablas/ajuste_med_S0008.json, ver [H-05]",
        "verificado": True,
    },
    "inicio_decaimiento": {
        "u_rms_superficie_m_s": 4.0e-3,
        "ic68_m_s": None,
        "que_es": "Velocidad al inicio de un decaimiento.",
        "procedencia": (
            "Heredado: figura en ESTADO.md sin script que lo produzca. NO se pudo "
            "rastrear a una salida de este árbol."
        ),
        "verificado": False,
    },
}

DECAIMIENTOS = ("med_S0003", "med_S0005", "med_S0006", "med_S0009")


if __name__ == "__main__":
    c = CAMPANA_02_06_25
    print("celda %s: h = %.1f mm, nu = %.2e m^2/s" % (c.nombre, c.h * 1e3, c.fluido.nu))
    print("  alfa (tope libre)   = %.6f 1/s   tau = %.2f s" %
          (c.alfa(), c.tau_alfa()))
    print("  alfa (canal)        = %.6f 1/s" % c.alfa("noslip-noslip"))
    print("  h^2/nu              = %.2f s" % c.tau_difusion())
    print("  olvido modal        = %.2f s" % c.t_olvido_modal())
    print("  k del forzado       = %.2f 1/m  (lambda = %.2f cm)" %
          (c.forzado.k_fundamental(), 2 * math.pi / c.forzado.k_fundamental() * 100))
    print("  duración del registro = %.1f s" % c.camara.duracion())
