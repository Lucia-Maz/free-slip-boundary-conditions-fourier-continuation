# Ronda 01 — worker 1

## Trabajo realizado

- Leídos `spec.json`, `state.json`, los reportes anteriores, la puerta completa y `laplace_z`; `inbox.jsonl` no existe. No había claims `refuted` o `unclear` previos.
- Completé `out/notes.tex`, `out/capa.py`, `out/checks.py` y `out/provenance.json` sin modificar SPECTER ni archivos fuera del job.
- `notes.tex` contiene ambos espectros completos y expansiones modales, la reducción promediada, las desigualdades de validez, el sesgo PIV, la rama Laplace DN estable, procedencia y salvedades no vacías.
- `capa.py` implementa las cinco funciones de la API congelada. El solver integra el sistema FD semidiscreto de orden dos mediante la acción de su exponencial; no inserta la solución continua.
- La rama DN usa `phi=C1 exp(k(z-Lz))+C2 exp(-kz)`, con `C1=(g1/k+r*g0)/(1+r^2)`, `C2=(g0-r*g1/k)/(1+r^2)`, `r=exp(-kLz)`; para `k=0`, los coeficientes lineales son pendiente `g1` e intercepto `g0`.
- Compilé `notes.tex` dos veces con `pdflatex`: PDF de 6 páginas, sin errores ni referencias indefinidas. Quedan sólo dos avisos `Overfull` tipográficos menores (4.073 pt y 0.546 pt).

## Resultados verificables

- Claim: `alpha_DN = pi^2 nu/(4 h^2)` y `alpha_DD = pi^2 nu/h^2`; el cociente es exactamente 4 en la implementación.
- Claim: para el modo DN lento, `u(h)/promedio = pi/2 = 1.570796326795`.
- Claim: la clausura con velocidad promediada tiene `beta=<F^2>=pi^2/8=1.233700550136`.
- Claim: checks propios finales: **5/5**, exit 0.
- Claim: espectros contra FD al pasar `nz=120 -> 240`: error relativo máximo DD `9.135e-04 -> 2.284e-04`, DN `6.995e-04 -> 1.749e-04`; `p=2.000` en ambos.
- Claim: solver al pasar `nz=48 -> 96`: error de tasa DD `3.523e-03 -> 8.808e-04`, DN `2.202e-04 -> 5.505e-05`; `p=2.000` en ambos.
- Claim: siete casos Laplace DN propios, incluidos `kLz=700`, `k=1e-10` y datos complejos, dieron residuo de borde máximo `4.441e-16`.
- Claim: se rechazaron correctamente 9 rutas de entrada inválida.
- Claim: checker de provenance final: **26 tags, 26 entries, 0 errors, 0 warnings**, exit 0.
- Claim: puerta oficial final: **6/6**, exit 0; errores espectrales `8.23e-05` DD y `6.30e-05` DN; orden del solver `2.00`; cinco casos DN bajo residuo `1e-10`.
- Claim: el SHA-256 de la puerta inspeccionada fue el contratado: `2db98ae5a286c64ff42e6bcbdf9e9d69ccbcae60cd631c16d28a93e9855c5527`.

## No resuelto / no verificado

- No se usaron datos experimentales ni se midieron viscosidad, espesor efectivo, rugosidad, deslizamiento o contaminación superficial; el valor para 6 mm sigue siendo paramétrico en `nu`.
- No se verificaron deformación, Marangoni, ondas, capa límite turbulenta, profundidad variable ni forzamiento con perfil vertical arbitrario.
- La clausura de modo lento es asintótica bajo las desigualdades declaradas, no una identidad para Navier–Stokes no lineal.
- La rama DN se verificó aisladamente en Python; no se implementó/compiló en Fortran, no se corrió MPI y no se probó la distribución espectral interna de SPECTER.
