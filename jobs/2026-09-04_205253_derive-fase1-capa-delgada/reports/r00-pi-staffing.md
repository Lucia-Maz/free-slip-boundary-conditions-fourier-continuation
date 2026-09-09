# Ronda 00 — staffing PI

## Lectura y estado inicial

- Leídos `spec.json`, `state.json`, los artefactos iniciales de `out/`, la puerta de aceptación completa y `laplace_z` en SPECTER (sólo lectura). No hay reportes anteriores.
- `inbox.jsonl` no existe en el job pese a figurar en el brief; por tanto no había direcciones humanas persistidas que incorporar.
- Ledger: 0 claims (0 verified, 0 refuted, 0 unclear), 0 backlog, sin plan anterior y sin resultados de checks/aceptación.
- La puerta conserva el SHA-256 contratado `2db98ae5a286c64ff42e6bcbdf9e9d69ccbcae60cd631c16d28a93e9855c5527`. Exige 6/6 chequeos y cinco funciones públicas: `tasas_decaimiento`, `cociente_superficie_promedio`, `alfa`, `solver_1d` y `coef_laplace_dn`.
- Estado de artefactos: `out/notes.tex` es un esqueleto de 4 líneas, `out/provenance.json` es `{}`, y faltan `out/capa.py` y `out/checks.py`.
- SPECTER representa los modos no nulos con la base acotada `C_1 exp(k(z-Lz)) + C_2 exp(-kz)`; la rama Dirichlet-abajo/Neumann-arriba debe seguir esa convención y tratar `k=0` aparte.

## Plan en vigor

Se conserva **1 worker**. Un único dueño evita divergencias entre fórmulas, código, pruebas y claves de provenance, y respeta la prohibición de ejecutar scripts en paralelo. Se mantienen los defaults ya contratados: worker `high`; verifier independiente en Claude `xhigh`; writer al final.

1. **Worker:** producir un candidato completo y autoconsistente en `out/notes.tex`, `out/capa.py`, `out/checks.py` y `out/provenance.json`: derivar ambos espectros completos y la clausura promediada con desigualdades explícitas; obtener el sesgo PIV; derivar e implementar DN en la base estable de SPECTER con caso `k=0`; implementar exactamente la API fijada (incluido un solver FD realmente integrado y de orden 2); generar por código todo número citado; documentar qué se derivó, qué proviene de literatura y salvedades no vacías; finalmente correr, **en serie**, checks propios, checker de provenance y puerta 6/6, dejando salidas y fallos exactos en su reporte. No modificar ni compilar SPECTER ni escribir fuera del job.

Después, el verifier deberá auditar de forma independiente ecuaciones, límites, estabilidad, API, provenance y resultados ejecutables antes de que el writer haga la reconciliación final. Ningún resultado queda verificado en esta ronda.

## No resuelto

- No se evaluó aún ninguna fórmula física ni numérica; hacerlo corresponde al worker y confirmarlo al verifier.
- No se conoce todavía el resultado de los checks propios, del checker de provenance ni de la puerta: faltan los dos módulos Python.
- No hay claims `refuted` o `unclear` previos que resolver.

## WHAT MORE COULD BE DONE

1. **Auditoría independiente y puerta completa** — imprescindible para convertir claims del worker en verified; cuesta una ronda de verifier más ejecuciones pequeñas en serie.
2. **Reparación dirigida de hallazgos** — corregir cualquier residuo, orden, desigualdad o clave de provenance que refute el verifier; cuesta aproximadamente una ronda de worker según el número de fallos.
3. **Edición final cruzada** — unificar notación, trazabilidad y salvedades sin perder detalle técnico; cuesta una ronda corta de writer.
4. **Contraste bibliográfico más amplio** — fortalecer la separación entre derivación propia y resultados recordados; cuesta una ronda de búsqueda y revisión de fuentes primarias, sin cambiar la aceptación numérica.
