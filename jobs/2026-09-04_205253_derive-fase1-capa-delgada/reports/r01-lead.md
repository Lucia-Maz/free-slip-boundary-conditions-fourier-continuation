# Ronda 01 — plan del PI

## Estado observado

- Leídos `spec.json`, `state.json`, el inbox (ausente), `reports/r00-pi-staffing.md`, `out/` y la puerta de aceptación.
- No hay direcciones humanas persistidas ni claims previos `refuted` o `unclear`.
- La entrega no está hecha: `notes.tex` tiene 5 líneas de esqueleto, `provenance.json` tiene sólo `{}`, y faltan `capa.py` y `checks.py`.
- La puerta conserva el SHA-256 contratado y exige 6 chequeos: dos espectros, cociente de fricciones, sesgo PIV, convergencia de orden 2 y Laplace DN estable. No se ejecutó porque falta `capa.py`; por tanto el estado efectivo es no evaluable, no 6/6.

## Plan en vigor

1. **Worker:** construir la entrega completa y autoconsistente en `out/notes.tex`, `out/capa.py`, `out/checks.py` y `out/provenance.json`, usando exactamente la API y cadenas de borde de la puerta: derivar los dos espectros completos, la clausura 2D y sus desigualdades de validez, el sesgo PIV y los coeficientes DN en la base exponencial acotada de SPECTER (incluido `k=0`); implementar un solver FD integrado de orden 2; producir por código cada número citado y resolver cada `\\src{}`; incluir procedencia y salvedades no vacías; finalmente ejecutar **en serie** `out/checks.py`, el checker de provenance y la puerta fijada, y registrar resultados exactos. No tocar SPECTER, `/mnt`, `~/Desktop/PhD` ni escribir fuera del job.

## No resuelto

- Ninguna fórmula, implementación ni afirmación cuantitativa está todavía disponible para verificación.
- No puede declararse terminado hasta que el verifier confirme los resultados y los controles ejecutables pasen.

## WHAT MORE COULD BE DONE

1. **Verificación independiente integral** — necesaria para elevar los claims del worker a verificados; cuesta una ronda de auditoría y tres ejecuciones breves en serie.
2. **Correcciones dirigidas** — resolver cualquier fórmula, residuo, interfaz o fuente refutada; cuesta una ronda según los hallazgos.
3. **Reconciliación editorial final** — asegurar consistencia entre texto, código, números y salvedades; cuesta una ronda corta de writer.
4. **Contraste bibliográfico primario ampliado** — mejora la atribución entre derivación propia y literatura, aunque no sustituye los controles; cuesta una ronda de revisión.
