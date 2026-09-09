# Ronda 01 — writer

## Trabajo realizado

- Leídos spec.json, state.json, todos los reportes previos, out/ y los scripts
  independientes pertinentes. inbox.jsonl no existe; no había direcciones humanas.
- Reconcilié out/notes.tex y out/provenance.json con los 16 claims verificados,
  el claim refutado y los 4 claims unclear del ledger. No modifiqué capa.py ni
  checks.py.
- Eliminé la promesa refutada de estabilidad sin cancelación para todo \(k>0\).
  La nota ahora acota la forma regular a la rama implementada \(kL_z<0.5\) y
  documenta el contraejemplo \(k=5\times10^{-11}\), \(L_z=10^{10}\),
  \(g_0=1\), \(g_1=7.3\), residuo \(7.629\times10^{-6}\), reproducible con
  work/verif/indep3.py.
- Añadí la consecuencia práctica de la desigualdad de validez,
  \(U\ll\nu\ell/h^2\), y dejé explícito que faltan \(U\) y \(\ell\)
  experimentales para decidir si la clausura aplica a la celda.
- Separé la procedencia de espectros y expansiones modales: las fórmulas completas
  apuntan a sus derivaciones, y ortogonalidad/coeficientes apuntan a
  work/verif/indep4.py. Esto corrige el puntero que antes prometía más que
  test_spectra_against_fd.
- Añadí salvedades explícitas para los puntos FC--Gram con \(z>L_z\) y para el
  margen del control C1: \(8.23\times10^{-5}\) frente a \(10^{-4}\), factor
  aproximado \(1.2\). El registro pasó de 26 a 31 entradas; las cinco nuevas
  claves están citadas en la nota y sus punteros existen.

## Puntos que siguen abiertos

1. Aplicabilidad experimental: sigue sin resolverse porque no hay valores de
   \(U\) ni de la menor escala horizontal \(\ell\).
2. Continuación FC--Gram: sigue sin resolverse porque no se implementó ni probó
   la rama DN en Fortran para los puntos con \(z>L_z\).
3. Margen de C1 fuera de la puerta congelada: el caso contratado pasa, pero no se
   verificó una variante con índices modales mayores.

No reejecuté checks, puerta, compilación LaTeX ni verificaciones físicas, conforme
al rol de writer. Sólo hice inspección estática de sintaxis JSON y existencia de
los nuevos punteros. Los resultados ejecutables vigentes continúan siendo los
registrados por el verifier: checks 5/5 y puerta 6/6.
