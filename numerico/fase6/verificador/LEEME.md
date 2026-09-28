# Revisión independiente de la Fase 6 (rol «verifier» del equipo de agentes)

Pedida por Lucía para las implementaciones: después de convencerse uno mismo, una revisión
adversarial con el rol `~/GW_AI/agent-team/roles/verifier.md`. La hizo un agente separado el
2026-09-28, sin correr SPECTER ni las puertas: sólo Python liviano, en secuencia, con
`ulimit -v 2500000` y un pico de 155 MB. No modificó el repositorio. Sus archivos se copiaron
acá tal cual.

- `informe_verificador_fase6.md`: el informe, en inglés, con el bloque `claims` al final.
- `informe_verificador_fase6_ronda2.md` y `r2_*.py`: la segunda ronda, sobre las correcciones.
- `c1_*.py` … `c9_*.py`: su código, uno por afirmación. Las rutas apuntan al repositorio en la
  laptop (`/home/lucia/GW_AI/proyecto-final-piv`); en otra máquina hay que ajustarlas.
- `*.out`: salidas de texto de sus corridas de la réplica.

Qué se hizo con cada hallazgo: `numerico/fase6/REGISTRO_fase6.md`, P14, y [V-22] en
`DECISIONES.md`.
