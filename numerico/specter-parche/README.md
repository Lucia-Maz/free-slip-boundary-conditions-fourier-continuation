# El aporte propio sobre SPECTER, como parche

El árbol de SPECTER **no está en este repositorio**: es código de terceros y son 5,4 MB
entre la copia intacta de upstream y la de trabajo. Lo que sí está es el parche, que es
exactamente lo que [D-26] define como el aporte de este proyecto. Desde el 2026-09-28 incluye
también la superficie **deformable** lineal ([D-49]) y la corrección de `freeslip_z` ([D-50]), y
desde la Fase 6 la superficie deformable **de orden N** ([D-51], [D-52]). En total son **1272
líneas agregadas y 20 borradas** en 5 archivos. Aparte van cuatro archivos nuevos:
- `fsorder.f90`, el camino de orden N, con 1063 líneas;
- tres tests: `laplace_dirneu.f90`, `laplace_neudir.f90` y `fs_orden.f90`.

## Cómo reconstruir el árbol de trabajo, en cualquier máquina

```bash
git clone https://github.com/mfontanaar/SPECTER.git SPECTER-upstream
cd SPECTER-upstream && git checkout 0ad1edb && cd ..

cp -r SPECTER-upstream SPECTER-trabajo
cd SPECTER-trabajo
git apply ../specter-parche/superficie-libre.patch
cp ../specter-parche/laplace_dirneu.f90 ../specter-parche/laplace_neudir.f90 \
   ../specter-parche/fs_orden.f90 src/tests/
cp ../specter-parche/fsorder.f90 src/boundary/
```

El commit `0ad1edb` es «Fixed bug preventing ny to be set to 1», del 2023-11-08. La copia
de upstream se deja intacta a propósito: es la referencia contra la que se compara, y de
ahí salen V5 y V6 de la puerta de la Fase 2.

Suma de verificación del parche, para que se pueda comprobar que es el mismo:

```
sha256  a1d31136b1543cec…  superficie-libre.patch   (2026-09-28, Fase 6; antes 5aa77f8ab65035da…, Fase 5)
sha256  8178a7e3b0df5f19…  fsorder.f90
sha256  e5a75527f2951174…  fs_orden.f90
```

## Qué toca el parche

| archivo | qué |
|---|---|
| `src/boundary/vboundary.f90` | `v_parsebc` acepta `freeslip`; `v_setup` carga las tablas de Neumann del plan en z; el guardián de `v_imposebc_and_project` acepta la clase 1 y despacha; subrutina nueva `freeslip_z`; `vdiagnostic` escribe además `freeslip_diagnostic.txt` con el residuo de **tensión** |
| `src/boundary/boundary_mod.fpp` | rama Dirichlet(z=0)–Neumann(z=Lz) de `laplace_z` |
| `src/specter.fpp` | engancha el test unitario nuevo |
| `src/SOLVERS_AND_BOUNDARY_CONDITIONS.md`, `bin/parameter.inp` | documentación de la cadena nueva |
| `src/tests/laplace_dirneu.f90` | **archivo nuevo**, 57 líneas: test unitario de la rama nueva |
| `src/boundary/vboundary.f90` ([D-50]) | `freeslip_z` recibe la traza de v*_z en la cara (`trace_z`) en vez del campo: vz ya no hace la ida y vuelta por las transformadas, que dejaba una disipación espuria por paso ([H-20]) |
| `src/boundary/vboundary.f90` (Fase 5) | cadena `freesurface` (sólo z = Lz), namelist `&freesurface`, `freesurface_imposebc` y sus auxiliares (acople exacto de la tensión en dos pasadas, trazas, `v_fs_init`, E/S de η, `fs_diagnostic`) |
| `src/boundary/boundary_mod.fpp` (Fase 5) | estado `fs_*`, rama Neumann(z=0)–Dirichlet(z=Lz) de `laplace_z`, `sol_project_fs` |
| `src/specter.fpp` (Fase 5) | `v_fs_init` después de la condición inicial, `v_fs_output` con los binarios, test nuevo |
| `src/tests/laplace_neudir.f90` | **archivo nuevo**: test unitario de la rama Neumann–Dirichlet |
| `src/boundary/fsorder.f90` (Fase 6) | **archivo nuevo**, 1063 líneas: el camino de orden N. Contiene trazas espectrales y la sonda β, las funcionales de orden N, la pasada, el subpaso con el precondicionador y Anderson real, las guardias, la presión inicial y el diagnóstico de orden |
| `src/boundary/boundary_mod.fpp` (Fase 6) | estado de orden N (`fs_order`, `fs_beta`, historia de Anderson, contadores), `INCLUDE 'fsorder.f90'`, y la salida opcional `dder` de `sol_project_fs` (derivadas normales del potencial) |
| `src/boundary/vboundary.f90` (Fase 6) | claves `fsorder`, `fsmean`, `fstol`, `fsmaxit` y `fsgen` de `&freesurface`; despacho al camino general; bloque de `v_fs_init` para el orden N; la llamada al diagnóstico de orden; y la guardia de Nyquist en y con ny = 1 |
| `src/specter.fpp`, `src/tests/fs_orden.f90` (Fase 6) | test unitario nuevo de trazas y funcionales, 129 líneas, que usan H2 y H3 de la puerta |

El camino no-deslizante queda intacto a propósito: `vz` se lleva al dominio (kx,ky,z)
sólo cuando alguna cara es `freeslip`. Eso importa porque V5 y V6 de la puerta comparan
contra él.

## Cómo verificar que quedó bien

```bash
/ruta/al/python verificacion/test_aceptacion_fase2.py
```

Tiene que dar **6/6**, y `verificacion/test_aceptacion_fase5.py` tiene que dar **8/8**.
`verificacion/test_aceptacion_fase6.py` da **4/10**, y eso es lo esperado: los seis fallos están
explicados en [V-21]. Sólo es regresión si cambian sus números.
(unos 8 minutos, 173 MB). La puerta compila SPECTER en un scratch temporal, corre unas quince
simulaciones chicas y tarda del orden de diez minutos. Está fijada por hash
(`6f0b4da4ba20617d483dfff0a20f8e2e11e82795b20351c259ffcf740ca4b1cc`): se puede leer, no se
edita para que pase.

## Licencia del parche

El repositorio está bajo MIT (`LICENSE` en la raíz), y eso cubre las líneas propias de
este parche y el test unitario. Pero el parche es **obra derivada de SPECTER**, cuyos
autores no publicaron una licencia (su README pide citar a Fontana, Bruno, Mininni &
Dmitruk 2020). Se distribuye acá como diff para que este trabajo se pueda reproducir; su
reutilización en otro contexto queda sujeta a lo que dispongan los autores de SPECTER.
