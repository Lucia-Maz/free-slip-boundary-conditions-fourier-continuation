# El aporte propio sobre SPECTER, como parche

El árbol de SPECTER **no está en este repositorio**: es código de terceros y son 5,4 MB
entre la copia intacta de upstream y la de trabajo. Lo que sí está es el parche, que es
exactamente lo que [D-26] define como el aporte de este proyecto: **237 líneas agregadas
y 22 borradas** en 5 archivos, más un archivo nuevo.

## Cómo reconstruir el árbol de trabajo, en cualquier máquina

```bash
git clone https://github.com/mfontanaar/SPECTER.git SPECTER-upstream
cd SPECTER-upstream && git checkout 0ad1edb && cd ..

cp -r SPECTER-upstream SPECTER-trabajo
cd SPECTER-trabajo
git apply ../specter-parche/superficie-libre.patch
cp ../specter-parche/laplace_dirneu.f90 src/tests/
```

El commit `0ad1edb` es «Fixed bug preventing ny to be set to 1», del 2023-11-08. La copia
de upstream se deja intacta a propósito: es la referencia contra la que se compara, y de
ahí salen V5 y V6 de la puerta de la Fase 2.

Suma de verificación del parche, para que se pueda comprobar que es el mismo:

```
sha256  bd06b93ce598998a…  superficie-libre.patch
```

## Qué toca el parche

| archivo | qué |
|---|---|
| `src/boundary/vboundary.f90` | `v_parsebc` acepta `freeslip`; `v_setup` carga las tablas de Neumann del plan en z; el guardián de `v_imposebc_and_project` acepta la clase 1 y despacha; subrutina nueva `freeslip_z`; `vdiagnostic` escribe además `freeslip_diagnostic.txt` con el residuo de **tensión** |
| `src/boundary/boundary_mod.fpp` | rama Dirichlet(z=0)–Neumann(z=Lz) de `laplace_z` |
| `src/specter.fpp` | engancha el test unitario nuevo |
| `src/SOLVERS_AND_BOUNDARY_CONDITIONS.md`, `bin/parameter.inp` | documentación de la cadena nueva |
| `src/tests/laplace_dirneu.f90` | **archivo nuevo**, 57 líneas: test unitario de la rama nueva |

El camino no-deslizante queda intacto a propósito: `vz` se lleva al dominio (kx,ky,z)
sólo cuando alguna cara es `freeslip`. Eso importa porque V5 y V6 de la puerta comparan
contra él.

## Cómo verificar que quedó bien

```bash
/ruta/al/python verificacion/test_aceptacion_fase2.py
```

Tiene que dar **6/6**. La puerta compila SPECTER en un scratch temporal, corre unas quince
simulaciones chicas y tarda del orden de diez minutos. Está fijada por hash
(`6f0b4da4ba20617d483dfff0a20f8e2e11e82795b20351c259ffcf740ca4b1cc`): se puede leer, no se
edita para que pase.

## Licencia del parche

El repositorio está bajo MIT (`LICENSE` en la raíz), y eso cubre las líneas propias de
este parche y el test unitario. Pero el parche es **obra derivada de SPECTER**, cuyos
autores no publicaron una licencia (su README pide citar a Fontana, Bruno, Mininni &
Dmitruk 2020). Se distribuye acá como diff para que este trabajo se pueda reproducir; su
reutilización en otro contexto queda sujeta a lo que dispongan los autores de SPECTER.
