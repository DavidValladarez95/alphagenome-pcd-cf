# alphagenome-pcd-cf

Pipeline reproducible que evalúa el puntaje de impacto de variantes (AVI) de **AlphaGenome Atlas** (Google DeepMind) en genes de **discinesia ciliar primaria (DCP)** y **fibrosis quística (FQ)**, con datos de ClinVar, CFTR2, gnomAD, CADD y REVEL.

Código de análisis del artículo *Predictive Genomic AI for Early Respiratory Risk Stratification in PCD and CF*.

> **Uso exclusivo para investigación.** El AVI no es una herramienta de diagnóstico clínico y no debe usarse por sí solo para clasificar variantes ni para tomar decisiones de tratamiento.

## Qué hace

En una sola ejecución, `pipeline_atlas.py`:

1. Lee las variantes de ClinVar de siete genes: `CFTR`, `CCDC39`, `CCDC40`, `CCNO`, `DNAH5`, `DNAH11` y `ODAD1`.
2. Prepara el conjunto de análisis: variantes de un solo nucleótido (SNV) patogénicas o probablemente patogénicas frente a benignas o probablemente benignas, con al menos una estrella de revisión y sin clasificaciones conflictivas.
3. Empareja las variantes de CFTR2 (Sosnay et al., 2013) con sus datos clínicos: cloruro en sudor, insuficiencia pancreática y función pulmonar (FEV1).
4. Consulta AlphaGenome Atlas (AVI, componentes del AVI y su importancia), con caché y reintentos.
5. Obtiene datos públicos externos: frecuencias de gnomAD v4.1, CADD v1.7 y REVEL.
6. Calcula los resultados: AUC con intervalos de confianza por bootstrap, pruebas de Mann-Whitney, correlaciones de Spearman, comparación con otros predictores, umbrales de especificidad y su aplicación a variantes de significado incierto (VSI).
7. Genera las figuras en inglés y en español, tablas CSV, un Excel con todos los resultados, un resumen y los metadatos de la ejecución.

## Requisitos

- Python 3.10 o superior (probado con Python 3.14).
- Una clave de la API de AlphaGenome.
- Conexión a Internet para la primera ejecución (Atlas, gnomAD, CADD y REVEL).

## Instalación

```bash
git clone https://github.com/DavidValladarez95/alphagenome-pcd-cf.git
cd alphagenome-pcd-cf
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

Crea un archivo `.env` en la carpeta del proyecto con tu clave (este archivo está excluido del repositorio):

```
ALPHA_GENOME_API_KEY=tu_clave
```

## Datos de entrada

| Archivo | Origen | ¿Incluido? |
|---|---|---|
| `clinvar_download_20261027_1.zip` | Descargas de [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) de los siete genes (27 de octubre de 2026) | Sí |
| `NIHMS512446-supplement-12.xlsx` | Tabla Suplementaria 2 de Sosnay PR et al. *Nat Genet.* 2013;45(10):1160-7 ([doi:10.1038/ng.2745](https://doi.org/10.1038/ng.2745)); disponible en PubMed Central como manuscrito NIHMS512446 | No: descárgala y colócala en la carpeta del proyecto |

## Uso

```bash
python pipeline_atlas.py --clinvar clinvar_download_20261027_1.zip --cftr2 NIHMS512446-supplement-12.xlsx
```

Opciones útiles:

| Opción | Descripción |
|---|---|
| `--salida CARPETA` | Carpeta de resultados (por defecto `resultados_AAAAMMDD_HHMMSS`) |
| `--prueba N` | Limita el análisis a unas N variantes para una prueba rápida |
| `--sin-consulta` | No consulta Atlas; usa solo la caché |
| `--sin-cadd` | Omite CADD |
| `--sin-vsi` | No consulta ni analiza las variantes de significado incierto |
| `--hilos N` / `--lote N` | Consultas en paralelo a Atlas (por defecto 4 hilos, lotes de 100) |

Las respuestas de Atlas se guardan en `cache_atlas.csv` y los datos externos en `cache_externos/`. Las ejecuciones posteriores reutilizan esa caché y solo consultan lo que falta.

## Resultados

Cada ejecución crea una carpeta nueva con:

- `01_` a `09_*.csv`: datos intermedios (variantes de ClinVar, conjunto de análisis, AVI por variante, CFTR2, datos externos y VSI).
- `tabla_*.csv`: tablas de resultados del artículo.
- `Figura1` a `Figura4` y `FiguraS1`, en inglés (`_en`) y en español (`_es`).
- `resultados_completos.xlsx`: todas las tablas en un solo Excel.
- `RESUMEN.md`: resumen de los resultados principales.
- `pipeline.log` y `metadatos.json`: registro completo, versiones de los paquetes y huellas SHA-256 de los archivos de entrada.

## Reproducibilidad

- Los intervalos de confianza se calculan por bootstrap (1000 remuestreos para el AUC y 2000 para Spearman) con semilla fija (`1`).
- `metadatos.json` registra las versiones de Python y de cada paquete, y la huella de los archivos de entrada.
- Las predicciones de AlphaGenome Atlas pueden cambiar entre versiones del modelo. El análisis original se hizo con `alphagenome` 0.10.0.

## Fuentes de datos

- [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) (NCBI)
- [CFTR2](https://cftr2.org/) mediante Sosnay et al., 2013
- [AlphaGenome](https://deepmind.google/) (Google DeepMind)
- [gnomAD v4.1](https://gnomad.broadinstitute.org/)
- [CADD v1.7](https://cadd.gs.washington.edu/)
- REVEL mediante [MyVariant.info](https://myvariant.info/) (dbNSFP)
