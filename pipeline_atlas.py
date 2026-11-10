#!/usr/bin/env python3
"""
Pipeline de análisis del artículo "Predictive Genomic AI for Early Respiratory Risk
Stratification in PCD and CF".

Evalúa el puntaje de impacto de variantes (AVI) de AlphaGenome Atlas en siete genes de
discinesia ciliar primaria (DCP) y fibrosis quística (FQ). Etapas:

  1. Lectura y filtrado de variantes de ClinVar.
  2. Emparejamiento con los datos clínicos de CFTR2 (Sosnay et al., 2013).
  3. Consulta de AlphaGenome Atlas, con caché en disco.
  4. Obtención de datos externos: gnomAD v4.1, CADD v1.7 y REVEL.
  5. Análisis estadístico principal, análisis adicionales y de robustez.
  6. Figuras en inglés y en español, tablas CSV, Excel, resumen y metadatos.

Uso:
  python pipeline_atlas.py --clinvar clinvar_download_20261027_1.zip --cftr2 NIHMS512446-supplement-12.xlsx

La clave de la API se lee de la variable ALPHA_GENOME_API_KEY o de un archivo .env.
"""
import argparse
import concurrent.futures as cf
import csv
import glob
import hashlib
import io
import json
import logging
import os
import platform
import re
import sys
import time
import traceback
import zipfile
from datetime import datetime

import requests

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------------------
# Constantes del estudio
# --------------------------------------------------------------------------------------
GENES = ["CFTR", "CCDC39", "CCDC40", "CCNO", "DNAH5", "DNAH11", "ODAD1"]
SCORERS = ["AVI_SCORE", "AVI_SCORE_FEATURE_IMPORTANCE", "AVI_SCORE_MODEL_FEATURES"]
FILAS_ESPERADAS = 37  # filas por variante en la caché: 1 AVI + 18 importancias + 18 componentes
# Estados de revisión de ClinVar con al menos una estrella y sin conflicto de clasificación
REVISION_ACEPTADA = [
    "criteria provided, single submitter",
    "criteria provided, multiple submitters, no conflicts",
    "reviewed by expert panel",
    "practice guideline",
]
ESTRELLAS = {
    "criteria provided, single submitter": 1,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting classifications": 1,
    "reviewed by expert panel": 3,
    "practice guideline": 4,
    "no assertion criteria provided": 0,
    "no classification provided": 0,
}
# Agrupación de los 18 componentes del AVI en mecanismos biológicos (agrupación propia del estudio)
MECANISMO = {
    "MERGED_SPLICING": "Splicing",
    "MAX_ABS_ATAC": "Regulación/expresión", "MAX_ABS_DNASE": "Regulación/expresión",
    "MAX_ABS_CHIP_TF": "Regulación/expresión", "MAX_ABS_CHIP_HISTONE": "Regulación/expresión",
    "MAX_ABS_CAGE": "Regulación/expresión", "MAX_ABS_PROCAP": "Regulación/expresión",
    "MAX_ABS_RNA_SEQ": "Regulación/expresión", "MAX_ABS_CONTACT_MAPS": "Regulación/expresión",
    "MAX_ABS_POLYADENYLATION": "Regulación/expresión",
    "ALPHAMISSENSE": "Proteína", "PROTEIN_TERMINATION": "Proteína", "START_LOST": "Proteína", "STOP_LOST": "Proteína",
    "CACTUS_241_WAY": "Conservación", "PHASTCONS_470_WAY": "Conservación",
    "IS_INSERTION": "Tipo de variante", "IS_DELETION": "Tipo de variante",
}
# Genes de DCP agrupados según la gravedad de la afectación pulmonar descrita en la literatura
GRUPO_DCP = {
    "CCDC39": "Peor función pulmonar (CCDC39, CCDC40, CCNO)",
    "CCDC40": "Peor función pulmonar (CCDC39, CCDC40, CCNO)",
    "CCNO": "Peor función pulmonar (CCDC39, CCDC40, CCNO)",
    "DNAH11": "Mejor función pulmonar (DNAH11, ODAD1)",
    "ODAD1": "Mejor función pulmonar (DNAH11, ODAD1)",
}
TIPOS = ["Missense", "Intrónica/UTR", "Sinónima", "Sin sentido (stop)", "Sitio de splicing"]
GNOMAD_API = "https://gnomad.broadinstitute.org/api"
GNOMAD_DATASET = "gnomad_r4"
CADD_URL = "https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz"
MYVARIANT_URL = "https://myvariant.info/v1/variant"
SEMILLA = 1
B_AUC = 2000  # remuestreos bootstrap para AUC y PR-AUC
UMBRAL_RARA = 0.001  # frecuencia máxima de una variante rara, criterio usado en el entrenamiento del AVI
# Grupos de ascendencia que gnomAD v4 usa para grpmax/FAF (excluye ami, asj, fin, mid y remaining)
GRUPOS_FAF = ["afr", "amr", "eas", "nfe", "sas"]
VENTANA_CLUSTER = 5000  # tamaño en pb de la ventana genómica para el bootstrap por clusters
B_RHO = 2000  # remuestreos bootstrap para la correlación de Spearman

log = logging.getLogger("pipeline")


# --------------------------------------------------------------------------------------
# Utilidades
# --------------------------------------------------------------------------------------
def configurar_logs(carpeta):
    """Envía el registro a la consola y a pipeline.log dentro de la carpeta de resultados."""
    log.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s | %(levelname)-7s | %(message)s", "%Y-%m-%d %H:%M:%S")
    fh = logging.FileHandler(os.path.join(carpeta, "pipeline.log"), encoding="utf-8")
    fh.setFormatter(fmt)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    log.handlers = [fh, sh]


def sha256(ruta):
    """Huella SHA-256 de un archivo, leída por bloques de 1 MB."""
    h = hashlib.sha256()
    with open(ruta, "rb") as f:
        for bloque in iter(lambda: f.read(1 << 20), b""):
            h.update(bloque)
    return h.hexdigest()


def version(paquete):
    """Versión instalada de un paquete, o 'no instalado'."""
    try:
        from importlib.metadata import version as v
        return v(paquete)
    except Exception:  # noqa: BLE001
        return "no instalado"


def paso(titulo):
    """Escribe en el registro el encabezado de una etapa del pipeline."""
    log.info("=" * 80)
    log.info(titulo)
    log.info("=" * 80)


# --------------------------------------------------------------------------------------
# 1. ClinVar
# --------------------------------------------------------------------------------------
def leer_clinvar(ruta):
    """Lee las descargas de ClinVar desde una carpeta, un zip o un único archivo .tsv y las concatena."""
    tablas = []
    if os.path.isdir(ruta):
        archivos = sorted(glob.glob(os.path.join(ruta, "*.tsv")) + glob.glob(os.path.join(ruta, "*.txt")))
        for a in archivos:
            tablas.append(pd.read_csv(a, sep="\t", dtype=str, index_col=False))
            log.info(f"ClinVar: leído {os.path.basename(a)} ({len(tablas[-1])} filas)")
    elif zipfile.is_zipfile(ruta):
        with zipfile.ZipFile(ruta) as z:
            for n in sorted(z.namelist()):
                if n.endswith((".tsv", ".txt")) and not os.path.basename(n).startswith("."):
                    with z.open(n) as f:
                        tablas.append(pd.read_csv(io.TextIOWrapper(f, encoding="utf-8"), sep="\t", dtype=str, index_col=False))
                    log.info(f"ClinVar: leído {n} ({len(tablas[-1])} filas)")
    else:
        tablas.append(pd.read_csv(ruta, sep="\t", dtype=str, index_col=False))
    if not tablas:
        raise SystemExit("No se encontraron archivos .tsv de ClinVar.")
    d = pd.concat(tablas, ignore_index=True)
    log.info(f"ClinVar: {len(d)} registros en total")
    return d


def categoria(c):
    """Agrupa la clasificación germinal de ClinVar en P/PP, B/PB, VSI, Conflictiva u Otra."""
    c = str(c).lower()
    if "conflicting" in c:
        return "Conflictiva"
    if "pathogenic" in c and "benign" not in c:
        return "P/PP"
    if "benign" in c and "pathogenic" not in c:
        return "B/PB"
    if "uncertain" in c:
        return "VSI"
    return "Otra"


def tipo_variante(c):
    """Asigna un tipo de variante a partir de la consecuencia molecular de ClinVar.

    El orden de las comprobaciones define la prioridad cuando una variante tiene varias consecuencias.
    """
    c = str(c)
    if "nonsense" in c:
        return "Sin sentido (stop)"
    if "splice donor" in c or "splice acceptor" in c:
        return "Sitio de splicing"
    if "initiator" in c or "stop lost" in c:
        return "Pérdida de inicio/parada"
    if "missense" in c:
        return "Missense"
    if "synonymous" in c:
        return "Sinónima"
    if "intron" in c or "UTR" in c or "non-coding" in c:
        return "Intrónica/UTR"
    return "Sin anotación"


def preparar_clinvar(d, out):
    """Filtra las SNV de los genes del estudio y construye el conjunto de análisis y el de VSI.

    Devuelve todas las SNV válidas, el conjunto de análisis (P/PP y B/PB con al menos una estrella)
    y las VSI con al menos una estrella. Las coordenadas se toman del SPDI canónico (base 0)
    y se convierten a base 1.
    """
    d = d.copy()
    d["Gen"] = d["Gene(s)"].map(lambda s: next((g for g in GENES if g in str(s).split("|")), None))
    sin_gen = d["Gen"].isna().sum()
    if sin_gen:
        log.warning(f"{sin_gen} registros sin uno de los 7 genes del estudio: se descartan")
    d = d[d["Gen"].notna()]
    d = d.drop_duplicates("VariationID")
    snv = d[d["Variant type"] == "single nucleotide variant"].copy()
    log.info(f"SNV: {len(snv)} (descartados {len(d) - len(snv)} registros que no son SNV, p. ej., haplotipos)")
    sp = snv["Canonical SPDI"].fillna("").str.split(":", expand=True).reindex(columns=range(4))
    snv["REF"], snv["ALT"] = sp[2], sp[3]
    snv["Posición GRCh38"] = pd.to_numeric(sp[1], errors="coerce") + 1
    valido = (snv["REF"].str.len() == 1) & (snv["ALT"].str.len() == 1) & snv["Posición GRCh38"].notna() & snv["GRCh38Chromosome"].notna()
    log.info(f"SNV sin REF/ALT de una letra o sin coordenada: {(~valido).sum()} (se descartan)")
    snv = snv[valido].copy()
    snv["Posición GRCh38"] = snv["Posición GRCh38"].astype(int)
    misma = snv["REF"] == snv["ALT"]
    if misma.any():
        log.info(f"Registros con REF = ALT (alelo de referencia, sin cambio): {misma.sum()} -> se excluyen: "
                 + ", ".join(snv.loc[misma, "Name"].astype(str)))
    snv = snv[~misma].copy()
    snv["Consulta"] = "chr" + snv["GRCh38Chromosome"].astype(str) + ":" + snv["Posición GRCh38"].astype(str) + ":" + snv["REF"] + ">" + snv["ALT"]
    snv["Categoría"] = snv["Germline classification"].map(categoria)
    snv["Estrellas"] = snv["Germline review status"].map(ESTRELLAS).fillna(0).astype(int)
    snv["Tipo de variante"] = snv["Molecular consequence"].map(tipo_variante)
    snv["c_HGVS"] = snv["Name"].str.extract(r"\):(c\.[^ ]+)")[0]
    snv = snv.sort_values(["Gen", "Posición GRCh38"])
    snv.to_csv(os.path.join(out, "01_clinvar_snv_todas.csv"), index=False)
    an = snv[snv["Categoría"].isin(["P/PP", "B/PB"]) & snv["Germline review status"].isin(REVISION_ACEPTADA)].copy()
    an = an.sort_values(["Gen", "Categoría", "Posición GRCh38"])
    an.to_csv(os.path.join(out, "02_conjunto_analisis.csv"), index=False)
    log.info(f"Conjunto de análisis: {len(an)} SNV (P/PP = {(an['Categoría'] == 'P/PP').sum()}, B/PB = {(an['Categoría'] == 'B/PB').sum()})")
    vsi = snv[(snv["Categoría"] == "VSI") & snv["Germline review status"].isin(REVISION_ACEPTADA)].copy()
    vsi.to_csv(os.path.join(out, "02b_vsi.csv"), index=False)
    log.info(f"Variantes de significado incierto (VSI) con >=1 estrella: {len(vsi)}")
    log.info("Conteos por gen y categoría (todas las SNV):\n" + pd.crosstab(snv["Gen"], snv["Categoría"], margins=True).to_string())
    return snv, an, vsi


# --------------------------------------------------------------------------------------
# 2. CFTR2 (Sosnay et al. 2013, Tabla Suplementaria 2)
# --------------------------------------------------------------------------------------
def leer_cftr2(ruta, snv, out):
    """Lee la tabla de fenotipos de CFTR2 y empareja cada variante con su coordenada de ClinVar.

    El emparejamiento se hace por nombre c.HGVS y, si falla, por identificador rsID.
    """
    # Las dos primeras filas son encabezados; las columnas se seleccionan por posición
    t = pd.read_excel(ruta, sheet_name="PhenotypeSummaryTable", header=None).iloc[2:].reset_index(drop=True)
    t = t[[0, 1, 2, 3, 7, 8, 9, 10, 11, 12, 15, 28]]
    t.columns = ["hgvs", "proteina", "nombre_clasico", "rsID", "cloruro_sudor", "n_sudor", "prop_IP", "n_IP",
                 "fev1", "n_fev1", "tipo_cftr2", "clasificacion_cftr2"]
    t = t[t["hgvs"].notna()].copy()
    for c in ["cloruro_sudor", "n_sudor", "prop_IP", "n_IP", "fev1", "n_fev1"]:
        t[c] = pd.to_numeric(t[c], errors="coerce")
    t["fev1"] = t["fev1"] * 100  # la tabla lo informa como fracción
    t["hgvs"] = t["hgvs"].astype(str).str.strip()
    log.info(f"CFTR2: {len(t)} variantes; clasificación: {t['clasificacion_cftr2'].value_counts().to_dict()}")
    cv = snv[snv["Gen"] == "CFTR"][["c_HGVS", "Consulta", "dbSNP ID"]]
    m = t.merge(cv.rename(columns={"c_HGVS": "hgvs"}), on="hgvs", how="left")
    falta = m["Consulta"].isna() & m["rsID"].notna()
    rs = cv.dropna(subset=["dbSNP ID"]).drop_duplicates("dbSNP ID").set_index("dbSNP ID")["Consulta"]
    m.loc[falta, "Consulta"] = m.loc[falta, "rsID"].astype(str).str.strip().map(rs)
    m = m.drop_duplicates("hgvs")
    # Solo las sustituciones simples (p. ej., c.1521G>A) pueden consultarse en Atlas
    es_snv = m["hgvs"].str.contains(r"^c\.[-*\d+]+[ACGT]>[ACGT]$", regex=True)
    m["es_SNV"] = es_snv
    # CFTR está en la hebra positiva: los alelos del nombre c. deben coincidir con REF>ALT genómicos.
    al_c = m["hgvs"].str.extract(r"([ACGT])>([ACGT])$")
    al_g = m["Consulta"].astype(str).str.extract(r":([ACGT])>([ACGT])$")
    distinto = es_snv & m["Consulta"].notna() & ((al_c[0] != al_g[0]) | (al_c[1] != al_g[1]))
    if distinto.any():
        log.info("CFTR2: alelos que no coinciden con la coordenada de ClinVar (p. ej., el genoma de referencia ya porta el alelo), "
                 f"se excluyen: {m.loc[distinto, 'hgvs'].tolist()}")
        m.loc[distinto, "Consulta"] = np.nan
    log.info(f"CFTR2: {es_snv.sum()} son SNV; {m.loc[es_snv, 'Consulta'].notna().sum()} emparejadas con ClinVar")
    sin = m[es_snv & m["Consulta"].isna()]["hgvs"].tolist()
    if sin:
        log.info(f"CFTR2: SNV sin coordenada en la descarga de ClinVar: {sin}")
    m.to_csv(os.path.join(out, "04_cftr2_emparejado.csv"), index=False)
    return m


# --------------------------------------------------------------------------------------
# 3. AlphaGenome Atlas (consulta con caché)
# --------------------------------------------------------------------------------------
# Columnas de la caché de Atlas: una fila por variante, scorer y componente
CAMPOS = ["consulta", "scorer", "gen_objetivo", "track", "score", "cuantil"]


def crear_cliente():
    """Crea el cliente de AlphaGenome Atlas con la clave de ALPHA_GENOME_API_KEY."""
    try:
        from dotenv import load_dotenv
        load_dotenv(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".env"))
        load_dotenv()
    except ImportError:
        log.warning("python-dotenv no instalado; se usa solo la variable de entorno")
    clave = (os.environ.get("ALPHA_GENOME_API_KEY") or "").strip()
    if not clave:
        raise SystemExit("No encuentro ALPHA_GENOME_API_KEY (archivo .env o variable de entorno).")
    from alphagenome.atlas import atlas
    return atlas.create(clave, timeout=60)


def aplanar(consulta, resultados):
    """Convierte la respuesta de Atlas (un AnnData por scorer) en filas con las columnas de CAMPOS."""
    filas = []
    for scorer, ad in (resultados or {}).items():
        x = np.asarray(ad.X)
        x = x.reshape(1, -1) if x.ndim == 1 else x
        q = np.asarray(ad.layers["quantiles"]) if "quantiles" in ad.layers else None
        if q is not None and q.ndim == 1:
            q = q.reshape(1, -1)
        nombres = list(ad.var["name"]) if "name" in ad.var.columns else list(ad.var_names)
        obs = ad.obs if ad.obs is not None else pd.DataFrame()
        for i in range(x.shape[0]):
            gen = str(obs["gene_name"].iloc[i]) if "gene_name" in obs.columns and i < len(obs) else ""
            for j, tr in enumerate(nombres):
                filas.append({"consulta": consulta, "scorer": scorer, "gen_objetivo": gen, "track": str(tr),
                              "score": float(x[i, j]), "cuantil": float(q[i, j]) if q is not None else ""})
    return filas


def consultar_una(cliente, consulta, reintentos=4):
    """Consulta una variante en formato chr:pos:REF>ALT con reintentos y espera exponencial.

    Devuelve (consulta, filas, error); si la consulta falla, filas está vacía y error describe la causa.
    """
    from alphagenome.data import genome
    crom, pos, cambio = consulta.split(":")
    ref, alt = cambio.split(">")
    v = genome.Variant(chromosome=crom, position=int(pos), reference_bases=ref, alternate_bases=alt)
    error = ""
    for intento in range(1, reintentos + 1):
        try:
            return consulta, aplanar(consulta, cliente.query_variant(v, requested_scorers=SCORERS)), ""
        except Exception as e:  # noqa: BLE001
            error = f"{type(e).__name__}: {e}"
            if "invalid argument" in str(e).lower():
                break  # error de la variante, no de la conexión: reintentar no cambia el resultado
            time.sleep(2 ** intento)
    return consulta, [], error


def importar_resultados_previos(cache):
    """Crea la caché a partir de resultados_atlas_largo.csv y resultados_extra.csv si la caché no existe."""
    previos = [f for f in ["resultados_atlas_largo.csv", "resultados_extra.csv"] if os.path.exists(f)]
    if os.path.exists(cache) or not previos:
        return
    partes = [pd.read_csv(f) for f in previos]
    pd.concat(partes, ignore_index=True)[CAMPOS].to_csv(cache, index=False)
    log.info(f"Caché inicial creada a partir de resultados previos: {previos}")


def consultas_completas(cache):
    """Variantes de la caché que tienen todas las filas esperadas (FILAS_ESPERADAS)."""
    if not os.path.exists(cache):
        return set()
    c = pd.read_csv(cache, usecols=["consulta", "scorer"])
    n = c.groupby("consulta")["scorer"].size()
    return set(n[n == FILAS_ESPERADAS].index)


def consultar_atlas(consultas, cache, out, hilos, lote, sin_consulta):
    """Consulta en Atlas las variantes que faltan en la caché y devuelve la caché completa.

    Las consultas se hacen en paralelo por lotes y cada lote se escribe en disco al terminar,
    de modo que una ejecución interrumpida se reanuda sin repetir lo ya consultado.
    """
    importar_resultados_previos(cache)
    hechas = consultas_completas(cache)
    pendientes = [c for c in dict.fromkeys(consultas) if c not in hechas]
    log.info(f"Atlas: {len(set(consultas))} consultas necesarias | en caché completas: {len(set(consultas) & hechas)} | pendientes: {len(pendientes)}")
    errores = []
    if pendientes and sin_consulta:
        log.warning("Modo --sin-consulta: las consultas pendientes no se harán y quedarán sin resultado.")
    elif pendientes:
        # Se eliminan las consultas incompletas para que no queden filas duplicadas al repetirlas
        if os.path.exists(cache):
            c = pd.read_csv(cache)
            c = c[c["consulta"].isin(hechas)]
            c.to_csv(cache, index=False)
        cliente = crear_cliente()
        nuevo = not os.path.exists(cache)
        t0 = time.time()
        with open(cache, "a", newline="", encoding="utf-8") as fs:
            w = csv.DictWriter(fs, fieldnames=CAMPOS)
            if nuevo:
                w.writeheader()
            for k in range(0, len(pendientes), lote):
                bloque = pendientes[k:k + lote]
                with cf.ThreadPoolExecutor(max_workers=hilos) as ex:
                    for consulta, filas, err in ex.map(lambda q: consultar_una(cliente, q), bloque):
                        if err:
                            errores.append({"consulta": consulta, "error": err})
                            log.warning(f"Error en {consulta}: {err}")
                        else:
                            w.writerows(filas)
                fs.flush()
                hecho = min(k + lote, len(pendientes))
                log.info(f"Atlas: {hecho}/{len(pendientes)} consultadas ({time.time() - t0:.0f} s)")
    pd.DataFrame(errores, columns=["consulta", "error"]).to_csv(os.path.join(out, "errores_atlas.csv"), index=False)
    log.info(f"Atlas: errores = {len(errores)}")
    return pd.read_csv(cache)


def resumir_atlas(largo, out):
    """Resume la caché en una fila por variante: AVI, cuantil, componente principal y componentes.

    El componente principal es el de mayor importancia en valor absoluto.
    """
    avi = largo[largo["scorer"] == "AVI_SCORE"].drop_duplicates("consulta").set_index("consulta")
    r = pd.DataFrame({"AVI": avi["score"], "AVI_cuantil": pd.to_numeric(avi["cuantil"], errors="coerce")})
    imp = largo[largo["scorer"] == "AVI_SCORE_FEATURE_IMPORTANCE"]
    piv_imp = imp.pivot_table(index="consulta", columns="track", values="score", aggfunc="first")
    piv_imp.columns = ["imp_" + c for c in piv_imp.columns]
    feat = largo[largo["scorer"] == "AVI_SCORE_MODEL_FEATURES"]
    piv_feat = feat.pivot_table(index="consulta", columns="track", values="score", aggfunc="first", dropna=False)
    piv_feat.columns = ["feat_" + c for c in piv_feat.columns]
    top = imp.assign(a=imp["score"].abs()).sort_values("a", ascending=False).drop_duplicates("consulta").set_index("consulta")
    r["Componente_principal"] = top["track"]
    r["Importancia_componente_principal"] = top["score"]
    r["Mecanismo_principal"] = top["track"].map(MECANISMO)
    r = r.join(piv_imp).join(piv_feat)
    r["AlphaMissense"] = r.get("feat_ALPHAMISSENSE")
    r.index.name = "Consulta"
    r.to_csv(os.path.join(out, "03_atlas_por_variante.csv"))
    log.info(f"Atlas: {len(r)} variantes con AVI resumidas")
    return r


# --------------------------------------------------------------------------------------
# 3b. Datos externos públicos: gnomAD, CADD, REVEL
# --------------------------------------------------------------------------------------
def _post(url, reintentos=5, **kw):
    """Petición POST con reintentos y espera creciente; trata el código 429 como error recuperable."""
    for intento in range(1, reintentos + 1):
        try:
            r = requests.post(url, timeout=300, **kw)
            if r.status_code == 429:
                raise RuntimeError("límite de peticiones (429)")
            r.raise_for_status()
            return r
        except Exception as e:  # noqa: BLE001
            log.warning(f"Petición a {url} falló (intento {intento}/{reintentos}): {e}")
            if intento == reintentos:
                raise
            time.sleep(10 * intento)


def obtener_gnomad(genes, carpeta_cache):
    """Frecuencias de gnomAD v4.1 por gen (exomas + genomas), incluidas la población 'amr' y la faf95.

    Cada respuesta se guarda en carpeta_cache y se reutiliza en ejecuciones posteriores.
    """
    os.makedirs(carpeta_cache, exist_ok=True)
    base = ("query($g: String!) { gene(gene_symbol: $g, reference_genome: GRCh38) { "
            f"variants(dataset: {GNOMAD_DATASET}) {{ variant_id "
            "exome { ac an %s populations { id ac an } } genome { ac an %s populations { id ac an } } } } }")
    # Primero se pide la frecuencia alélica filtrada (faf95, grupo de ascendencia máximo). Si la API
    # no acepta ese campo, se usa la consulta básica y se registra la advertencia.
    faf = "faf95 { popmax popmax_population }"
    consultas = [("faf", base % (faf, faf)), ("basica", base % ("", ""))]
    filas = []
    for g in genes:
        datos = None
        for etiqueta, consulta in consultas:
            ruta = os.path.join(carpeta_cache, f"gnomad_{GNOMAD_DATASET}_{g}_{etiqueta}.json")
            if os.path.exists(ruta):
                with open(ruta, encoding="utf-8") as fh:
                    datos = json.load(fh)
                log.info(f"gnomAD: {g} leído de caché ({etiqueta})")
                break
            try:  # un error de esquema llega como HTTP 400; en ese caso se pasa a la consulta básica
                r = _post(GNOMAD_API, reintentos=2 if etiqueta == "faf" else 5, json={"query": consulta, "variables": {"g": g}})
                datos = r.json()
            except Exception as e:  # noqa: BLE001
                log.warning(f"gnomAD: la consulta '{etiqueta}' falló para {g}: {e}")
                datos = None
                continue
            time.sleep(6)  # límite de peticiones de la API de gnomAD
            if datos.get("errors"):
                log.warning(f"gnomAD: la consulta '{etiqueta}' devolvió errores para {g}: {datos['errors']}")
                datos = None
                continue
            with open(ruta, "w", encoding="utf-8") as fh:
                json.dump(datos, fh)
            log.info(f"gnomAD: {g} descargado ({etiqueta})")
            break
        if datos is None:
            raise RuntimeError(f"gnomAD: ninguna consulta funcionó para {g}")
        for v in datos["data"]["gene"]["variants"]:
            ac = an = ac_amr = an_amr = 0
            fafs = []
            grupos = {g_: [0, 0] for g_ in GRUPOS_FAF}  # AC y AN conjuntos (exomas + genomas) por grupo
            for parte in ("exome", "genome"):
                p = v.get(parte)
                if p:
                    ac += p.get("ac") or 0
                    an += p.get("an") or 0
                    f95 = (p.get("faf95") or {}).get("popmax")
                    if f95 is not None:
                        fafs.append(float(f95))
                    for pop in p.get("populations") or []:
                        if pop.get("id") in grupos:
                            grupos[pop["id"]][0] += pop.get("ac") or 0
                            grupos[pop["id"]][1] += pop.get("an") or 0
                        if pop.get("id") == "amr":
                            ac_amr += pop.get("ac") or 0
                            an_amr += pop.get("an") or 0
            c, pos, ref, alt = v["variant_id"].split("-")
            filas.append({"Consulta": f"chr{c}:{pos}:{ref}>{alt}", "gnomad_AC": ac, "gnomad_AN": an,
                          "gnomad_AF": ac / an if an else 0.0, "gnomad_AC_amr": ac_amr, "gnomad_AN_amr": an_amr,
                          "gnomad_faf95_api": max(fafs) if fafs else np.nan,
                          **{f"_ac_{k}": x[0] for k, x in grupos.items()}, **{f"_an_{k}": x[1] for k, x in grupos.items()}})
    df = pd.DataFrame(filas).drop_duplicates("Consulta").set_index("Consulta")
    # La API de gnomAD v4 devuelve faf95 = null en la consulta por gen. Se calcula localmente como el
    # límite inferior exacto (Poisson, unilateral 95 %) de la frecuencia por grupo, y se toma el máximo
    # (grpmax), siguiendo Whiffin et al. 2017. Si la API trae el valor, se usa el de la API.
    from scipy.stats import chi2
    calc = []
    for g_ in GRUPOS_FAF:
        ac_, an_ = df.pop(f"_ac_{g_}").to_numpy(float), df.pop(f"_an_{g_}").to_numpy(float)
        with np.errstate(divide="ignore", invalid="ignore"):
            f = np.where((ac_ > 0) & (an_ > 0), chi2.ppf(0.05, 2 * np.maximum(ac_, 1)) / 2 / np.where(an_ > 0, an_, 1), 0.0)
        calc.append(f)
    df["gnomad_faf95_calculada"] = np.max(np.vstack(calc), axis=0)
    df["gnomad_faf95_max"] = df["gnomad_faf95_api"].fillna(df["gnomad_faf95_calculada"])
    df["gnomad_faf95_origen"] = np.where(df["gnomad_faf95_api"].notna(), "API", "calculada")
    log.info(f"gnomAD: {len(df)} variantes con datos en {len(genes)} genes; faf95 de la API en "
             f"{int(df['gnomad_faf95_api'].notna().sum())}, calculada en {int(df['gnomad_faf95_api'].isna().sum())}")
    return df


def obtener_cadd(variantes, url, carpeta_cache):
    """PHRED de CADD por lectura remota (tabix) de las regiones de los genes.

    Solo se lee el intervalo entre la primera y la última variante de cada gen. La caché se
    reutiliza únicamente si contiene todas las variantes pedidas.
    """
    ruta = os.path.join(carpeta_cache, "cadd.csv")
    if os.path.exists(ruta):
        c = pd.read_csv(ruta).set_index("Consulta")
        if set(variantes["Consulta"]) <= set(c.index):
            log.info("CADD: leído de caché")
            return c
    import pysam
    tbx = pysam.TabixFile(url)
    necesarios = set(variantes["Consulta"])
    filas = []
    for (crom, gen), sub in variantes.groupby(["GRCh38Chromosome", "Gen"]):
        ini, fin = int(sub["Posición GRCh38"].min()), int(sub["Posición GRCh38"].max())
        n = 0
        for linea in tbx.fetch(str(crom), ini - 1, fin):
            c, pos, ref, alt, raw, phred = linea.split("\t")[:6]
            q = f"chr{c}:{pos}:{ref}>{alt}"
            if q in necesarios:
                filas.append({"Consulta": q, "CADD_raw": float(raw), "CADD_PHRED": float(phred)})
            n += 1
        log.info(f"CADD: {gen} región {crom}:{ini}-{fin} ({n} líneas leídas)")
    c = pd.DataFrame(filas).drop_duplicates("Consulta")
    c.to_csv(ruta, index=False)
    log.info(f"CADD: {len(c)} de {len(necesarios)} variantes con puntaje")
    return c.set_index("Consulta")


def _revel(nodo):
    """Extrae el puntaje REVEL (máximo si hay varios) de la respuesta de MyVariant.info."""
    valores = []

    def recorrer(x):
        if isinstance(x, list):
            for y in x:
                recorrer(y)
        elif isinstance(x, dict):
            if "revel" in x:
                sc = x["revel"].get("score") if isinstance(x["revel"], dict) else x["revel"]
                for v in (sc if isinstance(sc, list) else [sc]):
                    try:
                        valores.append(float(v))
                    except (TypeError, ValueError):
                        pass
            else:
                for y in x.values():
                    recorrer(y)
    recorrer(nodo)
    return max(valores) if valores else None


def obtener_revel(variantes, carpeta_cache):
    """REVEL (solo missense) mediante MyVariant.info (dbNSFP), consultas por lotes de 1.000."""
    ruta = os.path.join(carpeta_cache, "revel.csv")
    previo = pd.read_csv(ruta) if os.path.exists(ruta) else pd.DataFrame(columns=["Consulta", "REVEL"])
    hechas = set(previo["Consulta"])
    pend = variantes[~variantes["Consulta"].isin(hechas)]
    filas = []
    for k in range(0, len(pend), 1000):
        bloque = pend.iloc[k:k + 1000]
        ids = {f"chr{r['GRCh38Chromosome']}:g.{r['Posición GRCh38']}{r['REF']}>{r['ALT']}": r["Consulta"] for _, r in bloque.iterrows()}
        resp = _post(MYVARIANT_URL, params={"assembly": "hg38"},
                     data={"ids": ",".join(ids), "fields": "dbnsfp.revel.score"}).json()
        for item in resp:
            q = ids.get(item.get("query"))
            if q is None:
                continue
            filas.append({"Consulta": q, "REVEL": _revel(item.get("dbnsfp"))})
        log.info(f"REVEL: {min(k + 1000, len(pend))}/{len(pend)} consultadas")
        time.sleep(1)
    r = pd.concat([previo, pd.DataFrame(filas, columns=["Consulta", "REVEL"])], ignore_index=True).drop_duplicates("Consulta")
    r.to_csv(ruta, index=False)
    r["REVEL"] = pd.to_numeric(r["REVEL"], errors="coerce")
    r = r[r["Consulta"].isin(set(variantes["Consulta"]))]  # la caché puede contener variantes de otras ejecuciones
    log.info(f"REVEL: {r['REVEL'].notna().sum()} de {len(variantes)} variantes missense con puntaje")
    return r.set_index("Consulta")


def datos_externos(todas, a, out):  # noqa: C901
    """Une gnomAD, CADD y REVEL en una tabla por variante.

    Si una fuente falla, la ejecución se detiene, salvo que se use --permitir-fuentes-faltantes.
    """
    cache = a.cache_externos
    ext = pd.DataFrame(index=pd.Index(todas["Consulta"].unique(), name="Consulta"))
    fuentes = {}
    try:
        g = obtener_gnomad(GENES, cache)
        ext = ext.join(g)
        ext[["gnomad_AC", "gnomad_AC_amr"]] = ext[["gnomad_AC", "gnomad_AC_amr"]].fillna(0)
        ext["en_gnomAD"] = ext["gnomad_AC"] > 0
        fuentes["gnomAD"] = f"API GraphQL, dataset {GNOMAD_DATASET}"
    except Exception as e:  # noqa: BLE001
        log.error(f"gnomAD no disponible; se omiten los análisis que dependen de él: {e}")
        _fuente_faltante(a)
    if not a.sin_cadd:
        try:
            ext = ext.join(obtener_cadd(todas, a.cadd_url, cache))
            fuentes["CADD"] = a.cadd_url
        except Exception as e:  # noqa: BLE001
            log.error(f"CADD no disponible; se omite: {e}")
            _fuente_faltante(a)
    try:
        mis = todas[todas["Tipo de variante"] == "Missense"]
        ext = ext.join(obtener_revel(mis, cache))
        fuentes["REVEL"] = "MyVariant.info (dbNSFP), assembly hg38"
    except Exception as e:  # noqa: BLE001
        log.error(f"REVEL no disponible; se omite: {e}")
        _fuente_faltante(a)
    ext.to_csv(os.path.join(out, "07_datos_externos.csv"))
    log.info(f"Datos externos: columnas = {list(ext.columns)}")
    return ext, fuentes


def _fuente_faltante(a):
    """Detiene la ejecución cuando falla una fuente externa, para no producir resultados parciales."""
    if not getattr(a, "permitir_fuentes_faltantes", False):
        raise SystemExit("Falló una fuente externa (ver el log). Repite la ejecución o usa "
                         "--permitir-fuentes-faltantes si de verdad quieres continuar sin ella.")


def anotar_frecuencia(atl):
    """Marca como 'Rara' cada variante con el criterio de entrenamiento del AVI.

    Se usa faf95 < 0,1 % cuando está disponible y AF < 0,1 % en caso contrario. Las variantes
    ausentes de gnomAD se consideran raras.
    """
    atl = atl.copy()
    if "gnomad_AC" not in atl.columns:
        return atl
    atl["gnomad_AC"] = atl["gnomad_AC"].fillna(0)
    atl["gnomad_AF"] = atl["gnomad_AF"].fillna(0.0)
    atl["en_gnomAD"] = atl["gnomad_AC"] > 0
    if "gnomad_faf95_max" in atl.columns and atl["gnomad_faf95_max"].notna().any():
        f = atl["gnomad_faf95_max"].where(atl["en_gnomAD"], 0.0)
        atl["Rara"] = f.fillna(atl["gnomad_AF"]) < UMBRAL_RARA
        atl["Rara_criterio"] = np.where(f.notna(), "faf95", "AF")
        if "gnomad_faf95_origen" in atl.columns:
            atl.loc[atl["en_gnomAD"] & (atl["gnomad_faf95_origen"] == "calculada"), "Rara_criterio"] = "faf95 calculada"
    else:
        atl["Rara"] = atl["gnomad_AF"] < UMBRAL_RARA
        atl["Rara_criterio"] = "AF"
        log.warning("gnomAD sin faf95: 'Rara' se define con AF total (exomas + genomas). Decláralo en Métodos.")
    return atl


# --------------------------------------------------------------------------------------
# 4. Estadística
# --------------------------------------------------------------------------------------
def auc_ic(y, s, rng, B=B_AUC):
    """AUC con IC 95 % por bootstrap de percentiles; omite los remuestreos con una sola clase."""
    from sklearn.metrics import roc_auc_score
    y, s = np.asarray(y), np.asarray(s)
    if len(set(y)) < 2:
        return np.nan, np.nan, np.nan
    a = roc_auc_score(y, s)
    bs = []
    for _ in range(B):
        i = rng.integers(0, len(y), len(y))
        if y[i].min() != y[i].max():
            bs.append(roc_auc_score(y[i], s[i]))
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return a, lo, hi


def rho_ic(x, y, rng, B=B_RHO):
    """Rho de Spearman con IC 95 % por bootstrap de percentiles y valor p."""
    from scipy.stats import spearmanr
    x, y = np.asarray(x, float), np.asarray(y, float)
    r, p = spearmanr(x, y)
    bs = []
    for _ in range(B):
        i = rng.integers(0, len(x), len(x))
        bs.append(spearmanr(x[i], y[i])[0])
    lo, hi = np.nanpercentile(bs, [2.5, 97.5])
    return r, lo, hi, p


def med_iqr(s):
    """Texto 'mediana (P25 a P75)' con dos decimales."""
    return f"{s.median():.2f} ({s.quantile(.25):.2f} a {s.quantile(.75):.2f})" if len(s) else "—"


def analizar(an, atl, cftr2, out):
    """Análisis principal del artículo (tablas A a F).

    A: AVI por gen. B: AVI por tipo de variante. C: sensibilidad y comparación con AlphaMissense.
    D: mecanismo principal. E: genes de DCP según gravedad. F: correlaciones con fenotipos de CFTR2.
    Devuelve el conjunto de análisis con AVI, las variantes de CFTR2 con AVI y las tablas.
    """
    from scipy.stats import mannwhitneyu
    res = {}
    d = an.merge(atl, left_on="Consulta", right_index=True, how="left")
    sin_avi = d["AVI"].isna().sum()
    if sin_avi:
        log.warning(f"Variantes del conjunto de análisis sin AVI: {sin_avi}")
    d = d[d["AVI"].notna()].copy()
    d["y"] = (d["Categoría"] == "P/PP").astype(int)  # 1 = patogénica o probablemente patogénica
    d.to_csv(os.path.join(out, "05_analisis_con_avi.csv"), index=False)
    log.info(f"Análisis: {len(d)} variantes con AVI (P/PP = {d['y'].sum()}, B/PB = {(1 - d['y']).sum()})")

    rng = np.random.default_rng(SEMILLA)
    # Tabla A: por gen
    filas = []
    for g in GENES + ["Total"]:
        s = d if g == "Total" else d[d["Gen"] == g]
        p, b = s[s["y"] == 1]["AVI"], s[s["y"] == 0]["AVI"]
        a, lo, hi = auc_ic(s["y"], s["AVI"], rng)
        filas.append({"Gen": g, "n_P/PP": len(p), "AVI_P/PP_mediana(RIC)": med_iqr(p), "n_B/PB": len(b),
                      "AVI_B/PB_mediana(RIC)": med_iqr(b), "AUC": round(a, 3), "IC95_inf": round(lo, 3), "IC95_sup": round(hi, 3),
                      "p_MannWhitney": mannwhitneyu(p, b).pvalue if len(p) and len(b) else np.nan})
    res["A_AVI_por_gen"] = pd.DataFrame(filas)
    # Tabla B: por tipo de variante; la AUC solo se calcula con al menos 5 variantes por clase
    filas = []
    for t in TIPOS:
        s = d[d["Tipo de variante"] == t]
        p, b = s[s["y"] == 1]["AVI"], s[s["y"] == 0]["AVI"]
        if len(p) >= 5 and len(b) >= 5:
            a, lo, hi = auc_ic(s["y"], s["AVI"], rng)
        else:
            a = lo = hi = np.nan
        filas.append({"Tipo": t, "n_P/PP": len(p), "AVI_P/PP_mediana(RIC)": med_iqr(p), "n_B/PB": len(b),
                      "AVI_B/PB_mediana(RIC)": med_iqr(b), "AUC": round(a, 3) if a == a else "No calculable (<5 en un grupo)",
                      "IC95_inf": round(lo, 3) if lo == lo else "", "IC95_sup": round(hi, 3) if hi == hi else ""})
    res["B_AVI_por_tipo"] = pd.DataFrame(filas)
    # Tabla C: sensibilidad y comparador
    mis = d[(d["Tipo de variante"] == "Missense") & d["AlphaMissense"].notna()]
    e2 = d[d["Estrellas"] >= 2]
    filas = []
    for nombre, s, col in [("AVI, solo missense", mis, "AVI"), ("AlphaMissense, solo missense", mis, "AlphaMissense"),
                           ("AVI, solo variantes con >=2 estrellas", e2, "AVI")]:
        a, lo, hi = auc_ic(s["y"], s[col], rng)
        filas.append({"Análisis": nombre, "n": len(s), "n_P/PP": int(s["y"].sum()), "AUC": round(a, 3), "IC95_inf": round(lo, 3), "IC95_sup": round(hi, 3)})
    res["C_sensibilidad_comparador"] = pd.DataFrame(filas)
    # Tabla D: mecanismo
    ct = pd.crosstab(d["Mecanismo_principal"], d["Categoría"])
    ctp = (pd.crosstab(d["Mecanismo_principal"], d["Categoría"], normalize="columns") * 100).round(1)
    res["D_mecanismo"] = ct.join(ctp.add_suffix("_%")).reset_index()
    # Tabla E: P/PP de genes de DCP con peor frente a mejor función pulmonar
    p = d[d["y"] == 1].assign(Grupo=lambda x: x["Gen"].map(GRUPO_DCP))
    q = p[p["Grupo"].notna()]
    grupos = sorted(q["Grupo"].unique())
    a_, b_ = [q[q["Grupo"] == g]["AVI"] for g in grupos]
    pv = mannwhitneyu(a_, b_).pvalue
    filas = [{"Grupo": g, "n": int((q["Grupo"] == g).sum()), "AVI_mediana(RIC)": med_iqr(q[q["Grupo"] == g]["AVI"])} for g in grupos]
    filas.append({"Grupo": "Prueba de Mann-Whitney (p)", "n": "", "AVI_mediana(RIC)": round(pv, 3)})
    res["E_DCP_gravedad"] = pd.DataFrame(filas)
    res["E2_DCP_composicion_tipo"] = (pd.crosstab(q["Grupo"], q["Tipo de variante"], normalize="index") * 100).round(1).reset_index()

    # Tabla F: correlación del AVI con los fenotipos de CFTR2. Generador propio para que los IC
    # no dependan del número de remuestreos usados en las tablas anteriores.
    rng2 = np.random.default_rng(SEMILLA)
    m = cftr2.merge(atl[["AVI", "AVI_cuantil", "Componente_principal", "Mecanismo_principal"]], left_on="Consulta", right_index=True, how="left")
    m = m[m["AVI"].notna()].copy()
    m.to_csv(os.path.join(out, "06_cftr2_con_avi.csv"), index=False)
    log.info(f"CFTR2: {len(m)} variantes con AVI; clasificación: {m['clasificacion_cftr2'].value_counts().to_dict()}")
    # Missense: cambio de un aminoácido por otro en la nomenclatura proteica (p. ej., p.Gly551Asp)
    mis_c = m[m["proteina"].astype(str).str.match(r"p\.[A-Z][a-z]{2}\d+[A-Z][a-z]{2}$")]
    filas = []
    nombres = {"cloruro_sudor": "Cloruro en sudor", "prop_IP": "Insuficiencia pancreática", "fev1": "FEV1 % predicho"}
    for grupo, s, vars_ in [("Todas", m, ["cloruro_sudor", "prop_IP", "fev1"]), ("Missense", mis_c, ["cloruro_sudor", "prop_IP", "fev1"]),
                            ("FEV1 con n>=10 pacientes", m[m["n_fev1"] >= 10], ["fev1"])]:
        for v in vars_:
            t = s[s[v].notna()]
            r, lo, hi, pval = rho_ic(t["AVI"], t[v], rng2)
            filas.append({"Grupo": grupo, "Variable": nombres[v], "n": len(t), "rho_Spearman": round(r, 2),
                          "IC95_inf": round(lo, 2), "IC95_sup": round(hi, 2), "p": pval})
    res["F_CFTR2_correlaciones"] = pd.DataFrame(filas)
    for k, v in res.items():
        v.to_csv(os.path.join(out, f"tabla_{k}.csv"), index=False)
        log.info(f"Resultado {k}:\n{v.to_string(index=False)}")
    return d, m, res


def _subconjuntos(d):
    """Subconjuntos del análisis adicional; los de gnomAD solo si hay datos de frecuencia."""
    subs = [("Todas", d), ("Missense", d[d["Tipo de variante"] == "Missense"])]
    if "en_gnomAD" in d.columns:
        aus = d[d["en_gnomAD"] == False]  # noqa: E712
        subs.insert(1, ("Ausentes de gnomAD", aus))
        subs.append(("Missense ausentes de gnomAD", aus[aus["Tipo de variante"] == "Missense"]))
    return subs


def _predictores(sub_nombre, d):
    """Predictores a comparar en un subconjunto; AlphaMissense y REVEL solo aplican a missense."""
    preds = ["AVI"]
    if "CADD_PHRED" in d.columns and d["CADD_PHRED"].notna().any():
        preds.append("CADD_PHRED")
    if sub_nombre.startswith("Missense"):
        preds.append("AlphaMissense")
        if "REVEL" in d.columns and d["REVEL"].notna().any():
            preds.append("REVEL")
    return preds


def analisis_adicional(d, vsi_av, out):
    """Análisis adicionales (tablas G a I).

    G: comparación pareada de predictores. H: umbrales de especificidad y su aplicación a las VSI.
    I: representación de las variantes P/PP en la población latinoamericana de gnomAD.
    Devuelve las tablas y los datos necesarios para la Figura 3 y la Figura S1.
    """
    from sklearn.metrics import roc_auc_score
    res, extra = {}, {}
    rng = np.random.default_rng(SEMILLA)
    # G. Comparación pareada de predictores, con y sin variantes presentes en gnomAD
    filas = []
    for nombre, s in _subconjuntos(d):
        preds = _predictores(nombre, d)
        c = s.dropna(subset=preds)
        if c["y"].nunique() < 2 or c["y"].sum() < 5 or (1 - c["y"]).sum() < 5:
            log.warning(f"Comparadores: subconjunto '{nombre}' con muy pocas variantes por clase; se omite")
            continue
        y = c["y"].to_numpy()
        X = {p: c[p].to_numpy(float) for p in preds}
        base = {p: roc_auc_score(y, X[p]) for p in preds}
        # Los mismos índices se aplican a todos los predictores para obtener la diferencia pareada
        boot = {p: [] for p in preds}
        for _ in range(B_AUC):
            i = rng.integers(0, len(y), len(y))
            if y[i].min() == y[i].max():
                continue
            for p in preds:
                boot[p].append(roc_auc_score(y[i], X[p][i]))
        for p in preds:
            lo, hi = np.percentile(boot[p], [2.5, 97.5])
            fila = {"Subconjunto": nombre, "n": len(c), "n_P/PP": int(y.sum()), "n_B/PB": int(len(y) - y.sum()),
                    "Predictor": p, "AUC": round(base[p], 3), "IC95_inf": round(lo, 3), "IC95_sup": round(hi, 3)}
            if p != "AVI":
                dif = np.array(boot["AVI"]) - np.array(boot[p])
                dl, dh = np.percentile(dif, [2.5, 97.5])
                fila.update({"Dif_AUC_AVI_menos_predictor": round(base["AVI"] - base[p], 3),
                             "Dif_IC95_inf": round(dl, 3), "Dif_IC95_sup": round(dh, 3)})
            filas.append(fila)
        extra[nombre] = (y, X)
    res["G_comparadores"] = pd.DataFrame(filas)

    # H. Umbrales con especificidad del 95 % y 99 % (definidos con las variantes B/PB)
    p_, b_ = d[d["y"] == 1]["AVI"], d[d["y"] == 0]["AVI"]
    umbrales = {"Especificidad 95 %": float(np.percentile(b_, 95)), "Especificidad 99 %": float(np.percentile(b_, 99))}
    filas = []
    for nombre, s in _subconjuntos(d):
        for et, u in umbrales.items():
            ps, bs = s[s["y"] == 1]["AVI"], s[s["y"] == 0]["AVI"]
            if len(ps) == 0 or len(bs) == 0:
                continue
            filas.append({"Subconjunto": nombre, "Umbral": et, "AVI_umbral": round(u, 3),
                          "AVI_cuantil_en_umbral(mediana)": round(float(d.loc[(d["AVI"] - u).abs().nsmallest(25).index, "AVI_cuantil"].median()), 4),
                          "Sensibilidad_%": round(100 * (ps >= u).mean(), 1), "Especificidad_%": round(100 * (bs < u).mean(), 1),
                          "n_P/PP": len(ps), "n_B/PB": len(bs)})
    res["H_umbrales"] = pd.DataFrame(filas)
    extra["umbrales"] = umbrales

    # H1b. Umbrales específicos para missense (las VSI son mayoritariamente missense)
    mis = d[d["Tipo de variante"] == "Missense"]
    umbrales_mis = {"Especificidad 95 %": float(np.percentile(mis[mis["y"] == 0]["AVI"], 95)),
                    "Especificidad 99 %": float(np.percentile(mis[mis["y"] == 0]["AVI"], 99))}
    extra["umbrales_missense"] = umbrales_mis
    vm = vsi_av[vsi_av["Tipo de variante"] == "Missense"] if len(vsi_av) else vsi_av
    filas = []
    for et, u in umbrales_mis.items():
        fila = {"Umbral (definido en missense B/PB)": et, "AVI_umbral": round(u, 3),
                "Sensibilidad_missense_%": round(100 * (mis[mis["y"] == 1]["AVI"] >= u).mean(), 1),
                "Especificidad_missense_%": round(100 * (mis[mis["y"] == 0]["AVI"] < u).mean(), 1),
                "VSI_missense": len(vm), "VSI_missense_AVI>=umbral": int((vm["AVI"] >= u).sum()) if len(vm) else 0}
        fila["%_VSI_missense"] = round(100 * fila["VSI_missense_AVI>=umbral"] / len(vm), 1) if len(vm) else np.nan
        if len(vm) and "en_gnomAD" in vm.columns:
            aus = vm[vm["en_gnomAD"] == False]  # noqa: E712
            fila["VSI_missense_ausentes_gnomAD"] = len(aus)
            fila["VSI_missense_ausentes_gnomAD_AVI>=umbral"] = int((aus["AVI"] >= u).sum())
        filas.append(fila)
    res["H4_umbrales_missense_VSI"] = pd.DataFrame(filas)

    # H2. Aplicación a las VSI
    if len(vsi_av):
        filas = []
        for g in GENES + ["Total"]:
            s = vsi_av if g == "Total" else vsi_av[vsi_av["Gen"] == g]
            fila = {"Gen": g, "n_VSI_con_AVI": len(s)}
            for et, u in umbrales.items():
                k = int((s["AVI"] >= u).sum())
                fila[f"VSI_AVI>=umbral_{et}"] = k
                fila[f"%_{et}"] = round(100 * k / len(s), 1) if len(s) else np.nan
            if "en_gnomAD" in s.columns:
                aus = s[s["en_gnomAD"] == False]  # noqa: E712
                fila["VSI_ausentes_gnomAD"] = len(aus)
                fila["VSI_ausentes_gnomAD_AVI>=umbral_95"] = int((aus["AVI"] >= umbrales["Especificidad 95 %"]).sum())
            filas.append(fila)
        res["H2_VSI_por_gen"] = pd.DataFrame(filas)
        alto = vsi_av[vsi_av["AVI"] >= umbrales["Especificidad 95 %"]]
        res["H3_VSI_altas_mecanismo_tipo"] = pd.crosstab(alto["Tipo de variante"], alto["Mecanismo_principal"], margins=True).reset_index()
        alto.sort_values("AVI", ascending=False).to_csv(os.path.join(out, "08_VSI_AVI_alto.csv"), index=False)
        log.info(f"VSI con AVI >= umbral de especificidad 95 %: {len(alto)} de {len(vsi_av)}")
    else:
        log.warning("No hay VSI con AVI: se omite la aplicación a VSI (¿se ejecutó con --sin-consulta?)")

    # I. Representación en la población latinoamericana de gnomAD ('amr')
    if "en_gnomAD" in d.columns:
        p = d[d["y"] == 1]
        filas = []
        for g in GENES + ["Total"]:
            s = p if g == "Total" else p[p["Gen"] == g]
            obs = s[s["en_gnomAD"]]
            amr = s[s["gnomad_AC_amr"] > 0]
            filas.append({"Gen": g, "P/PP": len(s), "P/PP_observadas_en_gnomAD": len(obs),
                          "P/PP_observadas_en_amr": len(amr),
                          "%_de_observadas_que_aparecen_en_amr": round(100 * len(amr) / len(obs), 1) if len(obs) else np.nan,
                          "AN_amr_mediana": float(obs["gnomad_AN_amr"].median()) if len(obs) else np.nan,
                          "AN_total_mediana": float(obs["gnomad_AN"].median()) if len(obs) else np.nan})
        res["I_poblacion_amr"] = pd.DataFrame(filas)
    for k in [x for x in res]:
        res[k].to_csv(os.path.join(out, f"tabla_{k}.csv"), index=False)
        log.info(f"Resultado {k}:\n{res[k].to_string(index=False)}")
    return res, extra


# --------------------------------------------------------------------------------------
# 4b. Análisis de robustez
# --------------------------------------------------------------------------------------
def _midrank(x):
    """Rangos medios (1-based) con empates promediados, usados por la prueba de DeLong."""
    orden = np.argsort(x, kind="mergesort")
    z = x[orden]
    n = len(x)
    t = np.zeros(n)
    i = 0
    while i < n:
        j = i
        while j < n and z[j] == z[i]:
            j += 1
        t[i:j] = 0.5 * (i + j - 1) + 1
        i = j
    salida = np.empty(n)
    salida[orden] = t
    return salida


def delong_pareado(y, s1, s2):
    """Prueba de DeLong para dos AUC correlacionadas (Sun y Xu, 2014). Devuelve AUC1, AUC2, dif, IC95, p."""
    from scipy.stats import norm
    y = np.asarray(y).astype(bool)
    preds = np.vstack([np.asarray(s1, float), np.asarray(s2, float)])
    m, n = int(y.sum()), int((~y).sum())
    pos, neg = preds[:, y], preds[:, ~y]
    tx = np.array([_midrank(p) for p in pos])
    ty = np.array([_midrank(p) for p in neg])
    tz = np.array([_midrank(p) for p in np.hstack([pos, neg])])
    aucs = tz[:, :m].sum(axis=1) / m / n - (m + 1.0) / (2.0 * n)
    v01 = (tz[:, :m] - tx) / n
    v10 = 1.0 - (tz[:, m:] - ty) / m
    s = np.cov(v01) / m + np.cov(v10) / n
    dif = aucs[0] - aucs[1]
    se = np.sqrt(max(s[0, 0] + s[1, 1] - 2 * s[0, 1], 1e-12))
    p = 2 * norm.sf(abs(dif / se))
    return aucs[0], aucs[1], dif, dif - 1.96 * se, dif + 1.96 * se, p


def ap_ic(y, s, rng, B=B_AUC):
    """Área bajo la curva precisión-exhaustividad (PR-AUC) con IC 95 % por bootstrap."""
    from sklearn.metrics import average_precision_score
    y, s = np.asarray(y), np.asarray(s, float)
    a = average_precision_score(y, s)
    bs = []
    for _ in range(B):
        i = rng.integers(0, len(y), len(y))
        if y[i].min() != y[i].max():
            bs.append(average_precision_score(y[i], s[i]))
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return a, lo, hi


def spearman_ponderado(x, y, w):
    """Correlación de Pearson ponderada sobre rangos (ponderación por número de pacientes)."""
    rx, ry = pd.Series(x).rank().to_numpy(), pd.Series(y).rank().to_numpy()
    w = np.asarray(w, float) / np.sum(w)
    mx, my = np.sum(w * rx), np.sum(w * ry)
    return np.sum(w * (rx - mx) * (ry - my)) / np.sqrt(np.sum(w * (rx - mx) ** 2) * np.sum(w * (ry - my) ** 2))


def holm(p):
    """Valores p ajustados por el método de Holm-Bonferroni, en el orden original."""
    p = np.asarray(p, float)
    orden = np.argsort(p)
    k = len(p)
    ajustado = np.empty(k)
    acumulado = 0.0
    for r, i in enumerate(orden):
        acumulado = max(acumulado, min(1.0, (k - r) * p[i]))
        ajustado[i] = acumulado
    return ajustado


def _preds_disponibles(d, missense):
    """Predictores con al menos un valor en d; AlphaMissense y REVEL solo para missense."""
    preds = ["AVI"]
    for p in (["CADD_PHRED", "AlphaMissense", "REVEL"] if missense else ["CADD_PHRED"]):
        if p in d.columns and d[p].notna().any():
            preds.append(p)
    return preds


def robustez(d, m, atl, out):  # noqa: C901
    """Análisis de robustez (tablas J a L).

    J: subconjuntos con menor solapamiento con el entrenamiento del AVI, PR-AUC, DeLong, AUC por
    estrato, bootstrap por clusters, umbrales con validación cruzada y razones de verosimilitud.
    K: comparadores pareados, techo, ponderación e influencia en CFTR2. L: corrección de Holm.
    """
    from scipy.stats import mannwhitneyu, spearmanr
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import StratifiedKFold
    rng = np.random.default_rng(SEMILLA)
    res, pruebas = {}, []  # pruebas: valores p secundarios que se corrigen por Holm en la tabla L
    mis = d[d["Tipo de variante"] == "Missense"]
    rara = d["Rara"] if "Rara" in d.columns else pd.Series(True, index=d.index)

    # J1. Subconjuntos que reducen la dependencia del entrenamiento: AUC, PR-AUC y DeLong (mismas variantes)
    subs = [("Todas", d, False), ("Raras", d[rara], False), ("Missense", mis, True),
            ("Missense raras", mis[rara.loc[mis.index]], True)]
    subs.append(("Missense raras, >=2 estrellas", mis[rara.loc[mis.index] & (mis["Estrellas"] >= 2)], True))
    if "en_gnomAD" in d.columns:
        subs.append(("Ausentes de gnomAD", d[d["en_gnomAD"] == False], False))  # noqa: E712
        subs.append(("Missense ausentes de gnomAD", mis[mis["en_gnomAD"] == False], True))  # noqa: E712
    filas = []
    for nombre, s, es_mis in subs:
        preds = _preds_disponibles(s, es_mis)
        c = s.dropna(subset=preds)
        if c["y"].sum() < 5 or (1 - c["y"]).sum() < 5:
            continue
        y = c["y"].to_numpy()
        for p in preds:
            a, lo, hi = auc_ic(y, c[p], rng)
            ap, alo, ahi = ap_ic(y, c[p], rng)
            fila = {"Subconjunto": nombre, "n": len(c), "n_P/PP": int(y.sum()), "Prevalencia": round(y.mean(), 3),
                    "Predictor": p, "AUC": round(a, 3), "AUC_IC95": f"{lo:.3f}–{hi:.3f}",
                    "PR_AUC": round(ap, 3), "PR_AUC_IC95": f"{alo:.3f}–{ahi:.3f}"}
            if p != "AVI":
                _, _, dif, dlo, dhi, pv = delong_pareado(y, c["AVI"], c[p])
                fila.update({"Dif_AUC_AVI_menos_pred": round(dif, 3), "Dif_IC95_DeLong": f"{dlo:.3f} a {dhi:.3f}", "p_DeLong": pv})
                pruebas.append((f"DeLong AVI vs {p}, {nombre}", pv))
            filas.append(fila)
    res["J1_subconjuntos_PR_DeLong"] = pd.DataFrame(filas)

    # J2. AUC dentro de cada tipo de variante y media ponderada (elimina la separación por consecuencia)
    filas, acum = [], {}
    for t, g in d.groupby("Tipo de variante"):
        nP, nB = int(g["y"].sum()), int((1 - g["y"]).sum())
        if nP < 5 or nB < 5:
            continue
        fila = {"Tipo": t, "n_P/PP": nP, "n_B/PB": nB}
        for p in _preds_disponibles(g, False):
            gg = g.dropna(subset=[p])
            fila[f"AUC_{p}"] = round(roc_auc_score(gg["y"], gg[p]), 3)
            acum.setdefault(p, []).append((fila[f"AUC_{p}"], nP * nB))
        filas.append(fila)
    fila = {"Tipo": "Media ponderada (nP x nB)"}
    for p, v in acum.items():
        fila[f"AUC_{p}"] = round(sum(a * w for a, w in v) / sum(w for _, w in v), 3)
    filas.append(fila)
    res["J2_AUC_intra_estrato"] = pd.DataFrame(filas)

    # J3. Missense por gen
    filas = []
    for g, s in mis.groupby("Gen"):
        if s["y"].sum() >= 5 and (1 - s["y"]).sum() >= 5:
            fila = {"Gen": g, "n_P/PP": int(s["y"].sum()), "n_B/PB": int((1 - s["y"]).sum())}
            for p in _preds_disponibles(s, True):
                ss = s.dropna(subset=[p])
                fila[f"AUC_{p}"] = round(roc_auc_score(ss["y"], ss[p]), 3)
            filas.append(fila)
    res["J3_missense_por_gen"] = pd.DataFrame(filas)

    # J4. Bootstrap por clusters (gen x ventana genómica) para la AUC missense
    c = mis.assign(cl=mis["Gen"] + "_" + (mis["Posición GRCh38"] // VENTANA_CLUSTER).astype(str))
    grupos = {k: v for k, v in c.groupby("cl")}
    claves = np.array(list(grupos))
    bs = []
    for _ in range(B_AUC):
        s = pd.concat([grupos[k] for k in rng.choice(claves, len(claves))])
        if s["y"].nunique() == 2:
            bs.append(roc_auc_score(s["y"], s["AVI"]))
    res["J4_bootstrap_cluster"] = pd.DataFrame([{
        "Análisis": "AVI, missense", "AUC": round(roc_auc_score(mis["y"], mis["AVI"]), 3), "Clusters": len(claves),
        "Ventana_pb": VENTANA_CLUSTER, "IC95_cluster": f"{np.percentile(bs, 2.5):.3f}–{np.percentile(bs, 97.5):.3f}"}])

    # J5. Umbrales con validación cruzada estratificada (5 particiones x 5 repeticiones)
    X, y = mis["AVI"].to_numpy(), mis["y"].to_numpy()
    filas = []
    for espec in (0.95, 0.99):
        sens, esp = [], []
        for r in range(5):
            for tr, te in StratifiedKFold(5, shuffle=True, random_state=SEMILLA + r).split(X, y):
                u = np.quantile(X[tr][y[tr] == 0], espec)
                sens.append((X[te][y[te] == 1] >= u).mean())
                esp.append((X[te][y[te] == 0] < u).mean())
        filas.append({"Especificidad objetivo": espec, "Sensibilidad_validada_%": round(100 * np.mean(sens), 1),
                      "Sensibilidad_DE": round(100 * np.std(sens), 1), "Especificidad_validada_%": round(100 * np.mean(esp), 1)})
    res["J5_umbrales_validacion_cruzada"] = pd.DataFrame(filas)

    # J6. Razones de verosimilitud por tramos del AVI (missense). Exploratorio: etiquetas de ClinVar.
    cortes = [-np.inf, 0.5, 0.8, 1.0, 1.2, 1.4, np.inf]
    t = mis.assign(tramo=pd.cut(mis["AVI"], cortes)).groupby("tramo", observed=False)["y"].agg(["sum", "count"])
    t["B"] = t["count"] - t["sum"]
    t["LR"] = (t["sum"] / mis["y"].sum()) / (t["B"] / (1 - mis["y"]).sum()).replace(0, np.nan)
    res["J6_LR_tramos_missense"] = t.rename(columns={"sum": "P/PP", "count": "n"}).round(2).reset_index().astype({"tramo": str})

    # J7. DCP: comparación peor frente a mejor función pulmonar dentro de cada tipo de variante, y
    # proporción de P/PP con cuantil genómico >= 0,99 (efecto techo del puntaje)
    p = d[(d["y"] == 1) & d["Gen"].isin(GRUPO_DCP)].assign(Grupo=lambda x: np.where(x["Gen"].isin(["CCDC39", "CCDC40", "CCNO"]), "Peor", "Mejor"))
    filas = []
    for t_, g in p.groupby("Tipo de variante"):
        a_, b_ = g[g["Grupo"] == "Peor"]["AVI"], g[g["Grupo"] == "Mejor"]["AVI"]
        if len(a_) >= 5 and len(b_) >= 5:
            pv = mannwhitneyu(a_, b_).pvalue
            filas.append({"Tipo": t_, "n_peor": len(a_), "mediana_peor": round(a_.median(), 2), "n_mejor": len(b_),
                          "mediana_mejor": round(b_.median(), 2), "p_MannWhitney": pv})
            pruebas.append((f"DCP peor vs mejor, {t_}", pv))
    filas.append({"Tipo": "P/PP de todos los genes con cuantil genómico >= 0,99 (techo)",
                  "n_peor": int(d[d["y"] == 1]["AVI_cuantil"].notna().sum()),
                  "mediana_peor": round(100 * (d[d["y"] == 1]["AVI_cuantil"] >= 0.99).mean(), 1)})
    res["J7_DCP_por_estrato"] = pd.DataFrame(filas)

    # K. CFTR2: comparadores en las mismas variantes, techo, ponderación, influencia y etiquetas independientes
    cols = [c for c in ["AlphaMissense", "CADD_PHRED", "REVEL", "en_gnomAD"] if c in atl.columns and c not in m.columns]
    mm = m.merge(atl[cols], left_on="Consulta", right_index=True, how="left")
    salidas = {"cloruro_sudor": "Cloruro en sudor", "prop_IP": "Insuficiencia pancreática", "fev1": "FEV1 % predicho"}
    filas = []
    for nombre, preds in [("AVI vs CADD", ["AVI", "CADD_PHRED"]), ("Missense: AVI vs CADD/AlphaMissense/REVEL",
                                                                    ["AVI", "CADD_PHRED", "AlphaMissense", "REVEL"])]:
        preds = [q for q in preds if q in mm.columns]
        for v, et in salidas.items():
            s = mm.dropna(subset=preds + [v])
            if len(s) < 10:
                continue
            for q in preds:
                r, lo, hi, pv = rho_ic(s[q], s[v], rng)
                filas.append({"Comparación (mismas variantes)": nombre, "Variable": et, "n": len(s), "Predictor": q,
                              "rho": round(r, 2), "IC95": f"{lo:.2f} a {hi:.2f}"})
    res["K1_CFTR2_comparadores_pareados"] = pd.DataFrame(filas)

    # Techo: la correlación AVI-FEV1 se compara con la de un marcador clínico directo (cloruro en sudor)
    filas = []
    s = mm.dropna(subset=["AVI", "fev1", "cloruro_sudor"]).reset_index(drop=True)
    r_avi, r_sud = spearmanr(s["AVI"], s["fev1"])[0], spearmanr(s["cloruro_sudor"], s["fev1"])[0]
    dif = []
    for _ in range(B_RHO):
        i = rng.integers(0, len(s), len(s))
        dif.append(abs(spearmanr(s["AVI"].iloc[i], s["fev1"].iloc[i])[0]) - abs(spearmanr(s["cloruro_sudor"].iloc[i], s["fev1"].iloc[i])[0]))
    filas.append({"Análisis": "Techo: |rho| AVI-FEV1 menos |rho| cloruro-FEV1 (mismas variantes)", "n": len(s),
                  "Valor": round(abs(r_avi) - abs(r_sud), 2), "IC95": f"{np.percentile(dif, 2.5):.2f} a {np.percentile(dif, 97.5):.2f}",
                  "Detalle": f"rho AVI = {r_avi:.2f}; rho cloruro = {r_sud:.2f}"})
    for v, nv in [("fev1", "n_fev1"), ("cloruro_sudor", "n_sudor")]:
        s = mm.dropna(subset=["AVI", v, nv])
        filas.append({"Análisis": f"rho ponderado por n de pacientes, AVI-{salidas[v]}", "n": len(s),
                      "Valor": round(spearman_ponderado(s["AVI"], s[v], s[nv]), 2), "IC95": "",
                      "Detalle": f"sin ponderar = {spearmanr(s['AVI'], s[v])[0]:.2f}; mediana de pacientes = {s[nv].median():.0f}"})
    # Influencia: rho recalculado excluyendo una variante cada vez
    s = mm.dropna(subset=["AVI", "fev1"]).reset_index(drop=True)
    base = spearmanr(s["AVI"], s["fev1"])[0]
    loo = np.array([spearmanr(s["AVI"].drop(i), s["fev1"].drop(i))[0] for i in s.index])
    top = np.argsort(-np.abs(loo - base))[:5]
    filas.append({"Análisis": "Influencia (dejar una fuera): rango de rho AVI-FEV1", "n": len(s), "Valor": round(base, 2),
                  "IC95": f"{loo.min():.2f} a {loo.max():.2f}",
                  "Detalle": "Más influyentes: " + "; ".join(f"{s.loc[i, 'hgvs']} (rho sin ella = {loo[i]:.2f})" for i in top)})
    # Etiqueta independiente: la clasificación de CFTR2 no procede de ClinVar
    et = mm[mm["clasificacion_cftr2"].isin(["CF-causing", "Not CF-causing"])].dropna(subset=["AVI"])
    if et["clasificacion_cftr2"].nunique() == 2:
        a, lo, hi = auc_ic((et["clasificacion_cftr2"] == "CF-causing").astype(int), et["AVI"], rng)
        filas.append({"Análisis": "Etiqueta independiente de ClinVar: CFTR2 causante vs no causante de FQ (AUC)", "n": len(et),
                      "Valor": round(a, 3), "IC95": f"{lo:.3f}–{hi:.3f}",
                      "Detalle": f"causantes = {(et['clasificacion_cftr2'] == 'CF-causing').sum()}; no causantes = {(et['clasificacion_cftr2'] == 'Not CF-causing').sum()}"})
    res["K2_CFTR2_techo_ponderacion_influencia"] = pd.DataFrame(filas)
    for v, et_ in salidas.items():
        s = mm.dropna(subset=["AVI", v])
        pruebas.append((f"Spearman AVI-{et_} (CFTR2)", spearmanr(s["AVI"], s[v])[1]))

    # L. Corrección de Holm sobre todas las pruebas secundarias
    if pruebas:
        t = pd.DataFrame(pruebas, columns=["Prueba", "p"])
        t["p_Holm"] = holm(t["p"])
        res["L_Holm"] = t
    for k, v in res.items():
        v.to_csv(os.path.join(out, f"tabla_{k}.csv"), index=False)
        log.info(f"Resultado {k}:\n{v.to_string(index=False)}")
    mm.to_csv(os.path.join(out, "06b_cftr2_con_predictores.csv"), index=False)
    return res


def manifiesto(d, m, vsi_av, atl, out):
    """Guarda en manifiesto_n.json el número de variantes con datos de cada fuente.

    Permite detectar ejecuciones incompletas comparando los conteos entre ejecuciones.
    """
    man = {"conjunto_analisis": len(d), "P/PP": int(d["y"].sum()), "B/PB": int((1 - d["y"]).sum()),
           "VSI_con_AVI": len(vsi_av), "CFTR2_con_AVI": len(m)}
    for c in ["AVI", "CADD_PHRED", "AlphaMissense", "REVEL", "gnomad_AF", "gnomad_faf95_max"]:
        if c in d.columns:
            man[f"conjunto_con_{c}"] = int(d[c].notna().sum())
    for c in ["CADD_PHRED", "REVEL"]:
        if c in atl.columns:
            man[f"CFTR2_con_{c}"] = int(atl.reindex(m["Consulta"])[c].notna().sum())
    if "Rara_criterio" in d.columns:
        man["criterio_rara"] = d["Rara_criterio"].value_counts().to_dict()
    with open(os.path.join(out, "manifiesto_n.json"), "w", encoding="utf-8") as fh:
        json.dump(man, fh, ensure_ascii=False, indent=2)
    log.info(f"Manifiesto: {man}")
    return man


# --------------------------------------------------------------------------------------
# 5. Figuras
# --------------------------------------------------------------------------------------
def figuras_adicionales(d, vsi_av, extra, out):
    """Figura 3 (curvas ROC de los predictores) y Figura S1 (distribución del AVI en las VSI)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_auc_score, roc_curve
    nombres = {"en": {"AVI": "AVI", "CADD_PHRED": "CADD", "AlphaMissense": "AlphaMissense", "REVEL": "REVEL",
                      "Todas": "A. All variants", "Ausentes de gnomAD": "B. Variants absent from gnomAD",
                      "Missense": "Missense variants", "Missense ausentes de gnomAD": "C. Missense variants absent from gnomAD",
                      "x": "1 − specificity", "y": "Sensitivity", "h": "AVI score", "dens": "Density",
                      "pp": "Pathogenic/likely pathogenic", "bb": "Benign/likely benign", "vus": "Uncertain significance",
                      "u95": "95% specificity threshold", "u99": "99% specificity threshold",
                      "pA": "A. All variant types", "pB": "B. Missense variants (missense-specific thresholds)"},
               "es": {"AVI": "AVI", "CADD_PHRED": "CADD", "AlphaMissense": "AlphaMissense", "REVEL": "REVEL",
                      "Todas": "A. Todas las variantes", "Ausentes de gnomAD": "B. Variantes ausentes de gnomAD",
                      "Missense": "Variantes missense", "Missense ausentes de gnomAD": "C. Missense ausentes de gnomAD",
                      "x": "1 − especificidad", "y": "Sensibilidad", "h": "Puntaje AVI", "dens": "Densidad",
                      "pp": "Patogénica/probablemente patogénica", "bb": "Benigna/probablemente benigna", "vus": "Significado incierto",
                      "u95": "Umbral de especificidad 95 %", "u99": "Umbral de especificidad 99 %",
                      "pA": "A. Todos los tipos de variante", "pB": "B. Variantes missense (umbrales específicos)"}}
    colores = {"AVI": "#C0392B", "CADD_PHRED": "#555555", "AlphaMissense": "#2E86C1", "REVEL": "#27AE60"}
    paneles = [n for n in ["Todas", "Ausentes de gnomAD", "Missense ausentes de gnomAD"] if n in extra]
    for lang, L in nombres.items():
        if paneles:
            fig, ax = plt.subplots(1, len(paneles), figsize=(5.2 * len(paneles), 5))
            ax = np.atleast_1d(ax)
            for a, n in zip(ax, paneles):
                y, X = extra[n]
                for p, x in X.items():
                    fpr, tpr, _ = roc_curve(y, x)
                    a.plot(fpr, tpr, color=colores[p], lw=1.6, label=f"{L[p]} (AUC {roc_auc_score(y, x):.3f})".replace(".", "," if lang == "es" else "."))
                a.plot([0, 1], [0, 1], ls="--", color="#bbbbbb", lw=1)
                a.set_title(f"{L[n]}\n(n = {len(y)}; P/LP = {int(y.sum())})" if lang == "en" else f"{L[n]}\n(n = {len(y)}; P/PP = {int(y.sum())})", fontsize=10)
                a.set_xlabel(L["x"]); a.set_ylabel(L["y"]); a.legend(fontsize=8, loc="lower right", frameon=False)
                a.spines[["top", "right"]].set_visible(False)
            plt.tight_layout(); plt.savefig(os.path.join(out, f"Figura3_ROC_predictores_{lang}.png"), dpi=300); plt.close()
        if len(vsi_av):
            fig, ax2 = plt.subplots(1, 2, figsize=(14, 4.8))
            for k, (a, sub_d, sub_v, umb, tit) in enumerate([
                    (ax2[0], d, vsi_av, extra["umbrales"], L["pA"]),
                    (ax2[1], d[d["Tipo de variante"] == "Missense"], vsi_av[vsi_av["Tipo de variante"] == "Missense"], extra["umbrales_missense"], L["pB"])]):
                bins = np.linspace(min(sub_d["AVI"].min(), sub_v["AVI"].min()), max(sub_d["AVI"].max(), sub_v["AVI"].max()), 60)
                a.hist(sub_d[sub_d["y"] == 0]["AVI"], bins=bins, density=True, histtype="step", lw=1.6, color="#2E86C1", label=L["bb"])
                a.hist(sub_d[sub_d["y"] == 1]["AVI"], bins=bins, density=True, histtype="step", lw=1.6, color="#C0392B", label=L["pp"])
                a.hist(sub_v["AVI"], bins=bins, density=True, histtype="stepfilled", alpha=.35, color="#E6A23C", label=L["vus"])
                a.axvline(umb["Especificidad 95 %"], color="#333", ls="--", lw=1, label=L["u95"])
                a.axvline(umb["Especificidad 99 %"], color="#333", ls=":", lw=1.2, label=L["u99"])
                a.set_title(f"{tit} (n VUS = {len(sub_v)})" if lang == "en" else f"{tit} (n VSI = {len(sub_v)})", fontsize=10)
                a.set_xlabel(L["h"]); a.set_ylabel(L["dens"]); a.spines[["top", "right"]].set_visible(False)
                if k == 0:
                    a.legend(fontsize=8, frameon=False)
            plt.tight_layout(); plt.savefig(os.path.join(out, f"FiguraS1_VSI_{lang}.png"), dpi=300); plt.close()
    log.info("Figura 3 y Figura suplementaria e1 guardadas en inglés y en español")


def figuras(d, m, res, out):
    """Figura 1 (AVI por gen y tipo), Figura 2 (AVI frente a fenotipos de CFTR2) y Figura 4 (modelo conceptual)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch, Patch
    T = {"en": dict(a="A. AVI by gene and ClinVar class", b="B. AVI by variant type (groups with n ≥ 5)", y="AVI score",
                    p="Pathogenic/likely pathogenic", bn="Benign/likely benign",
                    tipos=["Missense", "Intronic/UTR", "Synonymous", "Nonsense", "Splice site"],
                    f2=[("cloruro_sudor", "Mean sweat chloride (mmol/L)", "A. Sweat chloride"),
                        ("prop_IP", "Pancreatic insufficiency (proportion)", "B. Pancreatic insufficiency"),
                        ("fev1", "Mean FEV1 (% predicted)", "C. Lung function")],
                    clas={"CF-causing": "CF-causing", "Indeterminate": "Indeterminate", "Not CF-causing": "Not CF-causing"},
                    rho="Spearman ρ = {r} (95% CI {l} to {h}), n = {n}", dec="."),
         "es": dict(a="A. AVI por gen y clasificación de ClinVar", b="B. AVI por tipo de variante (grupos con n ≥ 5)", y="Puntaje AVI",
                    p="Patogénica/probablemente patogénica", bn="Benigna/probablemente benigna",
                    tipos=["Missense", "Intrónica/UTR", "Sinónima", "Sin sentido", "Sitio de splicing"],
                    f2=[("cloruro_sudor", "Cloruro en sudor medio (mmol/L)", "A. Cloruro en sudor"),
                        ("prop_IP", "Insuficiencia pancreática (proporción)", "B. Insuficiencia pancreática"),
                        ("fev1", "FEV1 medio (% del predicho)", "C. Función pulmonar")],
                    clas={"CF-causing": "Causante de FQ", "Indeterminate": "Indeterminada", "Not CF-causing": "No causante de FQ"},
                    rho="ρ de Spearman = {r} (IC 95 %: {l} a {h}), n = {n}", dec=",")}
    colP, colB = "#C0392B", "#2E86C1"
    corr = res["F_CFTR2_correlaciones"]
    for lang, L in T.items():
        fmt = lambda v: f"{v:.2f}".replace(".", L["dec"])  # noqa: E731
        # Figura 1: diagramas de caja del AVI, P/PP y B/PB lado a lado por gen y por tipo de variante
        fig, ax = plt.subplots(1, 2, figsize=(13, 5), gridspec_kw={"width_ratios": [1.4, 1]})
        datos, pos, cols = [], [], []
        for i, g in enumerate(GENES):
            for j, (y, c) in enumerate([(1, colP), (0, colB)]):
                datos.append(d[(d["Gen"] == g) & (d["y"] == y)]["AVI"]); pos.append(i * 3 + j); cols.append(c)
        bp = ax[0].boxplot(datos, positions=pos, widths=0.8, patch_artist=True, showfliers=False)
        for b_, c in zip(bp["boxes"], cols):
            b_.set_facecolor(c); b_.set_alpha(.75)
        ax[0].set_xticks([i * 3 + .5 for i in range(len(GENES))]); ax[0].set_xticklabels(GENES)
        ax[0].set_ylabel(L["y"]); ax[0].set_title(L["a"])
        datos, pos, cols, etiq = [], [], [], []
        for i, t in enumerate(TIPOS):
            for j, (y, c) in enumerate([(1, colP), (0, colB)]):
                s = d[(d["Tipo de variante"] == t) & (d["y"] == y)]["AVI"]
                if len(s) >= 5:
                    datos.append(s); pos.append(i * 3 + j); cols.append(c)
        bp = ax[1].boxplot(datos, positions=pos, widths=0.8, patch_artist=True, showfliers=False)
        for b_, c in zip(bp["boxes"], cols):
            b_.set_facecolor(c); b_.set_alpha(.75)
        ax[1].set_xticks([i * 3 + .5 for i in range(len(TIPOS))]); ax[1].set_xticklabels(L["tipos"], rotation=20)
        ax[1].set_title(L["b"])
        ax[0].legend(handles=[Patch(color=colP, alpha=.75, label=L["p"]), Patch(color=colB, alpha=.75, label=L["bn"])], loc="upper right", fontsize=8)
        for a in ax:
            a.axhline(0, color="grey", lw=.5); a.spines[["top", "right"]].set_visible(False)
        plt.tight_layout(); plt.savefig(os.path.join(out, f"Figura1_AVI_{lang}.png"), dpi=300); plt.close()
        # Figura 2: AVI frente a cloruro en sudor, insuficiencia pancreática y FEV1, coloreado por clasificación CFTR2
        colc = {"CF-causing": "#C0392B", "Indeterminate": "#E6A23C", "Not CF-causing": "#2E86C1"}
        fig, ax = plt.subplots(1, 3, figsize=(15, 4.8))
        for a, (v, yl, tit) in zip(ax, L["f2"]):
            for k, c in colc.items():
                s = m[m["clasificacion_cftr2"] == k]
                a.scatter(s["AVI"], s[v], s=22, c=c, alpha=.8, label=L["clas"][k], edgecolors="none")
            r = corr[(corr["Grupo"] == "Todas") & (corr["Variable"] == {"cloruro_sudor": "Cloruro en sudor", "prop_IP": "Insuficiencia pancreática", "fev1": "FEV1 % predicho"}[v])].iloc[0]
            a.set_title(tit + "\n" + L["rho"].format(r=fmt(r["rho_Spearman"]), l=fmt(r["IC95_inf"]), h=fmt(r["IC95_sup"]), n=r["n"]), fontsize=10)
            a.set_xlabel(L["y"]); a.set_ylabel(yl); a.spines[["top", "right"]].set_visible(False)
        ax[0].legend(fontsize=8, frameon=False)
        plt.tight_layout(); plt.savefig(os.path.join(out, f"Figura2_AVI_CFTR2_{lang}.png"), dpi=300); plt.close()
    # Figura 4: modelo conceptual. La IA genómica es un elemento exploratorio de priorización de variantes,
    # no de pronóstico. Las coordenadas de cajas y flechas están en un lienzo de 100 x 90 unidades.
    textos_f4 = {'en':dict(b1=('Suspected or diagnosed PCD or CF',['Clinical assessment and genetic testing']),
                 b2=('AlphaGenome Atlas prediction',['Raw AVI score and main mechanism (research use)']),
                 sd=('Variant not in Atlas',['Unobserved indels or structural','variants: registry and clinical','data only']),
                 b3=('Variant prioritization (exploratory)',['AVI as one input alongside ACMG/AMP evidence,','CFTR2 data and published cohorts; never used','alone (a low score does not rule out pathogenicity)']),
                 b4=('Multidisciplinary diagnostic confirmation',['Full genotype (all variants, both alleles), clinical','findings, specialized tests, genetic counseling']),
                 pt=('Clinical care and physical therapy',['Start at presentation; do not wait for genetic results','Guided by guidelines, clinical status and monitoring:','lung function, symptoms and exacerbations']),
                 no=('Not supported by current evidence',['Using AVI to predict lung','trajectory or to set physical','therapy intensity']),
                 par='In parallel, without waiting for genetic results',
                 rq='Research question: does genomic information add prognostic value in longitudinal individual data?'),
       'es':dict(b1=('DCP o FQ sospechada o diagnosticada',['Valoración clínica y estudio genético']),
                 b2=('Predicción de AlphaGenome Atlas',['Puntaje AVI bruto y mecanismo principal (investigación)']),
                 sd=('Variante no incluida en Atlas',['Indels no observados o variantes','estructurales: solo datos de','registros y clínicos']),
                 b3=('Priorización de la variante (exploratoria)',['El AVI es un elemento más junto con la evidencia','ACMG/AMP, CFTR2 y las cohortes publicadas; nunca se','usa solo (un puntaje bajo no descarta patogenicidad)']),
                 b4=('Confirmación diagnóstica multidisciplinaria',['Genotipo completo (todas las variantes, ambos alelos),','hallazgos clínicos, pruebas especializadas, asesoramiento']),
                 pt=('Atención clínica y fisioterapia',['Inicio en la presentación, sin esperar los resultados genéticos','Guiadas por las guías, el estado clínico y la monitorización:','función pulmonar, síntomas y exacerbaciones']),
                 no=('Sin respaldo en la evidencia actual',['Usar el AVI para predecir la','trayectoria pulmonar o fijar la','intensidad de la fisioterapia']),
                 par='En paralelo, sin esperar los resultados genéticos',
                 rq='Pregunta de investigación: ¿aporta la información genómica valor pronóstico en datos individuales longitudinales?')}
    for lang, T4 in textos_f4.items():
        fig,ax=plt.subplots(figsize=(10,9));ax.set_xlim(0,100);ax.set_ylim(0,90);ax.axis('off')
        G='#6b6b6b';A='#2E6FB7';Rd='#B03A2E'
        def box(x,y,w,h,tl,hl=False,dash=False,fs=9):
            t,lines=tl
            ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.4,rounding_size=1.2',fc='#e6eef8' if hl else ('#fbeeee' if dash else 'white'),
                         ec=A if hl else (Rd if dash else G),lw=2 if hl else 1.4,ls='--' if dash else '-'))
            top=y+h/2+(3.4+(len(lines)-1)*3.0)/2
            ax.text(x+w/2,top,t,ha='center',va='center',fontsize=10 if dash or len(t)>34 else 10.5,fontweight='bold',color=Rd if dash else 'black')
            for i,l in enumerate(lines): ax.text(x+w/2,top-3.4-i*3.0,l,ha='center',va='center',fontsize=fs,color='#333')
        def arr(x1,y1,x2,y2): ax.annotate('',xy=(x2,y2),xytext=(x1,y1),arrowprops=dict(arrowstyle='-|>',color=G,lw=1.2))
        C=37
        box(12,79,50,8,T4['b1']); box(12,66,50,8,T4['b2']); box(67,62,32,13,T4['sd'])
        box(12,47,50,14,T4['b3'],fs=8.6); box(12,30,50,11,T4['b4'],fs=8.6)
        box(1,7,57,16,T4['pt'],hl=True,fs=8.8); box(60,7,39,16,T4['no'],dash=True)
        arr(C,78.6,C,74.4); arr(C,65.6,C,61.4); arr(62.4,70,66.6,70); arr(C,46.6,C,41.4)
        ax.plot([83.5,83.5],[61.6,35.5],color=G,lw=1.2); arr(83.5,35.5,62.4,35.5)
        arr(C,29.6,C,23.4)
        ax.plot([11.6,5],[83,83],color=A,lw=1.6); arr(5,83,5,23.4)
        ax.annotate('',xy=(5,23.4),xytext=(5,83),arrowprops=dict(arrowstyle='-|>',color=A,lw=1.6))
        ax.text(3.0,53,T4['par'],rotation=90,ha='center',va='center',fontsize=8.5,color=A,style='italic')
        ax.add_patch(FancyBboxPatch((1,0.8),98,3.8,boxstyle='round,pad=0.3,rounding_size=1',fc='#efefef',ec=G,lw=1))
        ax.text(50,2.7,T4['rq'],ha='center',va='center',fontsize=9)
        plt.savefig(os.path.join(out, f"Figura4_modelo_{lang}.png"), dpi=300, bbox_inches="tight"); plt.close()
    log.info("Figuras 1, 2 y 4 guardadas en inglés y en español")


# --------------------------------------------------------------------------------------
# 6. Excel y resumen
# --------------------------------------------------------------------------------------
def guardar_excel(d, m, res, out):
    """Guarda todas las tablas, las variantes analizadas y las de CFTR2 en resultados_completos.xlsx."""
    ruta = os.path.join(out, "resultados_completos.xlsx")
    with pd.ExcelWriter(ruta) as w:
        pd.DataFrame({"Contenido": [
            "Resultados generados por pipeline_atlas.py",
            "Hojas A-F: tablas de resultados del artículo. Hoja 'Variantes': conjunto de análisis con AVI. Hoja 'CFTR2': variantes de Sosnay 2013 con AVI.",
            "AUC con IC 95 % por bootstrap (1.000 remuestreos); rho de Spearman con IC 95 % por bootstrap (2.000 remuestreos); semilla = 1.",
            "Valores numéricos calculados en Python (no son fórmulas de Excel)."]}).to_excel(w, sheet_name="Léeme", index=False)
        for k, v in res.items():
            v.to_excel(w, sheet_name=k[:31], index=False)  # Excel limita el nombre de hoja a 31 caracteres
        cols = ["Gen", "Name", "Categoría", "Germline classification", "Germline review status", "Estrellas", "Tipo de variante",
                "Molecular consequence", "Consulta", "VariationID", "AVI", "AVI_cuantil", "Componente_principal", "Mecanismo_principal", "AlphaMissense",
                "CADD_PHRED", "REVEL", "en_gnomAD", "gnomad_AF", "gnomad_AC_amr"]
        d[[c for c in cols if c in d.columns]].to_excel(w, sheet_name="Variantes", index=False)
        m.to_excel(w, sheet_name="CFTR2", index=False)
    log.info(f"Excel guardado: {ruta}")


def resumen(d, m, res, out):
    """Escribe RESUMEN.md con los resultados principales y lo muestra en el registro."""
    A, C, E, F = res["A_AVI_por_gen"], res["C_sensibilidad_comparador"], res["E_DCP_gravedad"], res["F_CFTR2_correlaciones"]
    tot = A[A["Gen"] == "Total"].iloc[0]
    f = lambda g, v: F[(F["Grupo"] == g) & (F["Variable"] == v)].iloc[0]  # noqa: E731
    lineas = [
        "# Resumen de resultados", "",
        f"- Variantes analizadas: {len(d)} (P/PP = {int(d['y'].sum())}; B/PB = {int((1 - d['y']).sum())}).",
        f"- AUC total: {tot['AUC']} (IC 95 % {tot['IC95_inf']}–{tot['IC95_sup']}).",
        f"- AUC missense (AVI): {C.iloc[0]['AUC']} ({C.iloc[0]['IC95_inf']}–{C.iloc[0]['IC95_sup']}); AlphaMissense: {C.iloc[1]['AUC']} ({C.iloc[1]['IC95_inf']}–{C.iloc[1]['IC95_sup']}).",
        f"- AUC con >=2 estrellas: {C.iloc[2]['AUC']} ({C.iloc[2]['IC95_inf']}–{C.iloc[2]['IC95_sup']}); n = {C.iloc[2]['n']}.",
        f"- DCP grave vs leve: {E.iloc[0]['Grupo']} {E.iloc[0]['AVI_mediana(RIC)']} (n={E.iloc[0]['n']}) vs {E.iloc[1]['Grupo']} {E.iloc[1]['AVI_mediana(RIC)']} (n={E.iloc[1]['n']}); p = {E.iloc[2]['AVI_mediana(RIC)']}.",
        f"- CFTR2 ({len(m)} variantes): cloruro en sudor rho = {f('Todas', 'Cloruro en sudor')['rho_Spearman']} ({f('Todas', 'Cloruro en sudor')['IC95_inf']} a {f('Todas', 'Cloruro en sudor')['IC95_sup']}); "
        f"insuficiencia pancreática rho = {f('Todas', 'Insuficiencia pancreática')['rho_Spearman']} ({f('Todas', 'Insuficiencia pancreática')['IC95_inf']} a {f('Todas', 'Insuficiencia pancreática')['IC95_sup']}); "
        f"FEV1 rho = {f('Todas', 'FEV1 % predicho')['rho_Spearman']} ({f('Todas', 'FEV1 % predicho')['IC95_inf']} a {f('Todas', 'FEV1 % predicho')['IC95_sup']}).",
    ]
    if "G_comparadores" in res and len(res["G_comparadores"]):
        G_ = res["G_comparadores"]
        for sub in G_["Subconjunto"].unique():
            t = G_[G_["Subconjunto"] == sub]
            lineas.append(f"- Comparadores, {sub} (n = {t['n'].iloc[0]}): " + "; ".join(
                f"{r['Predictor']} AUC {r['AUC']} ({r['IC95_inf']}–{r['IC95_sup']})" for _, r in t.iterrows()))
    if "H_umbrales" in res:
        for _, r in res["H_umbrales"][res["H_umbrales"]["Subconjunto"] == "Todas"].iterrows():
            lineas.append(f"- Umbral {r['Umbral']}: AVI >= {r['AVI_umbral']}; sensibilidad {r['Sensibilidad_%']} %; especificidad {r['Especificidad_%']} %.")
    if "H2_VSI_por_gen" in res:
        t = res["H2_VSI_por_gen"].set_index("Gen").loc["Total"]
        lineas.append(f"- VSI con AVI: {int(t['n_VSI_con_AVI'])}; con AVI >= umbral 95 %: {int(t['VSI_AVI>=umbral_Especificidad 95 %'])} ({t['%_Especificidad 95 %']} %); "
                      f">= umbral 99 %: {int(t['VSI_AVI>=umbral_Especificidad 99 %'])} ({t['%_Especificidad 99 %']} %).")
    if "I_poblacion_amr" in res:
        t = res["I_poblacion_amr"].set_index("Gen").loc["Total"]
        lineas.append(f"- P/PP observadas en gnomAD: {int(t['P/PP_observadas_en_gnomAD'])} de {int(t['P/PP'])}; observadas en población amr: {int(t['P/PP_observadas_en_amr'])} ({t['%_de_observadas_que_aparecen_en_amr']} % de las observadas).")
    if "J1_subconjuntos_PR_DeLong" in res:
        for _, r in res["J1_subconjuntos_PR_DeLong"].iterrows():
            if r["Subconjunto"] in ("Missense raras", "Missense raras, >=2 estrellas"):
                lineas.append(f"- {r['Subconjunto']} (n = {r['n']}): {r['Predictor']} AUC {r['AUC']} ({r['AUC_IC95']}); PR-AUC {r['PR_AUC']} ({r['PR_AUC_IC95']})"
                              + (f"; dif. vs AVI (DeLong) {r['Dif_AUC_AVI_menos_pred']}, p = {r['p_DeLong']:.3g}" if r["Predictor"] != "AVI" else ""))
    if "K2_CFTR2_techo_ponderacion_influencia" in res:
        for _, r in res["K2_CFTR2_techo_ponderacion_influencia"].iterrows():
            lineas.append(f"- CFTR2, {r['Análisis']} (n = {r['n']}): {r['Valor']} {('[' + r['IC95'] + ']') if r['IC95'] else ''} {r['Detalle']}".rstrip())
    with open(os.path.join(out, "RESUMEN.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lineas) + "\n")
    log.info("\n".join(lineas))


# --------------------------------------------------------------------------------------
# Principal
# --------------------------------------------------------------------------------------
def _pip_freeze(out):
    """Guarda la salida de pip freeze para documentar el entorno exacto de la ejecución."""
    import subprocess
    try:
        txt = subprocess.run([sys.executable, "-m", "pip", "freeze"], capture_output=True, text=True, timeout=120).stdout
        with open(os.path.join(out, "entorno_pip_freeze.txt"), "w", encoding="utf-8") as fh:
            fh.write(txt)
    except Exception as e:  # noqa: BLE001
        log.warning(f"No se pudo registrar pip freeze: {e}")


def main():
    """Ejecuta todas las etapas y guarda metadatos.json aunque la ejecución falle."""
    p = argparse.ArgumentParser(description="Pipeline completo ClinVar + CFTR2 + AlphaGenome Atlas")
    p.add_argument("--clinvar", required=True, help="zip o carpeta con los .tsv de ClinVar de los 7 genes")
    p.add_argument("--cftr2", required=True, help="NIHMS512446-supplement-12.xlsx (Sosnay 2013, Tabla Suplementaria 2)")
    p.add_argument("--salida", default=None, help="carpeta de resultados (por defecto resultados_AAAAMMDD_HHMMSS)")
    p.add_argument("--cache", default="cache_atlas.csv", help="archivo de caché de Atlas (se reutiliza entre ejecuciones)")
    p.add_argument("--hilos", type=int, default=4)
    p.add_argument("--lote", type=int, default=100)
    p.add_argument("--sin-consulta", action="store_true", help="no consultar Atlas; usar solo la caché")
    p.add_argument("--prueba", type=int, default=0, help="limitar a N variantes del conjunto de análisis (prueba rápida)")
    p.add_argument("--cache-externos", default="cache_externos", help="carpeta de caché para gnomAD, CADD y REVEL")
    p.add_argument("--cadd-url", default=CADD_URL, help="archivo CADD (tabix) remoto o local")
    p.add_argument("--sin-cadd", action="store_true", help="omitir CADD")
    p.add_argument("--sin-vsi", action="store_true", help="no consultar ni analizar las VSI")
    p.add_argument("--desde-resultados", default=None,
                   help="carpeta con 03_atlas_por_variante.csv y 07_datos_externos.csv: rehace la estadística sin API ni red")
    p.add_argument("--permitir-fuentes-faltantes", action="store_true",
                   help="continuar aunque falle gnomAD, CADD o REVEL (no recomendado para el artículo)")
    a = p.parse_args()

    out = a.salida or datetime.now().strftime("resultados_%Y%m%d_%H%M%S")
    os.makedirs(out, exist_ok=True)
    configurar_logs(out)
    inicio = time.time()
    meta = {"inicio": datetime.now().isoformat(timespec="seconds"), "argumentos": vars(a),
            "python": sys.version.split()[0], "sistema": platform.platform(),
            "paquetes": {k: version(k) for k in ["alphagenome", "pandas", "numpy", "scipy", "scikit-learn", "matplotlib", "openpyxl", "python-dotenv", "requests", "pysam"]},
            "entradas": {}}
    try:
        faltan = [r for r in (a.clinvar, a.cftr2) if not os.path.exists(r)]
        if faltan:
            raise SystemExit(f"No encuentro estos archivos de entrada en {os.getcwd()}: {faltan}. "
                             "Cópialos a esta carpeta o indica la ruta completa en --clinvar / --cftr2.")
        for clave, ruta in [("clinvar", a.clinvar), ("cftr2", a.cftr2)]:
            if os.path.isfile(ruta):
                meta["entradas"][clave] = {"archivo": os.path.basename(ruta), "sha256": sha256(ruta)}
        # Solo nombres de archivo, sin rutas absolutas, para no exponer directorios personales
        meta["argumentos"] = {k: (os.path.basename(v) if isinstance(v, str) and os.sep in v else v) for k, v in vars(a).items()}
        _pip_freeze(out)
        log.info(f"Inicio. Carpeta de salida: {os.path.basename(os.path.normpath(out))}")
        log.info(f"Versiones: {meta['paquetes']}")

        paso("1. ClinVar")
        snv, an, vsi = preparar_clinvar(leer_clinvar(a.clinvar), out)
        if a.sin_vsi:
            vsi = vsi.iloc[0:0]
        if a.prueba:
            # Muestra equilibrada: el mismo número de variantes por gen y categoría
            an = an.groupby(["Gen", "Categoría"], group_keys=False).head(max(1, a.prueba // (2 * len(GENES))))
            vsi = vsi.groupby("Gen", group_keys=False).head(max(1, a.prueba // len(GENES)))
            log.warning(f"MODO PRUEBA: análisis limitado a {len(an)} variantes y {len(vsi)} VSI")

        paso("2. CFTR2 (Sosnay 2013)")
        cftr2 = leer_cftr2(a.cftr2, snv, out)

        consultas = list(an["Consulta"]) + list(cftr2.loc[cftr2["es_SNV"] & cftr2["Consulta"].notna(), "Consulta"]) + list(vsi["Consulta"])
        if a.desde_resultados:
            paso("3-4. Reanálisis desde resultados derivados (sin API ni red)")
            atl = pd.read_csv(os.path.join(a.desde_resultados, "03_atlas_por_variante.csv")).set_index("Consulta")
            ext = pd.read_csv(os.path.join(a.desde_resultados, "07_datos_externos.csv")).set_index("Consulta")
            ext = ext.drop(columns=[c for c in ext.columns if c in atl.columns])
            atl = atl[atl.index.isin(set(consultas))].join(ext, how="left")
            meta["fuentes_externas"] = {"desde_resultados": os.path.basename(os.path.normpath(a.desde_resultados))}
            log.info(f"Reanálisis: {len(atl)} variantes con AVI leídas de resultados previos")
        else:
            paso("3. AlphaGenome Atlas")
            largo = consultar_atlas(consultas, a.cache, out, a.hilos, a.lote, a.sin_consulta)
            atl = resumir_atlas(largo[largo["consulta"].isin(set(consultas))], out)

            paso("4. Datos externos públicos (gnomAD, CADD, REVEL)")
            # Incluye las variantes de CFTR2 para comparar los predictores sobre las mismas variantes
            ext, fuentes = datos_externos(snv[snv["Consulta"].isin(set(consultas))].drop_duplicates("Consulta"), a, out)
            meta["fuentes_externas"] = fuentes
            atl = atl.join(ext, how="left")
        atl = anotar_frecuencia(atl)

        paso("5. Análisis estadístico principal")
        d, m, res = analizar(an, atl, cftr2, out)

        paso("6. Análisis adicionales (gnomAD, comparadores, umbrales, VSI, población amr)")
        vsi_av = vsi.merge(atl, left_on="Consulta", right_index=True, how="left")
        vsi_av = vsi_av[vsi_av["AVI"].notna()].copy()
        vsi_av.to_csv(os.path.join(out, "09_VSI_con_AVI.csv"), index=False)
        res2, extra = analisis_adicional(d, vsi_av, out)
        res.update(res2)

        paso("6b. Robustez: PR-AUC, DeLong, estratos, clusters, validación cruzada, CFTR2, Holm")
        res.update(robustez(d, m, atl, out))
        meta["manifiesto"] = manifiesto(d, m, vsi_av, atl, out)

        paso("7. Figuras")
        figuras(d, m, res, out)
        figuras_adicionales(d, vsi_av, extra, out)

        paso("8. Excel y resumen")
        guardar_excel(d, m, res, out)
        resumen(d, m, res, out)
        meta["estado"] = "completado"
    except SystemExit as e:
        meta["estado"] = f"detenido: {e}"
        log.error(str(e))
        raise
    except Exception as e:  # noqa: BLE001
        meta["estado"] = f"error: {type(e).__name__}: {e}"
        log.error("Error no controlado:\n" + traceback.format_exc())
        raise
    finally:
        meta["fin"] = datetime.now().isoformat(timespec="seconds")
        meta["duracion_s"] = round(time.time() - inicio, 1)
        with open(os.path.join(out, "metadatos.json"), "w", encoding="utf-8") as fh:
            json.dump(meta, fh, ensure_ascii=False, indent=2)
        log.info(f"Fin ({meta.get('estado')}). Duración: {meta['duracion_s']} s. Resultados en {os.path.basename(os.path.normpath(out))}")


if __name__ == "__main__":
    main()