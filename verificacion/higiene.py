#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Chequeo de higiene: que no entren credenciales ni datos personales al repositorio.

Corre antes de cada push (verificacion/hooks/pre-push) y a mano. Sólo biblioteca estándar,
así que anda con cualquier python3, en la laptop y en el clúster. Existe porque el
repositorio es público y una vez pusheado algo, sacarlo cuesta reescribir la historia
([D-47]).

Uso:
    python3 verificacion/higiene.py               # el árbol trackeado (lo que se va a pushear)
    python3 verificacion/higiene.py --historia    # además, todos los blobs de toda la historia
    python3 verificacion/higiene.py --autotest    # prueba que cada detector dispara y que
                                                  #   un texto limpio no dispara

Qué busca, y cómo:

  1. Patrones genéricos (regex, abajo): claves privadas, claves públicas SSH, tokens de
     GitHub / Anthropic / AWS / Slack, contraseñas o passphrases con valor asignado,
     direcciones IPv4 que no sean loopback, y nombres de archivo típicos de credenciales.
  2. Palabras prohibidas propias del proyecto —usuarios del clúster, hostnames con
     dominio— leídas de `verificacion/higiene.local`, un archivo NO trackeado (está en
     .gitignore) con una palabra por línea. Se guarda aparte justamente para que las
     palabras que no deben aparecer en el repo no aparezcan en el repo. Si el archivo no
     existe el chequeo avisa y sigue: la parte genérica corre igual.

Sale con 1 si encuentra algo, y lo imprime enmascarado: nunca escribe el secreto entero
en la terminal ni en un log.
"""

import argparse
import os
import re
import subprocess
import sys

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOCAL = os.path.join(RAIZ, "verificacion", "higiene.local")

PATRONES = [
    ("clave privada", re.compile(r"-----BEGIN [A-Z ]*PRIVATE KEY-----")),
    ("clave pública SSH", re.compile(r"\bssh-(?:rsa|ed25519|dss|ecdsa-sha2-nistp\d+) AAAA[0-9A-Za-z+/]{20,}")),
    ("token de GitHub", re.compile(r"\b(?:ghp|gho|ghu|ghs|ghr)_[A-Za-z0-9]{30,}\b|\bgithub_pat_[A-Za-z0-9_]{20,}\b")),
    ("clave de API de Anthropic", re.compile(r"\bsk-ant-[A-Za-z0-9_-]{20,}")),
    ("clave de AWS", re.compile(r"\bAKIA[0-9A-Z]{16}\b")),
    ("token de Slack", re.compile(r"\bxox[baprs]-[A-Za-z0-9-]{10,}")),
    ("contraseña con valor", re.compile(
        r"(?i)\b(?:password|passwd|contraseña|passphrase|clave_secreta)\s*[:=]\s*['\"]?[^\s'\"]{3,}")),
    ("línea de .netrc", re.compile(r"(?m)^\s*machine\s+\S+\s+login\s+\S+\s+password\s+\S+")),
    ("dirección IPv4", re.compile(r"(?<![\d.])(?!127\.|0\.0\.0\.0)(?:\d{1,3}\.){3}\d{1,3}(?![\d.])")),
]

ARCHIVOS_SOSPECHOSOS = re.compile(
    r"(?:^|/)(?:id_rsa|id_ed25519|id_ecdsa|id_dsa)(?:\.pub)?$|\.(?:pem|key|p12|pfx)$|(?:^|/)(?:\.netrc|\.env|credentials)$")

# Números de versión con cuatro campos parecen IPs sin serlo: se descartan cuando lo que
# precede al número es una palabra de versión o el nombre de un módulo.
FALSA_IP = re.compile(r"(?i)(?:version|versión|fftw|python|openmpi|gnu|cuda|rocm|v)\s*/?\s*$")


def es_binario(datos):
    return b"\x00" in datos[:8000]


def cargar_prohibidas():
    if not os.path.exists(LOCAL):
        return None
    with open(LOCAL, encoding="utf-8") as f:
        return [l.strip() for l in f if l.strip() and not l.startswith("#")]


def enmascarar(s):
    s = s.strip()
    return s if len(s) <= 12 else s[:4] + "…" + s[-3:]


def revisar_texto(texto, nombre, prohibidas):
    """Devuelve una lista de (nombre, línea, detector, muestra enmascarada)."""
    hallazgos = []
    for i, linea in enumerate(texto.splitlines(), 1):
        for etiqueta, rx in PATRONES:
            m = rx.search(linea)
            if not m:
                continue
            if etiqueta == "dirección IPv4":
                antes = linea[:m.start()]
                if FALSA_IP.search(antes) or "doi.org" in linea:
                    continue
            hallazgos.append((nombre, i, etiqueta, enmascarar(m.group(0))))
        if prohibidas:
            for p in prohibidas:
                if p in linea:
                    hallazgos.append((nombre, i, "palabra prohibida", enmascarar(p)))
    return hallazgos


def revisar_arbol(prohibidas):
    salida = subprocess.run(["git", "ls-files", "-z"], cwd=RAIZ, capture_output=True, check=True)
    rutas = [r for r in salida.stdout.decode().split("\0") if r]
    hallazgos = []
    for ruta in rutas:
        if ARCHIVOS_SOSPECHOSOS.search(ruta):
            hallazgos.append((ruta, 0, "nombre de archivo de credenciales", ruta))
        completa = os.path.join(RAIZ, ruta)
        if not os.path.isfile(completa):
            continue
        with open(completa, "rb") as f:
            datos = f.read()
        if es_binario(datos):
            continue
        hallazgos += revisar_texto(datos.decode("utf-8", "replace"), ruta, prohibidas)
    return hallazgos, len(rutas)


def revisar_historia(prohibidas):
    """Todos los blobs alcanzables desde cualquier ref, una sola vez cada uno."""
    objetos = subprocess.run(["git", "rev-list", "--all", "--objects"], cwd=RAIZ,
                             capture_output=True, check=True).stdout.decode().splitlines()
    vistos, hallazgos, n = set(), [], 0
    for linea in objetos:
        partes = linea.split(" ", 1)
        if len(partes) < 2:
            continue                       # commits y árboles no traen ruta
        sha, ruta = partes
        if sha in vistos:
            continue
        vistos.add(sha)
        tipo = subprocess.run(["git", "cat-file", "-t", sha], cwd=RAIZ,
                              capture_output=True).stdout.decode().strip()
        if tipo != "blob":
            continue
        n += 1
        datos = subprocess.run(["git", "cat-file", "-p", sha], cwd=RAIZ, capture_output=True).stdout
        if es_binario(datos):
            continue
        hallazgos += revisar_texto(datos.decode("utf-8", "replace"),
                                   "historia:%s:%s" % (sha[:7], ruta), prohibidas)
    return hallazgos, n


def autotest():
    """Cada detector tiene que disparar sobre una muestra sintética, y ninguno sobre un
    texto limpio. Es la prueba de que el chequeo es capaz de fallar."""
    # Las muestras se arman concatenando pedazos: si estuvieran escritas enteras acá, el
    # propio escáner las encontraría en este archivo y bloquearía el push (pasó).
    j = "".join
    muestras = {
        "clave privada": j(["-----BEGIN ", "OPENSSH PRIVATE", " KEY-----"]),
        "clave pública SSH": j(["ssh-ed25519", " AAAAC3NzaC1lZDI1NTE5", "AAAAIGb0000000000000000000000000000 x@y"]),
        "token de GitHub": j(["token = gh", "p_", "a" * 36]),
        "clave de API de Anthropic": j(["sk-", "ant-api03-", "x" * 30]),
        "clave de AWS": j(["AK", "IA", "Q" * 16]),
        "token de Slack": j(["xox", "b-123456789-abcdefghij"]),
        "contraseña con valor": j(["pass", "word = hunter22"]),
        "línea de .netrc": j(["machine github.com login yo ", "password s3creto"]),
        "dirección IPv4": j(["ssh 192.168", ".4.17"]),
    }
    limpio = ("La puerta de la Fase 2 dio 6/6 en el trabajo 8168. ulimit -v 2500000.\n"
              "module load fftw/3.3.11 python/3.13.13. Sin passphrase para no tipearla.\n"
              "fetch('http://127.0.0.1:8757/inject'). doi 10.1016/j.cpc.2020.107482\n"
              "IdentityFile ~/.ssh/cluster\n")
    fallas = 0
    for etiqueta, texto in muestras.items():
        h = [x for x in revisar_texto(texto, "muestra", ["palabra-prohibida-de-prueba"]) if x[2] == etiqueta]
        print("  %-28s %s" % (etiqueta, "dispara" if h else "NO DISPARA"))
        fallas += not h
    h = revisar_texto("usuario palabra-prohibida-de-prueba en una ruta", "muestra", ["palabra-prohibida-de-prueba"])
    print("  %-28s %s" % ("palabra prohibida", "dispara" if h else "NO DISPARA")); fallas += not h
    h = revisar_texto(limpio, "limpio", ["palabra-prohibida-de-prueba"])
    print("  %-28s %s" % ("texto limpio", "no dispara" if not h else "DISPARA: %s" % h)); fallas += bool(h)
    return fallas


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--historia", action="store_true", help="revisar también todos los blobs de la historia")
    ap.add_argument("--autotest", action="store_true", help="probar los detectores sobre muestras sintéticas")
    args = ap.parse_args()

    if args.autotest:
        fallas = autotest()
        print("autotest: %s" % ("OK" if not fallas else "%d falla(s)" % fallas))
        return 1 if fallas else 0

    prohibidas = cargar_prohibidas()
    if prohibidas is None:
        print("[aviso] no existe verificacion/higiene.local: sólo corre la parte genérica.")
        print("        Crealo con una palabra por línea (usuarios del clúster, hosts con dominio).")

    hallazgos, n = revisar_arbol(prohibidas)
    ambito = "%d archivos trackeados" % n
    if args.historia:
        hh, nb = revisar_historia(prohibidas)
        hallazgos += hh
        ambito += " y %d blobs de la historia" % nb

    if hallazgos:
        print("HIGIENE: %d hallazgo(s) en %s" % (len(hallazgos), ambito))
        for nombre, linea, etiqueta, muestra in hallazgos:
            print("  %s:%s  %s  (%s)" % (nombre, linea, etiqueta, muestra))
        return 1
    print("HIGIENE: limpio (%s)" % ambito)
    return 0


if __name__ == "__main__":
    sys.exit(main())
