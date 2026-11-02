#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# One-time setup for the BioTransformer 3.0 batch Shiny app.
#
# The official BioTransformer 3.0 JAR is ~110 MB, which is over GitHub's 100 MB
# file limit, so it is NOT committed to this repo. This script downloads it from
# the official Wishart-lab Bitbucket repository and prepares two ready-to-run
# JARs next to this script:
#
#   BioTransformer3.0.jar          - default; keeps the online HMDB/PathBank
#                                    annotation but adds a short JDBC connect
#                                    timeout so it never hangs when those
#                                    servers are unreachable.
#   BioTransformer3.0.offline.jar  - skips the online annotation entirely
#                                    (fastest; metabolite prediction only).
#
# The btkb/ and supportfiles/ folders and config.json are already in the repo.
#
# Requirements: git, java (>= 8), python3, and R with the 'shiny', 'DT' and
# 'processx' packages.
# -----------------------------------------------------------------------------
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

JAR_REPO="${BT_JAR_REPO:-https://bitbucket.org/wishartlab/biotransformer3.0jar.git}"
VENDOR_JAR_GLOB="BioTransformer3.0_*.jar"

echo "==> Checking tools"
for t in git java python3; do
  command -v "$t" >/dev/null 2>&1 || { echo "ERROR: '$t' not found on PATH." >&2; exit 1; }
done

# 1. Obtain the raw vendor JAR ----------------------------------------------
# Preferred: reassemble from the split parts committed under vendor/ (works
# with no external network access — only a GitHub checkout is needed).
# Fallback: clone the official Bitbucket repo (needs outbound HTTPS).
raw_jar="$(ls -1 $VENDOR_JAR_GLOB 2>/dev/null | head -n1 || true)"
parts=(vendor/BioTransformer3.0_*.jar.*.part)
if [[ -z "${raw_jar}" && -e "${parts[0]}" ]]; then
  base="$(basename "${parts[0]}")"; base="${base%.*.part}"   # -> BioTransformer3.0_YYYYMMDD.jar
  echo "==> Reassembling $base from $(ls -1 vendor/*.part | wc -l) committed part(s)"
  cat vendor/BioTransformer3.0_*.jar.*.part > "$base"
  sumfile="$(ls -1 vendor/BioTransformer3.0_*.jar.sha256 2>/dev/null | head -n1 || true)"
  if [[ -n "$sumfile" ]] && command -v sha256sum >/dev/null 2>&1; then
    want="$(awk '{print $1}' "$sumfile")"; got="$(sha256sum "$base" | awk '{print $1}')"
    if [[ "$want" != "$got" ]]; then
      echo "ERROR: checksum mismatch after reassembly ($got != $want)" >&2; exit 1
    fi
    echo "   checksum OK"
  fi
  raw_jar="$base"
fi
if [[ -z "${raw_jar}" ]]; then
  echo "==> No vendor parts found; downloading BioTransformer JAR from ${JAR_REPO}"
  tmp="$(mktemp -d)"
  trap 'rm -rf "$tmp"' EXIT
  git clone --depth 1 "$JAR_REPO" "$tmp/btjar"
  raw_jar="$(ls -1 "$tmp"/btjar/BioTransformer3.0*.jar | head -n1)"
  cp "$raw_jar" .
  raw_jar="$(basename "$raw_jar")"
  for d in btkb supportfiles config.json; do
    if [[ ! -e "$d" && -e "$tmp/btjar/$d" ]]; then
      echo "   + copying $d from vendor repo"
      cp -r "$tmp/btjar/$d" .
    fi
  done
fi
echo "==> Vendor JAR: $raw_jar"

# 2. Produce the two patched JARs -------------------------------------------
echo "==> Building BioTransformer3.0.jar (online annotation, fail-fast)"
python3 tools/patch_jar.py "$raw_jar" BioTransformer3.0.jar --mode timeout --connect-ms 2000

echo "==> Building BioTransformer3.0.offline.jar (prediction only)"
python3 tools/patch_jar.py "$raw_jar" BioTransformer3.0.offline.jar --mode offline

echo
echo "Setup complete. Launch the app with:"
echo "    R -e 'shiny::runApp(\"$HERE\", host=\"0.0.0.0\", port=3838)'"
