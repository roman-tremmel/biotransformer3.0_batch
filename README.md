# BioTransformer 3.0 — batch Shiny app

A small R/Shiny front-end for [BioTransformer 3.0](https://biotransformer.ca)
that lets you predict small-molecule metabolism for **many compounds at once**
by pasting a list of SMILES or uploading a CSV.

The official BioTransformer web interface ([biotransformer.ca/new](https://biotransformer.ca/new))
only accepts a single compound at a time. This app wraps the official
BioTransformer 3.0 JAR and runs a whole batch, keeping every input in the
results — including ones that cannot be parsed or predicted — so a single bad
SMILES never aborts the run.

## Features

- **Batch input** — paste SMILES (one per line, optional `name, SMILES`) or
  upload a CSV with a `SMILES` column.
- **Robust** — each molecule runs independently; parse errors, failures and
  timeouts are reported per compound and the batch continues.
- **All BioTransformer modes** — `allHuman`, `cyp450` (with CYP450 mode 1/2/3),
  `phaseII`, `hgut`, `ecbased`, `superbio`, `env`; 1–3 prediction steps.
- **Parallel** — runs several molecules concurrently.
- **Structure previews** — 2D depictions of query and metabolite structures
  (client-side via [smiles-drawer](https://github.com/reymond-group/smilesDrawer),
  bundled in `www/`, no internet needed).
- **Downloads** — full metabolite table and per-compound summary as CSV.

## Requirements

- **Java** 8 or newer on `PATH` (OpenJDK is fine).
- **R** with packages `shiny`, `DT`, `processx` (and base `parallel`).
- **git** + **python3** (used once by `setup.sh` to fetch and prepare the JAR).

```bash
R -e 'install.packages(c("shiny","DT","processx"))'
```

## Setup

The BioTransformer JAR is ~110 MB — over GitHub's 100 MB file limit — so it is
committed **split into parts** under `vendor/` (each < 95 MB). `setup.sh`
reassembles them locally (verifying a SHA-256 checksum) and prepares two
ready-to-run JARs — **no external download is needed**, a plain `git` checkout
is enough:

```bash
bash setup.sh
```

(If the `vendor/*.part` files are ever missing, `setup.sh` falls back to
cloning the official Wishart-lab Bitbucket repository, which needs outbound
HTTPS access.)

This produces, next to `app.R`:

| File | Online annotation | Notes |
|------|-------------------|-------|
| `BioTransformer3.0.jar`         | yes (HMDB / PathBank → PubChem CID & synonyms) | default; adds a short JDBC connect timeout so it never hangs if those servers are unreachable |
| `BioTransformer3.0.offline.jar` | no | fastest; metabolite prediction only, fully self-contained |

The knowledge base (`btkb/`), models (`supportfiles/`) and `config.json` are
already in the repo and must sit next to the JAR (they do).

> **Note on the online annotation.** BioTransformer's JAR connects to two
> Wishart-lab MySQL servers (`162.243.63.130`, `192.241.161.227`) only to add
> PubChem CIDs and synonyms. The metabolism prediction itself is fully local.
> If your server cannot reach those hosts, either tick **Offline mode** in the
> app or let the default JAR fall back automatically (it gives up after ~2 s per
> attempt instead of hanging). See `tools/patch_jar.py` for exactly what the
> two prepared JARs change — nothing else in the vendor JAR is modified.

## Run

```bash
R -e 'shiny::runApp(".", host="0.0.0.0", port=3838)'
```

Then open <http://SERVER:3838>.

### Posit Connect / shinyapps.io

Publish the directory as usual. The app sets a writable `user.home` for the JVM
per run, which avoids the common Posit Connect error where the JNI-InChI native
library cannot unpack into a read-only home directory. Make sure Java is
available in the deployment environment and run `setup.sh` there (or commit the
prepared JARs via Git LFS) so the JAR is present.

## Usage

1. Paste SMILES (or upload a CSV) — click *load example set* to try it.
2. Pick the biotransformer type, CYP450 mode and number of steps.
3. (Optional) tick *Offline mode* for the fastest, fully local run.
4. Click **Run prediction**. Progress is shown per molecule.
5. Browse the **Metabolites** and **Per-compound summary** tabs; download CSVs.

### Status values

| Status | Meaning |
|--------|---------|
| `predicted`      | metabolites were found |
| `no metabolites` | SMILES parsed, but no metabolites predicted |
| `parse error`    | SMILES could not be read |
| `error`          | BioTransformer failed for this molecule (batch continued) |
| `timeout`        | exceeded the per-molecule time limit (batch continued) |

## Project layout

```
app.R                      Shiny UI + server
R/biotransformer.R         batch engine (wraps the JAR, one run per molecule)
tools/patch_jar.py         prepares the timeout / offline JARs from the vendor JAR
setup.sh                   downloads the vendor JAR and runs patch_jar.py
www/smiles-drawer.min.js   client-side structure depiction
btkb/, supportfiles/       BioTransformer knowledge base + models (required)
config.json                BioTransformer configuration (required)
examples/                  example input files
```

## Citation & license

BioTransformer is developed by the Wishart Lab and distributed under GNU LGPL
v3. If you use this app, cite:

> Djoumbou Feunang Y, Fiamoncini J, de la Fuente AG, Manach C, Greiner R,
> Wishart DS. *BioTransformer: A Comprehensive Computational Tool for Small
> Molecule Metabolism Prediction and Metabolite Identification.* Journal of
> Cheminformatics (2019) 11:2. doi:10.1186/s13321-018-0324-5

The environmental microbial module uses EAWAG-BBD data licensed by
[enviPath](https://envipath.com/license/) (CC BY-NC-SA 4.0); commercial use of
that module requires a licence from enviPath.
