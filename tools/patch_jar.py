#!/usr/bin/env python3
"""
Patch the official BioTransformer 3.0 JAR so it never *hangs* on its two
hard-coded, external MySQL annotation servers.

Why this is needed
------------------
BioTransformer's downloadable JAR opens JDBC connections to two Wishart-lab
MySQL servers during prediction, purely to annotate predicted metabolites with
PubChem CIDs / synonyms (HMDB + PathBank):

    jdbc:mysql://162.243.63.130:3306/hmdb
    jdbc:mysql://192.241.161.227/pathbank

The metabolism prediction itself (metabolites, reactions, enzymes, formulas,
masses) is fully local and does **not** need these servers.

The MySQL Connector/J default connect timeout is 0 (= wait forever). On a
machine that *silently drops* outbound packets to those hosts (common on
locked-down servers / containers), the JVM blocks indefinitely and the whole
batch stalls. When the connection fails *fast* instead, BioTransformer catches
the exception and finishes the prediction normally.

This script rewrites only those two JDBC URL string constants inside the three
`biotransformer/dbrelevant/*.class` entries. Nothing else in the vendor JAR is
touched. It is byte-safe: Java class files reference the constant pool by index
(not byte offset), so a CONSTANT_Utf8 entry can be grown as long as its 2-byte
length prefix is updated, which is exactly what this does.

Modes
-----
  timeout  (default) keep the real hosts, but append
           ?connectTimeout=<ms>&socketTimeout=<ms>
           -> annotation still works when the servers are reachable;
              otherwise each attempt gives up after <connect ms> instead of hanging.
  offline  repoint both hosts to 127.0.0.1:1 so the connection is refused
           instantly -> fastest; metabolite prediction only, no online
           PubChem/synonym annotation. Best for fully air-gapped servers or
           very large batches.

Usage
-----
  python3 tools/patch_jar.py INPUT.jar OUTPUT.jar [--mode timeout|offline]
                             [--connect-ms 2000] [--socket-ms 15000]
"""
import argparse, io, os, shutil, struct, sys, tempfile, zipfile

DB_URLS = [
    "jdbc:mysql://162.243.63.130:3306/hmdb",
    "jdbc:mysql://192.241.161.227/pathbank",
]
DB_CLASSES_DIR = "biotransformer/dbrelevant/"


def _replace_utf8(buf: bytearray, old: bytes, new: bytes) -> int:
    """Replace every CONSTANT_Utf8 occurrence of `old` with `new`, fixing the
    2-byte length prefix. Returns the number of replacements."""
    n = 0
    i = 0
    while True:
        j = buf.find(old, i)
        if j < 0:
            break
        lenpos = j - 2
        # Validate this really is a Utf8 constant: tag 0x01, then u2 length == len(old)
        if (lenpos >= 1 and buf[lenpos - 1] == 0x01
                and struct.unpack(">H", bytes(buf[lenpos:lenpos + 2]))[0] == len(old)):
            buf[lenpos:lenpos + 2] = struct.pack(">H", len(new))
            buf[j:j + len(old)] = new
            n += 1
            i = lenpos + 2 + len(new)
        else:
            i = j + len(old)
    return n


def _new_url(old: str, mode: str, connect_ms: int, socket_ms: int) -> str:
    if mode == "offline":
        # keep the scheme + db name, point host:port at a guaranteed-refusing address
        db = old.rsplit("/", 1)[-1]
        return "jdbc:mysql://127.0.0.1:1/%s?connectTimeout=1000" % db
    # timeout mode: keep host, add fast connect + bounded socket timeout
    sep = "&" if "?" in old else "?"
    return "%s%sconnectTimeout=%d&socketTimeout=%d" % (old, sep, connect_ms, socket_ms)


def patch_jar(src: str, dst: str, mode: str, connect_ms: int, socket_ms: int) -> int:
    edits = [(u, _new_url(u, mode, connect_ms, socket_ms)) for u in DB_URLS]
    total = 0
    tmp_fd, tmp_path = tempfile.mkstemp(suffix=".jar", dir=os.path.dirname(os.path.abspath(dst)) or ".")
    os.close(tmp_fd)
    try:
        with zipfile.ZipFile(src, "r") as zin, \
             zipfile.ZipFile(tmp_path, "w", zipfile.ZIP_DEFLATED) as zout:
            for item in zin.infolist():
                data = zin.read(item.filename)
                if item.filename.startswith(DB_CLASSES_DIR) and item.filename.endswith(".class"):
                    buf = bytearray(data)
                    changed = 0
                    for old, new in edits:
                        changed += _replace_utf8(buf, old.encode("utf-8"), new.encode("utf-8"))
                    if changed:
                        data = bytes(buf)
                        total += changed
                # Preserve original compression type where possible
                zi = zipfile.ZipInfo(item.filename, date_time=item.date_time)
                zi.compress_type = item.compress_type
                zi.external_attr = item.external_attr
                zout.writestr(zi, data)
        shutil.move(tmp_path, dst)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
    return total


def main():
    ap = argparse.ArgumentParser(description="Patch BioTransformer JAR DB timeouts.")
    ap.add_argument("src")
    ap.add_argument("dst")
    ap.add_argument("--mode", choices=["timeout", "offline"], default="timeout")
    ap.add_argument("--connect-ms", type=int, default=2000)
    ap.add_argument("--socket-ms", type=int, default=15000)
    a = ap.parse_args()
    if not os.path.isfile(a.src):
        sys.exit("input jar not found: %s" % a.src)
    n = patch_jar(a.src, a.dst, a.mode, a.connect_ms, a.socket_ms)
    if n == 0:
        sys.exit("ERROR: no JDBC URL constants were patched -- is this the "
                 "official BioTransformer3.0 JAR? (%s)" % a.src)
    print("Patched %d JDBC URL constant(s) [mode=%s] -> %s" % (n, a.mode, a.dst))


if __name__ == "__main__":
    main()
