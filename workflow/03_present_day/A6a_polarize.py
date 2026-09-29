"""A6a: polarise 1KGP-CHM13 SNPs with the human-Pan ancestor (Anc4) from the 8-way T2T Cactus MAF.

inputs : MAF chunk (hal2maf, ref hs1, targets Anc4, chimp, bonobo, gorilla)
         SNP table (bcftools query): CHROM POS REF ALT AN AC AN_AFR AC_AFR   (POS 1-based)
output : TSV: chrom pos anc der cls cpg pars dac_all n_all dac_afr n_afr
Site filters: 11-bp window around the site gap-free in hs1, Anc4, chimp, bonobo (as A1);
anc base must be REF or ALT. pars = chimp == bonobo == gorilla == anc4 (parsimony support).
cpg = CpG context in Anc4 (C followed by G, or G preceded by C) or in hs1 at the site.
usage: python A6a_polarize.py chunk.maf snps.tsv out.tsv
"""
import sys

import numpy as np
import pandas as pd

W = 5
LUT = np.full(256, 4, dtype=np.uint8)
for i, c in enumerate("ACGT"):
    LUT[ord(c)] = i
    LUT[ord(c.lower())] = i
BASES = "ACGTN"
WEAK = set("AT")
IDS = {"hs1": "hs1", "anc4": "Anc4", "chimp": "GCA_028858775.2", "bonobo": "GCA_029289425.2", "gorilla": "GCA_029281585.2"}
CORE = ["hs1", "anc4", "chimp", "bonobo"]


def win_ok(bad):
    n = len(bad)
    ok = np.zeros(n, bool)
    if n < 2 * W + 1:
        return ok
    cs = np.concatenate(([0], np.cumsum(bad.astype(np.int32))))
    ok[W:n - W] = (cs[2 * W + 1:] - cs[:n - 2 * W]) == 0
    return ok


maf, snpf, out = sys.argv[1:4]
snps = pd.read_csv(snpf, sep="\t", header=None, names=["chrom", "pos", "ref", "alt", "an", "ac", "an_afr", "ac_afr"])
snps = snps[snps.ref.str.len().eq(1) & snps.alt.str.len().eq(1)]
want = dict(zip(snps.pos - 1, range(len(snps))))  # 0-based position -> row
res = {}
blk = []


def process(blk):
    byid = {}
    for s in blk:
        name = s[1]
        gid = ".".join(name.split(".")[:2]) if name.startswith("GCA_") else name.split(".", 1)[0]
        if gid not in byid:
            byid[gid] = s
    rows = {k: byid[v] for k, v in IDS.items() if v in byid}
    if not all(k in rows for k in CORE):
        return
    start = int(rows["hs1"][2])
    raw = {k: np.frombuffer(v[6].encode(), np.uint8) for k, v in rows.items()}
    L = len(raw["hs1"])
    if any(len(v) != L for v in raw.values()):
        return
    code = {k: LUT[v] for k, v in raw.items()}
    hpos = start + np.cumsum(raw["hs1"] != ord("-")) - 1
    bad = np.zeros(L, bool)
    for k in CORE:
        bad |= code[k] == 4
    ok = win_ok(bad)
    for i in np.nonzero(ok)[0]:
        j = want.get(int(hpos[i]))
        if j is None:
            continue
        a = code["anc4"][i]
        cpg = (a == 1 and code["anc4"][i + 1] == 2) or (a == 2 and code["anc4"][i - 1] == 1) or \
              (code["hs1"][i] == 1 and code["hs1"][i + 1] == 2) or (code["hs1"][i] == 2 and code["hs1"][i - 1] == 1)
        pars = all(k in code and code[k][i] == a for k in ("chimp", "bonobo", "gorilla"))
        res[j] = (BASES[a], int(cpg), int(pars))


with open(maf) as fh:
    for line in fh:
        if line.startswith("a"):
            if blk:
                process(blk)
            blk = []
        elif line.startswith("s"):
            blk.append(line.split())
    if blk:
        process(blk)

rows = []
for j, (anc, cpg, pars) in res.items():
    r = snps.iloc[j]
    if anc == r.ref:
        der, dac, dafr = r.alt, r.ac, r.ac_afr
    elif anc == r.alt:
        der, dac, dafr = r.ref, r.an - r.ac, r.an_afr - r.ac_afr
    else:
        continue
    cls = ("W" if anc in WEAK else "S") + ("W" if der in WEAK else "S")
    rows.append((r.chrom, r.pos, anc, der, cls, cpg, pars, dac, r.an, dafr, r.an_afr))
pd.DataFrame(rows, columns=["chrom", "pos", "anc", "der", "cls", "cpg", "pars", "dac_all", "n_all", "dac_afr", "n_afr"]
             ).to_csv(out, sep="\t", index=False)
