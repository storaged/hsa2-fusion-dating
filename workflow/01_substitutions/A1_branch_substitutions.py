"""A1: branch-specific substitution catalogue from the T2T-apes 8-way Cactus alignment.

Input : orthology-only MAF chunk (hal2maf --noDupes, reference hs1, with ancestors).
Output: <out>.subs.tsv.gz  one row per substitution on a branch
        <out>.callable.tsv  per 100-kb hs1 window x branch: callable sites and parent W/S counts

Branches (Cactus ancestral states):
  human   Anc4 -> hs1        chimp   Anc5 -> chimp      bonobo  Anc5 -> bonobo
  panstem Anc4 -> Anc5       gorilla Anc3 -> gorilla
Parsimony flag (scheme b): the branch change is also supported by extant outgroups
  human: chimp == bonobo == gorilla != human (orangutan, if aligned, agrees with chimp)
  chimp: bonobo == human == gorilla != chimp;  bonobo: chimp == human == gorilla != bonobo
  panstem: chimp == bonobo != human == gorilla; gorilla: human == chimp == bonobo != gorilla (orang agrees)
Quality filter (after Dreszer et al. 2007 / tytus): 11-bp window centred on the site with
  no gaps/N in any core genome (hs1, Anc4, Anc5, chimp, bonobo) and <=2 parent/child differences
  in the window (site included). Gorilla branch additionally requires gorilla and Anc3 gap-free.
Site annotations: parent CpG context, hs1 soft-mask (repeat) flag, species coordinate.
"""
import gzip
import sys
from collections import defaultdict

import numpy as np

W = 5  # half window
WIN = 100_000
# Roles -> HAL genome ids. "hs1" is the coordinate (MAF reference) genome; "human" is the tip of the
# human branch (defaults to the same genome; can be e.g. an HG002 haplotype in the 16-way HAL).
G = dict(hs1="hs1", human="hs1", anc4="Anc4", anc5="Anc5", anc3="Anc3", chimp="GCA_028858775.2",
         bonobo="GCA_029289425.2", gorilla="GCA_029281585.2", sorang="GCA_028885655.2", borang="GCA_028885625.2")
CORE = ["hs1", "human", "anc4", "anc5", "chimp", "bonobo"]
BR = {  # branch: (parent, child, extra genomes required gap-free)
    "human": ("anc4", "human", []),
    "chimp": ("anc5", "chimp", []),
    "bonobo": ("anc5", "bonobo", []),
    "panstem": ("anc4", "anc5", []),
    "gorilla": ("anc3", "gorilla", ["anc3", "gorilla"]),
}
LUT = np.full(256, 4, dtype=np.uint8)  # A C G T -> 0 1 2 3, else 4 (gap/N)
for i, c in enumerate("ACGT"):
    LUT[ord(c)] = i
    LUT[ord(c.lower())] = i
BASES = np.array(list("ACGTN"))
WEAK = np.array([1, 0, 0, 1, 0], bool)  # A,T weak


def win_sum(x):
    """sum over centred 11-bp windows; edges (incomplete window) -> large."""
    n = len(x)
    c = np.full(n, 99, np.int32)
    if n < 2 * W + 1:
        return c
    cs = np.concatenate(([0], np.cumsum(x.astype(np.int32))))
    c[W:n - W] = cs[2 * W + 1:] - cs[:n - 2 * W]
    return c


def blocks(fh):
    blk = []
    for line in fh:
        if line.startswith("a"):
            if blk:
                yield blk
            blk = []
        elif line.startswith("s"):
            blk.append(line.split())
    if blk:
        yield blk


def main(maf, out):
    op = gzip.open if maf.endswith(".gz") else open
    fout = gzip.open(out + ".subs.tsv.gz", "wt")
    fout.write("chrom\tpos\tbranch\tanc\tder\tcls\tcpg\trepeat\tpars\tsp_chrom\tsp_pos\tsp_len\n")
    callable_ = defaultdict(lambda: np.zeros(3, np.int64))  # (win, branch) -> [n, parentW, parentS]
    # ILS site patterns (A8): per window [n_valid, HP|G, HG|P, PG|H], orangutan as outgroup, Pan = chimp == bonobo
    ils = defaultdict(lambda: np.zeros(4, np.int64))
    with op(maf, "rt") as fh:
        for blk in blocks(fh):
            byid = {}
            for s in blk:
                gid, chrom = s[1].split(".", 1) if not s[1].startswith("GCA_") else (".".join(s[1].split(".")[:2]), ".".join(s[1].split(".")[2:]))
                if gid not in byid:
                    byid[gid] = (chrom, int(s[2]), int(s[3]), s[4], int(s[5]), s[6])
            rows = {k: byid[v] for k, v in G.items() if v in byid}
            if "hs1" not in rows or not all(k in rows for k in CORE):
                continue
            hchrom, hstart = rows["hs1"][0], rows["hs1"][1]
            raw = {k: np.frombuffer(v[5].encode(), np.uint8) for k, v in rows.items()}
            L = len(raw["hs1"])
            if L < 2 * W + 1 or any(len(v) != L for v in raw.values()):
                continue
            code = {k: LUT[v] for k, v in raw.items()}
            hs_nongap = code["hs1"] < 4
            # hs1 coordinate for each column (ref normally ungapped; handle gaps defensively)
            hpos = hstart + np.cumsum(raw["hs1"] != ord("-")) - 1
            rep = (raw["hs1"] >= ord("a")) & (raw["hs1"] <= ord("z"))
            bad_core = np.zeros(L, bool)
            for k in CORE:
                bad_core |= code[k] == 4
            core_ok = win_sum(bad_core) == 0
            og = "sorang" if "sorang" in code else ("borang" if "borang" in code else None)
            if og and "gorilla" in code:
                bad5 = bad_core | (code["gorilla"] == 4) | (code[og] == 4)
                v = (win_sum(bad5) == 0) & hs_nongap & (code["chimp"] == code["bonobo"])
                h, pn, gg, oo = code["human"], code["chimp"], code["gorilla"], code[og]
                pats = [v, v & (h == pn) & (gg == oo) & (h != gg), v & (h == gg) & (pn == oo) & (h != pn),
                        v & (pn == gg) & (h == oo) & (pn != h)]
                wins_all = hpos // WIN
                for j, m in enumerate(pats):
                    if m.any():
                        for wv, cnt in zip(*np.unique(wins_all[m], return_counts=True)):
                            ils[(hchrom, int(wv))][j] += cnt
            for br, (par, chi, extra) in BR.items():
                if par not in rows or chi not in rows:
                    continue
                ok = core_ok.copy()
                if extra:
                    bad = np.zeros(L, bool)
                    for k in extra:
                        bad |= code[k] == 4
                    ok &= win_sum(bad) == 0
                p, c = code[par], code[chi]
                diff = (p != c) & (p < 4) & (c < 4)
                ok &= win_sum(diff) <= 2
                ok &= hs_nongap
                if not ok.any():
                    continue
                # callable counts per 100-kb window
                idx = np.nonzero(ok)[0]
                wins = hpos[idx] // WIN
                pw = WEAK[p[idx]]
                for wv in np.unique(wins):
                    m = wins == wv
                    a = callable_[(hchrom, int(wv), br)]
                    a[0] += m.sum()
                    a[1] += pw[m].sum()
                    a[2] += (~pw[m]).sum()
                sites = np.nonzero(ok & diff)[0]
                if not len(sites):
                    continue
                # parsimony support
                def eq(a, b):
                    return (code[a] == code[b]) & (code[a] < 4) if (a in code and b in code) else np.zeros(L, bool)
                has = lambda k: (code[k] < 4) if k in code else np.zeros(L, bool)
                orang = "sorang" if "sorang" in code else ("borang" if "borang" in code else None)
                if br == "human":
                    pars = eq("chimp", "bonobo") & eq("chimp", "gorilla") & ~eq("human", "chimp")
                    if orang:
                        pars &= ~has(orang) | eq(orang, "chimp")
                elif br in ("chimp", "bonobo"):
                    other = "bonobo" if br == "chimp" else "chimp"
                    pars = eq(other, "human") & eq(other, "gorilla") & ~eq(br, other)
                elif br == "panstem":
                    pars = eq("chimp", "bonobo") & eq("human", "gorilla") & ~eq("chimp", "human")
                else:
                    pars = eq("human", "chimp") & eq("human", "bonobo") & ~eq("gorilla", "human")
                    if orang:
                        pars &= ~has(orang) | eq(orang, "human")
                # CpG context in parent (no gaps within window, so neighbours are adjacent columns)
                cpg = np.zeros(L, bool)
                cpg[:-1] |= (p[:-1] == 1) & (p[1:] == 2)   # C followed by G
                cpg[1:] |= (p[1:] == 2) & (p[:-1] == 1)    # G preceded by C
                # species coordinate of the child (or hs1 for internal branches)
                sp = chi if chi in ("chimp", "bonobo", "gorilla") else "hs1"
                schrom, sstart, _, sstrand, ssize, stext = rows[sp]
                sp_nongap = np.frombuffer(stext.encode(), np.uint8) != ord("-")
                sp_off = np.cumsum(sp_nongap) - 1
                for i in sites:
                    a, d = p[i], c[i]
                    cls = ("W" if WEAK[a] else "S") + ("W" if WEAK[d] else "S")
                    spos = sstart + sp_off[i]
                    if sstrand == "-":
                        spos = ssize - 1 - spos
                    fout.write(f"{hchrom}\t{hpos[i]}\t{br}\t{BASES[a]}\t{BASES[d]}\t{cls}\t{int(cpg[i])}\t{int(rep[i])}"
                               f"\t{int(pars[i])}\t{schrom}\t{spos}\t{ssize}\n")
    fout.close()
    with open(out + ".callable.tsv", "w") as f:
        f.write("chrom\twin_start\tbranch\tn_callable\tparent_W\tparent_S\n")
        for (ch, w, br), a in sorted(callable_.items()):
            f.write(f"{ch}\t{w * WIN}\t{br}\t{a[0]}\t{a[1]}\t{a[2]}\n")
    with open(out + ".ils.tsv", "w") as f:
        f.write("chrom\twin_start\tn_valid\tHP_G\tHG_P\tPG_H\n")
        for (ch, w), a in sorted(ils.items()):
            f.write(f"{ch}\t{w * WIN}\t{a[0]}\t{a[1]}\t{a[2]}\t{a[3]}\n")


if __name__ == "__main__":
    # optional 3rd arg: either a HAL genome id used as coordinate+human genome (e.g. hg38),
    # or a JSON file mapping roles -> genome ids (e.g. for the 16-way diploid HAL)
    if len(sys.argv) > 3:
        if sys.argv[3].endswith(".json"):
            import json
            G.update(json.load(open(sys.argv[3])))
        elif sys.argv[3] != "hs1":
            G["hs1"] = G["human"] = sys.argv[3]
    main(sys.argv[1], sys.argv[2])
