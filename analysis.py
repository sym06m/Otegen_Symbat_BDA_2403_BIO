"""Reference-based analysis of all assemblies:
 1. reference regions NOT covered by each assembly (and overlap with rRNA / tRNA / IS elements, GC)
 2. indels / SNPs of each assembly vs reference, and fraction of indels in homopolymers >= 4
Usage: python scripts/analysis.py   (run from the repository root, after assemblies exist)
"""
import glob, os, re, subprocess
from pathlib import Path
import pandas as pd

REF = "data/ref/ref.fasta"
GFF = "data/ref/ref.gff"
T = 2

def sh(cmd):
    subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True)

def read_fasta(p):
    seq = []
    for l in open(p):
        if l.startswith(">"):
            if seq:
                break
            continue
        seq.append(l.strip().upper())
    return "".join(seq)

def list_assemblies():
    pats = ["asm/short_*x/contigs.fasta", "asm/long_*x*/assembly.fasta",
            "asm/hybridspades_*/contigs.fasta", "asm/polish_*/polished.fasta"]
    return [(Path(p).parent.name, p) for pat in pats for p in sorted(glob.glob(pat)) if os.path.getsize(p) > 0]

REFSEQ = read_fasta(REF)

def gc(s):
    return (s.count("G") + s.count("C")) / max(1, len(s))

def ov(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

# ---------- 1. uncovered reference regions ----------
if not Path(GFF).exists():
    url = open("config.yaml").read().split("reference_url:")[1].split('"')[1].replace("_genomic.fna.gz", "_genomic.gff.gz")
    sh(f"wget -q -O data/ref/ref.gff.gz {url} && gunzip -c data/ref/ref.gff.gz > {GFF}")
FEAT = {"rRNA": [], "tRNA": [], "mobile_genetic_element": []}
for l in open(GFF):
    if l.startswith("#"):
        continue
    c = l.rstrip("\n").split("\t")
    if len(c) >= 9 and c[2] in FEAT:
        FEAT[c[2]].append((int(c[3]) - 1, int(c[4])))

def uncovered(fasta, label, min_len=100):
    paf = f"results/{label}.paf"
    sh(f"mkdir -p results && minimap2 -x asm10 -t {T} {REF} {fasta} > {paf}")
    iv = sorted((int(c[7]), int(c[8])) for c in (l.split("\t") for l in open(paf) if l.strip()))
    merged = []
    for s, e in iv:
        if merged and s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    gaps, prev = [], 0
    for s, e in merged:
        if s - prev >= min_len:
            gaps.append((prev, s))
        prev = max(prev, e)
    if len(REFSEQ) - prev >= min_len:
        gaps.append((prev, len(REFSEQ)))
    return gaps

rows = []
for lab, fa in list_assemblies():
    gaps = uncovered(fa, lab)
    tot = sum(e - s for s, e in gaps)
    row = {"assembly": lab, "n_gaps": len(gaps), "missing_bp": tot, "missing_pct": round(100 * tot / len(REFSEQ), 2)}
    for k, fl in FEAT.items():
        row[f"bp_in_{k}"] = sum(ov(g, f) for g in gaps for f in fl)
    row["GC_missing"] = round(gc("".join(REFSEQ[s:e] for s, e in gaps)), 3) if gaps else None
    rows.append(row)
pd.DataFrame(rows).to_csv("results/missing_regions.csv", index=False)

# ---------- 2. indels in homopolymers ----------
def run_len(seq, pos, base):
    i = pos
    while i > 0 and seq[i - 1] == base:
        i -= 1
    j = pos
    while j < len(seq) and seq[j] == base:
        j += 1
    return j - i

def ctx(st, b):
    best = 1
    for p in (st - 1, st):
        if 0 <= p < len(REFSEQ) and REFSEQ[p] == b:
            best = max(best, run_len(REFSEQ, p, b))
    return best

pos_hp4, i = 0, 0
while i < len(REFSEQ):
    j = i
    while j < len(REFSEQ) and REFSEQ[j] == REFSEQ[i]:
        j += 1
    if j - i >= 4:
        pos_hp4 += j - i
    i = j
base_hp4 = 100 * pos_hp4 / len(REFSEQ)

amap = dict(list_assemblies())
targets = [l for l in ["short_30x", "long_20x", "long_30x-rep43", "polish_s30_l20", "hybridspades_s30_l20"] if l in amap]
rows = []
for lab in targets:
    calls = f"results/{lab}.calls.vcf"
    sh(f"minimap2 -cx asm10 --cs -t {T} {REF} {amap[lab]} | sort -k6,6 -k8,8n | paftools.js call -l 1000 -L 1000 -f {REF} - > {calls}")
    snp = ins = dele = hp4 = 0
    bins = {k: 0 for k in range(1, 7)}
    for l in open(calls):
        if l.startswith("#"):
            continue
        c = l.rstrip("\n").split("\t")
        try:
            pos, r, a = int(c[1]), c[3], c[4]
        except Exception:
            continue
        if a.startswith("<") or "," in a:
            continue
        if len(r) == len(a):
            snp += 1
            continue
        if len(a) > len(r):
            ins += 1; b = a[1]
        else:
            dele += 1; b = r[1]
        h = ctx(pos, b)
        bins[min(h, 6)] += 1
        hp4 += h >= 4
    n = ins + dele
    rows.append({"assembly": lab, "snps": snp, "insertions": ins, "deletions": dele,
                 "indels_in_hp>=4_pct": round(100 * hp4 / n, 1) if n else None,
                 "baseline_ref_pos_in_hp>=4_pct": round(base_hp4, 1),
                 "hp1": bins[1], "hp2": bins[2], "hp3": bins[3], "hp4": bins[4], "hp5": bins[5], "hp6plus": bins[6]})
pd.DataFrame(rows).to_csv("results/indel_context.csv", index=False)
