"""Download all FASTQ files of one ENA run, verify md5, append to data/PROVENANCE.tsv."""
import sys, os, subprocess, hashlib, time, urllib.request

acc, outdir = sys.argv[1], sys.argv[2]
fields = "run_accession,fastq_ftp,fastq_md5"
url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={acc}&result=read_run&fields={fields}&format=tsv"
rows = urllib.request.urlopen(url).read().decode().strip().split("\n")
meta = dict(zip(rows[0].split("\t"), rows[1].split("\t")))
os.makedirs(outdir, exist_ok=True)

def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

for u, m in zip(meta["fastq_ftp"].split(";"), meta["fastq_md5"].split(";")):
    dest = os.path.join(outdir, os.path.basename(u))
    if not (os.path.exists(dest) and md5(dest) == m):
        subprocess.run(["wget", "-c", "-q", "-P", outdir, "https://" + u], check=True)
    assert md5(dest) == m, f"md5 mismatch for {dest}"
    with open("data/PROVENANCE.tsv", "a") as f:
        f.write(f"{acc}\t{os.path.basename(dest)}\t{u}\t{m}\t{time.strftime('%Y-%m-%d')}\n")
