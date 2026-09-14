# AI Use Log
- Tool/model & version:
- What I asked for:
- Snippet of prompt(s):
- What I changed before committing:
- How I verified correctness (tests, sample data):


- Tool/model & version: Gemini 2.5 Flash
- What I asked for: Store the content, close the file, and join ( ) the list
- Snippet of prompt(s):with open("practice.txt","r") as data:

  list_of_lines=data.readlines()

print (list_of_lines)
- What I changed before committing: with open("practice.txt","r") as data:

  list_of_lines=data.readlines()

print (list_of_lines)
- How I verified correctness (tests, sample data): The result was correct, no error given


Tool/model & version (9/13/26): Gemini 2.5 Flash
- What I asked for:Run the cell below to scan for TFBS
- Snippet of prompt: The error KeyError: 'passes_threshold' is happening because no promoters were successfully loaded. This is likely due to a mismatch between the chromosome identifiers in your FASTA file (e.g., 'NC_000022.11') and your GFF3 file (e.g., '22'). I'll modify the code to normalize the chromosome IDs when loading the genome, so they match and promoters can be correctly identified.

- Fixed Code:
- import io, re, math, sys
from typing import Dict, List, Tuple
from collections import defaultdict

import pandas as pd
from google.colab import files
from Bio import SeqIO, motifs
from Bio.Seq import Seq

# ---------------------------
# Settings for the exercise
# ---------------------------
UPSTREAM_BP   = 1000     # length of upstream region to scan (5' of TSS)
PSEUDOCOUNTS  = 0.5      # for PFM -> PPM adjustment for low read count
PCT_OF_MAX    = None     # threshold as % of the max possible PWM score
ABS_THRESHOLD = 14     # or set to e.g. 6.0 (log-odds) to use absolute threshold


# ---------------------------
# Helpers
# ---------------------------
def compute_bg_from_sequences(seq_dict):
    counts = defaultdict(int); total = 0
    for s in seq_dict.values():
        for ch in s.upper():
            if ch in "ACGT":
                counts[ch]+=1; total+=1
    if total == 0:
        return {"A":0.25,"C":0.25,"G":0.25,"T":0.25}
    return {b: counts.get(b,0)/total for b in "ACGT"}

def parse_gff3(path):
    feats = []
    with open(path, "r") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = parts
            if ftype not in ("gene", "mRNA"):
                continue
            try:
                start = int(start); end = int(end)
            except:
                continue
            attrd = {}
            for kv in attrs.split(";"):
                if "=" in kv:
                    k,v = kv.split("=",1); attrd[k]=v
            gene_id = attrd.get("Name") or attrd.get("gene") or attrd.get("gene_name") or attrd.get("locus_tag") or attrd.get("ID") or f"{seqid}:{start}-{end}"
            feats.append({"seqid": seqid, "start": start, "end": end, "strand": strand, "id": gene_id})
    return feats

def clamp(a, lo, hi): return max(lo, min(a, hi))

# ---------------------------
# Load genome & promoters
# ---------------------------

# Normalize sequence IDs from FASTA to match GFF3 (e.g., 'NC_000022.11' to '22')
parsed_genome = SeqIO.parse(genome_fa, "fasta")
genome = {}
for rec in parsed_genome:
    chrom_id = rec.id
    if chrom_id.startswith('NC_0000'):
        chrom_id = chrom_id.split('.')[0].lstrip('NC_0000')
    elif chrom_id.startswith('chr'):
        chrom_id = chrom_id.lstrip('chr')
    genome[chrom_id] = str(rec.seq).upper()

bg = compute_bg_from_sequences(genome)

features = parse_gff3(gff3_file)
promoters = {}
for f in features:
    chrom = genome.get(f["seqid"])
    if chrom is None:
        continue
    L = len(chrom)
    if f["strand"] == "+":
        tss = f["start"]
        up_end   = clamp(tss-1, 0, L)
        up_start = clamp(up_end - UPSTREAM_BP, 0, L)
        prom = chrom[up_start:up_end]
    else:
        tss = f["end"]
        down_start = clamp(tss, 0, L)
        down_end   = clamp(tss + UPSTREAM_BP, 0, L)
        prom = chrom[down_start:down_end]
        prom = str(Seq(prom).reverse_complement())
    if len(prom) >= 2:   # must be >= motif length; we’ll check later
        promoters[f["id"]] = prom.upper()

print(f"Loaded {len(promoters)} promoters.")
print("Background (A,C,G,T):", bg)

# ---------------------------
# Motif: PFM -> PPM -> PSSM (log-odds)
# ---------------------------
with open(motif_file, "r") as fh:
    m = motifs.read(fh, "jaspar")

ppm = m.counts.normalize(pseudocounts=PSEUDOCOUNTS)         # probabilities
pssm = ppm.log_odds(background=bg)                          # log2-odds

# Biopython’s matrices are row-major (rows = alphabet, cols = positions).
alphabet = getattr(pssm, "alphabet", "ACGT")
row_index = {b:i for i,b in enumerate(alphabet)}
motif_len = len(pssm[0])

PWM = []
for i in range(motif_len):
    PWM.append({
        'A': float(pssm[row_index['A']][i]),
        'C': float(pssm[row_index['C']][i]),
        'G': float(pssm[row_index['G']][i]),
        'T': float(pssm[row_index['T']][i]),
    })

def max_possible_score(pwm):
    return sum(max(col.values()) for col in pwm)

MAX = max_possible_score(PWM)
if ABS_THRESHOLD is None:
    cutoff = PCT_OF_MAX * MAX
    cutoff_desc = f"{PCT_OF_MAX:.0%} of max (={cutoff:.3f})"
else:
    cutoff = float(ABS_THRESHOLD)
    cutoff_desc = f"absolute log-odds ≥ {cutoff:.3f}"

print(f"Motif length: {motif_len}; threshold: {cutoff_desc}")

# ---------------------------
# Scanner (both strands)
# ---------------------------
def score_kmer(kmer, pwm):
    if len(kmer) != len(pwm): return None
    s = 0.0
    for i,b in enumerate(kmer.upper()):
        if b not in "ACGT": return None
        s += pwm[i][b]
    return s

def best_hit(seq, pwm):
    L = len(seq); k = len(pwm)
    best = (-1e18, -1, '+')
    for i in range(L-k+1):
        s = score_kmer(seq[i:i+k], pwm)
        if s is not None and s > best[0]:
            best = (s, i, '+')
    rc = str(Seq(seq).reverse_complement())
    for i in range(L-k+1):
        s = score_kmer(rc[i:i+k], pwm)
        if s is not None and s > best[0]:
            orig_i = L - k - i
            best = (s, orig_i, '-')
    return best

# ---------------------------
# Scan promoters & report genes
# ---------------------------
rows = []
for gene_id, prom in promoters.items():
    if len(prom) < motif_len:
        rows.append(dict(gene_id=gene_id, promoter_len=len(prom), best_score=None, best_pos0=None, best_strand=None, passes_threshold=False))
        continue
    s,i,strand = best_hit(prom, PWM)
    rows.append(dict(gene_id=gene_id, promoter_len=len(prom), best_score=round(s,3), best_pos0=i, best_strand=strand, passes_threshold=bool(s>=cutoff)))

df = pd.DataFrame(rows).sort_values(["passes_threshold","best_score"], ascending=[False, False])
df.to_csv("tfbs_hits.csv", index=False)

print("\nTop hits:")
display(df[df.passes_threshold].head(10))
print("\nSaved: tfbs_hits.csv")
files.download("tfbs_hits.csv")
