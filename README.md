# Tracing the origin of an ancestral F420-dependent dehydrogenase

Ana Ledesma · 1002083048 · BIOL 5340-001 Bioinformatics, Fall 2026

**Question.** Where does the reconstructed ancestral enzyme AncD1 sit within the
F420-dependent dehydrogenase family, and which residues distinguish the
substrate-specific (FGD) lineage from the promiscuous (FSD) one?

AncD1 was reconstructed by Mascotti et al. (2018) as the common ancestor of the
FGD and FSD lineages and has been expressed and kinetically characterised in the
Johnson-Winters laboratory. This repository holds the computational analysis:
dataset assembly, phylogeny, ancestral sequence reconstruction, a
specificity-determining position (SDP) predictor written for this project, and
structural mapping of the candidate residues.

---

## Contents

| File | What it does |
|---|---|
| `Stage1_bioinformatics_final.ipynb` | Dataset assembly: multi-seed BLASTP against UniProtKB, filtering, group-wise redundancy reduction. **Documentation of provenance — see the note below.** |
| `Stage1b_alignment_tree_finalfinal.ipynb` | Alignment (MAFFT + Clustal Omega), trimming, maximum-likelihood tree, rooting tests, clade definition, ancestral reconstruction, SDP scoring, structural mapping, figures |
| `Stage2_bioinformatics_final.ipynb` | Earlier combined version of the above, retained for reference |
| `environment.yml` | Conda specification: Python packages and the four external tools, with versions |
| `ai_usage.md` | What AI assistance was used for, how outputs were verified, and the errors that verification caught |

### Sequences

| File | What it is |
|---|---|
| `ancd1_set_final.fasta` | **The dataset of record** — 239 sequences, md5 `98084ed542c339e64eb5449550336ff1` |
| `AncD1_mature.fasta` | The 333-residue mature AncD1 sequence, His6–SUMO tag removed. All reported positions are numbered against this. |
| `merged_hits.json` | Merged output of the five BLASTP searches, 23 October 2026 — the frozen record of what the database returned |
| `homologs_nr_balanced.fasta` | Hits after filtering and group-wise redundancy reduction |
| `input_safe.fasta` | The dataset with tip names sanitised for the phylogenetics tools |
| `name_map.json` | Maps the sanitised names back to accessions, organisms and clade labels |
| `TableS1_sequences.xlsx` | All 239 sequences with accessions, organisms and clade assignments |

### Alignment and phylogeny

| File | What it is |
|---|---|
| `aln_mafft.fasta` | MAFFT L-INS-i alignment, 653 columns — the primary alignment |
| `aln_clustalo.fasta` | Clustal Omega alignment, 571 columns — the independent control |
| `aln_trimmed.fasta` | trimAl `-gappyout` output, 313 columns — the tree inference input |
| `ancd1_tree.treefile` | Maximum-likelihood tree with SH-aLRT and ultrafast bootstrap at every node |
| `ancd1_tree.iqtree` | IQ-TREE report: model selection, likelihood, parameter estimates |
| `ancd1_tree.log` | Full IQ-TREE run log |
| `tree_rooted_outgroup.nwk`, `tree_rooted_midpoint.nwk`, `tree_rooted_adh.nwk` | The three rootings tested; clade memberships are identical under all three |
| `tree_full_detailed.pdf` | The complete phylogeny — all 239 tips labelled, all 236 support values shown. Vector, zoomable. |
| `tree_groups.png` | The collapsed-clade summary figure |
| `itol_labels.txt`, `itol_styles.txt` | Annotation files for viewing the tree in iTOL |

### Ancestral reconstruction

| File | What it is |
|---|---|
| `ancd1_asr.treefile` | The tree with internal nodes labelled, as used for reconstruction |
| `ancd1_asr.iqtree`, `ancd1_asr.log` | IQ-TREE reconstruction report and log |
| `ancestral_seqs.fasta` | Reconstructed sequences at the three target nodes |
| `asr_nodes.json` | Node identities: C = Node84 (common ancestor), F = Node92 (FGD), S = Node15 (FSD) |
| `asr_summary.json` | Per-node mean posterior probabilities and sequence lengths |
| `substitutions_C_to_F_renum.csv` | Substitutions from the common ancestor to the FGD clade ancestor |
| `substitutions_C_to_S_renum.csv` | Substitutions from the common ancestor to the FSD clade ancestor |
| `substitutions_F_vs_S_renum.csv` | Positions differing between the two clade ancestors — 108 total, 86 confident |

### Specificity-determining positions

| File | What it is |
|---|---|
| `stage3_groups.json` | Clade membership used for the comparison, derived from support values |
| `sdp_ranked.csv` | All 313 columns scored, with bootstrap confidence intervals |
| `sdp_with_structure.csv` | The 63 robust positions with ligand distances and structural band |
| `ancd1_fgd_like_positions_renum.csv` | Positions where AncD1 carries the FGD-type residue despite its FSD-clade placement |
| `trimcol_to_ancd1.json` | Trimmed-alignment column → mature AncD1 residue mapping |

All `.csv` files whose name ends `_renum` use mature AncD1 numbering. Earlier
versions of these files numbered positions by trimmed-alignment column, which does
not correspond to residue number; those versions are superseded and are not
included here. See `ai_usage.md`.

---

## Reproducibility — read this first

**Stage 1 queries UniProtKB live and will not reproduce.** The BLASTP searches go
out to the database as it stands on the day they run. UniProt is updated
continuously, and re-running the searches two weeks after the original analysis
returned a different sequence set, a different tree, and support values across
the backbone that fell from 95–100 to 50–80. No random seed fixes this, because
the variation is in the database rather than in the code.

The dataset is therefore treated as a deposited artifact rather than something
regenerated on each run:

- Homologues were retrieved from UniProtKB on **23 October 2026**.
- The resulting 239-sequence dataset is committed here as
  `ancd1_set_final.fasta` (md5 `98084ed542c339e64eb5449550336ff1`).
- **All reported results were produced from that file**, by
  `Stage1b_alignment_tree_finalfinal.ipynb`.
- `Stage1_bioinformatics_final.ipynb` documents how the dataset was assembled.
  It is not re-run as part of the analysis.

Within that fixed input, the pipeline is deterministic: MAFFT, Clustal Omega,
trimAl and CD-HIT are deterministic on fixed input, and the SDP bootstrap uses a
fixed seed. IQ-TREE's maximum-likelihood search and ultrafast bootstrap are
stochastic and should be run with `-seed` for an exact reproduction of the
reported support values.

### To reproduce

```bash
conda env create -f environment.yml
conda activate ancd1
```

Then open `Stage1b_alignment_tree_finalfinal.ipynb`, upload
`ancd1_set_final.fasta`, and run all cells. The notebook installs the external
tools if they are not already present.

---

## Pipeline summary

**Dataset.** Five seeded BLASTP searches (AncD1 plus one representative of each
characterised lineage) against UniProtKB, E ≤ 1e-10. A single query saturates its
hit quota inside the densely sampled FGD cluster and never reaches the divergent
FSDs, so multi-seed searching is necessary rather than cosmetic. Hits filtered to
≥70% query coverage and 250–450 residues, then reduced for redundancy
**group-wise** — 80% identity for the densely sampled FGD and ADH lineages, 95%
for the sparse FSD lineage. A single global threshold collapses the sparse group
along with the dense one; group-wise reduction raised FSD representation from 18
to 76 sequences. Final set: 239 sequences, 123 genera.

**Alignment and tree.** MAFFT v7.505 L-INS-i (653 columns) and, as an independent
control, Clustal Omega v1.2.4 (571 columns). trimAl `-gappyout` retains 313
columns. IQ-TREE v2.0.7, ModelFinder restricted to LG/WAG/JTT/Blosum62 with Γ,
I+Γ and R4 rate models; **LG+F+I+G4** selected by BIC. Support from 1000
ultrafast bootstrap and 1000 SH-aLRT replicates. Three rootings tested — outgroup,
ADH-clade and midpoint — all giving identical clade memberships.

**Clades.** Defined by topology and support (SH-aLRT ≥ 80, bootstrap ≥ 95), not by
database annotation. FGD clade: 28 sequences, 99.3/100. FSD clade: 22 sequences
plus AncD1, 92.6/99. AncD1 falls **inside** the FSD clade.

**Ancestral reconstruction.** States reconstructed on the fixed topology at the
FGD/FSD common ancestor (Node84), the FGD clade ancestor (Node92) and the FSD
clade ancestor (Node15). Mean posterior 0.96–0.97.

**SDP prediction.** See below.

---

## The SDP scoring method

A specificity determinant should be **conserved within each clade and different
between them**. Jensen–Shannon divergence measures only the second half, so a
column that is scattered in one clade and fixed in the other scores highly
without meaning anything. The score used here is JSD weighted by the **weaker**
of the two within-clade conservation values, so both conditions must hold:

| test case | raw JSD | weighted score |
|---|---|---|
| all Ala vs all Ser | 0.99 | 0.96 |
| same residue in both | 0.00 | 0.00 |
| **variable in FGD, fixed in FSD** | **0.85** | **0.10** |
| both variable | 0.00 | 0.00 |
| mostly Ala vs mostly Ser, some noise | 0.99 | 0.83 |

The function was validated against these constructed cases before being applied
to real data. Confidence intervals come from 500 bootstrap replicates resampling
sequences within each clade — necessary because the FSD clade has only 22 members
and a point estimate there can rest on two or three sequences.

Three filters are applied before anything is reported: a bootstrap CI lower bound
above 0.20, concordance between the two alignments (they disagree on 12.6% of
columns, so a misplaced gap would otherwise produce a false signal), and the
structural context, measured as minimum distance to the bound cofactor (F420) and
to the citrate ion marking the anion site in PDB 3B4Y.

**Result:** 63 of 313 positions pass the statistical filters; 55 of those also
appear in the confident ancestral substitution set. Structurally, 7 are in direct
contact with a ligand, 10 in the active site, 14 in the second shell and 32
distal. Two positions — A108 and H257 — are robust SDPs where AncD1 retains the
FGD-type residue despite sitting in the FSD clade. H257 contacts both the
cofactor (3.08 Å) and the anion site (3.01 Å).

---

## All positions are mature AncD1 numbering

Positions refer to mature AncD1, M1 = start of `MLRFGYKAS…`, excluding the
99-residue N-terminal His6–SUMO tag of the expression construct. Add 99 for
construct numbering. The trimmed alignment has its own column numbering which
does **not** match residue numbers; `trimcol_to_ancd1.json` records the mapping,
and every reported position is verified against the AncD1 sequence before being
printed.

---

## Limitations

- Sampling is dense around the FGD and FSD lineages and sparse elsewhere in the
  superfamily. Root-to-tip distances average 6.1 substitutions per site, so the
  deep backbone is near saturation and is not interpreted.
- The two ADH-type references do not form a clade here; this reflects the absence
  of intervening ADH sequences from the dataset, not evidence against the group.
- Structural context comes from an extant FGD crystal structure used as a proxy
  for an enzyme from the other clade. Predicted models of AncD1 and the ancestors
  are in progress.
- Clade membership is defined by branch support, so the FGD and FSD groups are
  deliberately conservative cores rather than every sequence carrying the
  corresponding database annotation. 154 further sequences are annotated as
  F420-dependent glucose-6-phosphate dehydrogenases but fall outside the
  well-supported clade.

---

## References

Bashiri, G., Squire, C. J., Moreland, N. J., & Baker, E. N. (2008). Crystal
structures of F420-dependent glucose-6-phosphate dehydrogenase FGD1. *Journal of
Biological Chemistry*, 283(25), 17531–17541.

Capella-Gutiérrez, S., Silla-Martínez, J. M., & Gabaldón, T. (2009). trimAl.
*Bioinformatics*, 25(15), 1972–1973.

Katoh, K., & Standley, D. M. (2013). MAFFT version 7. *Molecular Biology and
Evolution*, 30(4), 772–780.

Mascotti, M. L., Kumar, H., Nguyen, Q.-T., Ayub, M. J., & Fraaije, M. W. (2018).
Reconstructing the evolutionary history of F420-dependent dehydrogenases.
*Scientific Reports*, 8, 17571.

Minh, B. Q., et al. (2020). IQ-TREE 2. *Molecular Biology and Evolution*, 37(5),
1530–1534.

Sievers, F., et al. (2011). Clustal Omega. *Molecular Systems Biology*, 7, 539.
