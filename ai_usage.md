# AI use

## Project authorship

The research question, the enzyme, the experimental data and the analytical
decisions in this project are mine. AncD1 is the subject of my dissertation work
in the Johnson-Winters laboratory; I expressed and purified the protein and
measured its steady-state kinetics. The computational question addressed
here — where AncD1 sits within the F420-dependent dehydrogenase family, and which
residues separate the specific from the promiscuous lineage — is one I posed, and
the design of the analysis is mine.

AI assistance was used as a programming tool. This file records what it was used
for and how its output was checked, as required by the course.

## Tool

**Claude Opus 5** (Anthropic), accessed through Claude Code, between 22 October
and 1 December 2026. The session was configured for the model identifier
`claude-opus-5`; a runtime may substitute a different model for an individual
request, which is not visible from within the session, so this is the configured
identifier rather than a verified per-request record.

## Division of work

**My contributions.** Choosing the question and the analytical strategy; all
decisions about scope, thresholds and interpretation; the biochemical and
taxonomic knowledge used to judge whether outputs were plausible; identifying
errors (below); deciding what the results mean and how to present them.

**Tool contributions.** Writing and debugging Python for the BLASTP queries, XML
parsing, CD-HIT orchestration, IQ-TREE invocation, ancestral-state parsing, the
SDP scoring module, structural mapping and the figure scripts; diagnosing
infrastructure failures (a hung NCBI BLAST job, a silent apt failure, IQ-TREE
flag changes between versions); drafting prose that I then revised.

**Not used for.** Any experimental data. Every kinetic measurement is laboratory
work.

## Errors I identified

These were caught by domain knowledge, not by code, and each changed the
analysis.

**Reference enzymes mislabelled.** The initial reference panel followed the clade
assignments in Mascotti et al. (2018) and listed the *Nocardioides* and
*Cryptosporangium arvum* enzymes as FGDs. Later kinetic characterisation
identifies both as FSDs. I corrected the assignments before the clade analysis
ran; left uncorrected, the two comparison groups would have been defined wrongly
and every downstream result would have been invalid.

**Wrong accession for the *Nocardioides* enzyme.** An incorrect accession was
initially proposed for FSD-Noca. I supplied the correct one (E9V3D1).

**Premature abandonment of model selection.** When the full ModelFinder search ran
long, the suggestion was to fix the substitution model instead. I questioned
whether the search was worth waiting for; the compromise adopted — restricting the
candidate set to matrices appropriate for cellular proteins — preserved
BIC-based model selection rather than discarding it.

**Inadequate figures.** Three successive versions of the phylogeny figure were
rejected as unreadable or as hiding information. I required that support values
be shown at every node rather than filtered for tidiness, which is the convention
in this literature and the reason the supplementary tree is presented as it is.

**Over-interpretation of tip placement.** I identified that AncD1, placed as a tip
among extant sequences, cannot be shown by that placement to be ancestral to any
of them. This distinction is now stated explicitly in the results.

## Errors the code checks caught

Verification was built into the pipeline rather than done by inspection, because
the dangerous errors here produce plausible-looking wrong numbers.

**Trimmed-alignment numbering.** Candidate residues were first reported using
positions counted within the trimmed alignment (313 columns) rather than the
full-length protein (333 residues). Every position was wrong by 3–12 residues:
the corrected set is 54, 105, 108, 180, 218, 245, 252, 257, 279, 298, 319 against
an original 50, 100, 103, 174, 209, 236, 243, 247, 267, 286, 307. Caught by a
check that verifies each reported position against the AncD1 sequence — if
position 257 is reported as histidine, residue 257 must be histidine. Without
this, the mutagenesis targets would have been the wrong residues.

**A monophyly test that measured nothing.** Clade membership was initially tested
by asking whether the reference enzymes alone formed a monophyletic group. With
~230 homologs interspersed, that is false regardless of whether the tree is
correct. Clade definitions were re-derived from the most recent common ancestor
of the references.

**Ancestral nodes resolving to the tree root.** IQ-TREE writes its own rooting;
because the FSD clade spanned that root, two target nodes resolved to the root
itself, producing an impossible "zero substitutions" result. Fixed by re-rooting
on the outgroup before locating nodes.

**Non-reproducible dataset assembly.** The pipeline as originally written
re-queried UniProtKB on every run. Re-running the searches two weeks after the
original analysis returned a different sequence set: the tree changed, the FSD
clade went from 22 to 25 sequences, and backbone support fell from 95-100 to
50-80. Nothing in the output looked wrong - the second run produced a perfectly
plausible tree with different numbers. I found this by re-running the pipeline
end to end, which is the only way it could have surfaced. The analysis is now
run from a deposited dataset (`ancd1_set_final.fasta`, md5
98084ed542c339e64eb5449550336ff1, retrieved 23 October 2026) rather than from a
live database query, and the retrieval step is retained as documentation of
provenance only.

## Verification built into the pipeline

1. **Input assertions** — the construct is checked as 432 aa with a 99-aa
   His6–SUMO tag and a 333-aa mature protein beginning `MLRFGYKAS`, so a wrong
   sequence fails loudly rather than shifting every position by 99.
2. **Position verification** — every column-to-residue conversion is checked
   against the AncD1 sequence before any position is reported.
3. **Known-answer tests** — the SDP scoring function was validated on constructed
   cases with known correct answers before use on real data.
4. **Independent replication** — all position-specific results must hold under two
   independently computed alignments, which disagree on 12.6% of columns.
5. **Robustness checks** — three rooting schemes for the topology; bootstrap
   confidence intervals (500 replicates) for every SDP score.

## Ethics and data

All sequence data come from public repositories (UniProtKB, PDB). No human
subject data, no personally identifying information, no restricted data. The
AncD1 sequence originates from Mascotti et al. (2018) and is cited throughout.
Database annotations are treated as hypotheses rather than measurements: FGD2
from *M. tuberculosis* is annotated as a glucose-6-phosphate dehydrogenase but is
reported to act on glucose, and this analysis places it outside the FGD clade
accordingly.

The AI tool is not an author, and is disclosed in the methods of any resulting
manuscript.
