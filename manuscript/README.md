# Proof of Notik's Singer-ruler conjecture with sparse certificates

Haobo Ma (ChronoAI Pte Ltd, Singapore), followed by Wenlin Zhang
(National University of Singapore; ChronoAI Pte Ltd, Singapore).
License: Creative Commons Attribution 4.0 International, CC BY 4.0.

## Mathematical contents

`note.tex` proves the largest-gap chamber conjecture in Notik's
*Golomb rulers with real marks*, Section 7, page 7:
Zenodo version 2, 31 October 2026, DOI
[10.5281/zenodo.23844781](https://doi.org/10.5281/zenodo.23844781),
and [arXiv:2611.22656v1](https://arxiv.org/abs/2611.22656v1).
It applies to every cyclic perfect difference set, not just the Singer
construction. The complete-prefix lemma also gives gapwise domination,
sharp coordinate stability and uniqueness for every positive weighted
gap objective. A unit-triangular certificate uses exactly one tight row
per adjacent gap. Largest-gap examples answer the Bose and Ruzsa family
questions negatively, with shorter chamber lengths 28 < 29 and 8 < 9.
The families are defined as largest-gap cuts of the cyclic Bose and Ruzsa
Sidon constructions, allowing unit multipliers. The Bose example also
belongs to the construction without multipliers:
`5 B(theta) = 8 + B(theta^11)` modulo 48. Perfect difference sets have
distinct cyclic gaps, hence a unique largest gap.

The bibliography is embedded in `note.tex`. No additional figures or
bibliography files are required.

## Build

Use Tectonic 0.17.0:

```sh
tectonic note.tex
```

Alternatively, run `pdflatex note.tex` twice; the embedded bibliography
does not require BibTeX. The archive contains the complete scientific
sources and a dependency-free Python checker.

## Exact finite checks

```sh
python3 verify.py > checks-new.json
cmp checks.json checks-new.json
python3 -O verify.py > checks-optimized.json
cmp checks.json checks-optimized.json
```

The checker enumerates all 94 cyclic perfect difference sets containing
zero with cardinalities 2 through 6, modulo 3, 7, 13, 21 and 31, and each
largest-gap cut. It checks the dense and sparse integer dual identities,
triangularity and exact objective values. It also checks every normalized
integer Golomb ruler with 2 through 5 marks and length at most 18 that
has a complete distance prefix: 27 rulers in all. The field and CRT
constructions, labelled distance orders and strict shortenings of the
Bose and Ruzsa examples are checked exactly, as is Notik's twelve-mark
non-largest-cut example. Invalid cyclic differences and an invalid cut
are rejected. `checks.json` contains the full finite ranges and witnesses.
The checker also verifies distinct cyclic gaps, strict inequality between
each retained gap and the removed gap, and the Bose translate identity.
These finite statements accompany the all-orders written proofs.

## Formal statements

The pinned
[Singer chamber source](https://github.com/the-omega-institute/trureturing/blob/a1ecee6d9de8eff7f246e1569bbd6211f2c471cd/D5/S3/Geometry/FiniteGeometry/NotikSingerChamber.lean)
contains `D5.S3.Geometry.FiniteGeometry.NotikSingerChamber.result`.
For normalized strictly increasing natural marks, a positive cut gap
bounding the adjacent gaps, and perfect cyclic differences modulo
length plus cut gap, it proves ordinary-distance uniqueness, real
chamber feasibility, and equality of all marks for any chamber point
of no greater length.

The pinned
[Bose and Ruzsa source](https://github.com/the-omega-institute/trureturing/blob/1a9d63013ef1b813c385f84da25351bc88cc7a98/D5/S3/Geometry/FiniteGeometry/NotikBoseRuzsaChamber.lean)
contains `D5.S3.Geometry.FiniteGeometry.NotikBoseRuzsaChamber.result`.
It proves the negation of the disjunction of the two universal family
statements using the same field and CRT witnesses displayed in the note.

The general complete-prefix lemma, sparse and dense certificate formulas,
stability and weighted-objective assertions, and the finite-field
construction are written proofs. The checker establishes the finite
instances specified above.
