# Third-party notices

Arachnia embeds or derives from the following third-party works. The original
Arachnia source code is licensed under the MIT license; these notices apply to
the listed components and model data.

## Rule Set 3 model (`assets/rs3_model.bin`)

- Source: <https://github.com/gpp-rnd/rs3>
- File used: `rs3/RuleSet3.pkl` (LightGBM regressor)
- License: Apache License 2.0
- Copyright: Broad Institute / GPP (see upstream repository)
- Modification: the 2,279-tree LightGBM ensemble was flattened into a compact
  binary (`scripts/export_models.py`). Tree structure, thresholds and leaf
  values are unchanged. Rule Set 3 is evaluated by original Rust code.

The full Apache License 2.0 text is in [`licenses/Apache-2.0.txt`](licenses/Apache-2.0.txt).

The 632-feature extraction is a Rust reimplementation of
`sglearn.featurize_guides` from <https://github.com/gpp-rnd/sglearn>
(MIT license, Copyright (c) 2020 Peter C. DeWeirdt). See
[`licenses/MIT-sglearn.txt`](licenses/MIT-sglearn.txt).

## SeqDeepCpf1 (`assets/seqdeepcpf1_weights.bin`)

- Source: <https://github.com/alisakumbara/SeqDeepCpf1-Py3>
- File used: `weights/Seq_deepCpf1_weights.h5`
- Original model: Kim, H. et al. *Deep learning improves prediction of
  CRISPR–Cpf1 guide RNA activity.* Nature Biotechnology 36, 239–241 (2018).
  <https://doi.org/10.1038/nbt.4061>
- License: MIT, Copyright (c) 2025 Alisa Kumbara
- Modification: Keras weights were converted to a compact little-endian
  binary, using the same axis permutations and convolution flip as the
  upstream PyTorch conversion script. The network is evaluated by original
  Rust code.

The full MIT license is in
[`licenses/MIT-SeqDeepCpf1.txt`](licenses/MIT-SeqDeepCpf1.txt).

## SeqFold

- Source: <https://github.com/Lattice-Automation/seqfold>
- Version: v0.10.2, used as a Rust dependency for RNA secondary-structure ΔG.
- License: MIT
- Copyright: Lattice Automation

The full MIT license is in [`licenses/MIT-seqfold.txt`](licenses/MIT-seqfold.txt).

## Biopython thermodynamic tables

The Rust melting-temperature implementation in `src/tm.rs` reproduces
Biopython's `Bio.SeqUtils.MeltingTemp.Tm_NN` with the `DNA_NN3` and
`R_DNA_NN1` tables and salt correction method 5. Thermodynamic constants are
from the published nearest-neighbour parameter sets (Allawi & SantaLucia 1997;
Sugimoto et al. 1995) and were cross-checked against Biopython 1.88.

- Biopython: <https://biopython.org>
- License: Biopython License Agreement (BSD-3-Clause-like)
- Full text: [`licenses/Biopython-License.txt`](licenses/Biopython-License.txt)

## Cas-OFFinder

- Source: <https://github.com/snugel/cas-offinder>
- Reference used: Cas-OFFinder's input/output format and mismatch-search model
  (Bae, Park & Kim, *Bioinformatics* 2014).
- License: BSD-3-Clause, Copyright (c) 2013 Jeongbin Park and Sangsu Bae
- Relationship: `src/offtarget.rs` is an independent Rust reimplementation of
  the exact-PAM Hamming search and Cas-OFFinder-style TSV output. No
  Cas-OFFinder source code is compiled into or vendored in Arachnia. DNA/RNA
  bulges are not implemented.

The full BSD-3-Clause license is in
[`licenses/BSD-Cas-OFFinder.txt`](licenses/BSD-Cas-OFFinder.txt).

## hg38 reference data

- Source: UCSC Genome Browser, `hg38.fa.gz`
- URL: <https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz>
- Relationship: hg38 is **not** bundled in the repository or binaries.
  `scripts/fetch_hg38.sh` downloads and decompresses it on demand into the
  user's data directory.
- Usage is subject to the
  [UCSC Genome Browser data usage policy](https://genome.ucsc.edu/license/).

## GENCODE annotation data

- Source: GENCODE human (GRCh38/hg38) GTF/GFF3 annotation
- URL: <https://www.gencodegenes.org/human/>
- Relationship: the annotation is **not** bundled in the repository or
  binaries. `scripts/fetch_annotation.sh` downloads and decompresses it on
  demand into the user's data directory.
- Citation: Frankish et al., “GENCODE reference annotation for the human and
  mouse genomes”, *Nucleic Acids Research* (2019),
  <https://doi.org/10.1093/nar/gky955>.
- Usage is subject to the
  [EMBL-EBI Terms of Use](https://www.ebi.ac.uk/about/terms-of-use/).

## GUI and runtime dependencies

Arachnia uses the following crates under their respective permissive licenses:

- `eframe`, `egui`, `egui_extras` — MIT OR Apache-2.0
- `rayon` — MIT OR Apache-2.0
- `seqfold` — MIT (see above)
- `serde_json` (test-only) — MIT OR Apache-2.0

Cargo's dependency graph contains additional transitive dependencies; their
license metadata is available via `cargo metadata` / `cargo deny`.

## Citation

If you use Arachnia in academic work, please cite the underlying methods:

- Rule Set 3: <https://github.com/gpp-rnd/rs3> and the CRISPR/Cas9 design
  resources from the Broad Institute GPP.
- SeqDeepCpf1: Kim, H., Min, S., Song, M. et al. Deep learning improves
  prediction of CRISPR–Cpf1 guide RNA activity. *Nat. Biotechnol.* **36**,
  239–241 (2018). doi:10.1038/nbt.4061.
- Cas-OFFinder: Bae, S., Park, J. & Kim, J.-S. Cas-OFFinder: a fast and
  versatile algorithm that searches for potential off-target sites of Cas9
  RNA-guided endonucleases. *Bioinformatics* **30**, 1473–1475 (2014).
  doi:10.1093/bioinformatics/btu048.
- crisprScore / crisprVerse: <https://github.com/crisprVerse/crisprScore>.
