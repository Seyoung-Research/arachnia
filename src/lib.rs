//! # Arachnia
//!
//! Arachnia is a pure-Rust desktop application for designing CRISPR guide RNAs
//! and enumerating Cas-OFFinder-style off-target sites. It implements two
//! published on-target scoring models and a packed, parallel off-target search
//! engine, so no Python runtime or external aligner is required.
//!
//! The crate is both the engine behind the desktop UI and a reusable library.
//! The binary targets are:
//!
//! * `arachnia` — the `eframe`/`egui` desktop application.
//! * `cas-offinder` — a command-line front end with Cas-OFFinder-compatible
//!   input and TSV output formats.
//!
//! # Feature overview
//!
//! ## On-target design
//!
//! * **SpCas9**: 20 nt guide + `NGG` PAM, scored with Rule Set 3 sequence mode.
//!   The model is the LightGBM ensemble from
//!   [`gpp-rnd/rs3`](https://github.com/gpp-rnd/rs3); its 632 features and tree
//!   traversal are reimplemented in [`ruleset3`].
//! * **Cpf1 (Cas12a)**: 23 nt guide + `TTTV` PAM, scored with SeqDeepCpf1 from
//!   Kim et al. (2018). The 1-D CNN forward pass is implemented in [`cpf1`].
//! * Both strands are scanned by default. Candidates are filtered by GC
//!   content, minimum score and poly-T run length, then sorted by score.
//! * Each candidate records its guide and PAM in guide orientation, its
//!   1-based plus-strand coordinates, approximate cut site, GC percentage,
//!   score, full model context and poly-T run length.
//! * The desktop app can sort, search, filter, copy and export results as TSV.
//! * A built-in TP53 coding-sequence demo makes the application useful
//!   immediately.
//!
//! ## Off-target analysis
//!
//! [`offtarget`] implements the mismatch model of
//! [Cas-OFFinder](https://github.com/snugel/cas-offinder):
//!
//! * Exact IUPAC PAM matching (`NGG` for Cas9, `TTTV` for Cpf1, or any custom
//!   IUPAC pattern).
//! * Hamming mismatches over the guide only; the PAM is not counted.
//! * Both strands, with 0-based plus-strand locations using the same
//!   guide + PAM left-edge convention as Cas-OFFinder and Bowtie.
//! * Packed 2-bit guide comparison and Rayon-parallel scanning.
//! * Progress reporting, cancellation and a 200,000-site output cap.
//! * Streaming multi-record FASTA loading from a file or a directory.
//! * Optional local hg38 discovery via [`offtarget::find_hg38`] and the
//!   `hg38`/`GRCh38` shorthand accepted by the CLI.
//! * Optional chromosome-region annotation from GTF/GFF3: promoter, exon,
//!   coding exon, 5′/3′ UTR, intron, gene body or intergenic, with the
//!   overlapping gene. Promoters are 500 bp to 5 kb upstream of a TSS.
//! * Cas-OFFinder-compatible TSV export, with optional `Region`/`Gene`
//!   columns when annotation is loaded.
//!
//! DNA and RNA **bulges are not supported**. The search is intentionally
//! mismatch-only and exact-PAM.
//!
//! ## Desktop application
//!
//! The GUI provides:
//!
//! * A target-sequence editor with plain DNA, FASTA, multi-record FASTA
//!   handling and a one-click TP53 demo.
//! * SpCas9/Cpf1 selection, tracrRNA selection (Hsu2013 or Chen2013), both-
//!   strand toggle, GC range, minimum score and poly-T filters.
//! * A sortable, searchable guide table with score heat colors and flags.
//! * A selected-guide panel with the exact model context, approximate cut
//!   site, copy actions and an **Off-targets…** action.
//! * An off-target window with PAM pattern, mismatch limit, both-strands
//!   toggle, current-target or reference-genome search, background loading and
//!   cancellation, ranked mismatch table, mismatch highlighting, optional
//!   GTF/GFF3 region annotation, and copy or export of Cas-OFFinder-compatible
//!   TSV.
//! * A local hg38 reference is detected and preselected when present. Install
//!   one with `scripts/fetch_hg38.sh`; the app loads it on demand in a
//!   background thread.
//! * A local GTF/GFF3 annotation is detected automatically. Install a GENCODE
//!   file with `scripts/fetch_annotation.sh`.
//!
//! ### Desktop workflow
//!
//! 1. Paste plain DNA or a FASTA record, use **Load file…**, or click
//!    **TP53 demo**.
//! 2. Choose SpCas9 or Cpf1. For Cas9, choose the Hsu2013 or Chen2013
//!    tracrRNA. Adjust the both-strands toggle, GC range, minimum score and
//!    poly-T filter.
//! 3. Press **Find gRNAs** (or `Ctrl+Enter`). Changing the enzyme, tracrRNA or
//!    strand setting rescans automatically when results are already shown.
//! 4. Sort or search the results table. The rank, guide, PAM, strand, start,
//!    end, cut, GC, score, flags and context columns can all be inspected.
//! 5. Click a row to open the selected-guide panel. It shows the full model
//!    context and approximate cut site, and offers **Copy guide**,
//!    **Copy context**, **Copy TSV**, **Export TSV…** and **Off-targets…**.
//! 6. In the off-target window, set the PAM pattern, mismatch limit and
//!    both-strands toggle, then search the current target sequence or a
//!    reference genome. Load a GTF/GFF3 annotation to add `Region`/`Gene`
//!    labels; genome loading and searching run in the background with progress
//!    and **Cancel**.
//! 7. Copy or export the ranked off-target table as Cas-OFFinder-compatible
//!    TSV.
//!
//! The window is styled with the “obsidian & crimson” palette defined in the
//! private `theme` module. Results are not uploaded anywhere; the application
//! is fully offline.
//!
//! ## Command-line tool
//!
//! The `cas-offinder` binary accepts Cas-OFFinder-style input files:
//!
//! ```text
//! cas-offinder <input.txt|-> [C|G|A[device_ids]] [output.tsv|-]
//! ```
//!
//! The first line is a genome path or the `hg38`/`GRCh38` shorthand, followed
//! by a pattern line and one query line per guide:
//!
//! ```text
//! hg38
//! NNNNNNNNNNNNNNNNNNNNNGG
//! GAGCGCTGCTCAGATAGCGANNN 3 tp53
//! ```
//!
//! Cas9 queries put the PAM wildcards at the 3′ end (`guideNNN`). Cpf1 queries
//! put them at the 5′ end (`NNNNguide`) with a pattern such as `TTTV` followed
//! by 23 `N`s. Guide positions in the pattern must be `N`; only the PAM end
//! may be constrained. Optional bulge sizes must be `0 0`.
//!
//! The device argument (`C`, `G`, `A`, or a device list) is accepted for
//! compatibility and ignored; the Rust implementation always uses the CPU via
//! Rayon. Progress goes to standard error, so TSV on standard output is clean.
//!
//! The output columns are `Id`, `Bulge Type`, `crRNA`, `DNA`, `Chromosome`,
//! `Location`, `Direction`, `Mismatches` and `Bulge Size`. `Location` is
//! **0-based** and refers to the plus-strand left edge of the guide + PAM
//! window. `Bulge Type` is always `X` and `Bulge Size` is always `0`.
//!
//! See the `cas_offinder` binary documentation for the full format, examples,
//! reference resolution and exit status. After `cargo doc`, that page is
//! generated at `target/doc/cas_offinder/index.html`.
//!
//! ## hg38 reference setup
//!
//! Arachnia does not bundle hg38. Install a verified UCSC copy with:
//!
//! ```console
//! ./scripts/fetch_hg38.sh
//! ```
//!
//! The script downloads `hg38.fa.gz`, verifies its MD5 when available,
//! decompresses it, and writes
//! `${XDG_DATA_HOME:-$HOME/.local/share}/arachnia/hg38.fa`. It is idempotent
//! and honours `ARACHNIA_HG38_DIR`.
//!
//! The GUI and CLI look for hg38 in this order:
//!
//! 1. `$ARACHNIA_HG38` (file or directory);
//! 2. `$ARACHNIA_HG38_DIR/hg38.fa`;
//! 3. `$XDG_DATA_HOME/arachnia/hg38.fa`;
//! 4. `~/.local/share/arachnia/hg38.fa`;
//! 5. `~/genomes/hg38.fa`, `~/genomes/hg38.fasta`, `~/hg38.fa`;
//! 6. `/data/genomes/hg38.fa`, `/srv/genomes/hg38.fa`.
//!
//! The GUI preselects a discovered hg38 but loads it only when **Load** or
//! **Use hg38** is clicked, because the FASTA occupies about 3.1 GiB of RAM.
//!
//! ## Gene annotation setup
//!
//! Arachnia does not bundle a gene annotation. Install a GENCODE human
//! annotation with:
//!
//! ```console
//! ./scripts/fetch_annotation.sh
//! ```
//!
//! The script finds the latest GENCODE release, verifies the annotation
//! against the release's `MD5SUMS`, decompresses it, and writes
//! `${XDG_DATA_HOME:-$HOME/.local/share}/arachnia/gencode.annotation.gtf`. The
//! GUI detects this file automatically. The script accepts `gtf` or `gff3`
//! output via `ARACHNIA_ANNOTATION_FORMAT`, an `all`/`basic`/`primary_assembly`
//! annotation via `ARACHNIA_ANNOTATION_KIND`, and a pinned release via
//! `ARACHNIA_GENCODE_RELEASE`. `ARACHNIA_ANNOTATION_DIR` changes the
//! destination and is also searched by annotation discovery. The CLI still
//! requires `--annotation` or `ARACHNIA_ANNOTATION`.
//!
//! ## On-target scoring and filters
//!
//! Both enzymes are scored from a fixed context in guide orientation:
//!
//! * **SpCas9 / Rule Set 3** — 30 nt:
//!   `[4 nt upstream][20 nt guide][3 nt PAM][3 nt downstream]`. The 632
//!   features reproduce `sglearn.featurize_guides` plus the tracrRNA one-hot
//!   columns: GC/mononucleotide/dinucleotide fractions, position-specific
//!   nucleotides and dinucleotides, PAM-interaction features, five
//!   thermodynamic features, homopolymer runs and the tracrRNA choice. The
//!   2,279-tree LightGBM ensemble is traversed directly in Rust. The output is
//!   the raw model score; higher is more active.
//! * **Cpf1 / SeqDeepCpf1** — 34 nt:
//!   `[4 nt upstream][4 nt TTTV PAM][23 nt guide][3 nt downstream]`. The
//!   network is `Conv1d(4→80, k=5) → ReLU → AvgPool1d(2) → flatten → dense
//!   80 → 40 → 40 → 1`. The output is the paper’s 0–100 predicted activity;
//!   higher is better.
//!
//! [`Filters`] are applied after scoring: inclusive GC bounds, a minimum
//! score, an optional strict poly-T limit, and both-strand scanning. The
//! library default rejects T runs of four or more; the GUI keeps poly-T as a
//! flag/filter control and does not filter by default. Results are sorted by
//! descending score, then by plus-strand start and strand.
//!
//! Candidate coordinates are 1-based on the input strand. The approximate cut
//! site is 3 bp upstream of the Cas9 PAM, and a staggered position around
//! 18/23 of the Cpf1 protospacer. Scores are design aids, not measured editing
//! efficiencies.
//!
//! # Quick start (library)
//!
//! Design SpCas9 guides:
//!
//! ```
//! use arachnia::{design, Enzyme, Filters, Tracr};
//!
//! let sequence = "AAAAACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTGGCCC";
//! let filters = Filters {
//!     max_poly_t: None,
//!     both_strands: false,
//!     ..Filters::default()
//! };
//! let (prepared, result) =
//!     design(sequence, Enzyme::Cas9, Tracr::Hsu2013, &filters).unwrap();
//!
//! println!("{} bp, {} candidates", prepared.bases.len(), result.candidates.len());
//! for candidate in &result.candidates {
//!     println!("{} {} {} {}", candidate.guide, candidate.pam, candidate.strand.symbol(), candidate.score);
//! }
//! ```
//!
//! Search a small reference for off-targets:
//!
//! ```
//! use arachnia::offtarget::{
//!     search_offtargets, Genome, OffTargetParams, SearchControl,
//! };
//!
//! let genome = Genome::from_fasta_str(
//!     ">chr1\nGGGGACGTACGTACGTACGTACGTTGGAAAA\n",
//! )
//! .unwrap();
//! let params = OffTargetParams::default();
//! let result = search_offtargets(
//!     "ACGTACGTACGTACGTACGT",
//!     &genome,
//!     &params,
//!     &SearchControl::new(),
//! )
//! .unwrap();
//! assert_eq!(result.perfect_matches(), 1);
//! ```
//!
//! # Output formats
//!
//! ## On-target TSV
//!
//! [`results_tsv`] writes the visible candidates in display order with these
//! columns:
//!
//! ```text
//! rank  guide  pam  strand  start  end  cut  gc_percent  score  context  poly_t_run
//! ```
//!
//! `guide` and `pam` are in guide orientation; `start`, `end` and `cut` are
//! 1-based plus-strand coordinates. The GUI uses the same format for clipboard
//! copy and file export.
//!
//! ## Off-target TSV
//!
//! [`offtarget::results_tsv`] writes Cas-OFFinder-compatible rows:
//!
//! ```text
//! #Id  Bulge Type  crRNA  DNA  Chromosome  Location  Direction  Mismatches  Bulge Size
//! ```
//!
//! `Location` is 0-based and refers to the plus-strand left edge of the
//! guide + PAM window. `crRNA` uses `N`s for PAM positions; `DNA` contains the
//! actual guide-oriented reference sequence and PAM. `Bulge Type` is always
//! `X`, and `Bulge Size` is always `0`. When sites carry an
//! [`annotation::SiteAnnotation`], `Region` and `Gene` columns are appended
//! automatically.
//!
//! # Documentation map
//!
//! * [`mod@design`] — sequence parsing, enzyme/strand models, candidate discovery,
//!   filtering and on-target TSV export.
//! * [`ruleset3`] — Rule Set 3 featurization and LightGBM inference for SpCas9.
//! * [`cpf1`] — SeqDeepCpf1 CNN inference for Cpf1/Cas12a.
//! * [`tm`] — Biopython-compatible nearest-neighbour melting temperatures used
//!   by Rule Set 3.
//! * [`annotation`] — GTF/GFF3 parsing, TSS promoter intervals and
//!   off-target region labels.
//! * [`offtarget`] — genome loading, hg38 discovery, packed off-target search,
//!   progress/cancellation and Cas-OFFinder TSV output.
//!
//! # Coordinate conventions
//!
//! * On-target [`GuideCandidate`] coordinates are **1-based** on the input
//!   (plus) strand. `start` is the lowest coordinate and `end` the highest,
//!   regardless of guide strand. The guide and PAM themselves are always
//!   presented in guide orientation (5′→3′).
//! * Off-target [`offtarget::OffTargetSite::location`] values are **0-based**
//!   and refer to the left edge of the guide + PAM window on the plus strand,
//!   matching Cas-OFFinder and Bowtie. The guide-oriented DNA and PAM are
//!   reported separately.
//! * [`annotation::Annotation`] stores 0-based half-open intervals converted
//!   from 1-based GTF/GFF3 records. Off-target labels use the full guide + PAM
//!   window. A promoter spans 500 bp to 5 kb upstream of a TSS.
//!
//! # Performance and memory
//!
//! * Off-target scanning is parallelized with Rayon. Guide windows are packed
//!   into 2-bit integers, so mismatch counting is a handful of bit operations
//!   per candidate.
//! * The reference is held in memory. hg38 has 455 contigs and
//!   3,209,286,105 bases, so loading it needs roughly 3.1 GiB plus allocator
//!   overhead. The GUI and CLI load it once per process.
//! * A three-mismatch hg38 search for a single 20 nt guide typically takes a
//!   few seconds on a multi-core CPU after the reference is loaded.
//! * Results are capped at 200,000 sites per search. The cap is applied after
//!   sorting, and [`offtarget::OffTargetResult::truncated`] reports whether it
//!   was hit.
//! * On-target scoring is fast enough for live re-filtering; changing the
//!   enzyme, tracrRNA or strand setting triggers a rescan.
//!
//! # Validation and development
//!
//! `cargo test` checks Rule Set 3 featurization and inference against 84
//! regression cases, SeqDeepCpf1 against upstream examples, Biopython Tm
//! values, candidate coordinates and strand handling, off-target enumeration
//! against a naive implementation, and end-to-end CLI input/output.
//!
//! The optional real-hg38 smoke test validates a TP53 guide:
//!
//! ```console
//! cargo test --release --test hg38_smoke -- --ignored --nocapture
//! ```
//!
//! It expects a local hg38 installation. The model assets can be regenerated
//! with `scripts/export_models.py`; see the repository `README.md` for the
//! Python environment and upstream repositories.
//!
//! # Licenses and citations
//!
//! Arachnia’s original code is MIT-licensed. The embedded Rule Set 3 model is
//! Apache-2.0; SeqDeepCpf1 code and weights are MIT; SeqFold is MIT; the
//! Biopython thermodynamic tables follow the Biopython License Agreement; and
//! Cas-OFFinder is BSD-3-Clause. See `THIRD_PARTY_LICENSES.md` for details.
//!
//! Key references:
//!
//! * Kim, H. et al. Deep learning improves prediction of CRISPR–Cpf1 guide RNA
//!   activity. *Nature Biotechnology* **36**, 239–241 (2018).
//!   <https://doi.org/10.1038/nbt.4061>
//! * Bae, S., Park, J. & Kim, J.-S. Cas-OFFinder: a fast and versatile
//!   algorithm that searches for potential off-target sites of Cas9
//!   RNA-guided endonucleases. *Bioinformatics* **30**, 1473–1475 (2014).
//!   <https://doi.org/10.1093/bioinformatics/btu048>
//! * Rule Set 3 / GPP: <https://github.com/gpp-rnd/rs3>
//! * SeqDeepCpf1-Py3: <https://github.com/alisakumbara/SeqDeepCpf1-Py3>
//! * Cas-OFFinder: <https://github.com/snugel/cas-offinder>
//!
//! # Limitations
//!
//! * Off-target search is mismatch-only: PAMs must match exactly and DNA/RNA
//!   bulges are not modelled.
//! * Region annotation uses GTF/GFF3 transcript models. Compressed annotation
//!   files are not supported, and regulatory elements other than the defined
//!   TSS promoters are not annotated.
//! * Genome-wide searches load the reference into memory; compressed FASTA is
//!   not supported. hg38 needs roughly 3.1 GiB of RAM.
//! * Rule Set 3 is used in sequence mode; the upstream target-context model is
//!   not included.
//! * Cut positions are approximate (Cas9: 3 bp upstream of the PAM; Cpf1:
//!   staggered around positions 18/23 of the protospacer).
//! * Scores are design aids, not experimental measurements.
//!
//! See the repository `README.md` for installation, screenshots, validation
//! details, licenses and citations.

#![warn(missing_docs)]

pub mod annotation;
pub mod cpf1;
pub mod design;
pub mod offtarget;
pub mod ruleset3;
pub mod tm;

pub use annotation::{Annotation, AnnotationError, Region, SiteAnnotation};
pub use design::{
    design, design_bases, prepare_sequence, results_tsv, reverse_complement, DesignResult,
    DnaError, Enzyme, Filters, GuideCandidate, PreparedSequence, Strand,
};
pub use ruleset3::Tracr;
