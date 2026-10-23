//! Cas-OFFinder-compatible command-line off-target search.
//!
//! This binary is a pure-Rust, CPU-parallel front end to
//! [`arachnia::offtarget`](../arachnia/offtarget/index.html). It reads the
//! familiar Cas-OFFinder text format, searches one or more guide queries
//! against a local FASTA reference, and writes Cas-OFFinder-style TSV.
//!
//! # Usage
//!
//! ```text
//! cas-offinder [--annotation <GTF|GFF3>] <input.txt|-> [C|G|A[device_ids]] [output.tsv|-]
//! ```
//!
//! * `<input.txt|->` — input file, or `-` to read from standard input.
//! * `--annotation` — optional GTF/GFF3 gene annotation. Adds `Region` and
//!   `Gene` columns. The path can also be set with `ARACHNIA_ANNOTATION`.
//!   `scripts/fetch_annotation.sh` downloads and installs a GENCODE file.
//! * `[C|G|A[device_ids]]` — accepted for Cas-OFFinder compatibility and
//!   ignored. The Rust implementation always uses the CPU via Rayon.
//! * `[output.tsv|-]` — output path, or `-`/omitted for standard output.
//!
//! Progress and diagnostics go to standard error, so stdout remains clean TSV.
//!
//! # Input format
//!
//! Blank lines and lines beginning with `#` are ignored. The first content line
//! is the reference, the second is the search pattern, and every subsequent
//! line is a query:
//!
//! ```text
//! hg38
//! NNNNNNNNNNNNNNNNNNNNNGG
//! GAGCGCTGCTCAGATAGCGANNN 3 tp53
//! ```
//!
//! * **Reference** — `hg38`/`GRCh38` shorthand, a FASTA file, or a directory
//!   of FASTA files. The shorthand is resolved by
//!   [`resolve_reference_path`](../arachnia/offtarget/fn.resolve_reference_path.html).
//! * **Pattern** — a string whose length equals each query. Guide positions
//!   must be `N`; only the PAM end may contain IUPAC constraints. An optional
//!   trailing `dna_bulge rna_bulge` pair is accepted, but both values must be
//!   `0`.
//! * **Query** — `query max_mismatches [id]`. `query` is the guide plus PAM
//!   wildcards at one end. `id` defaults to the zero-based query number.
//!
//! ## PAM orientation
//!
//! * **Cas9 (3′ PAM)** — wildcards at the end: `GUIDE` + `NNN`, with a pattern
//!   such as 20 `N`s followed by `NGG`. Example:
//!   `GAGCGCTGCTCAGATAGCGANNN 3 tp53`.
//! * **Cpf1/Cas12a (5′ PAM)** — wildcards at the start: `NNNN` + `GUIDE`, with
//!   a pattern such as `TTTV` followed by 23 `N`s.
//!
//! A query with wildcards at both ends is rejected. The CLI determines the PAM
//! side from the wildcard position, so the same input format works for both
//! enzymes.
//!
//! # Reference resolution
//!
//! `hg38` and `GRCh38` are case-insensitive. The search order is:
//!
//! 1. `$ARACHNIA_HG38` (FASTA file or directory);
//! 2. `$ARACHNIA_HG38_DIR/hg38.fa`;
//! 3. `$XDG_DATA_HOME/arachnia/hg38.fa`;
//! 4. `~/.local/share/arachnia/hg38.fa`;
//! 5. `~/genomes/hg38.fa`, `~/genomes/hg38.fasta`, `~/hg38.fa`;
//! 6. `/data/genomes/hg38.fa`, `/srv/genomes/hg38.fa`.
//!
//! Run `scripts/fetch_hg38.sh` to install a verified UCSC hg38 assembly. Other
//! references can be passed as a path; `~` is expanded. Compressed FASTA is
//! not supported.
//!
//! # Output
//!
//! The TSV header is:
//!
//! ```text
//! #Id  Bulge Type  crRNA  DNA  Chromosome  Location  Direction  Mismatches  Bulge Size
//! ```
//!
//! Columns:
//!
//! * `Id` — query id from the input, or the query index.
//! * `Bulge Type` — always `X`.
//! * `crRNA` — guide plus `N`s for the PAM positions.
//! * `DNA` — guide-oriented reference guide + actual PAM.
//! * `Chromosome` — contig name.
//! * `Location` — **0-based** plus-strand left edge of the guide + PAM window
//!   (Cas-OFFinder/Bowtie convention).
//! * `Direction` — `+` or `-`.
//! * `Mismatches` — Hamming mismatches over the guide only.
//! * `Bulge Size` — always `0`.
//! * `Region`, `Gene` — appended only when `--annotation` is used. Region is
//!   one of promoter, exon, coding exon, 5′ UTR, 3′ UTR, intron, gene body or
//!   intergenic. Promoters span 500 bp to 5 kb upstream of a TSS.
//!
//! # Examples
//!
//! Cas9 against hg38:
//!
//! ```text
//! hg38
//! NNNNNNNNNNNNNNNNNNNNNGG
//! GAGCGCTGCTCAGATAGCGANNN 3 tp53
//! ```
//!
//! ```console
//! $ cas-offinder queries.txt C results.tsv
//! ```
//!
//! Cpf1 against a local FASTA, reading from stdin:
//!
//! ```text
//! genome.fa
//! TTTVNNNNNNNNNNNNNNNNNNNNNNN
//! NNNNTGACTTTGAATGGAGTCGTGAGCGCAAGAACGCT 2 cpf1_demo
//! ```
//!
//! # Differences from the original Cas-OFFinder
//!
//! * The device argument is ignored; the implementation is CPU/Rayon only.
//! * DNA and RNA bulges are not supported; use `0 0`.
//! * The PAM is matched exactly against the IUPAC pattern; mismatches are
//!   counted over the guide only.
//! * Results are capped at 200,000 sites per query; a warning is printed if
//!   the cap is reached.
//! * Reference FASTA is loaded fully into memory (hg38 ≈ 3.1 GiB).
//!
//! # Exit status
//!
//! `0` on success; `1` for invalid input, a missing reference, an I/O error,
//! or an invalid query. Errors are printed to standard error.
//!
//! For library-level details, see [`arachnia::offtarget`](../arachnia/offtarget/index.html)
//! and the rest of the [`arachnia`](../arachnia/index.html) documentation.

use std::env;
use std::fs;
use std::io::Read;
use std::process::ExitCode;
use std::time::Instant;

use arachnia::annotation::{Annotation, ANNOTATION_ENV};
use arachnia::offtarget::{
    expand_tilde, resolve_reference_path, results_tsv_header, results_tsv_header_with_annotation,
    results_tsv_rows, search_offtargets, validate_search, Genome, OffTargetParams, PamSide,
    SearchControl,
};

/// Run the CLI and convert errors into a process exit status.
fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(message) => {
            eprintln!("Error: {message}");
            ExitCode::FAILURE
        }
    }
}

/// Print usage and input-format help to standard error.
fn usage() {
    eprintln!(
        "Cas-OFFinder-compatible off-target search (Rust)\n\
         \n\
         Usage: cas-offinder [--annotation <GTF|GFF3>] <input.txt|-> [C|G|A[device_ids]] [output.tsv|-]\n\
         \n\
         Input format:\n\
           <hg38|genome.fa|directory>\n\
           <pattern> [dna_bulge] [rna_bulge]\n\
           <query> <max_mismatches> [id]\n\
         \n\
         --annotation adds Region and Gene columns. The path can also be set\n\
         with {ANNOTATION_ENV}.\n\
         \n\
         hg38/GRCh38 is resolved from ARACHNIA_HG38, ARACHNIA_HG38_DIR, or\n\
         standard data directories. Cas9 queries end with PAM wildcards\n\
         (guideNNN). Cpf1 queries start with them (NNNNguide). Guide positions\n\
         in the pattern must be N. Bulge sizes must be 0; the device argument\n\
         is accepted for compatibility."
    );
}

/// Parse arguments, load the reference, execute every query and write TSV.
///
/// Returns a human-readable error suitable for the `Error:` prefix printed by
/// [`main`]. Progress and diagnostics are written to standard error.
fn run() -> Result<(), String> {
    let arguments: Vec<String> = env::args().skip(1).collect();
    let mut annotation_path: Option<String> = env::var(ANNOTATION_ENV)
        .ok()
        .filter(|path| !path.trim().is_empty());
    let mut positional = Vec::new();
    let mut index = 0;
    while index < arguments.len() {
        let argument = arguments[index].as_str();
        match argument {
            "--help" | "-h" => {
                usage();
                return Ok(());
            }
            "--annotation" | "-a" => {
                index += 1;
                let path = arguments
                    .get(index)
                    .ok_or_else(|| "--annotation requires a GTF/GFF3 path".to_string())?;
                annotation_path = Some(path.clone());
            }
            _ if argument.starts_with("--annotation=") => {
                annotation_path = Some(argument["--annotation=".len()..].to_string());
            }
            _ => positional.push(arguments[index].clone()),
        }
        index += 1;
    }
    let mut args = positional;
    if args.is_empty() {
        usage();
        return Err("missing input file".into());
    }
    let input_path = args.remove(0);
    if !args.is_empty() && is_device_arg(&args[0]) {
        args.remove(0);
    }
    let output_path = args.first().cloned().unwrap_or_else(|| "-".to_string());
    if args.len() > 1 {
        usage();
        return Err("too many arguments".into());
    }

    let input = if input_path == "-" {
        let mut text = String::new();
        std::io::stdin()
            .read_to_string(&mut text)
            .map_err(|err| format!("could not read stdin: {err}"))?;
        text
    } else {
        fs::read_to_string(&input_path)
            .map_err(|err| format!("could not read {input_path}: {err}"))?
    };

    let mut lines = input
        .lines()
        .map(|line| line.trim_end_matches('\r').trim())
        .filter(|line| !line.is_empty() && !line.starts_with('#'));
    let genome_path = lines
        .next()
        .ok_or_else(|| "input is missing the genome path line".to_string())?
        .to_string();
    let pattern_line = lines
        .next()
        .ok_or_else(|| "input is missing the pattern line".to_string())?;
    let pattern_tokens: Vec<&str> = pattern_line.split_whitespace().collect();
    if pattern_tokens.is_empty() || !(pattern_tokens.len() == 1 || pattern_tokens.len() == 3) {
        return Err("pattern line must be: <pattern> [dna_bulge] [rna_bulge]".into());
    }
    if pattern_tokens.len() == 3 {
        let dna_bulge: usize = pattern_tokens[1]
            .parse()
            .map_err(|_| "DNA bulge size must be a number".to_string())?;
        let rna_bulge: usize = pattern_tokens[2]
            .parse()
            .map_err(|_| "RNA bulge size must be a number".to_string())?;
        if dna_bulge != 0 || rna_bulge != 0 {
            return Err(
                "DNA/RNA bulges are not supported by this Rust implementation; use 0 0".into(),
            );
        }
    }
    let pattern = pattern_tokens[0].to_ascii_uppercase();
    if !pattern.is_ascii() {
        return Err("pattern must be ASCII".into());
    }

    let genome_path = resolve_reference_path(&genome_path)
        .map_err(|err| format!("could not resolve reference {genome_path:?}: {err}"))?;
    let genome_path = expand_tilde(genome_path);
    eprintln!("Loading reference {} …", genome_path.display());
    let genome = Genome::from_path(&genome_path)
        .map_err(|err| format!("could not load {}: {err}", genome_path.display()))?;
    eprintln!(
        "Loaded {} contig(s), {} bp",
        genome.contig_count(),
        genome.total_bases()
    );

    let annotation = match annotation_path {
        Some(path) => {
            let path = expand_tilde(std::path::PathBuf::from(path));
            eprintln!("Loading annotation {} …", path.display());
            let annotation = Annotation::from_path(&path)
                .map_err(|err| format!("could not load {}: {err}", path.display()))?;
            eprintln!(
                "Loaded {} gene(s), {} transcript(s), {} feature(s)",
                annotation.gene_count(),
                annotation.transcript_count(),
                annotation.feature_count()
            );
            Some(annotation)
        }
        None => None,
    };

    let mut output = String::from(if annotation.is_some() {
        results_tsv_header_with_annotation()
    } else {
        results_tsv_header()
    });
    let mut query_count = 0usize;
    for (line_index, line) in lines.enumerate() {
        let tokens: Vec<&str> = line.split_whitespace().collect();
        if tokens.len() != 2 && tokens.len() != 3 {
            return Err(format!(
                "query line {} must be: <query> <max_mismatches> [id]",
                line_index + 1
            ));
        }
        let query = tokens[0].to_ascii_uppercase();
        if !query.is_ascii() {
            return Err("query must be ASCII".into());
        }
        let max_mismatches: usize = tokens[1]
            .parse()
            .map_err(|_| format!("invalid mismatch count on query line {}", line_index + 1))?;
        let id = tokens
            .get(2)
            .map(|id| id.to_string())
            .unwrap_or_else(|| query_count.to_string());

        if pattern.len() != query.len() {
            return Err(format!(
                "query line {}: query length ({}) must match pattern length ({}) when bulges are 0",
                line_index + 1,
                query.len(),
                pattern.len()
            ));
        }
        let leading_pam = query.bytes().take_while(|byte| *byte == b'N').count();
        let trailing_pam = query.bytes().rev().take_while(|byte| *byte == b'N').count();
        let (pam_side, pam_len, guide) = match (leading_pam, trailing_pam) {
            (0, 0) => {
                return Err(format!(
                    "query line {}: query must have PAM wildcards at one end (for example NNN)",
                    line_index + 1
                ))
            }
            (0, pam_len) => {
                if pam_len >= query.len() {
                    return Err(format!("query line {}: guide is empty", line_index + 1));
                }
                let guide_start = pattern.len() - pam_len;
                if !pattern[..guide_start].bytes().all(|byte| byte == b'N') {
                    return Err(format!(
                        "query line {}: only the 3' PAM region of the pattern may be constrained; guide positions must be N",
                        line_index + 1
                    ));
                }
                (PamSide::ThreePrime, pam_len, &query[..guide_start])
            }
            (pam_len, 0) => {
                if pam_len >= query.len() {
                    return Err(format!("query line {}: guide is empty", line_index + 1));
                }
                if !pattern[pam_len..].bytes().all(|byte| byte == b'N') {
                    return Err(format!(
                        "query line {}: only the 5' PAM region of the pattern may be constrained; guide positions must be N",
                        line_index + 1
                    ));
                }
                (PamSide::FivePrime, pam_len, &query[pam_len..])
            }
            _ => {
                return Err(format!(
                    "query line {}: query has PAM wildcards at both ends; use either guideNNN (3' PAM) or NNNNguide (5' PAM)",
                    line_index + 1
                ))
            }
        };
        let pam_pattern = match pam_side {
            PamSide::ThreePrime => &pattern[pattern.len() - pam_len..],
            PamSide::FivePrime => &pattern[..pam_len],
        };
        let params = OffTargetParams {
            max_mismatches,
            pam_pattern: pam_pattern.to_string(),
            pam_side,
            both_strands: true,
        };
        validate_search(guide, &params)
            .map_err(|err| format!("query line {}: {err}", line_index + 1))?;

        eprintln!(
            "[{}] {} … (PAM {} at {}, up to {} mismatches)",
            id,
            guide,
            pam_pattern,
            match pam_side {
                PamSide::ThreePrime => "3'",
                PamSide::FivePrime => "5'",
            },
            max_mismatches
        );
        let started = Instant::now();
        let mut result = search_offtargets(guide, &genome, &params, &SearchControl::new())
            .map_err(|err| format!("query line {}: {err}", line_index + 1))?;
        if let Some(annotation) = &annotation {
            result.annotate(annotation);
        }
        eprintln!(
            "[{}] {} site(s), {} perfect, {:.2} s",
            id,
            result.sites.len(),
            result.perfect_matches(),
            started.elapsed().as_secs_f64()
        );
        if result.truncated {
            eprintln!("[{}] warning: result cap reached; output was truncated", id);
        }
        output.push_str(&results_tsv_rows(
            &result.sites,
            guide,
            pam_len,
            pam_side,
            &id,
        ));
        query_count += 1;
    }

    if query_count == 0 {
        return Err("no queries found in the input file".into());
    }

    if output_path == "-" {
        print!("{output}");
    } else {
        fs::write(&output_path, output)
            .map_err(|err| format!("could not write {output_path}: {err}"))?;
        eprintln!("Wrote {output_path}");
    }
    Ok(())
}

/// Whether `value` looks like a Cas-OFFinder device argument.
///
/// Accepts a leading `C`, `G` or `A` followed only by digits, commas and
/// colons, e.g. `C`, `G`, `A0`, or `C0,1`. The value is accepted for
/// compatibility and otherwise ignored.
fn is_device_arg(value: &str) -> bool {
    let mut chars = value.chars();
    matches!(chars.next(), Some('C' | 'G' | 'A'))
        && chars.all(|ch| ch.is_ascii_digit() || ch == ',' || ch == ':')
}
