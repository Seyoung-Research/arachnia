//! Pure-Rust off-target enumeration for guide RNAs.
//!
//! This is an independent Rust implementation of the search model used by
//! [Cas-OFFinder](https://github.com/snugel/cas-offinder) (Bae, Park & Kim,
//! *Bioinformatics* 2014): find every genomic site whose guide region is
//! within a Hamming distance of the query and whose PAM matches an IUPAC
//! pattern exactly. Cas-OFFinder uses OpenCL; this implementation uses packed
//! 2-bit comparison plus Rayon and keeps the same location convention
//! (0-based, left edge of the guide + PAM window on the plus strand).
//!
//! # Search model
//!
//! * The guide is 8–31 nt; the PAM is 1–8 IUPAC positions.
//! * Mismatches are counted over the guide only. The PAM must match exactly.
//! * `NGG` (Cas9, 3′ PAM) and `TTTV` (Cpf1, 5′ PAM) are supported through
//!   [`PamSide`] and an arbitrary IUPAC `pam_pattern`.
//! * Both strands are searched by default. For a minus-strand hit, the
//!   reported guide and PAM are in guide orientation; the location is still
//!   the 0-based plus-strand left edge of the guide + PAM window.
//! * Ambiguous reference bases (`N` or IUPAC codes) never match and cause the
//!   window to be skipped, matching Cas-OFFinder’s treatment of unknown bases.
//! * Results are sorted by mismatch count, then contig, location and strand,
//!   and capped at 200,000 sites per search.
//!
//! **DNA and RNA bulges are not considered.** The PAM is always matched
//! exactly against the supplied IUPAC pattern; mismatches are counted over the
//! guide only, matching the usual Cas-OFFinder invocation where the query’s
//! PAM positions are `N`s.
//!
//! # Genome input
//!
//! [`Genome`] holds one or more [`GenomeRecord`]s loaded from memory or from
//! FASTA. [`Genome::from_path`] accepts a single FASTA file or a directory of
//! `.fa`, `.fasta`, `.fna`, `.fas`, `.ffn` and `.seq` files. Compressed FASTA
//! is not supported. A whole hg38 reference is about 3.1 GiB in memory and
//! should be loaded in a background thread in a GUI.
//!
//! # Local hg38 discovery
//!
//! [`find_hg38`] checks, in order:
//!
//! 1. `$ARACHNIA_HG38` (a FASTA file or directory);
//! 2. `$ARACHNIA_HG38_DIR/hg38.fa`;
//! 3. `$XDG_DATA_HOME/arachnia/hg38.fa`;
//! 4. `~/.local/share/arachnia/hg38.fa`;
//! 5. `~/genomes/hg38.fa`, `~/genomes/hg38.fasta`, `~/hg38.fa`;
//! 6. `/data/genomes/hg38.fa` and `/srv/genomes/hg38.fa`.
//!
//! The repository’s `scripts/fetch_hg38.sh` installs a verified UCSC hg38
//! assembly at the XDG location by default. [`is_hg38_alias`] and
//! [`resolve_reference_path`] implement the `hg38`/`GRCh38` shorthand used by
//! the command-line tool.
//!
//! # Progress and cancellation
//!
//! [`search_offtargets`] updates the atomics in [`SearchControl`] as it scans.
//! A GUI can poll [`SearchControl::progress`] from another thread and call
//! [`SearchControl::cancel`] to stop early; cancellation returns
//! [`SearchError::Cancelled`].
//!
//! # TSV output
//!
//! [`results_tsv`] and [`results_tsv_rows`] emit Cas-OFFinder-compatible rows:
//! `Id`, `Bulge Type`, `crRNA`, `DNA`, `Chromosome`, `Location`, `Direction`,
//! `Mismatches`, and `Bulge Size`. `Bulge Type` is always `X` and `Bulge Size`
//! is always `0` because bulges are unsupported. The PAM is written as `N`s in
//! the crRNA column, as Cas-OFFinder does.

use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::time::{Duration, Instant};

use rayon::prelude::*;

use crate::design::{reverse_complement, Strand};

/// Longest guide that fits a packed `u64` (2 bits per base).
const MAX_GUIDE_LEN: usize = 31;
const MAX_PAM_LEN: usize = 8;
/// Output cap; very large searches are truncated after sorting.
const MAX_RESULTS: usize = 200_000;
/// Number of start positions scanned per parallel work item.
const CHUNK_STARTS: usize = 1 << 22;
const PROGRESS_MASK: usize = (1 << 20) - 1;

/// A single FASTA record.
///
/// Bases are stored one byte per nucleotide after normalization: uppercase
/// A/C/G/T/N, with `U` converted to `T` and all other characters converted to
/// `N`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GenomeRecord {
    /// Record name. FASTA parsing stores the first whitespace-delimited token
    /// of the header; [`GenomeRecord::new`] stores the supplied name as-is.
    pub name: String,
    /// Normalized nucleotide sequence.
    pub bases: Vec<u8>,
}

impl GenomeRecord {
    /// Build a record from any byte sequence. Whitespace is ignored, `U` is
    /// treated as `T`, and anything outside A/C/G/T becomes `N`.
    pub fn new(name: impl Into<String>, bases: &[u8]) -> Self {
        let bases = bases
            .iter()
            .filter(|b| !b.is_ascii_whitespace())
            .map(|&b| normalize_base(b))
            .collect();
        Self {
            name: name.into(),
            bases,
        }
    }
}

/// One or more reference contigs.
///
/// A genome is a collection of independently named [`GenomeRecord`]s. The
/// off-target scanner treats each record separately, so contig boundaries are
/// never crossed.
#[derive(Clone, Debug, Default)]
pub struct Genome {
    /// Reference contigs in input order.
    pub records: Vec<GenomeRecord>,
}

impl Genome {
    /// Build a genome from records, discarding empty records.
    ///
    /// # Errors
    ///
    /// Returns [`GenomeError`] if no non-empty record remains.
    pub fn from_records(records: Vec<GenomeRecord>) -> Result<Self, GenomeError> {
        let records: Vec<GenomeRecord> = records
            .into_iter()
            .filter(|record| !record.bases.is_empty())
            .collect();
        if records.is_empty() {
            return Err(GenomeError("no sequence found in reference".into()));
        }
        Ok(Self { records })
    }

    /// Build a single-contig genome from a name and sequence.
    ///
    /// The sequence is normalized exactly like [`GenomeRecord::new`].
    ///
    /// # Errors
    ///
    /// Returns [`GenomeError`] if the normalized sequence is empty.
    pub fn from_sequence(name: impl Into<String>, bases: &[u8]) -> Result<Self, GenomeError> {
        Self::from_records(vec![GenomeRecord::new(name, bases)])
    }

    /// Parse FASTA text (used by tests and callers that already have the
    /// sequence in memory).
    pub fn from_fasta_str(text: &str) -> Result<Self, GenomeError> {
        let mut records = Vec::new();
        parse_fasta(text.as_bytes(), "sequence", &mut records)
            .map_err(|err| GenomeError(err.to_string()))?;
        Self::from_records(records)
    }

    /// Load a FASTA file, or every FASTA file in a directory.
    ///
    /// Directories are scanned for `.fa`, `.fasta`, `.fna`, `.fas`, `.ffn`
    /// and `.seq` files (case-insensitive), in sorted order. Compressed
    /// FASTA files are not supported.
    pub fn from_path(path: &Path) -> Result<Self, GenomeError> {
        let files = if path.is_dir() {
            let mut files = Vec::new();
            let mut saw_gzip = false;
            let entries = std::fs::read_dir(path).map_err(|err| {
                GenomeError(format!("cannot read directory {}: {err}", path.display()))
            })?;
            for entry in entries {
                let entry = entry.map_err(|err| GenomeError(err.to_string()))?;
                let file_path = entry.path();
                if !file_path.is_file() {
                    continue;
                }
                let ext = file_path
                    .extension()
                    .and_then(|ext| ext.to_str())
                    .unwrap_or("")
                    .to_ascii_lowercase();
                if ext == "gz" {
                    saw_gzip = true;
                } else if matches!(ext.as_str(), "fa" | "fasta" | "fna" | "fas" | "ffn" | "seq") {
                    files.push(file_path);
                }
            }
            files.sort();
            if files.is_empty() {
                if saw_gzip {
                    return Err(GenomeError(
                        "compressed FASTA is not supported; decompress the genome first".into(),
                    ));
                }
                return Err(GenomeError(format!(
                    "no FASTA files (.fa, .fasta, .fna, .fas, .ffn, .seq) in {}",
                    path.display()
                )));
            }
            files
        } else {
            if path
                .extension()
                .and_then(|ext| ext.to_str())
                .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"))
            {
                return Err(GenomeError(
                    "compressed FASTA is not supported; decompress the genome first".into(),
                ));
            }
            vec![path.to_path_buf()]
        };

        let mut records = Vec::new();
        for file_path in files {
            let fallback = file_path
                .file_stem()
                .and_then(|stem| stem.to_str())
                .unwrap_or("sequence")
                .to_string();
            let file = File::open(&file_path).map_err(|err| {
                GenomeError(format!("cannot open {}: {err}", file_path.display()))
            })?;
            let reader = BufReader::with_capacity(1 << 20, file);
            parse_fasta(reader, &fallback, &mut records)
                .map_err(|err| GenomeError(format!("{}: {err}", file_path.display())))?;
        }
        Self::from_records(records)
    }

    /// Total number of normalized bases across all contigs.
    pub fn total_bases(&self) -> usize {
        self.records.iter().map(|record| record.bases.len()).sum()
    }

    /// Number of reference contigs.
    pub fn contig_count(&self) -> usize {
        self.records.len()
    }
}

/// Environment variable that can point directly at an hg38 FASTA file or
/// directory.
pub const HG38_ENV: &str = "ARACHNIA_HG38";

/// Environment variable that can point at the directory containing
/// `hg38.fa` (used by `scripts/fetch_hg38.sh`).
pub const HG38_DIR_ENV: &str = "ARACHNIA_HG38_DIR";

/// Paths searched for a local hg38 reference, in priority order.
///
/// The first existing path wins. Set `ARACHNIA_HG38` to override the search.
pub fn hg38_candidates() -> Vec<PathBuf> {
    let mut candidates = Vec::new();
    if let Some(path) = std::env::var_os(HG38_ENV) {
        candidates.push(PathBuf::from(path));
    }
    if let Some(dir) = std::env::var_os(HG38_DIR_ENV) {
        candidates.push(PathBuf::from(dir).join("hg38.fa"));
    }
    if let Some(data_home) = std::env::var_os("XDG_DATA_HOME") {
        candidates.push(PathBuf::from(data_home).join("arachnia/hg38.fa"));
    }
    if let Some(home) = std::env::var_os("HOME") {
        let home = PathBuf::from(home);
        candidates.push(home.join(".local/share/arachnia/hg38.fa"));
        candidates.push(home.join("genomes/hg38.fa"));
        candidates.push(home.join("genomes/hg38.fasta"));
        candidates.push(home.join("hg38.fa"));
    }
    candidates.push(PathBuf::from("/data/genomes/hg38.fa"));
    candidates.push(PathBuf::from("/srv/genomes/hg38.fa"));
    candidates
}

/// Locate a local hg38 FASTA file or directory, if one is installed.
pub fn find_hg38() -> Option<PathBuf> {
    hg38_candidates()
        .into_iter()
        .find(|path| path.is_file() || path.is_dir())
}

/// Whether a reference argument is the `hg38` shorthand rather than a path.
pub fn is_hg38_alias(value: &str) -> bool {
    value.eq_ignore_ascii_case("hg38") || value.eq_ignore_ascii_case("grch38")
}

/// Resolve the `hg38`/`GRCh38` shorthand to a discovered local reference.
///
/// Any other value is returned as a path unchanged.
pub fn resolve_reference_path(value: &str) -> Result<PathBuf, String> {
    if is_hg38_alias(value) {
        find_hg38().ok_or_else(|| {
            format!(
                "hg38 was not found. Run scripts/fetch_hg38.sh, set {HG38_ENV}, or pass a FASTA path."
            )
        })
    } else {
        Ok(PathBuf::from(value))
    }
}

/// Expand a leading `~` using `$HOME`, leaving other paths unchanged.
pub fn expand_tilde(path: PathBuf) -> PathBuf {
    if let Ok(rest) = path.strip_prefix("~") {
        if let Some(home) = std::env::var_os("HOME") {
            return Path::new(&home).join(rest);
        }
    }
    path
}

/// Error returned while loading or validating a reference genome.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GenomeError(
    /// Human-readable error message.
    pub String,
);

impl fmt::Display for GenomeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for GenomeError {}

fn normalize_base(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'A',
        b'C' => b'C',
        b'G' => b'G',
        b'T' | b'U' => b'T',
        _ => b'N',
    }
}

fn parse_fasta<R: BufRead>(
    mut reader: R,
    fallback_name: &str,
    records: &mut Vec<GenomeRecord>,
) -> Result<(), GenomeError> {
    let mut current_name: Option<String> = None;
    let mut current_bases: Vec<u8> = Vec::new();
    let mut saw_header = false;
    let mut line = String::new();

    let flush =
        |name: &mut Option<String>, bases: &mut Vec<u8>, records: &mut Vec<GenomeRecord>| {
            if bases.is_empty() {
                *name = None;
                return;
            }
            let name = name.take().unwrap_or_else(|| fallback_name.to_string());
            records.push(GenomeRecord {
                name,
                bases: std::mem::take(bases),
            });
        };

    loop {
        line.clear();
        let bytes = reader
            .read_line(&mut line)
            .map_err(|err| GenomeError(err.to_string()))?;
        if bytes == 0 {
            break;
        }
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if let Some(header) = trimmed.strip_prefix('>') {
            flush(&mut current_name, &mut current_bases, records);
            let name = header
                .split_whitespace()
                .next()
                .filter(|name| !name.is_empty())
                .unwrap_or(fallback_name)
                .to_string();
            current_name = Some(name);
            saw_header = true;
        } else {
            for byte in trimmed.bytes() {
                if byte.is_ascii_whitespace() {
                    continue;
                }
                current_bases.push(normalize_base(byte));
            }
        }
    }
    if saw_header || !current_bases.is_empty() {
        flush(&mut current_name, &mut current_bases, records);
    }
    Ok(())
}

/// Shared progress and cancellation state for a running search.
///
/// Create one per search, pass it to [`search_offtargets`], and share it with a
/// UI thread. All fields are atomics and can be read while the search runs.
#[derive(Debug, Default)]
pub struct SearchControl {
    /// Set to `true` to request cancellation.
    pub cancel: AtomicBool,
    /// Number of candidate windows scanned so far.
    pub scanned: AtomicU64,
    /// Total number of candidate windows for this search (both strands).
    pub total: AtomicU64,
}

impl SearchControl {
    /// Create a fresh, uncancelled control block.
    pub fn new() -> Self {
        Self::default()
    }

    /// Request cancellation. The search returns [`SearchError::Cancelled`] at
    /// its next cancellation check.
    pub fn cancel(&self) {
        self.cancel.store(true, Ordering::Relaxed);
    }

    /// Whether cancellation has been requested.
    pub fn is_cancelled(&self) -> bool {
        self.cancel.load(Ordering::Relaxed)
    }

    /// Fraction of candidate windows scanned, clamped to `0.0..=1.0`.
    ///
    /// Returns `0.0` until [`search_offtargets`] has initialized
    /// [`SearchControl::total`].
    pub fn progress(&self) -> f32 {
        let total = self.total.load(Ordering::Relaxed);
        if total == 0 {
            return 0.0;
        }
        (self.scanned.load(Ordering::Relaxed) as f64 / total as f64).clamp(0.0, 1.0) as f32
    }
}

/// Error returned by off-target validation or search.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SearchError {
    /// The guide, PAM pattern, or mismatch limit is invalid.
    Invalid(String),
    /// The caller requested cancellation through [`SearchControl::cancel`].
    Cancelled,
}

impl fmt::Display for SearchError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SearchError::Invalid(message) => f.write_str(message),
            SearchError::Cancelled => f.write_str("search cancelled"),
        }
    }
}

impl std::error::Error for SearchError {}

/// Where the PAM sits relative to the guide in guide orientation.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PamSide {
    /// Cas9-style 3' PAM: `[guide][PAM]`.
    #[default]
    ThreePrime,
    /// Cpf1/Cas12a-style 5' PAM: `[PAM][guide]`.
    FivePrime,
}

/// Parameters controlling a Cas-OFFinder-style off-target search.
///
/// The defaults are Cas9-like: up to three guide mismatches, exact `NGG` 3′
/// PAM, and both strands. Use [`PamSide::FivePrime`] and a `TTTV` pattern for
/// Cpf1.
#[derive(Clone, Debug)]
pub struct OffTargetParams {
    /// Maximum Hamming mismatches over the guide (PAM excluded).
    pub max_mismatches: usize,
    /// IUPAC PAM pattern in guide orientation, e.g. `NGG` or `TTTV`.
    pub pam_pattern: String,
    /// Whether the PAM is 3' (Cas9) or 5' (Cpf1) of the guide.
    pub pam_side: PamSide,
    /// Search the reverse-complement strand as well as the plus strand.
    pub both_strands: bool,
}

impl Default for OffTargetParams {
    fn default() -> Self {
        Self {
            max_mismatches: 3,
            pam_pattern: "NGG".to_string(),
            pam_side: PamSide::ThreePrime,
            both_strands: true,
        }
    }
}

/// A single off-target site.
///
/// `location` uses the Cas-OFFinder/Bowtie convention: 0-based plus-strand
/// coordinate of the left edge of the guide + PAM window. For a minus-strand
/// hit, `guide_dna` and `pam` are reverse-complemented into guide orientation,
/// but `location` still refers to the plus-strand window.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OffTargetSite {
    /// Reference contig name.
    pub contig: String,
    /// 0-based plus-strand coordinate of the left edge of the guide + PAM
    /// window (Cas-OFFinder / Bowtie convention).
    pub location: usize,
    /// Strand on which the hit was found.
    pub strand: Strand,
    /// Number of guide mismatches (the PAM always matches exactly).
    pub mismatches: usize,
    /// Guide-oriented reference sequence; mismatched guide bases are shown
    /// in lower case.
    pub guide_dna: String,
    /// Actual PAM in guide orientation.
    pub pam: String,
}

/// Result of an off-target search.
#[derive(Clone, Debug)]
pub struct OffTargetResult {
    /// Reported sites, sorted by mismatch count, contig, location and strand.
    pub sites: Vec<OffTargetSite>,
    /// Number of reported sites at each mismatch count (index = mismatches).
    pub counts_by_mismatch: Vec<usize>,
    /// Total reference bases loaded for the search.
    pub genome_bases: usize,
    /// Candidate windows scanned (both strands if enabled).
    pub scanned_positions: u64,
    /// Wall-clock search time, including parallel scanning and sorting.
    pub elapsed: Duration,
    /// Whether the 200,000-site output cap was reached and results truncated.
    pub truncated: bool,
}

impl OffTargetResult {
    /// Number of perfect (zero-mismatch) sites in the reported result.
    ///
    /// This counts only sites that survived the output cap; use
    /// [`OffTargetResult::truncated`] to detect a capped search.
    pub fn perfect_matches(&self) -> usize {
        self.counts_by_mismatch.first().copied().unwrap_or(0)
    }
}

fn iupac_mask(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => 0b0001,
        b'C' => 0b0010,
        b'G' => 0b0100,
        b'T' | b'U' => 0b1000,
        b'R' => 0b0101,
        b'Y' => 0b1010,
        b'S' => 0b0110,
        b'W' => 0b1001,
        b'K' => 0b1100,
        b'M' => 0b0011,
        b'B' => 0b1110,
        b'D' => 0b1101,
        b'H' => 0b1011,
        b'V' => 0b0111,
        b'N' => 0b1111,
        _ => 0,
    }
}

fn complement_iupac(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' | b'U' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        _ => b'N',
    }
}

fn reverse_complement_iupac(pattern: &[u8]) -> Vec<u8> {
    pattern
        .iter()
        .rev()
        .map(|&base| complement_iupac(base))
        .collect()
}

fn normalize_guide(guide: &str) -> Result<Vec<u8>, SearchError> {
    let mut bases = Vec::new();
    for ch in guide.chars() {
        if ch.is_whitespace() {
            continue;
        }
        let base = match ch.to_ascii_uppercase() {
            'A' => b'A',
            'C' => b'C',
            'G' => b'G',
            'T' | 'U' => b'T',
            _ => {
                return Err(SearchError::Invalid(format!(
                    "guide contains invalid base '{ch}' (only A, C, G, T, U are accepted)"
                )))
            }
        };
        bases.push(base);
    }
    if bases.len() < 8 || bases.len() > MAX_GUIDE_LEN {
        return Err(SearchError::Invalid(format!(
            "guide length must be 8-{MAX_GUIDE_LEN} nt, got {}",
            bases.len()
        )));
    }
    Ok(bases)
}

fn normalize_pam(pattern: &str) -> Result<Vec<u8>, SearchError> {
    let mut bases = Vec::new();
    for ch in pattern.chars() {
        if ch.is_whitespace() {
            continue;
        }
        if iupac_mask(ch as u8) == 0 {
            return Err(SearchError::Invalid(format!(
                "PAM pattern contains invalid IUPAC base '{ch}'"
            )));
        }
        bases.push(ch.to_ascii_uppercase() as u8);
    }
    if bases.is_empty() || bases.len() > MAX_PAM_LEN {
        return Err(SearchError::Invalid(format!(
            "PAM pattern length must be 1-{MAX_PAM_LEN}, got {}",
            bases.len()
        )));
    }
    if bases.iter().all(|&base| base == b'N') {
        return Err(SearchError::Invalid(
            "PAM pattern must constrain at least one position (all-N patterns are not supported)"
                .into(),
        ));
    }
    Ok(bases)
}

/// Validate a query before a (potentially expensive) genome scan.
///
/// Checks that the guide contains only A/C/G/T/U and is 8–31 nt, that the PAM
/// is a valid IUPAC pattern of 1–8 positions with at least one constrained
/// position, and that `max_mismatches <= guide.len()`.
///
/// # Errors
///
/// Returns [`SearchError::Invalid`] with a human-readable explanation.
pub fn validate_search(guide: &str, params: &OffTargetParams) -> Result<(), SearchError> {
    let guide = normalize_guide(guide)?;
    normalize_pam(&params.pam_pattern)?;
    if params.max_mismatches > guide.len() {
        return Err(SearchError::Invalid(format!(
            "maximum mismatches ({}) cannot exceed the guide length ({})",
            params.max_mismatches,
            guide.len()
        )));
    }
    Ok(())
}

#[inline]
fn base_code(base: u8) -> u64 {
    match base {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 0,
    }
}

#[inline]
fn is_acgt(base: u8) -> bool {
    matches!(base, b'A' | b'C' | b'G' | b'T')
}

fn pack(seq: &[u8]) -> u64 {
    let mut packed = 0u64;
    for (index, &base) in seq.iter().enumerate() {
        packed |= base_code(base) << (2 * index);
    }
    packed
}

/// Count positions where two 2-bit packed sequences differ. Ambiguous bases
/// never reach this function because windows containing them are skipped.
#[inline]
fn packed_mismatches(a: u64, b: u64) -> usize {
    let diff = a ^ b;
    ((diff | (diff >> 1)) & 0x5555_5555_5555_5555).count_ones() as usize
}

#[inline]
fn pam_matches(seq: &[u8], masks: &[u8]) -> bool {
    for (index, &base) in seq.iter().enumerate() {
        if iupac_mask(base) & masks[index] == 0 {
            return false;
        }
    }
    true
}

fn mark_mismatches(guide: &[u8], dna: &[u8]) -> String {
    let mut marked = String::with_capacity(dna.len());
    for (&expected, &actual) in guide.iter().zip(dna.iter()) {
        let base = if expected == actual {
            actual.to_ascii_uppercase()
        } else {
            actual.to_ascii_lowercase()
        };
        marked.push(base as char);
    }
    marked
}

/// Search a genome for off-target sites of `guide`.
///
/// The guide may contain U (treated as T) and whitespace (ignored). The PAM
/// pattern is IUPAC and is matched exactly against the reference. Mismatches
/// are counted over the guide only. With [`OffTargetParams::both_strands`]
/// enabled, the reverse-complement strand is searched too.
///
/// The function updates `control` as it scans and checks for cancellation
/// between chunks. Results are sorted by mismatch count, contig, location and
/// strand, and capped at 200,000 sites.
///
/// # Errors
///
/// Returns [`SearchError::Invalid`] for an invalid guide, PAM or mismatch
/// limit, and [`SearchError::Cancelled`] if cancellation is requested.
///
/// # Example
///
/// ```
/// use arachnia::offtarget::{
///     search_offtargets, Genome, OffTargetParams, SearchControl,
/// };
///
/// let genome = Genome::from_fasta_str(">chr1\nGGGGACGTACGTACGTACGTACGTTGGAAAA\n").unwrap();
/// let result = search_offtargets(
///     "ACGTACGTACGTACGTACGT",
///     &genome,
///     &OffTargetParams::default(),
///     &SearchControl::new(),
/// )
/// .unwrap();
/// assert_eq!(result.perfect_matches(), 1);
/// ```
pub fn search_offtargets(
    guide: &str,
    genome: &Genome,
    params: &OffTargetParams,
    control: &SearchControl,
) -> Result<OffTargetResult, SearchError> {
    let started = Instant::now();
    let guide_bases = normalize_guide(guide)?;
    let pam_pattern = normalize_pam(&params.pam_pattern)?;
    if params.max_mismatches > guide_bases.len() {
        return Err(SearchError::Invalid(format!(
            "maximum mismatches ({}) cannot exceed the guide length ({})",
            params.max_mismatches,
            guide_bases.len()
        )));
    }

    let guide_len = guide_bases.len();
    let pam_len = pam_pattern.len();
    let window_len = guide_len + pam_len;
    let max_mismatches = params.max_mismatches;
    let pam_side = params.pam_side;
    let both_strands = params.both_strands;

    let guide_packed = pack(&guide_bases);
    let guide_rc = reverse_complement(&guide_bases);
    let guide_rc_packed = pack(&guide_rc);
    let pam_plus = pam_pattern
        .iter()
        .map(|&base| iupac_mask(base))
        .collect::<Vec<u8>>();
    let pam_minus = reverse_complement_iupac(&pam_pattern)
        .iter()
        .map(|&base| iupac_mask(base))
        .collect::<Vec<u8>>();

    let strand_count = if both_strands { 2 } else { 1 } as u64;
    let total: u64 = genome
        .records
        .iter()
        .map(|record| {
            if record.bases.len() < window_len {
                0
            } else {
                (record.bases.len() - window_len + 1) as u64
            }
        })
        .sum::<u64>()
        * strand_count;
    control.total.store(total, Ordering::Relaxed);
    control.scanned.store(0, Ordering::Relaxed);

    if control.is_cancelled() {
        return Err(SearchError::Cancelled);
    }

    let mut jobs: Vec<(usize, bool, usize, usize)> = Vec::new();
    for (record_index, record) in genome.records.iter().enumerate() {
        if record.bases.len() < window_len {
            continue;
        }
        let positions = record.bases.len() - window_len + 1;
        let mut start = 0;
        while start < positions {
            let end = (start + CHUNK_STARTS).min(positions);
            jobs.push((record_index, false, start, end));
            if both_strands {
                jobs.push((record_index, true, start, end));
            }
            start = end;
        }
    }

    let chunks: Vec<Vec<OffTargetSite>> = jobs
        .par_iter()
        .map(|&(record_index, minus, start, end)| {
            scan_range(
                &genome.records[record_index],
                minus,
                start,
                end,
                guide_len,
                pam_len,
                pam_side,
                &guide_bases,
                guide_packed,
                guide_rc_packed,
                &pam_plus,
                &pam_minus,
                max_mismatches,
                control,
            )
        })
        .collect();

    if control.is_cancelled() {
        return Err(SearchError::Cancelled);
    }

    let mut sites: Vec<OffTargetSite> = chunks.into_iter().flatten().collect();
    sites.sort_by(|a, b| {
        a.mismatches
            .cmp(&b.mismatches)
            .then_with(|| a.contig.cmp(&b.contig))
            .then_with(|| a.location.cmp(&b.location))
            .then_with(|| a.strand.symbol().cmp(&b.strand.symbol()))
    });

    let truncated = sites.len() > MAX_RESULTS;
    if truncated {
        sites.truncate(MAX_RESULTS);
    }

    let mut counts_by_mismatch = vec![0usize; max_mismatches + 1];
    for site in &sites {
        counts_by_mismatch[site.mismatches] += 1;
    }

    Ok(OffTargetResult {
        sites,
        counts_by_mismatch,
        genome_bases: genome.total_bases(),
        scanned_positions: control.scanned.load(Ordering::Relaxed),
        elapsed: started.elapsed(),
        truncated,
    })
}

#[allow(clippy::too_many_arguments)]
fn scan_range(
    record: &GenomeRecord,
    minus: bool,
    start: usize,
    end: usize,
    guide_len: usize,
    pam_len: usize,
    pam_side: PamSide,
    guide_bases: &[u8],
    guide_packed: u64,
    guide_rc_packed: u64,
    pam_plus: &[u8],
    pam_minus: &[u8],
    max_mismatches: usize,
    control: &SearchControl,
) -> Vec<OffTargetSite> {
    if start >= end {
        return Vec::new();
    }
    let bases = &record.bases;
    let (target_packed, pam_masks) = if minus {
        (guide_rc_packed, pam_minus)
    } else {
        (guide_packed, pam_plus)
    };
    // On the plus strand the window is in guide orientation. On the minus
    // strand it is the reverse complement of the guide-oriented window.
    let (guide_offset, pam_offset) = match (minus, pam_side) {
        (false, PamSide::ThreePrime) | (true, PamSide::FivePrime) => (0, guide_len),
        (false, PamSide::FivePrime) | (true, PamSide::ThreePrime) => (pam_len, 0),
    };

    let mut window = pack(&bases[start + guide_offset..start + guide_offset + guide_len]);
    let mut ambiguous = bases[start + guide_offset..start + guide_offset + guide_len]
        .iter()
        .filter(|&&base| !is_acgt(base))
        .count();
    let mut sites = Vec::new();

    for position in start..end {
        if (position - start) & PROGRESS_MASK == 0 && control.is_cancelled() {
            control
                .scanned
                .fetch_add((position - start) as u64, Ordering::Relaxed);
            return sites;
        }

        if pam_matches(
            &bases[position + pam_offset..position + pam_offset + pam_len],
            pam_masks,
        ) && ambiguous == 0
        {
            let mismatches = packed_mismatches(window, target_packed);
            if mismatches <= max_mismatches {
                let dna_slice =
                    &bases[position + guide_offset..position + guide_offset + guide_len];
                let pam_slice = &bases[position + pam_offset..position + pam_offset + pam_len];
                let (dna, pam) = if minus {
                    (reverse_complement(dna_slice), reverse_complement(pam_slice))
                } else {
                    (dna_slice.to_vec(), pam_slice.to_vec())
                };
                sites.push(OffTargetSite {
                    contig: record.name.clone(),
                    location: position,
                    strand: if minus { Strand::Minus } else { Strand::Plus },
                    mismatches,
                    guide_dna: mark_mismatches(guide_bases, &dna),
                    pam: String::from_utf8_lossy(&pam).into_owned(),
                });
            }
        }

        if position + 1 < end {
            let leaving = bases[position + guide_offset];
            let entering = bases[position + guide_offset + guide_len];
            window >>= 2;
            window |= base_code(entering) << (2 * (guide_len - 1));
            if !is_acgt(leaving) {
                ambiguous -= 1;
            }
            if !is_acgt(entering) {
                ambiguous += 1;
            }
        }
    }

    control
        .scanned
        .fetch_add((end - start) as u64, Ordering::Relaxed);
    sites
}

/// The Cas-OFFinder-compatible TSV header (including its two comment lines).
pub fn results_tsv_header() -> &'static str {
    "##Generated by Arachnia (Cas-OFFinder-compatible off-target search)\n\
     #Id\tBulge Type\tcrRNA\tDNA\tChromosome\tLocation\tDirection\tMismatches\tBulge Size\n"
}

/// Render sites in the tab-separated Cas-OFFinder output format.
///
/// `guide` is the target guide (guide orientation), `pam_len` is the number
/// of PAM positions shown as `N` in the crRNA column, `pam_side` determines
/// whether those wildcards precede or follow the guide, and `id` is written to
/// the `Id` column. The result includes [`results_tsv_header`].
///
/// Bulge columns are always `X` and `0` because bulges are unsupported.
pub fn results_tsv(
    sites: &[OffTargetSite],
    guide: &str,
    pam_len: usize,
    pam_side: PamSide,
    id: &str,
) -> String {
    let mut out = String::from(results_tsv_header());
    out.push_str(&results_tsv_rows(sites, guide, pam_len, pam_side, id));
    out
}

/// Rows only (no header), for callers that batch several queries into one
/// Cas-OFFinder-style output file.
///
/// The crRNA column uses `N` for every PAM position, while the DNA column
/// contains the actual guide-oriented reference sequence and PAM. Tabs and
/// newlines in `id` and contig names are replaced with spaces.
pub fn results_tsv_rows(
    sites: &[OffTargetSite],
    guide: &str,
    pam_len: usize,
    pam_side: PamSide,
    id: &str,
) -> String {
    let pam_wildcards = "N".repeat(pam_len);
    let crna = match pam_side {
        PamSide::ThreePrime => format!("{guide}{pam_wildcards}"),
        PamSide::FivePrime => format!("{pam_wildcards}{guide}"),
    };
    let id = id.replace(['\t', '\n', '\r'], " ");
    let mut out = String::new();
    for site in sites {
        let contig = site.contig.replace(['\t', '\n', '\r'], " ");
        let dna = match pam_side {
            PamSide::ThreePrime => format!("{}{}", site.guide_dna, site.pam),
            PamSide::FivePrime => format!("{}{}", site.pam, site.guide_dna),
        };
        out.push_str(&format!(
            "{id}\tX\t{crna}\t{dna}\t{contig}\t{}\t{}\t{}\t0\n",
            site.location,
            site.strand.symbol(),
            site.mismatches
        ));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hg38_shorthand_and_plain_paths_resolve() {
        assert!(is_hg38_alias("hg38"));
        assert!(is_hg38_alias("GRCh38"));
        assert!(!is_hg38_alias("hg38.fa"));
        assert!(!is_hg38_alias("hg19"));
        assert_eq!(
            resolve_reference_path("custom.fa").unwrap(),
            PathBuf::from("custom.fa")
        );
        assert!(!hg38_candidates().is_empty());
    }

    fn pattern_masks(pattern: &[u8]) -> Vec<u8> {
        pattern.iter().map(|&b| iupac_mask(b)).collect()
    }

    /// Straightforward reference implementation used to validate the packed
    /// scanner.
    fn naive_sites(
        genome: &Genome,
        guide: &[u8],
        pam: &[u8],
        pam_side: PamSide,
        max_mismatches: usize,
        both_strands: bool,
    ) -> Vec<OffTargetSite> {
        let mut sites = Vec::new();
        let masks = pattern_masks(pam);
        let window_len = guide.len() + pam.len();
        for record in &genome.records {
            if record.bases.len() < window_len {
                continue;
            }
            for position in 0..=record.bases.len() - window_len {
                let plus_window = &record.bases[position..position + window_len];
                for minus in [false, true] {
                    if minus && !both_strands {
                        continue;
                    }
                    let window = if minus {
                        reverse_complement(plus_window)
                    } else {
                        plus_window.to_vec()
                    };
                    let (guide_window, pam_window) = match pam_side {
                        PamSide::ThreePrime => (&window[..guide.len()], &window[guide.len()..]),
                        PamSide::FivePrime => (&window[pam.len()..], &window[..pam.len()]),
                    };
                    if !pam_matches(pam_window, &masks) {
                        continue;
                    }
                    if guide_window.iter().any(|&b| !is_acgt(b)) {
                        continue;
                    }
                    let mismatches = guide
                        .iter()
                        .zip(guide_window.iter())
                        .filter(|(a, b)| a != b)
                        .count();
                    if mismatches <= max_mismatches {
                        sites.push(OffTargetSite {
                            contig: record.name.clone(),
                            location: position,
                            strand: if minus { Strand::Minus } else { Strand::Plus },
                            mismatches,
                            guide_dna: mark_mismatches(guide, guide_window),
                            pam: String::from_utf8(pam_window.to_vec()).unwrap(),
                        });
                    }
                }
            }
        }
        sites.sort_by(|a, b| {
            a.mismatches
                .cmp(&b.mismatches)
                .then_with(|| a.contig.cmp(&b.contig))
                .then_with(|| a.location.cmp(&b.location))
                .then_with(|| a.strand.symbol().cmp(&b.strand.symbol()))
        });
        sites
    }

    struct Lcg(u64);

    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            self.0 >> 33
        }

        fn base(&mut self) -> u8 {
            b"ACGT"[(self.next() & 3) as usize]
        }

        fn sequence(&mut self, len: usize) -> Vec<u8> {
            (0..len).map(|_| self.base()).collect()
        }
    }

    #[test]
    fn fasta_parsing_handles_multiple_records() {
        let genome = Genome::from_fasta_str(">chr1 human\nACGTacguNN\n>chr2\nTTTT\n").unwrap();
        assert_eq!(genome.contig_count(), 2);
        assert_eq!(genome.records[0].name, "chr1");
        assert_eq!(genome.records[0].bases, b"ACGTACGTNN");
        assert_eq!(genome.records[1].name, "chr2");
        assert_eq!(genome.records[1].bases, b"TTTT");
        assert_eq!(genome.total_bases(), 14);
    }

    #[test]
    fn empty_reference_is_rejected() {
        assert!(Genome::from_fasta_str(">chr1\n\n").is_err());
        assert!(Genome::from_records(Vec::new()).is_err());
    }

    #[test]
    fn iupac_pam_matching() {
        let ngg = pattern_masks(b"NGG");
        assert!(pam_matches(b"AGG", &ngg));
        assert!(pam_matches(b"TGG", &ngg));
        assert!(!pam_matches(b"NAG", &ngg));
        assert!(!pam_matches(b"GAG", &ngg));

        let tttv = pattern_masks(b"TTTV");
        assert!(pam_matches(b"TTTA", &tttv));
        assert!(pam_matches(b"TTTC", &tttv));
        assert!(pam_matches(b"TTTG", &tttv));
        assert!(!pam_matches(b"TTTT", &tttv));
        assert!(!pam_matches(b"CTTG", &tttv));
    }

    #[test]
    fn plus_and_minus_sites_are_reported_in_guide_orientation() {
        let guide = b"GTCGCTTAGTGTCCCTGG";
        let mut guide_mut = guide.to_vec();
        guide_mut[5] = if guide_mut[5] == b'A' { b'C' } else { b'A' };

        let mut plus = b"GCGCGCGCGC".to_vec();
        let plus_site = plus.len();
        plus.extend_from_slice(guide);
        plus.extend_from_slice(b"TGG");
        plus.extend_from_slice(b"ATATATATAT");

        let minus_orientation = [guide_mut.clone(), b"CGG".to_vec()].concat();
        let minus_plus = reverse_complement(&minus_orientation);
        let minus_site = plus.len();
        plus.extend_from_slice(&minus_plus);
        plus.extend_from_slice(b"GAGAGAGAGA");

        let genome = Genome::from_sequence("chrT", &plus).unwrap();
        let params = OffTargetParams {
            max_mismatches: 1,
            pam_pattern: "NGG".into(),
            pam_side: PamSide::ThreePrime,
            both_strands: true,
        };
        let result = search_offtargets(
            std::str::from_utf8(guide).unwrap(),
            &genome,
            &params,
            &SearchControl::new(),
        )
        .unwrap();

        let plus_hit = result
            .sites
            .iter()
            .find(|site| site.location == plus_site && site.strand == Strand::Plus)
            .unwrap();
        assert_eq!(plus_hit.mismatches, 0);
        assert_eq!(plus_hit.guide_dna, std::str::from_utf8(guide).unwrap());
        assert_eq!(plus_hit.pam, "TGG");

        let mut expected = String::from_utf8(guide_mut.clone()).unwrap();
        let lower = expected[5..6].to_ascii_lowercase();
        expected.replace_range(5..6, &lower);
        let minus_hit = result
            .sites
            .iter()
            .find(|site| site.location == minus_site && site.strand == Strand::Minus)
            .unwrap();
        assert_eq!(minus_hit.mismatches, 1);
        assert_eq!(minus_hit.guide_dna, expected);
        assert_eq!(minus_hit.pam, "CGG");
        assert_eq!(result.counts_by_mismatch.len(), 2);
        assert_eq!(result.perfect_matches(), result.counts_by_mismatch[0]);
    }

    #[test]
    fn ambiguous_windows_are_skipped() {
        let guide = b"ACGTACGTACGTACGTACGT";
        let mut sequence = b"TTTTTTTTTT".to_vec();
        sequence.extend_from_slice(guide);
        sequence.extend_from_slice(b"TGG");
        sequence.extend_from_slice(b"AAAAAAAAAA");
        let genome = Genome::from_sequence("chrN", &sequence).unwrap();
        let params = OffTargetParams {
            max_mismatches: 0,
            pam_pattern: "NGG".into(),
            pam_side: PamSide::ThreePrime,
            both_strands: false,
        };
        let clean = search_offtargets(
            std::str::from_utf8(guide).unwrap(),
            &genome,
            &params,
            &SearchControl::new(),
        )
        .unwrap();
        assert!(clean.sites.iter().any(|site| site.location == 10));

        let mut with_n = sequence;
        with_n[12] = b'N';
        let genome = Genome::from_sequence("chrN", &with_n).unwrap();
        let ambiguous = search_offtargets(
            std::str::from_utf8(guide).unwrap(),
            &genome,
            &params,
            &SearchControl::new(),
        )
        .unwrap();
        assert!(!ambiguous.sites.iter().any(|site| site.location == 10));
    }

    #[test]
    fn packed_scanner_matches_naive_reference() {
        let mut rng = Lcg(0x5eed_1234_9876_abcd);
        let cases = [
            ("NGG", b"TGG".as_slice(), PamSide::ThreePrime),
            ("NAG", b"TAG".as_slice(), PamSide::ThreePrime),
            ("TTTV", b"TTTA".as_slice(), PamSide::FivePrime),
            ("TTTV", b"TTTG".as_slice(), PamSide::FivePrime),
        ];
        for case in 0..12 {
            let (pattern, pam, pam_side) = cases[case % cases.len()];
            let guide = rng.sequence(20);
            let window_len = guide.len() + pam.len();
            let mut mutated = guide.clone();
            mutated[3] = match mutated[3] {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                _ => b'A',
            };
            let oriented = |guide: &[u8], pam: &[u8]| match pam_side {
                PamSide::ThreePrime => [guide, pam].concat(),
                PamSide::FivePrime => [pam, guide].concat(),
            };
            let mut genome_bases = rng.sequence(520);

            // Plus-strand target, exact.
            let plus_site = 40;
            genome_bases[plus_site..plus_site + window_len].copy_from_slice(&oriented(&guide, pam));
            // Minus-strand target with one mismatch.
            let minus_window = reverse_complement(&oriented(&mutated, pam));
            let minus_site = 240;
            genome_bases[minus_site..minus_site + window_len].copy_from_slice(&minus_window);
            // Second plus-strand target, possibly with two mismatches.
            let mut mutated2 = mutated.clone();
            mutated2[11] = b'A';
            mutated2[15] = b'T';
            let mut pam2 = pam.to_vec();
            pam2[0] = b'T';
            let plus2_site = 400;
            genome_bases[plus2_site..plus2_site + window_len]
                .copy_from_slice(&oriented(&mutated2, &pam2));

            let genome = Genome::from_records(vec![
                GenomeRecord::new("chr1", &genome_bases),
                GenomeRecord::new("chr2", &rng.sequence(180)),
            ])
            .unwrap();

            for max_mismatches in 0..=3 {
                for both_strands in [false, true] {
                    let params = OffTargetParams {
                        max_mismatches,
                        pam_pattern: pattern.into(),
                        pam_side,
                        both_strands,
                    };
                    let result = search_offtargets(
                        std::str::from_utf8(&guide).unwrap(),
                        &genome,
                        &params,
                        &SearchControl::new(),
                    )
                    .unwrap();
                    let expected = naive_sites(
                        &genome,
                        &guide,
                        pattern.as_bytes(),
                        pam_side,
                        max_mismatches,
                        both_strands,
                    );
                    assert_eq!(
                        result.sites, expected,
                        "case {case}, pattern {pattern}, side {pam_side:?}, max {max_mismatches}, both {both_strands}"
                    );
                    assert_eq!(
                        result.sites.len(),
                        result.counts_by_mismatch.iter().sum::<usize>()
                    );
                }
            }
        }
    }

    #[test]
    fn tsv_uses_cas_offinder_columns() {
        let guide = "ACGTACGTACGTACGTACGT";
        let sequence = [
            b"TTTTTTTTTT".as_slice(),
            guide.as_bytes(),
            b"TGG",
            b"AAAAAAAAAA",
        ]
        .concat();
        let genome = Genome::from_sequence("chr1", &sequence).unwrap();
        let result = search_offtargets(
            guide,
            &genome,
            &OffTargetParams {
                max_mismatches: 0,
                pam_pattern: "NGG".into(),
                pam_side: PamSide::ThreePrime,
                both_strands: false,
            },
            &SearchControl::new(),
        )
        .unwrap();
        let tsv = results_tsv(&result.sites, guide, 3, PamSide::ThreePrime, "guide_1");
        let mut lines = tsv.lines();
        assert!(lines.next().unwrap().starts_with("##Generated by Arachnia"));
        assert_eq!(
            lines.next().unwrap(),
            "#Id\tBulge Type\tcrRNA\tDNA\tChromosome\tLocation\tDirection\tMismatches\tBulge Size"
        );
        let fields: Vec<&str> = lines.next().unwrap().split('\t').collect();
        assert_eq!(fields.len(), 9);
        assert_eq!(fields[0], "guide_1");
        assert_eq!(fields[1], "X");
        assert_eq!(fields[2], format!("{guide}NNN"));
        assert_eq!(fields[3], format!("{guide}TGG"));
        assert_eq!(fields[4], "chr1");
        assert_eq!(fields[5], "10");
        assert_eq!(fields[6], "+");
        assert_eq!(fields[7], "0");
        assert_eq!(fields[8], "0");
    }

    #[test]
    fn five_prime_pam_is_searched_on_both_strands() {
        let guide = "AATGGAGTCGTGAGCGCAAGAAC";
        let pam = "TTTG";
        let plus_window = format!("{pam}{guide}");
        let minus_window = reverse_complement(plus_window.as_bytes());
        let sequence = [
            b"AAAA".as_slice(),
            plus_window.as_bytes(),
            b"CCCCCCCC".as_slice(),
            minus_window.as_slice(),
            b"GGGG".as_slice(),
        ]
        .concat();
        let plus_site = 4;
        let minus_site = plus_site + plus_window.len() + 8;

        let genome = Genome::from_sequence("chrCpf1", &sequence).unwrap();
        let result = search_offtargets(
            guide,
            &genome,
            &OffTargetParams {
                max_mismatches: 0,
                pam_pattern: "TTTV".into(),
                pam_side: PamSide::FivePrime,
                both_strands: true,
            },
            &SearchControl::new(),
        )
        .unwrap();

        for (location, strand) in [(plus_site, Strand::Plus), (minus_site, Strand::Minus)] {
            let site = result
                .sites
                .iter()
                .find(|site| site.location == location && site.strand == strand)
                .unwrap();
            assert_eq!(site.mismatches, 0);
            assert_eq!(site.guide_dna, guide);
            assert_eq!(site.pam, pam);
        }

        let tsv = results_tsv(&result.sites, guide, 4, PamSide::FivePrime, "cpf1");
        let fields: Vec<&str> = tsv.lines().nth(2).unwrap().split('\t').collect();
        assert_eq!(fields[2], format!("NNNN{guide}"));
        assert_eq!(fields[3], format!("{pam}{guide}"));
        assert_eq!(fields[5], plus_site.to_string());
        assert_eq!(fields[6], "+");
    }

    #[test]
    fn cpf1_on_target_site_matches_design_coordinates() {
        use crate::design::{design, Filters};
        use crate::{Enzyme, Tracr};

        let sequence = "TGACTTTGAATGGAGTCGTGAGCGCAAGAACGCT";
        let (_, designed) = design(
            sequence,
            Enzyme::Cpf1,
            Tracr::Hsu2013,
            &Filters {
                both_strands: false,
                ..Default::default()
            },
        )
        .unwrap();
        let candidate = &designed.candidates[0];
        assert_eq!(candidate.pam, "TTTG");

        let genome = Genome::from_sequence("target", sequence.as_bytes()).unwrap();
        let result = search_offtargets(
            &candidate.guide,
            &genome,
            &OffTargetParams {
                max_mismatches: 0,
                pam_pattern: "TTTV".into(),
                pam_side: PamSide::FivePrime,
                both_strands: false,
            },
            &SearchControl::new(),
        )
        .unwrap();
        let site = result
            .sites
            .iter()
            .find(|site| site.strand == Strand::Plus)
            .unwrap();
        assert_eq!(site.mismatches, 0);
        // `candidate.start` is the 1-based guide start; the 5' PAM starts
        // `pam.len()` bases earlier on the plus strand.
        assert_eq!(site.location, candidate.start - 1 - candidate.pam.len());
        assert_eq!(site.pam, candidate.pam);
        assert_eq!(site.guide_dna, candidate.guide);
    }

    #[test]
    fn invalid_queries_are_rejected() {
        let genome = Genome::from_sequence("chr1", b"ACGTACGTACGTACGTACGTACGTACGT").unwrap();
        let control = SearchControl::new();
        let err =
            search_offtargets("ACGTN", &genome, &OffTargetParams::default(), &control).unwrap_err();
        assert!(matches!(err, SearchError::Invalid(_)));

        let err = search_offtargets(
            "ACGTACGTACGTACGTACGT",
            &genome,
            &OffTargetParams {
                pam_pattern: "NNN".into(),
                ..Default::default()
            },
            &control,
        )
        .unwrap_err();
        assert!(matches!(err, SearchError::Invalid(_)));

        let control = SearchControl::new();
        control.cancel();
        assert_eq!(
            search_offtargets(
                "ACGTACGTACGTACGTACGT",
                &genome,
                &OffTargetParams::default(),
                &control,
            )
            .unwrap_err(),
            SearchError::Cancelled
        );
    }
}
