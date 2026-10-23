//! Candidate discovery, filtering and on-target TSV export.
//!
//! This module turns a raw DNA/FASTA sequence into scored [`GuideCandidate`]s
//! for either SpCas9 or Cpf1/Cas12a.
//!
//! # Model contexts
//!
//! Both enzymes are scored from a fixed-length context in guide orientation:
//!
//! * **SpCas9** — 30 nt: `[4 nt upstream][20 nt guide][3 nt PAM][3 nt downstream]`.
//!   The PAM is 3′ of the guide (`NGG`).
//! * **Cpf1/Cas12a** — 34 nt:
//!   `[4 nt upstream][4 nt TTTV PAM][23 nt guide][3 nt downstream]`.
//!   The PAM is 5′ of the guide.
//!
//! Candidates whose context contains a base outside A/C/G/T are skipped and
//! counted in [`DesignResult::skipped_ambiguous`]. PAM sites are counted before
//! scoring, so [`DesignResult::pam_sites`] can be larger than the number of
//! returned candidates.
//!
//! # Strands and coordinates
//!
//! With [`Filters::both_strands`] enabled, the input is scanned on both
//! strands. A minus-strand candidate stores its guide and PAM in guide
//! orientation, but its [`GuideCandidate::start`], [`GuideCandidate::end`] and
//! [`GuideCandidate::cut`] are always **1-based coordinates on the input
//! (plus) strand**. `start` is the lowest coordinate and `end` the highest.
//!
//! # Filtering and ordering
//!
//! [`design_bases`] applies the inclusive GC bounds, minimum score and optional
//! poly-T limit after scoring, then sorts survivors by descending score, then
//! by ascending plus-strand start and strand symbol. The GUI can re-sort the
//! result for display without rescoring.
//!
//! # Example
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
//! assert_eq!(prepared.bases.len(), sequence.len());
//! assert!(!result.candidates.is_empty());
//! ```

use crate::cpf1::{Cpf1Model, CONTEXT_LEN as CPF1_CONTEXT_LEN, GUIDE_LEN as CPF1_GUIDE_LEN};
use crate::ruleset3::{
    Rs3Model, Tracr, CONTEXT_LEN as CAS9_CONTEXT_LEN, GUIDE_LEN as CAS9_GUIDE_LEN,
};

/// The supported CRISPR nuclease families.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Enzyme {
    /// *Streptococcus pyogenes* Cas9: 20 nt guide, 3′ `NGG` PAM, Rule Set 3.
    Cas9,
    /// Cpf1/Cas12a: 23 nt guide, 5′ `TTTV` PAM, SeqDeepCpf1.
    Cpf1,
}

impl Enzyme {
    /// Human-readable enzyme name used in the UI and reports.
    pub fn label(self) -> &'static str {
        match self {
            Enzyme::Cas9 => "SpCas9",
            Enzyme::Cpf1 => "Cpf1 (Cas12a)",
        }
    }

    /// Canonical PAM motif in guide orientation.
    pub fn pam_label(self) -> &'static str {
        match self {
            Enzyme::Cas9 => "NGG",
            Enzyme::Cpf1 => "TTTV",
        }
    }

    /// Guide length in nucleotides (20 for Cas9, 23 for Cpf1).
    pub fn guide_len(self) -> usize {
        match self {
            Enzyme::Cas9 => CAS9_GUIDE_LEN,
            Enzyme::Cpf1 => CPF1_GUIDE_LEN,
        }
    }

    /// Name of the on-target scoring model shown in the UI.
    pub fn score_label(self) -> &'static str {
        match self {
            Enzyme::Cas9 => "Rule Set 3 score",
            Enzyme::Cpf1 => "SeqDeepCpf1 score",
        }
    }
}

/// Which strand of the input sequence a candidate was found on.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Strand {
    /// The input (plus) strand.
    Plus,
    /// The reverse-complement (minus) strand.
    Minus,
}

impl Strand {
    /// `'+'` or `'-'`, used in tables and TSV output.
    pub fn symbol(self) -> char {
        match self {
            Strand::Plus => '+',
            Strand::Minus => '-',
        }
    }
}

/// Pre-scoring filters applied by [`design_bases`].
///
/// The GUI exposes these as live controls. They are applied after scoring, so
/// changing them does not rescore candidates.
#[derive(Clone, Debug)]
pub struct Filters {
    /// Inclusive minimum GC content in percent (0–100).
    pub gc_min: f64,
    /// Inclusive maximum GC content in percent (0–100).
    pub gc_max: f64,
    /// Minimum on-target score. Use `f64::NEG_INFINITY` to disable.
    pub min_score: f64,
    /// Reject guides containing a T run at least this long (`None` disables).
    pub max_poly_t: Option<usize>,
    /// Scan the reverse-complement strand in addition to the input strand.
    pub both_strands: bool,
}

impl Default for Filters {
    fn default() -> Self {
        Self {
            gc_min: 0.0,
            gc_max: 100.0,
            min_score: f64::NEG_INFINITY,
            max_poly_t: Some(4),
            both_strands: true,
        }
    }
}

/// A scored guide candidate.
///
/// `guide` and `pam` are always in guide orientation (5′→3′), while `start`,
/// `end` and `cut` are 1-based coordinates on the input (plus) strand.
#[derive(Clone, Debug)]
pub struct GuideCandidate {
    /// Guide in guide orientation, 5′→3′.
    pub guide: String,
    /// PAM in guide orientation (so it always matches the enzyme motif).
    pub pam: String,
    /// Strand on which the candidate was found.
    pub strand: Strand,
    /// Lowest 1-based coordinate on the input (plus) strand.
    pub start: usize,
    /// Highest 1-based coordinate on the input (plus) strand.
    pub end: usize,
    /// Approximate 1-based cut position on the plus strand.
    pub cut: usize,
    /// GC content of the guide in percent.
    pub gc: f64,
    /// On-target score (Rule Set 3 raw score or SeqDeepCpf1 0–100 activity).
    pub score: f64,
    /// Full model context in guide orientation.
    pub context: String,
    /// Length of the longest T run in the guide.
    pub poly_t_run: usize,
    /// Enzyme that produced this candidate.
    pub enzyme: Enzyme,
}

/// Result of scanning a sequence for guides.
#[derive(Clone, Debug)]
pub struct DesignResult {
    /// Scored candidates that passed all filters, sorted by descending score.
    pub candidates: Vec<GuideCandidate>,
    /// Number of PAM sites found before score/context/filter checks.
    pub pam_sites: usize,
    /// Sites skipped because the model context contained a non-ACGT base.
    pub skipped_ambiguous: usize,
    /// Candidates removed by the user's filters.
    pub filtered_out: usize,
}

/// A normalized sequence plus optional FASTA metadata.
#[derive(Clone, Debug)]
pub struct PreparedSequence {
    /// FASTA record name, if the input contained a header.
    pub name: Option<String>,
    /// Normalized bases: uppercase A/C/G/T/N, with U converted to T.
    pub bases: Vec<u8>,
    /// Non-fatal input warning, such as multiple FASTA records being ignored.
    pub warning: Option<String>,
}

/// Error returned when a sequence cannot be parsed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DnaError(
    /// Human-readable error message.
    pub String,
);

impl std::fmt::Display for DnaError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for DnaError {}

/// Parse a raw sequence or the first record of a FASTA file.
///
/// Whitespace is ignored, lower-case is accepted, U is treated as T, and N is
/// kept as an ambiguous base (candidates touching it are skipped later).
pub fn prepare_sequence(input: &str) -> Result<PreparedSequence, DnaError> {
    let mut name = None;
    let mut bases = Vec::new();
    let mut warning = None;
    let mut saw_bases = false;

    for (line_no, raw_line) in input.lines().enumerate() {
        let line = raw_line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            if saw_bases {
                warning = Some("multiple FASTA records found; only the first was used".into());
                break;
            }
            name = Some(header.trim().to_string());
            continue;
        }
        for ch in line.chars() {
            if ch.is_whitespace() {
                continue;
            }
            let base = match ch.to_ascii_uppercase() {
                'A' => b'A',
                'C' => b'C',
                'G' => b'G',
                'T' | 'U' => b'T',
                'N' => b'N',
                other => {
                    return Err(DnaError(format!(
                        "invalid base '{other}' on line {} (only A, C, G, T, U, N are accepted)",
                        line_no + 1
                    )));
                }
            };
            bases.push(base);
            saw_bases = true;
        }
    }

    if bases.is_empty() {
        return Err(DnaError("no sequence found".into()));
    }
    Ok(PreparedSequence {
        name,
        bases,
        warning,
    })
}

/// Parse a raw sequence and design guides in one call.
///
/// This is a convenience wrapper around [`prepare_sequence`] followed by
/// [`design_bases`]. It accepts the same input as [`prepare_sequence`] and
/// returns both the normalized sequence and the scored/filtered candidates.
///
/// # Errors
///
/// Returns [`DnaError`] if the input contains no sequence or an unsupported
/// base character.
///
/// # Example
///
/// ```no_run
/// use arachnia::{design, Enzyme, Filters, Tracr};
///
/// let (prepared, result) = design(
///     ">target\nAAAAACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTGGCCC\n",
///     Enzyme::Cas9,
///     Tracr::Hsu2013,
///     &Filters::default(),
/// )
/// .unwrap();
/// println!("{} bp, {} guides", prepared.bases.len(), result.candidates.len());
/// ```
pub fn design(
    input: &str,
    enzyme: Enzyme,
    tracr: Tracr,
    filters: &Filters,
) -> Result<(PreparedSequence, DesignResult), DnaError> {
    let prepared = prepare_sequence(input)?;
    let result = design_bases(&prepared.bases, enzyme, tracr, filters);
    Ok((prepared, result))
}

/// Design guides for an already-normalized base vector.
///
/// `bases` must contain uppercase A/C/G/T/N. Unlike [`prepare_sequence`], this
/// function performs no parsing or case conversion. It scans the requested
/// strand(s), scores every complete context with the selected model, applies
/// `filters`, and returns candidates sorted by descending score.
///
/// # Context and strand behavior
///
/// * SpCas9 contexts are 30 nt (`[4][20][3][3]`) with a 3′ `NGG` PAM.
/// * Cpf1 contexts are 34 nt (`[4][4][23][3]`) with a 5′ `TTTV` PAM.
/// * When [`Filters::both_strands`] is true, the reverse complement is scanned
///   and minus-strand candidates are reported in guide orientation with
///   plus-strand coordinates.
///
/// # Filtering
///
/// GC bounds are inclusive. A candidate passes the poly-T filter when its
/// longest T run is strictly shorter than the configured limit. Sites with
/// ambiguous context bases are skipped and counted in
/// [`DesignResult::skipped_ambiguous`].
pub fn design_bases(bases: &[u8], enzyme: Enzyme, tracr: Tracr, filters: &Filters) -> DesignResult {
    let reverse = if filters.both_strands {
        Some(reverse_complement(bases))
    } else {
        None
    };

    let mut candidates = Vec::new();
    let mut pam_sites = 0;
    let mut skipped_ambiguous = 0;

    scan_orientation(
        bases,
        false,
        bases.len(),
        enzyme,
        tracr,
        &mut candidates,
        &mut pam_sites,
        &mut skipped_ambiguous,
    );
    if let Some(ref rc) = reverse {
        scan_orientation(
            rc,
            true,
            bases.len(),
            enzyme,
            tracr,
            &mut candidates,
            &mut pam_sites,
            &mut skipped_ambiguous,
        );
    }

    let raw = candidates.len();
    candidates.retain(|c| {
        c.gc >= filters.gc_min
            && c.gc <= filters.gc_max
            && c.score >= filters.min_score
            && filters.max_poly_t.is_none_or(|max| c.poly_t_run < max)
    });
    let filtered_out = raw - candidates.len();
    candidates.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.start.cmp(&b.start))
            .then_with(|| a.strand.symbol().cmp(&b.strand.symbol()))
    });

    DesignResult {
        candidates,
        pam_sites,
        skipped_ambiguous,
        filtered_out,
    }
}

/// Render the selected candidates as the on-target TSV table used by the GUI
/// export and clipboard actions.
///
/// `view` is a list of indices into `candidates`, in display order. Invalid
/// indices are skipped. The header columns are:
///
/// ```text
/// rank  guide  pam  strand  start  end  cut  gc_percent  score  context  poly_t_run
/// ```
///
/// Coordinates are 1-based plus-strand positions; the guide and PAM are in
/// guide orientation. `rank` follows the order of `view`, not the score.
pub fn results_tsv(candidates: &[GuideCandidate], view: &[usize]) -> String {
    let mut out = String::from(
        "rank\tguide\tpam\tstrand\tstart\tend\tcut\tgc_percent\tscore\tcontext\tpoly_t_run\n",
    );
    for (rank, &idx) in view.iter().enumerate() {
        let Some(c) = candidates.get(idx) else {
            continue;
        };
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.6}\t{}\t{}\n",
            rank + 1,
            c.guide,
            c.pam,
            c.strand.symbol(),
            c.start,
            c.end,
            c.cut,
            c.gc,
            c.score,
            c.context,
            c.poly_t_run
        ));
    }
    out
}

#[allow(clippy::too_many_arguments)]
fn scan_orientation(
    seq: &[u8],
    minus: bool,
    plus_len: usize,
    enzyme: Enzyme,
    tracr: Tracr,
    out: &mut Vec<GuideCandidate>,
    pam_sites: &mut usize,
    skipped: &mut usize,
) {
    let n = seq.len();
    match enzyme {
        Enzyme::Cas9 => {
            if n < CAS9_CONTEXT_LEN {
                return;
            }
            for p in 24..=n - 6 {
                if !matches!(seq[p], b'A' | b'C' | b'G' | b'T')
                    || seq[p + 1] != b'G'
                    || seq[p + 2] != b'G'
                {
                    continue;
                }
                *pam_sites += 1;
                let context = &seq[p - 24..p + 6];
                let Some(score) = Rs3Model::embedded().score(context, tracr) else {
                    *skipped += 1;
                    continue;
                };
                let guide = &context[4..24];
                let pam = &context[24..27];
                let (start, end, cut) = if minus {
                    let start0 = plus_len - p;
                    (start0 + 1, start0 + 20, start0 + 4)
                } else {
                    let start0 = p - 20;
                    (start0 + 1, start0 + 20, start0 + 17)
                };
                out.push(make_candidate(
                    guide, pam, context, minus, start, end, cut, score, enzyme,
                ));
            }
        }
        Enzyme::Cpf1 => {
            if n < CPF1_CONTEXT_LEN {
                return;
            }
            for p in 4..=n - 30 {
                if seq[p] != b'T'
                    || seq[p + 1] != b'T'
                    || seq[p + 2] != b'T'
                    || !matches!(seq[p + 3], b'A' | b'C' | b'G')
                {
                    continue;
                }
                *pam_sites += 1;
                let context = &seq[p - 4..p + 30];
                let Some(score) = Cpf1Model::embedded().score(context) else {
                    *skipped += 1;
                    continue;
                };
                let guide = &context[8..31];
                let pam = &context[4..8];
                // A Keras model trained on 23 nt guides: the cut is staggered
                // after position 18 from the PAM-proximal 5' end.
                let (start, end, cut) = if minus {
                    let start0 = plus_len - p - 27;
                    (start0 + 1, start0 + 23, start0 + 6)
                } else {
                    let start0 = p + 4;
                    (start0 + 1, start0 + 23, start0 + 18)
                };
                out.push(make_candidate(
                    guide, pam, context, minus, start, end, cut, score, enzyme,
                ));
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn make_candidate(
    guide: &[u8],
    pam: &[u8],
    context: &[u8],
    minus: bool,
    start: usize,
    end: usize,
    cut: usize,
    score: f64,
    enzyme: Enzyme,
) -> GuideCandidate {
    let gc = guide.iter().filter(|&&b| b == b'G' || b == b'C').count() as f64 * 100.0
        / guide.len() as f64;
    GuideCandidate {
        guide: String::from_utf8_lossy(guide).into_owned(),
        pam: String::from_utf8_lossy(pam).into_owned(),
        strand: if minus { Strand::Minus } else { Strand::Plus },
        start,
        end,
        cut,
        gc,
        score,
        context: String::from_utf8_lossy(context).into_owned(),
        poly_t_run: max_run(guide, b'T'),
        enzyme,
    }
}

fn max_run(seq: &[u8], base: u8) -> usize {
    let mut best = 0;
    let mut current = 0;
    for &b in seq {
        if b == base {
            current += 1;
            best = best.max(current);
        } else {
            current = 0;
        }
    }
    best
}

/// Return the reverse complement of a nucleotide sequence.
///
/// A/C/G/T are complemented; every other byte (including `N` and IUPAC codes)
/// becomes `N`. The returned vector has the same length as the input and is
/// suitable for scanning the minus strand.
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'N',
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cpf1_example_fields() {
        let seq = "TGACTTTGAATGGAGTCGTGAGCGCAAGAACGCT";
        let (_, result) = design(
            seq,
            Enzyme::Cpf1,
            Tracr::Hsu2013,
            &Filters {
                both_strands: false,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(result.candidates.len(), 1);
        let c = &result.candidates[0];
        assert_eq!(c.guide, "AATGGAGTCGTGAGCGCAAGAAC");
        assert_eq!(c.pam, "TTTG");
        assert_eq!((c.start, c.end, c.cut), (9, 31, 26));
        assert_eq!(c.strand, Strand::Plus);
        assert!((c.score - 55.69930648803711).abs() < 1e-3);
    }

    #[test]
    fn cas9_example_fields() {
        // 4 nt upstream + guide + TGG + 3 nt downstream
        let seq = "AAAAACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTGGCCC";
        let (_, result) = design(
            seq,
            Enzyme::Cas9,
            Tracr::Hsu2013,
            &Filters {
                both_strands: false,
                max_poly_t: None,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(result.candidates.len(), 1);
        let c = &result.candidates[0];
        assert_eq!(c.guide.len(), 20);
        assert_eq!(c.pam, "TGG");
        assert_eq!(c.strand, Strand::Plus);
        assert!(c.start > 1);
        assert_eq!(c.end, c.start + 19);
        assert_eq!(c.cut, c.start + 16);
    }

    #[test]
    fn reverse_strand_coordinates() {
        // A Cas9 site on the minus strand of `plus`: the guide and NGG PAM
        // are shown in guide orientation, and on the plus strand the PAM is
        // CCT.
        let guide = "ACGTACGTACGTACGTACGT";
        let target_minus = format!("GATC{guide}AGGTTA");
        let plus = String::from_utf8(reverse_complement(target_minus.as_bytes())).unwrap();
        assert_eq!(plus, format!("TAACCT{guide}GATC"));

        let filters = Filters {
            both_strands: false,
            max_poly_t: None,
            ..Default::default()
        };
        // Search the plus strand for a site on the minus strand.
        let (_, result) = design(&plus, Enzyme::Cas9, Tracr::Hsu2013, &filters).unwrap();
        // There is no NGG on the plus strand, so no plus-strand candidate.
        assert!(result.candidates.is_empty());

        let (_, both) = design(
            &plus,
            Enzyme::Cas9,
            Tracr::Hsu2013,
            &Filters {
                both_strands: true,
                max_poly_t: None,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(both.candidates.len(), 1);
        let c = &both.candidates[0];
        assert_eq!(c.strand, Strand::Minus);
        assert_eq!(c.guide, guide);
        assert_eq!(c.pam, "AGG");
        assert_eq!((c.start, c.end, c.cut), (7, 26, 10));

        // The same site looks like a plus-strand site in the reverse
        // complement, with mirrored coordinates.
        let rc = String::from_utf8(reverse_complement(plus.as_bytes())).unwrap();
        let (_, rc_result) = design(&rc, Enzyme::Cas9, Tracr::Hsu2013, &filters).unwrap();
        assert_eq!(rc_result.candidates.len(), 1);
        let m = &rc_result.candidates[0];
        assert_eq!((m.start, m.end), (5, 24));
        assert_eq!(c.start, plus.len() + 1 - m.end);
        assert_eq!(c.end, plus.len() + 1 - m.start);
    }

    #[test]
    fn n_skips_candidates() {
        let seq = "NNNNTTTGAATGGAGTCGTGAGCGCAAGAACGCT";
        let (_, result) = design(seq, Enzyme::Cpf1, Tracr::Hsu2013, &Filters::default()).unwrap();
        assert_eq!(result.candidates.len(), 0);
        assert_eq!(result.pam_sites, 1);
        assert_eq!(result.skipped_ambiguous, 1);
    }

    #[test]
    fn tsv_export_follows_view_order() {
        let seq = "AAAAACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTGGCCC";
        let (_, result) = design(
            seq,
            Enzyme::Cas9,
            Tracr::Hsu2013,
            &Filters {
                both_strands: false,
                max_poly_t: None,
                ..Default::default()
            },
        )
        .unwrap();
        let tsv = results_tsv(&result.candidates, &[0]);
        let mut lines = tsv.lines();
        assert_eq!(
            lines.next().unwrap(),
            "rank\tguide\tpam\tstrand\tstart\tend\tcut\tgc_percent\tscore\tcontext\tpoly_t_run"
        );
        let fields: Vec<&str> = lines.next().unwrap().split('\t').collect();
        assert_eq!(fields.len(), 11);
        assert_eq!(fields[0], "1");
        assert_eq!(fields[1], result.candidates[0].guide);
        assert_eq!(fields[2], "TGG");
        assert_eq!(fields[3], "+");
        assert!(lines.next().is_none());

        assert_eq!(results_tsv(&result.candidates, &[]).lines().count(), 1);
    }
}
