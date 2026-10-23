//! Pure-Rust Rule Set 3 scoring for SpCas9 guides.
//!
//! Rule Set 3 is the gradient-boosted tree model described by Xiang et al.
//! (2021) and distributed by the GPP RNA Design group in
//! [`gpp-rnd/rs3`](https://github.com/gpp-rnd/rs3). This module reimplements
//! both the 632-feature featurizer and the LightGBM tree ensemble in Rust, so
//! the application can score guides offline without Python or a pickle reader.
//!
//! # Input layout
//!
//! A 30 nt context in guide orientation:
//!
//! ```text
//! [-4..-1] [guide 1..20] [PAM P1..P3] [+1..+3]
//!    4 nt       20 nt          3 nt          3 nt
//! ```
//!
//! The guide is `context[4..24]`; the PAM is `context[24..27]`. Only A/C/G/T
//! contexts are accepted. The PAM sequence is part of the model input even
//! though Rule Set 3 is normally applied to `NGG` sites.
//!
//! # Features
//!
//! [`features`] produces the same columns, in the same order, as
//! `sglearn.featurize_guides` plus the tracrRNA one-hot columns:
//!
//! * guide GC content and A/C/T/G mononucleotide fractions;
//! * non-overlapping dinucleotide fractions over the guide;
//! * position-specific one-hot nucleotides for all 30 positions;
//! * position-specific one-hot dinucleotides;
//! * five thermodynamic features (two Biopython Tm values, `seqfold` folding
//!   energy, and two GC-based values) computed by [`crate::tm`] and `seqfold`;
//! * two tracrRNA one-hot columns ([`Tracr::Hsu2013`], [`Tracr::Chen2013`]).
//!
//! # Output
//!
//! [`Rs3Model::score`] returns the raw ensemble score. Higher is predicted to
//! be more active. The GUI displays the raw score with three decimal places;
//! no calibration to a 0–100 scale is applied.

use std::sync::OnceLock;

use crate::tm;

/// Number of nucleotides in a Rule Set 3 context (4 + 20 + 3 + 3).
pub const CONTEXT_LEN: usize = 30;
/// SpCas9 guide length in nucleotides.
pub const GUIDE_LEN: usize = 20;
/// Number of features expected by the embedded LightGBM model.
pub const N_FEATURES: usize = 632;

/// The two tracrRNA variants supported by Rule Set 3.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Tracr {
    /// The tracrRNA used in Hsu et al. (2013).
    Hsu2013,
    /// The tracrRNA used in Chen et al. (2013).
    Chen2013,
}

impl Tracr {
    /// Human-readable tracrRNA name used in the UI and TSV metadata.
    pub fn label(self) -> &'static str {
        match self {
            Tracr::Hsu2013 => "Hsu2013",
            Tracr::Chen2013 => "Chen2013",
        }
    }
}

#[derive(Clone, Copy)]
struct Node {
    feature: u16,
    left: i32,
    right: i32,
    threshold: f64,
    value: f64,
}

struct Tree {
    nodes: Vec<Node>,
}

/// The embedded Rule Set 3 LightGBM ensemble.
///
/// Obtain the shared instance with [`Rs3Model::embedded`], or load a compatible
/// binary model with [`Rs3Model::from_bytes`]. The model file starts with the
/// magic `RS3LGBM1`, followed by little-endian version/count fields and
/// per-tree node arrays.
pub struct Rs3Model {
    trees: Vec<Tree>,
}

impl Rs3Model {
    /// Load the model compiled into the binary.
    ///
    /// The model is initialized once and shared for the life of the process.
    /// This is the normal entry point used by the GUI and library.
    pub fn embedded() -> &'static Self {
        static MODEL: OnceLock<Rs3Model> = OnceLock::new();
        MODEL.get_or_init(|| {
            Self::from_bytes(include_bytes!("../assets/rs3_model.bin"))
                .expect("embedded Rule Set 3 model is valid")
        })
    }

    /// Load a Rule Set 3 model from the compact binary format.
    ///
    /// # Errors
    ///
    /// Returns a description of the problem if the magic, version, feature
    /// count, node count, or tree data is invalid or truncated.
    pub fn from_bytes(data: &[u8]) -> Result<Self, String> {
        if data.len() < 24 || &data[..8] != b"RS3LGBM1" {
            return Err("not a Rule Set 3 model file".into());
        }
        let read_u32 = |offset: usize| -> u32 {
            u32::from_le_bytes(data[offset..offset + 4].try_into().unwrap())
        };
        let version = read_u32(8);
        if version != 1 {
            return Err(format!("unsupported model version {version}"));
        }
        let n_trees = read_u32(12) as usize;
        let n_features = read_u32(16) as usize;
        let total_nodes = read_u32(20) as usize;
        if n_features != N_FEATURES {
            return Err(format!("expected {N_FEATURES} features, got {n_features}"));
        }

        let mut trees = Vec::with_capacity(n_trees);
        let mut offset = 24;
        let mut seen_nodes = 0usize;
        for _ in 0..n_trees {
            if offset + 4 > data.len() {
                return Err("truncated tree header".into());
            }
            let count = read_u32(offset) as usize;
            offset += 4;
            if offset + count * 28 > data.len() {
                return Err("truncated tree data".into());
            }
            let mut nodes = Vec::with_capacity(count);
            for _ in 0..count {
                let feature = u16::from_le_bytes(data[offset..offset + 2].try_into().unwrap());
                let left = i32::from_le_bytes(data[offset + 4..offset + 8].try_into().unwrap());
                let right = i32::from_le_bytes(data[offset + 8..offset + 12].try_into().unwrap());
                let threshold =
                    f64::from_le_bytes(data[offset + 12..offset + 20].try_into().unwrap());
                let value = f64::from_le_bytes(data[offset + 20..offset + 28].try_into().unwrap());
                nodes.push(Node {
                    feature,
                    left,
                    right,
                    threshold,
                    value,
                });
                offset += 28;
            }
            seen_nodes += count;
            trees.push(Tree { nodes });
        }
        if seen_nodes != total_nodes {
            return Err("model node count mismatch".into());
        }
        Ok(Self { trees })
    }

    /// Raw Rule Set 3 sequence score. Higher is more active.
    ///
    /// Returns `None` if `context` is not exactly [`CONTEXT_LEN`] nt or contains
    /// a base outside A/C/G/T. The score is the sum of the leaf values of all
    /// trees in the ensemble; it is not clamped to 0–100.
    pub fn score(&self, context: &[u8], tracr: Tracr) -> Option<f64> {
        let features = features(context, tracr)?;
        Some(self.score_features(&features))
    }

    /// Evaluate a precomputed feature vector.
    ///
    /// `features` must have exactly [`N_FEATURES`] entries and be in the order
    /// produced by [`features`]. This is useful for batch scoring or for
    /// comparing the Rust featurizer with the Python reference.
    ///
    /// # Panics
    ///
    /// In debug builds, panics if `features.len() != N_FEATURES`. In release
    /// builds the length is only a debug assertion, so passing a short vector
    /// can panic while indexing a tree feature.
    pub fn score_features(&self, features: &[f64]) -> f64 {
        debug_assert_eq!(features.len(), N_FEATURES);
        let mut total = 0.0;
        for tree in &self.trees {
            let mut index = 0usize;
            loop {
                let node = &tree.nodes[index];
                if node.feature == u16::MAX {
                    total += node.value;
                    break;
                }
                let value = features[node.feature as usize];
                index = if value <= node.threshold {
                    node.left as usize
                } else {
                    node.right as usize
                };
            }
        }
        total
    }

    /// Number of decision trees in the ensemble.
    ///
    /// Useful for diagnostics and for checking that a loaded model matches the
    /// expected export.
    pub fn tree_count(&self) -> usize {
        self.trees.len()
    }
}

/// Build the 632-feature vector in exactly the order used by the training
/// DataFrame (LightGBM feature names are in `assets/rs3_features.json`).
///
/// `context` must be a 30 nt A/C/G/T window in guide orientation. Returns
/// `None` for the wrong length or any ambiguous base. The returned vector is
/// suitable for [`Rs3Model::score_features`].
///
/// Feature blocks, in order:
///
/// 1. guide GC fraction;
/// 2. guide A/C/T/G fractions;
/// 3. 16 non-overlapping guide dinucleotide fractions;
/// 4. 30 × 4 position-specific one-hot nucleotides;
/// 5. 29 × 4 position-specific one-hot dinucleotides;
/// 6. five thermodynamic features;
/// 7. two tracrRNA one-hot columns.
pub fn features(context: &[u8], tracr: Tracr) -> Option<Vec<f64>> {
    if context.len() != CONTEXT_LEN
        || !context
            .iter()
            .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
    {
        return None;
    }
    let guide = &context[4..24];
    let mut f = Vec::with_capacity(N_FEATURES);

    // GC content and mononucleotide fractions over the guide.
    let gc = guide.iter().filter(|&&b| b == b'G' || b == b'C').count() as f64 / 20.0;
    f.push(gc);
    for base in *b"ACTG" {
        f.push(guide.iter().filter(|&&b| b == base).count() as f64 / 20.0);
    }

    // Dinucleotide fractions over the guide. Python's str.count() is
    // non-overlapping, so "AAA".count("AA") == 1.
    for a in *b"ACTG" {
        for b in *b"ACTG" {
            let count = count_non_overlapping(guide, &[a, b]);
            f.push(count as f64 / 19.0);
        }
    }

    // Position-dependent single nucleotides: -4..-1, 1..20, P1..P3, +1..+3.
    for &base in context {
        for nt in *b"ACTG" {
            f.push(if base == nt { 1.0 } else { 0.0 });
        }
    }

    // Position-dependent dinucleotides. sglearn uses the label of the first
    // base of each adjacent pair, hence the last pair is labelled +2.
    for i in 0..CONTEXT_LEN - 1 {
        let pair = [context[i], context[i + 1]];
        for a in *b"ACTG" {
            for b in *b"ACTG" {
                f.push(if pair == [a, b] { 1.0 } else { 0.0 });
            }
        }
    }

    // PAM interaction: base at P1 crossed with base at +1.
    for a in *b"ACTG" {
        for b in *b"ACTG" {
            f.push(if context[24] == a && context[27] == b {
                1.0
            } else {
                0.0
            });
        }
    }

    // Thermodynamics, identical to sglearn.get_thermo.
    let tm_dd = tm::tm_dna_nn3(guide);
    let tm_rd = tm::tm_rna_dna_nn1(guide);
    let dg = dg_rr(guide);
    let tm_rd_1 = tm::tm_rna_dna_nn1(&guide[..10]);
    let tm_rd_2 = tm::tm_rna_dna_nn1(&guide[10..]);
    for value in [tm_dd, tm_rd, dg, tm_rd_1, tm_rd_2] {
        if !value.is_finite() {
            return None;
        }
        f.push(value);
    }

    // Longest homopolymer run in the guide.
    for nt in *b"ACTG" {
        f.push(max_run(guide, nt) as f64);
    }

    // tracrRNA one-hot.
    f.push(if tracr == Tracr::Hsu2013 { 1.0 } else { 0.0 });
    f.push(if tracr == Tracr::Chen2013 { 1.0 } else { 0.0 });

    debug_assert_eq!(f.len(), N_FEATURES);
    Some(f)
}

/// `seqfold.dg(guide.replace('T','U'))` with the special values sglearn maps.
fn dg_rr(guide: &[u8]) -> f64 {
    let rna: String = guide
        .iter()
        .map(|&b| if b == b'T' { 'U' } else { b as char })
        .collect();
    match seqfold::core::fold::dg(&rna, 37.0) {
        Ok(v) if v == f64::NEG_INFINITY => -20.0,
        Ok(v) if v == f64::INFINITY => 8.0,
        Ok(1600.0) => 7.0,
        Ok(v) => v,
        Err(_) => f64::NAN,
    }
}

fn count_non_overlapping(haystack: &[u8], needle: &[u8; 2]) -> usize {
    let mut count = 0;
    let mut i = 0;
    while i + 2 <= haystack.len() {
        if haystack[i] == needle[0] && haystack[i + 1] == needle[1] {
            count += 1;
            i += 2;
        } else {
            i += 1;
        }
    }
    count
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn model_loads() {
        let model = Rs3Model::embedded();
        assert_eq!(model.tree_count(), 2279);
    }

    #[test]
    fn readme_examples_match() {
        let model = Rs3Model::embedded();
        let cases = [
            (
                b"GACGAAAGCGACAACGCGTTCATCCGGGCA".as_slice(),
                Tracr::Hsu2013,
                -0.9003094354582534,
            ),
            (
                b"AGAAAACACTAGCATCCCCACCCGCGGACT".as_slice(),
                Tracr::Hsu2013,
                1.1145162188484359,
            ),
        ];
        for (context, tracr, expected) in cases {
            let got = model.score(context, tracr).unwrap();
            assert!((got - expected).abs() < 1e-12, "{got} != {expected}");
        }
    }

    #[test]
    fn featurization_matches_sglearn_regression_set() {
        let model = Rs3Model::embedded();
        let cases: serde_json::Value =
            serde_json::from_str(include_str!("../assets/rs3_tests.json")).unwrap();
        let mut checked = 0;
        for case in cases.as_array().unwrap() {
            let context = case["context"].as_str().unwrap().as_bytes();
            let tracr = match case["tracr"].as_str().unwrap() {
                "Hsu2013" => Tracr::Hsu2013,
                "Chen2013" => Tracr::Chen2013,
                other => panic!("unexpected tracr {other}"),
            };
            let expected = case["score"].as_f64().unwrap();
            let got = model.score(context, tracr).unwrap();
            assert!(
                (got - expected).abs() < 1e-9,
                "{} {}: {got} != {expected}",
                String::from_utf8_lossy(context),
                tracr.label()
            );
            checked += 1;
        }
        assert_eq!(checked, 84);
    }

    #[test]
    fn invalid_context_rejected() {
        let model = Rs3Model::embedded();
        assert!(model
            .score(b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN", Tracr::Hsu2013)
            .is_none());
    }

    #[test]
    fn homopolymer_dg_mapping() {
        assert_eq!(dg_rr(b"AAAAAAAAAAAAAAAAAAAA"), -20.0);
        assert_eq!(dg_rr(b"TTTTTTTTTTTTTTTTTTTT"), -20.0);
        assert_eq!(dg_rr(b"AAAGCGACAACGCGTTCATCC"), -4.4);
    }
}
