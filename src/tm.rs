//! Biopython-compatible nearest-neighbour melting temperatures.
//!
//! `sglearn` (used to featurize Rule Set 3 inputs) computes five
//! thermodynamic features with Biopython's `MeltingTemp.Tm_NN` defaults and
//! with `seqfold.dg`. This module reproduces the two Tm calculations exactly:
//!
//! * `Tm_NN(seq)` with `DNA_NN3` (Allawi & SantaLucia, 1997)
//! * `Tm_NN(seq, nn_table=R_DNA_NN1)` (Sugimoto et al., 1995)
//!
//! Only perfect-match, blunt-ended duplexes with the default concentrations
//! are needed. Biopython's terminal-mismatch and internal-mismatch tables are
//! never consulted for those inputs (verified against all 16 perfect
//! dinucleotides), and the default `saltcorr=5` entropy correction is applied.

const R_GAS: f64 = 1.987;
const DNAC1: f64 = 25.0;
const DNAC2: f64 = 25.0;
const NA_MM: f64 = 50.0;

#[derive(Clone, Copy)]
struct NnEntry {
    seq: [u8; 2],
    dh: f64,
    ds: f64,
}

const fn e(a: u8, b: u8, dh: f64, ds: f64) -> NnEntry {
    NnEntry {
        seq: [a, b],
        dh,
        ds,
    }
}

/// Unique nearest-neighbour stacks of `DNA_NN3` (Biopython stores 10).
const DNA_NN3: &[NnEntry] = &[
    e(b'A', b'A', -7.9, -22.2),  // AA/TT
    e(b'A', b'T', -7.2, -20.4),  // AT/TA
    e(b'T', b'A', -7.2, -21.3),  // TA/AT
    e(b'C', b'A', -8.5, -22.7),  // CA/GT
    e(b'G', b'T', -8.4, -22.4),  // GT/CA
    e(b'C', b'T', -7.8, -21.0),  // CT/GA
    e(b'G', b'A', -8.2, -22.2),  // GA/CT
    e(b'C', b'G', -10.6, -27.2), // CG/GC
    e(b'G', b'C', -9.8, -24.4),  // GC/CG
    e(b'G', b'G', -8.0, -19.9),  // GG/CC
];

/// `R_DNA_NN1` (Sugimoto et al. 1995). This table is asymmetric, so direct
/// matches must be preferred over reverse-complement matches.
const R_DNA_NN1: &[NnEntry] = &[
    e(b'T', b'T', -11.5, -36.4),
    e(b'G', b'T', -7.8, -21.6),
    e(b'C', b'T', -7.0, -19.7),
    e(b'A', b'T', -8.3, -23.9),
    e(b'T', b'G', -10.4, -28.4),
    e(b'G', b'G', -12.8, -31.9),
    e(b'C', b'G', -16.3, -47.1),
    e(b'A', b'G', -9.1, -23.5),
    e(b'T', b'C', -8.6, -22.9),
    e(b'G', b'C', -8.0, -17.1),
    e(b'C', b'C', -9.3, -23.2),
    e(b'A', b'C', -5.9, -12.3),
    e(b'T', b'A', -7.8, -23.2),
    e(b'G', b'A', -5.5, -13.5),
    e(b'C', b'A', -9.0, -26.1),
    e(b'A', b'A', -7.8, -21.9),
];

#[inline]
fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        other => other,
    }
}

/// Look up a perfect-match stack. Biopython first tries `seq/comp(seq)` and
/// then the reversed key.
fn nn_lookup(entries: &[NnEntry], a: u8, b: u8) -> Option<(f64, f64)> {
    for entry in entries {
        let [x, y] = entry.seq;
        if a == x && b == y {
            return Some((entry.dh, entry.ds));
        }
    }
    for entry in entries {
        let [x, y] = entry.seq;
        // reverse("xy/comp(x)comp(y)") == "comp(y)comp(x)/yx"
        if a == complement(y) && b == complement(x) {
            return Some((entry.dh, entry.ds));
        }
    }
    None
}

/// Reproduce `Bio.SeqUtils.MeltingTemp.Tm_NN` for a perfect-match duplex.
///
/// `init` is `(dH, dS)` for the general initiation term, `init_at` and
/// `init_gc` are the per-terminal-base terms, and `init_5t` is the 5' T /
/// 3' A penalty. All parameters use Biopython's defaults (`dnac1=dnac2=25`,
/// `Na=50`, `saltcorr=5`).
fn tm_nn(
    seq: &[u8],
    init: (f64, f64),
    init_at: (f64, f64),
    init_gc: (f64, f64),
    init_5t: (f64, f64),
    entries: &[NnEntry],
) -> f64 {
    if seq.len() < 2 {
        return f64::NAN;
    }

    let (mut dh, mut ds) = init;
    let ends = [seq[0], seq[seq.len() - 1]];
    let at = ends.iter().filter(|&&b| b == b'A' || b == b'T').count() as f64;
    let gc = ends.iter().filter(|&&b| b == b'G' || b == b'C').count() as f64;
    dh += init_at.0 * at + init_gc.0 * gc;
    ds += init_at.1 * at + init_gc.1 * gc;
    if seq[0] == b'T' || seq[seq.len() - 1] == b'A' {
        dh += init_5t.0;
        ds += init_5t.1;
    }

    for pair in seq.windows(2) {
        if let Some((h, s)) = nn_lookup(entries, pair[0], pair[1]) {
            dh += h;
            ds += s;
        } else {
            return f64::NAN;
        }
    }

    // Biopython saltcorr=5: entropy correction, monovalent Na+ only.
    let mon = NA_MM * 1e-3;
    let corr = 0.368 * (seq.len() as f64 - 1.0) * mon.ln();
    ds += corr;

    // k = (dnac1 - dnac2 / 2) * 1e-9
    let k = (DNAC1 - DNAC2 / 2.0) * 1e-9;
    (1000.0 * dh) / (ds + R_GAS * k.ln()) - 273.15
}

/// `MeltingTemp.Tm_NN(seq)` for DNA/DNA with the default `DNA_NN3` table.
///
/// Reproduces Biopython's `Bio.SeqUtils.MeltingTemp.Tm_NN` defaults:
/// 25 µM strand concentrations, 50 mM Na⁺, and salt correction mode 5.
/// `seq` must contain uppercase A/C/G/T; behaviour for ambiguous bases is not
/// defined by the model and callers should validate first.
///
/// Returns the melting temperature in °C.
pub fn tm_dna_nn3(seq: &[u8]) -> f64 {
    tm_nn(
        seq,
        (0.0, 0.0),
        (2.3, 4.1),
        (0.1, -2.8),
        (0.0, 0.0),
        DNA_NN3,
    )
}

/// `MeltingTemp.Tm_NN(seq, nn_table=R_DNA_NN1)` for an RNA/DNA hybrid.
///
/// Reproduces Biopython's `Tm_NN` defaults with the Sugimoto et al. (1995)
/// `R_DNA_NN1` table, including the RNA initiation terms. `seq` must contain
/// uppercase A/C/G/T.
///
/// Returns the melting temperature in °C.
pub fn tm_rna_dna_nn1(seq: &[u8]) -> f64 {
    tm_nn(
        seq,
        (1.9, -3.9),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        R_DNA_NN1,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) {
        assert!((a - b).abs() < 1e-9, "{a} != {b}");
    }

    #[test]
    fn matches_biopython_dna_nn3() {
        close(tm_dna_nn3(b"AAAGCGACAACGCGTTCATC"), 54.53053045281678);
        close(tm_dna_nn3(b"AAAAAAAAAAAAAAAAAAAA"), 35.96102091833649);
        close(tm_dna_nn3(b"TTTTTTTTTTTTTTTTTTTT"), 35.96102091833649);
        close(tm_dna_nn3(b"ACGTACGTACGTACGTACGT"), 53.09666582369118);
        close(tm_dna_nn3(b"GGCCTTAAGGCCTTAAGGCC"), 54.868842816433926);
    }

    #[test]
    fn matches_biopython_rna_dna_nn1() {
        close(tm_rna_dna_nn1(b"AAAGCGACAACGCGTTCATC"), 49.81945535488734);
        close(tm_rna_dna_nn1(b"AAAAAAAAAAAAAAAAAAAA"), 33.49131014685321);
        close(tm_rna_dna_nn1(b"TTTTTTTTTTTTTTTTTTTT"), 14.65050091804062);
        close(tm_rna_dna_nn1(b"ACGTACGTACGTACGTACGT"), 47.71346224764068);
        close(tm_rna_dna_nn1(b"GGCCTTAAGGCCTTAAGGCC"), 58.53586755536401);
    }

    #[test]
    fn short_fragments_match_biopython() {
        // First/last ten bases of the guide used by sglearn.
        close(tm_rna_dna_nn1(b"AAAGCGACAA"), 21.819799287642695);
        close(tm_rna_dna_nn1(b"CGCGTTCATC"), 20.425888283149277);
    }
}
