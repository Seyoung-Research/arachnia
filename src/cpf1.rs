//! Pure-Rust inference for SeqDeepCpf1.
//!
//! SeqDeepCpf1 is the Cpf1/Cas12a activity model from Kim et al.,
//! *Nat. Biotechnol.* 2018. The reference implementation is a Keras 1-D CNN;
//! this module evaluates the same network with weights converted to a compact
//! embedded binary by `scripts/export_models.py`.
//!
//! # Input layout
//!
//! A 34 nt context in guide orientation:
//!
//! ```text
//! [-4..-1] [TTTV PAM] [guide 1..23] [+1..+3]
//!    4 nt      4 nt        23 nt         3 nt
//! ```
//!
//! The PAM is 5′ of the guide. Only A/C/G/T contexts are accepted; the PAM
//! itself is part of the model input.
//!
//! # Output
//!
//! [`Cpf1Model::score`] returns the predicted activity on the paper's 0–100
//! scale. Higher is predicted to be more active. The value is a model output,
//! not a measured editing efficiency.
//!
//! # Network
//!
//! `Conv1d(4→80, k=5)` → ReLU → `AvgPool1d(2)` → flatten → dense 80 →
//! dense 40 → dense 40 → scalar. All arithmetic is `f32`, matching the
//! exported weights and the original inference path closely enough for ranking
//! guides.

use std::sync::OnceLock;

/// Number of nucleotides in a SeqDeepCpf1 context (4 + 4 + 23 + 3).
pub const CONTEXT_LEN: usize = 34;
/// Cpf1/Cas12a guide length in nucleotides.
pub const GUIDE_LEN: usize = 23;
/// Length of the 5′ Cpf1 PAM (`TTTV`).
pub const PAM_LEN: usize = 4;

const CONV_OUT: usize = 80;
const CONV_IN: usize = 4;
const CONV_K: usize = 5;
const POOL_OUT: usize = 15;
const FC1_IN: usize = 1200;
const FC1_OUT: usize = 80;
const FC2_OUT: usize = 40;
const FC3_OUT: usize = 40;

/// The embedded SeqDeepCpf1 convolutional network.
///
/// Use [`Cpf1Model::embedded`] for the shared model compiled into the binary,
/// or [`Cpf1Model::from_bytes`] to load a compatible weight file. The weight
/// format starts with the magic `CPF1W001`, followed by convolution shape
/// fields and little-endian `f32` tensors.
pub struct Cpf1Model {
    conv_w: Vec<f32>, // [80][4][5]
    conv_b: Vec<f32>, // [80]
    fc1_w: Vec<f32>,  // [80][1200]
    fc1_b: Vec<f32>,
    fc2_w: Vec<f32>, // [40][80]
    fc2_b: Vec<f32>,
    fc3_w: Vec<f32>, // [40][40]
    fc3_b: Vec<f32>,
    out_w: Vec<f32>, // [1][40]
    out_b: Vec<f32>,
}

impl Cpf1Model {
    /// Load the SeqDeepCpf1 weights compiled into the binary.
    ///
    /// The model is initialized once and shared for the life of the process.
    /// This is the normal entry point used by the GUI and library.
    pub fn embedded() -> &'static Self {
        static MODEL: OnceLock<Cpf1Model> = OnceLock::new();
        MODEL.get_or_init(|| {
            Self::from_bytes(include_bytes!("../assets/seqdeepcpf1_weights.bin"))
                .expect("embedded SeqDeepCpf1 weights are valid")
        })
    }

    /// Load SeqDeepCpf1 weights from the compact binary format.
    ///
    /// # Errors
    ///
    /// Returns a description of the problem if the magic, convolution shape,
    /// total file length, or tensor data is invalid or truncated.
    pub fn from_bytes(data: &[u8]) -> Result<Self, String> {
        if data.len() < 20 || &data[..8] != b"CPF1W001" {
            return Err("not a SeqDeepCpf1 weights file".into());
        }
        let u32_at = |offset: usize| {
            u32::from_le_bytes(data[offset..offset + 4].try_into().unwrap()) as usize
        };
        if u32_at(8) != CONV_OUT || u32_at(12) != CONV_IN || u32_at(16) != CONV_K {
            return Err("unexpected convolution shape".into());
        }
        let mut offset = 20;
        let lengths = [
            CONV_OUT * CONV_IN * CONV_K,
            CONV_OUT,
            FC1_OUT * FC1_IN,
            FC1_OUT,
            FC2_OUT * FC1_OUT,
            FC2_OUT,
            FC3_OUT * FC3_OUT,
            FC3_OUT,
            FC3_OUT,
            1,
        ];
        let expected: usize = 20 + lengths.iter().sum::<usize>() * 4;
        if data.len() != expected {
            return Err(format!(
                "weights file has {} bytes, expected {expected}",
                data.len()
            ));
        }
        let mut take = |n: usize| -> Vec<f32> {
            let mut out = Vec::with_capacity(n);
            for _ in 0..n {
                out.push(f32::from_le_bytes(
                    data[offset..offset + 4].try_into().unwrap(),
                ));
                offset += 4;
            }
            out
        };
        Ok(Self {
            conv_w: take(lengths[0]),
            conv_b: take(lengths[1]),
            fc1_w: take(lengths[2]),
            fc1_b: take(lengths[3]),
            fc2_w: take(lengths[4]),
            fc2_b: take(lengths[5]),
            fc3_w: take(lengths[6]),
            fc3_b: take(lengths[7]),
            out_w: take(lengths[8]),
            out_b: take(lengths[9]),
        })
    }

    /// Predict activity from a 34 nt context containing only A/C/G/T.
    ///
    /// The context must be in guide orientation:
    /// `[4 nt upstream][4 nt PAM][23 nt guide][3 nt downstream]`. Returns
    /// `None` for the wrong length or an ambiguous base. The output is on the
    /// paper's 0–100 activity scale (higher is better).
    pub fn score(&self, context: &[u8]) -> Option<f64> {
        if context.len() != CONTEXT_LEN
            || !context
                .iter()
                .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
        {
            return None;
        }

        // One-hot, channel order A, C, G, T (as in SeqDeepCpf1.py).
        let mut input = [0.0f32; CONV_IN * CONTEXT_LEN];
        for (pos, &base) in context.iter().enumerate() {
            let channel = match base {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => unreachable!(),
            };
            input[channel * CONTEXT_LEN + pos] = 1.0;
        }

        // Conv1d(4 -> 80, k=5), no padding.
        let conv_len = CONTEXT_LEN - CONV_K + 1; // 30
        let mut conv = vec![0.0f32; CONV_OUT * conv_len];
        for out_ch in 0..CONV_OUT {
            let w_base = out_ch * CONV_IN * CONV_K;
            let b = self.conv_b[out_ch];
            for pos in 0..conv_len {
                let mut sum = b;
                for in_ch in 0..CONV_IN {
                    let x_base = in_ch * CONTEXT_LEN + pos;
                    let w = &self.conv_w[w_base + in_ch * CONV_K..][..CONV_K];
                    for k in 0..CONV_K {
                        sum += w[k] * input[x_base + k];
                    }
                }
                // ReLU immediately follows the convolution.
                conv[out_ch * conv_len + pos] = sum.max(0.0);
            }
        }

        // AvgPool1d(2), then position-major flatten (15 * 80 = 1200).
        let mut flat = vec![0.0f32; FC1_IN];
        for out_ch in 0..CONV_OUT {
            for pos in 0..POOL_OUT {
                let a = conv[out_ch * conv_len + pos * 2];
                let b = conv[out_ch * conv_len + pos * 2 + 1];
                flat[pos * CONV_OUT + out_ch] = (a + b) * 0.5;
            }
        }

        // FC1 1200 -> 80, ReLU.
        let mut h1 = [0.0f32; FC1_OUT];
        for (out, h1_out) in h1.iter_mut().enumerate() {
            let row = &self.fc1_w[out * FC1_IN..][..FC1_IN];
            let mut sum = self.fc1_b[out];
            for (w, x) in row.iter().zip(flat.iter()) {
                sum += w * x;
            }
            *h1_out = sum.max(0.0);
        }

        // FC2 80 -> 40, ReLU.
        let mut h2 = [0.0f32; FC2_OUT];
        for (out, h2_out) in h2.iter_mut().enumerate() {
            let row = &self.fc2_w[out * FC1_OUT..][..FC1_OUT];
            let mut sum = self.fc2_b[out];
            for (w, x) in row.iter().zip(h1.iter()) {
                sum += w * x;
            }
            *h2_out = sum.max(0.0);
        }

        // FC3 40 -> 40, ReLU.
        let mut h3 = [0.0f32; FC3_OUT];
        for (out, h3_out) in h3.iter_mut().enumerate() {
            let row = &self.fc3_w[out * FC3_OUT..][..FC3_OUT];
            let mut sum = self.fc3_b[out];
            for (w, x) in row.iter().zip(h2.iter()) {
                sum += w * x;
            }
            *h3_out = sum.max(0.0);
        }

        // Linear 40 -> 1.
        let mut out = self.out_b[0];
        for (w, x) in self.out_w.iter().zip(h3.iter()) {
            out += w * x;
        }
        Some(out as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn example_scores_match_pytorch() {
        let model = Cpf1Model::embedded();
        let cases = [
            (
                b"TGACTTTGAATGGAGTCGTGAGCGCAAGAACGCT".as_slice(),
                55.69930648803711,
            ),
            (
                b"GTTATTTGAGCAATGCCACTTAATAAACATGTAA".as_slice(),
                53.4698371887207,
            ),
        ];
        for (context, expected) in cases {
            let got = model.score(context).unwrap();
            assert!((got - expected).abs() < 1e-3, "{got} != {expected}");
        }
    }

    #[test]
    fn invalid_context_rejected() {
        let model = Cpf1Model::embedded();
        assert!(model.score(&[b'N'; CONTEXT_LEN]).is_none());
        assert!(model.score(b"ACGT").is_none());
    }
}
