#!/usr/bin/env python3
"""Export the RuleSet3 LightGBM model and SeqDeepCpf1 weights to compact
little-endian binaries embedded by the Rust application.

Run from the repository root inside a Python environment with joblib,
lightgbm, numpy, h5py and sglearn installed (see .refs/venv).

    python scripts/export_models.py
"""
from __future__ import annotations

import json
import struct
from pathlib import Path

import h5py
import joblib
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
REFS = ROOT / ".refs"
ASSETS = ROOT / "assets"
ASSETS.mkdir(exist_ok=True)

RS3_PICKLE = REFS / "rs3" / "rs3" / "RuleSet3.pkl"
CPF1_H5 = REFS / "seqdeepcpf1" / "weights" / "Seq_deepCpf1_weights.h5"


def export_rs3() -> None:
    model = joblib.load(RS3_PICKLE)
    booster = model.booster_
    dump = booster.dump_model()
    feature_names = dump["feature_names"]
    assert len(feature_names) == 632, len(feature_names)

    trees = []
    for tree in dump["tree_info"]:
        nodes = []

        def visit(node: dict) -> int:
            idx = len(nodes)
            if "leaf_value" in node:
                nodes.append((0xFFFF, -1, -1, 0.0, float(node["leaf_value"])))
            else:
                assert node["decision_type"] == "<=", node["decision_type"]
                nodes.append([None, None, None, float(node["threshold"]), 0.0])
                left = visit(node["left_child"])
                right = visit(node["right_child"])
                nodes[idx][0] = int(node["split_feature"])
                nodes[idx][1] = left
                nodes[idx][2] = right
            return idx

        visit(tree["tree_structure"])
        trees.append(nodes)

    out = ASSETS / "rs3_model.bin"
    with out.open("wb") as fh:
        fh.write(b"RS3LGBM1")
        fh.write(struct.pack("<IIII", 1, len(trees), len(feature_names),
                             sum(len(t) for t in trees)))
        for nodes in trees:
            fh.write(struct.pack("<I", len(nodes)))
            for feat, left, right, threshold, value in nodes:
                fh.write(struct.pack("<HHiidd", feat, 0, left, right,
                                     threshold, value))
    print(f"wrote {out} ({out.stat().st_size / 1e6:.1f} MB, "
          f"{len(trees)} trees, {sum(len(t) for t in trees)} nodes)")

    # Save the feature order next to the binary for documentation/tests.
    (ASSETS / "rs3_features.json").write_text(json.dumps(feature_names, indent=0))


def export_cpf1() -> None:
    with h5py.File(CPF1_H5, "r") as fh:
        conv_w = fh["convolution1d_157"]["convolution1d_157_W"][:]
        conv_b = fh["convolution1d_157"]["convolution1d_157_b"][:]
        fc1_w = fh["dense_490"]["dense_490_W"][:]
        fc1_b = fh["dense_490"]["dense_490_b"][:]
        fc2_w = fh["dense_491"]["dense_491_W"][:]
        fc2_b = fh["dense_491"]["dense_491_b"][:]
        fc3_w = fh["dense_492"]["dense_492_W"][:]
        fc3_b = fh["dense_492"]["dense_492_b"][:]
        out_w = fh["dense_493"]["dense_493_W"][:]
        out_b = fh["dense_493"]["dense_493_b"][:]

    # Match the PyTorch conversion in SeqDeepCpf1.py exactly.
    conv_w = np.squeeze(conv_w, axis=1).transpose(2, 1, 0)[..., ::-1].copy()
    conv_b = np.asarray(conv_b).copy()
    # Keras dense kernels are (in, out); PyTorch stores (out, in).
    fc1_w = fc1_w.T.copy()
    fc2_w = fc2_w.T.copy()
    fc3_w = fc3_w.T.copy()
    out_w = out_w.T.copy()

    arrays = [conv_w, conv_b, fc1_w, fc1_b, fc2_w, fc2_b,
              fc3_w, fc3_b, out_w, out_b]
    arrays = [np.ascontiguousarray(a, dtype="<f4") for a in arrays]

    out = ASSETS / "seqdeepcpf1_weights.bin"
    with out.open("wb") as fh:
        fh.write(b"CPF1W001")
        fh.write(struct.pack("<III", 80, 4, 5))
        for array in arrays:
            fh.write(array.tobytes())
    print(f"wrote {out} ({out.stat().st_size / 1e6:.2f} MB)")

    # Validate the exported tensors with a NumPy reimplementation, using the
    # two examples shipped with SeqDeepCpf1-Py3.
    x = np.zeros((2, 4, 34), dtype=np.float32)
    seqs = ["TGACTTTGAATGGAGTCGTGAGCGCAAGAACGCT",
            "GTTATTTGAGCAATGCCACTTAATAAACATGTAA"]
    bases = {"A": 0, "C": 1, "G": 2, "T": 3}
    for n, seq in enumerate(seqs):
        for i, base in enumerate(seq):
            x[n, bases[base], i] = 1.0

    conv = np.zeros((2, 80, 30), dtype=np.float32)
    for n in range(2):
        for o in range(80):
            for p in range(30):
                conv[n, o, p] = conv_b[o] + np.sum(
                    conv_w[o] * x[n, :, p:p + 5])
    relu = np.maximum(conv, 0)
    pooled = relu.reshape(2, 80, 15, 2).mean(axis=3)
    flat = pooled.transpose(0, 2, 1).reshape(2, 1200)
    h1 = np.maximum(flat @ fc1_w.T + fc1_b, 0)
    h2 = np.maximum(h1 @ fc2_w.T + fc2_b, 0)
    h3 = np.maximum(h2 @ fc3_w.T + fc3_b, 0)
    y = (h3 @ out_w.T + out_b).ravel()
    print("SeqDeepCpf1 numpy check:", y.tolist())


def export_rs3_tests() -> None:
    """Create regression vectors: contexts, tracr, expected RS3 score."""
    import sglearn

    model = joblib.load(RS3_PICKLE)
    booster = model.booster_

    contexts = [
        "GACGAAAGCGACAACGCGTTCATCCGGGCA",
        "AGAAAACACTAGCATCCCCACCCGCGGACT",
    ]
    rng = np.random.default_rng(12345)
    for _ in range(40):
        contexts.append("".join(rng.choice(list("ACGT"), size=30)))

    cases = []
    for tracr in ("Hsu2013", "Chen2013"):
        df = sglearn.featurize_guides(contexts, n_jobs=1)
        df["Hsu2013 tracr"] = int(tracr == "Hsu2013")
        df["Chen2013 tracr"] = int(tracr == "Chen2013")
        scores = booster.predict(df.to_numpy())
        for context, score in zip(contexts, scores):
            cases.append({"context": context, "tracr": tracr,
                          "score": float(score)})
    (ASSETS / "rs3_tests.json").write_text(json.dumps(cases, indent=1))
    print(f"wrote {len(cases)} RuleSet3 regression cases")


if __name__ == "__main__":
    export_rs3()
    export_cpf1()
    export_rs3_tests()
