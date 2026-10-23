//! End-to-end tests for the Cas-OFFinder-compatible command-line interface.

use std::fs;
use std::path::PathBuf;
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering};

static NEXT_DIR: AtomicUsize = AtomicUsize::new(0);

fn temp_dir(label: &str) -> PathBuf {
    let id = NEXT_DIR.fetch_add(1, Ordering::Relaxed);
    let dir =
        std::env::temp_dir().join(format!("arachnia-cli-{label}-{}-{id}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn run(input: &str, genome: &str, label: &str) -> Output {
    let dir = temp_dir(label);
    let genome_path = dir.join("genome.fa");
    let input_path = dir.join("input.txt");
    fs::write(&genome_path, genome).unwrap();
    fs::write(
        &input_path,
        input.replace("{genome}", &genome_path.display().to_string()),
    )
    .unwrap();
    Command::new(env!("CARGO_BIN_EXE_cas-offinder"))
        .arg(&input_path)
        .arg("-")
        .output()
        .unwrap()
}

fn rows(output: &Output) -> Vec<Vec<String>> {
    let stdout = String::from_utf8(output.stdout.clone()).unwrap();
    stdout
        .lines()
        .skip(2)
        .map(|line| line.split('\t').map(str::to_string).collect())
        .collect()
}

fn reverse_complement(sequence: &str) -> String {
    sequence
        .chars()
        .rev()
        .map(|base| match base {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            other => other,
        })
        .collect()
}

#[test]
fn cas9_query_reports_both_strands_and_lowercase_mismatches() {
    let guide = "ACGTACGTACGTACGTACGT";
    let mut mutated = guide.as_bytes().to_vec();
    mutated[5] = if mutated[5] == b'A' { b'C' } else { b'A' };
    let mutated = String::from_utf8(mutated).unwrap();

    let exact_window = format!("{guide}TGG");
    let mismatch_window = format!("{mutated}TGG");
    let minus_window = reverse_complement(&exact_window);
    let plus_site = 10;
    let minus_site = plus_site + exact_window.len() + 12;
    let mismatch_site = minus_site + minus_window.len() + 12;
    let mismatch_tail = format!("{mismatch_window}{}", "A".repeat(10));
    let genome = format!(
        ">chr1 test\n{}{}{}{}{}{}\n",
        "G".repeat(plus_site),
        exact_window,
        "C".repeat(12),
        minus_window,
        "C".repeat(12),
        mismatch_tail,
    );

    let input =
        "{genome}\nNNNNNNNNNNNNNNNNNNNNNGG\n".to_string() + &format!("{guide}NNN 1 cas9_query\n");
    let output = run(&input, &genome, "cas9");
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let rows = rows(&output);
    let find = |location: usize, strand: &str| {
        rows.iter()
            .find(|row| row[5] == location.to_string() && row[6] == strand)
            .unwrap_or_else(|| panic!("missing {strand} site at {location}: {rows:?}"))
    };

    let plus = find(plus_site, "+");
    assert_eq!(plus[0], "cas9_query");
    assert_eq!(plus[1], "X");
    assert_eq!(plus[2], format!("{guide}NNN"));
    assert_eq!(plus[3], exact_window);
    assert_eq!(plus[4], "chr1");
    assert_eq!(plus[7], "0");
    assert_eq!(plus[8], "0");

    let minus = find(minus_site, "-");
    assert_eq!(minus[3], exact_window);
    assert_eq!(minus[7], "0");

    let mismatch = find(mismatch_site, "+");
    assert_eq!(mismatch[7], "1");
    let mut expected = mutated.clone();
    let lower = expected[5..6].to_ascii_lowercase();
    expected.replace_range(5..6, &lower);
    assert_eq!(mismatch[3], format!("{expected}TGG"));
}

#[test]
fn cpf1_query_uses_a_five_prime_pam() {
    let guide = "AATGGAGTCGTGAGCGCAAGAAC";
    let pam = "TTTG";
    let plus_window = format!("{pam}{guide}");
    let minus_window = reverse_complement(&plus_window);
    let plus_site = 4;
    let minus_site = plus_site + plus_window.len() + 8;
    let genome = format!(
        ">chrCpf1 test\n{}{}{}{}{}\n",
        "A".repeat(plus_site),
        plus_window,
        "C".repeat(8),
        minus_window,
        "G".repeat(8),
    );

    let pattern = format!("TTTV{}", "N".repeat(guide.len()));
    let query = format!("NNNN{guide}");
    let input = format!("{{genome}}\n{pattern}\n{query} 0 cpf1_query\n");
    let output = run(&input, &genome, "cpf1");
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let rows = rows(&output);
    let find = |location: usize, strand: &str| {
        rows.iter()
            .find(|row| row[5] == location.to_string() && row[6] == strand)
            .unwrap_or_else(|| panic!("missing {strand} site at {location}: {rows:?}"))
    };

    for (location, strand) in [(plus_site, "+"), (minus_site, "-")] {
        let row = find(location, strand);
        assert_eq!(row[2], format!("NNNN{guide}"));
        assert_eq!(row[3], plus_window);
        assert_eq!(row[7], "0");
    }
}

#[test]
fn bulge_sizes_are_rejected() {
    let genome = ">chr1\nACGTACGTACGTACGTACGTACGTACGTACGT\n";
    let input = "{genome}\nNNNNNNNNNNNNNNNNNNNNNGG 0 1\nACGTACGTACGTACGTACGTNNN 0\n";
    let output = run(input, genome, "bulge");
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("bulges are not supported"), "{stderr}");
}

#[test]
fn hg38_alias_resolves_from_environment() {
    let dir = temp_dir("hg38");
    let genome_path = dir.join("tiny-hg38.fa");
    let guide = "ACGTACGTACGTACGTACGT";
    fs::write(&genome_path, format!(">chr1 test\nGGGG{guide}TGGAAAA\n")).unwrap();
    let input_path = dir.join("input.txt");
    fs::write(
        &input_path,
        format!("hg38\nNNNNNNNNNNNNNNNNNNNNNGG\n{guide}NNN 0 hg38_query\n"),
    )
    .unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_cas-offinder"))
        .arg(&input_path)
        .arg("-")
        .env("ARACHNIA_HG38", &genome_path)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("tiny-hg38.fa"), "{stderr}");
    let rows = rows(&output);
    assert_eq!(rows.len(), 1, "{rows:?}");
    assert_eq!(rows[0][0], "hg38_query");
    assert_eq!(rows[0][4], "chr1");
    assert_eq!(rows[0][7], "0");
}
