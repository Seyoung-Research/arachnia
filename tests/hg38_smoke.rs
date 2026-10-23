//! Optional hg38 smoke test.
//!
//! This is ignored by default because it requires a local hg38 FASTA and about
//! 3.1 GiB of RAM. Install one with `scripts/fetch_hg38.sh`, then run:
//!
//! ```text
//! cargo test --release --test hg38_smoke -- --ignored --nocapture
//! ```

use arachnia::offtarget::{
    find_hg38, search_offtargets, Genome, OffTargetParams, PamSide, SearchControl,
};
use arachnia::Strand;

#[test]
#[ignore = "requires a local hg38 and ~3.1 GiB of RAM"]
fn tp53_guide_is_found_in_hg38() {
    let path = find_hg38().expect("hg38 not found; run scripts/fetch_hg38.sh or set ARACHNIA_HG38");
    eprintln!("loading {}", path.display());
    let genome = Genome::from_path(&path).expect("load hg38");
    assert!(genome.total_bases() > 3_000_000_000);

    let params = OffTargetParams {
        max_mismatches: 0,
        pam_pattern: "NGG".into(),
        pam_side: PamSide::ThreePrime,
        both_strands: true,
    };
    let result = search_offtargets(
        "GAGCGCTGCTCAGATAGCGA",
        &genome,
        &params,
        &SearchControl::new(),
    )
    .expect("search hg38");

    let hit = result.sites.iter().find(|site| {
        site.contig == "chr17"
            && site.location == 7_675_051
            && site.strand == Strand::Minus
            && site.mismatches == 0
    });
    assert!(
        hit.is_some(),
        "expected TP53 on-target at chr17:7675051 (-); found {:?}",
        result.sites
    );
}
