use arachnia::ruleset3::{features, Tracr};

fn main() {
    let context = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "GACGAAAGCGACAACGCGTTCATCCGGGCA".into());
    let tracr = match std::env::args().nth(2).as_deref() {
        Some("Chen2013") => Tracr::Chen2013,
        _ => Tracr::Hsu2013,
    };
    let f = features(context.as_bytes(), tracr).expect("valid context");
    for (i, value) in f.iter().enumerate() {
        println!("{i}\t{value}");
    }
}
