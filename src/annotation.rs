//! GTF/GFF3 genome annotation for chromosome-region labels.
//!
//! [`Annotation`] parses gene, transcript, exon, CDS and UTR records from a
//! GTF or GFF3 file and answers region queries for off-target sites. Results
//! are classified as promoter, exon, coding exon, 5′/3′ UTR, intron, gene body
//! or intergenic.
//!
//! # Promoter definition
//!
//! A promoter is the interval from **5 kb to 500 bp upstream of a
//! transcription start site (TSS)**. Both edges are measured in the direction
//! of transcription, so a minus-strand promoter lies at higher plus-strand
//! coordinates than its TSS. The 500 bp immediately upstream of the TSS and
//! the region beyond 5 kb are not promoters.
//!
//! # Region precedence
//!
//! A query may overlap several transcripts. The most specific region wins:
//! coding exon, 5′ UTR, 3′ UTR, exon, promoter, intron, gene body, then
//! intergenic. Ties are broken by the shortest interval, the closest TSS and
//! finally the gene name.
//!
//! # Input
//!
//! GTF and GFF3 are detected from the attribute syntax. Both 1-based inclusive
//! formats are accepted; all coordinates exposed by this module are 0-based
//! half-open. GFF3 `Parent` chains are resolved, and 5′/3′ UTRs are inferred
//! from CDS and exon intervals when the file only contains `UTR` records or no
//! explicit UTR records at all. Compressed annotation files are not supported;
//! `scripts/fetch_annotation.sh` downloads and decompresses a GENCODE
//! annotation into the standard data directory.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Cursor};
use std::path::{Path, PathBuf};

use crate::design::Strand;

/// Distance from a transcription start site to the far edge of a promoter,
/// in bases (5 kb).
pub const PROMOTER_UPSTREAM_START: u64 = 5_000;

/// Distance from a transcription start site to the near edge of a promoter,
/// in bases (500 bp). The interval is exclusive at this edge.
pub const PROMOTER_UPSTREAM_END: u64 = 500;

/// Environment variable pointing at a GTF/GFF3 annotation file.
pub const ANNOTATION_ENV: &str = "ARACHNIA_ANNOTATION";

/// Environment variable pointing at a directory that may contain an
/// annotation. `scripts/fetch_annotation.sh` uses this as its destination.
pub const ANNOTATION_DIR_ENV: &str = "ARACHNIA_ANNOTATION_DIR";

/// A chromosome-region category assigned to a site.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Region {
    /// Outside every annotated gene and promoter.
    #[default]
    Intergenic,
    /// Inside a gene but not covered by an exon, intron or promoter.
    GeneBody,
    /// Between two exons of a transcript.
    Intron,
    /// 500 bp to 5 kb upstream of a transcription start site.
    Promoter,
    /// An exon that is neither coding nor annotated as a UTR.
    Exon,
    /// An exon interval annotated as, or inferred to be, a 3′ UTR.
    ThreePrimeUtr,
    /// An exon interval annotated as, or inferred to be, a 5′ UTR.
    FivePrimeUtr,
    /// A coding sequence interval inside an exon.
    CodingExon,
}

impl Region {
    /// Human-readable label used in the GUI and TSV export.
    pub fn label(self) -> &'static str {
        match self {
            Region::Intergenic => "intergenic",
            Region::GeneBody => "gene body",
            Region::Intron => "intron",
            Region::Promoter => "promoter",
            Region::Exon => "exon",
            Region::ThreePrimeUtr => "3′ UTR",
            Region::FivePrimeUtr => "5′ UTR",
            Region::CodingExon => "coding exon",
        }
    }

    fn priority(self) -> u8 {
        match self {
            Region::Intergenic => 0,
            Region::GeneBody => 1,
            Region::Intron => 2,
            Region::Promoter => 3,
            Region::Exon => 4,
            Region::ThreePrimeUtr => 5,
            Region::FivePrimeUtr => 6,
            Region::CodingExon => 7,
        }
    }
}

/// The region and gene assigned to one genomic interval or off-target site.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SiteAnnotation {
    /// Most specific region overlapping the query.
    pub region: Region,
    /// Gene whose transcript model produced the region, if any.
    pub gene: Option<String>,
    /// Transcript whose model produced the region, if any.
    pub transcript: Option<String>,
}

/// Error returned while reading or validating an annotation file.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AnnotationError(
    /// Human-readable error message.
    pub String,
);

impl std::fmt::Display for AnnotationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for AnnotationError {}

/// Parsed gene models for one or more reference contigs.
///
/// Build an annotation with [`Annotation::from_path`], [`Annotation::from_text`]
/// or [`Annotation::from_reader`], then call [`Annotation::annotate_interval`].
#[derive(Clone, Debug, Default)]
pub struct Annotation {
    contigs: HashMap<String, ContigAnnotation>,
    gene_names: Vec<String>,
    transcript_names: Vec<String>,
    feature_count: usize,
}

#[derive(Clone, Debug, Default)]
struct ContigAnnotation {
    features: Vec<Feature>,
    prefix_max_end: Vec<u64>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct Feature {
    start: u64,
    end: u64,
    region: Region,
    gene: usize,
    transcript: Option<usize>,
    tss: Option<u64>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum RawKind {
    Gene,
    Transcript,
    Exon,
    Cds,
    FivePrimeUtr,
    ThreePrimeUtr,
    Utr,
}

#[derive(Clone, Debug)]
struct RawFeature {
    contig: String,
    kind: RawKind,
    start: u64,
    end: u64,
    strand: Strand,
    gene_id: Option<String>,
    gene_name: Option<String>,
    transcript_id: Option<String>,
    transcript_name: Option<String>,
    tags: Vec<String>,
}

#[derive(Default)]
struct GeneBuilder {
    contig: String,
    name: String,
    strand: Option<Strand>,
    start: u64,
    end: u64,
    seen: bool,
}

#[derive(Default)]
struct TranscriptBuilder {
    contig: String,
    id: String,
    gene_id: String,
    name: Option<String>,
    strand: Option<Strand>,
    start: u64,
    end: u64,
    seen: bool,
    exons: Vec<(u64, u64)>,
    cds: Vec<(u64, u64)>,
    utr5: Vec<(u64, u64)>,
    utr3: Vec<(u64, u64)>,
    utrs: Vec<(u64, u64)>,
}

#[derive(Default)]
struct Attributes {
    gene_id: Option<String>,
    gene_name: Option<String>,
    transcript_id: Option<String>,
    transcript_name: Option<String>,
    id: Option<String>,
    parent: Option<String>,
    name: Option<String>,
    tags: Vec<String>,
}

impl Annotation {
    /// Load a GTF or GFF3 annotation file.
    ///
    /// # Errors
    ///
    /// Returns [`AnnotationError`] for unreadable files, compressed files,
    /// malformed records, or files without any usable gene/transcript/exon
    /// records.
    pub fn from_path(path: &Path) -> Result<Self, AnnotationError> {
        if path
            .extension()
            .and_then(|ext| ext.to_str())
            .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"))
        {
            return Err(AnnotationError(
                "compressed annotation is not supported; decompress the GTF/GFF3 first".into(),
            ));
        }
        let file = File::open(path)
            .map_err(|err| AnnotationError(format!("cannot open {}: {err}", path.display())))?;
        let reader = BufReader::with_capacity(1 << 20, file);
        Self::from_reader(reader)
            .map_err(|err| AnnotationError(format!("{}: {err}", path.display())))
    }

    /// Parse GTF or GFF3 text.
    pub fn from_text(text: &str) -> Result<Self, AnnotationError> {
        Self::from_reader(Cursor::new(text.as_bytes()))
    }

    /// Parse GTF or GFF3 records from a buffered reader.
    ///
    /// # Errors
    ///
    /// Returns [`AnnotationError`] if a relevant record is malformed or the
    /// input contains no usable annotation.
    pub fn from_reader<R: BufRead>(reader: R) -> Result<Self, AnnotationError> {
        let mut raw = Vec::new();
        for (index, line) in reader.lines().enumerate() {
            let line = line.map_err(|err| AnnotationError(err.to_string()))?;
            let trimmed = line.trim();
            if trimmed == "##FASTA" {
                break;
            }
            if trimmed.is_empty() || trimmed.starts_with('#') || trimmed.starts_with('>') {
                continue;
            }
            parse_record(trimmed, index + 1, &mut raw)?;
        }
        build_annotation(raw)
    }

    /// Number of contigs with at least one annotated feature.
    pub fn contig_count(&self) -> usize {
        self.contigs.len()
    }

    /// Number of genes in the annotation.
    pub fn gene_count(&self) -> usize {
        self.gene_names.len()
    }

    /// Number of transcripts in the annotation.
    pub fn transcript_count(&self) -> usize {
        self.transcript_names.len()
    }

    /// Number of interval features stored after merging and deduplication.
    pub fn feature_count(&self) -> usize {
        self.feature_count
    }

    /// Annotate a 0-based half-open interval.
    ///
    /// The interval is matched against exons, introns, promoters and gene
    /// bodies on the named contig. An empty interval or an unknown contig is
    /// reported as intergenic.
    pub fn annotate_interval(&self, contig: &str, start: u64, end: u64) -> SiteAnnotation {
        let Some(annotation) = self.contigs.get(contig) else {
            return SiteAnnotation::default();
        };
        if start >= end || annotation.features.is_empty() {
            return SiteAnnotation::default();
        }
        let midpoint = start + (end - start) / 2;
        let limit = annotation
            .features
            .partition_point(|feature| feature.start < end);
        let mut best: Option<&Feature> = None;
        for index in (0..limit).rev() {
            if annotation.prefix_max_end[index] <= start {
                break;
            }
            let feature = &annotation.features[index];
            if feature.end <= start || feature.start >= end {
                continue;
            }
            let better = match best {
                None => true,
                Some(current) => {
                    let candidate_priority = feature.region.priority();
                    let current_priority = current.region.priority();
                    if candidate_priority != current_priority {
                        candidate_priority > current_priority
                    } else {
                        let candidate_length = feature.end - feature.start;
                        let current_length = current.end - current.start;
                        if candidate_length != current_length {
                            candidate_length < current_length
                        } else {
                            let candidate_distance =
                                feature.tss.map_or(u64::MAX, |tss| midpoint.abs_diff(tss));
                            let current_distance =
                                current.tss.map_or(u64::MAX, |tss| midpoint.abs_diff(tss));
                            if candidate_distance != current_distance {
                                candidate_distance < current_distance
                            } else {
                                let candidate_gene = &self.gene_names[feature.gene];
                                let current_gene = &self.gene_names[current.gene];
                                if candidate_gene != current_gene {
                                    candidate_gene < current_gene
                                } else {
                                    let candidate_tx = feature
                                        .transcript
                                        .and_then(|index| self.transcript_names.get(index));
                                    let current_tx = current
                                        .transcript
                                        .and_then(|index| self.transcript_names.get(index));
                                    candidate_tx < current_tx
                                }
                            }
                        }
                    }
                }
            };
            if better {
                best = Some(feature);
            }
        }
        let Some(feature) = best else {
            return SiteAnnotation::default();
        };
        SiteAnnotation {
            region: feature.region,
            gene: self.gene_names.get(feature.gene).cloned(),
            transcript: feature
                .transcript
                .and_then(|index| self.transcript_names.get(index).cloned()),
        }
    }

    /// Annotate a single 0-based position.
    pub fn annotate_position(&self, contig: &str, position: u64) -> SiteAnnotation {
        self.annotate_interval(contig, position, position.saturating_add(1))
    }
}

fn parse_record(
    line: &str,
    line_number: usize,
    raw: &mut Vec<RawFeature>,
) -> Result<(), AnnotationError> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 9 {
        return Err(AnnotationError(format!(
            "line {line_number}: expected 9 tab-separated GTF/GFF3 fields, got {}",
            fields.len()
        )));
    }
    let Some(kind) = raw_kind(fields[2]) else {
        return Ok(());
    };
    let contig = fields[0].trim();
    if contig.is_empty() {
        return Err(AnnotationError(format!(
            "line {line_number}: missing sequence name"
        )));
    }
    let start: u64 = fields[3]
        .parse()
        .map_err(|_| AnnotationError(format!("line {line_number}: invalid start coordinate")))?;
    let end: u64 = fields[4]
        .parse()
        .map_err(|_| AnnotationError(format!("line {line_number}: invalid end coordinate")))?;
    if start == 0 || end < start {
        return Err(AnnotationError(format!(
            "line {line_number}: invalid 1-based interval {start}-{end}"
        )));
    }
    let strand = match fields[6].trim() {
        "+" => Strand::Plus,
        "-" => Strand::Minus,
        _ => return Ok(()),
    };
    let attributes = parse_attributes(fields[8]);
    let (gene_id, gene_name, transcript_id, transcript_name) = match kind {
        RawKind::Gene => (
            attributes.gene_id.clone().or_else(|| attributes.id.clone()),
            attributes
                .gene_name
                .clone()
                .or_else(|| attributes.name.clone()),
            None,
            None,
        ),
        RawKind::Transcript => (
            attributes
                .gene_id
                .clone()
                .or_else(|| attributes.parent.clone()),
            attributes.gene_name.clone(),
            attributes
                .transcript_id
                .clone()
                .or_else(|| attributes.id.clone()),
            attributes
                .transcript_name
                .clone()
                .or_else(|| attributes.name.clone()),
        ),
        _ => (
            attributes.gene_id.clone(),
            attributes.gene_name.clone(),
            attributes
                .transcript_id
                .clone()
                .or_else(|| attributes.parent.clone()),
            attributes.transcript_name.clone(),
        ),
    };
    if gene_id.is_none() && gene_name.is_none() && transcript_id.is_none() {
        return Ok(());
    }
    let transcript_ids: Vec<Option<String>> = transcript_id
        .as_deref()
        .map(|ids| {
            let ids: Vec<Option<String>> = ids
                .split(',')
                .map(str::trim)
                .filter(|id| !id.is_empty())
                .map(|id| Some(id.to_string()))
                .collect();
            if ids.is_empty() {
                vec![None]
            } else {
                ids
            }
        })
        .unwrap_or_else(|| vec![None]);
    for transcript_id in transcript_ids {
        raw.push(RawFeature {
            contig: contig.to_string(),
            kind,
            start: start - 1,
            end,
            strand,
            gene_id: gene_id.clone(),
            gene_name: gene_name.clone(),
            transcript_id,
            transcript_name: transcript_name.clone(),
            tags: attributes.tags.clone(),
        });
    }
    Ok(())
}

fn raw_kind(feature: &str) -> Option<RawKind> {
    match feature.trim().to_ascii_lowercase().as_str() {
        "gene" => Some(RawKind::Gene),
        "transcript" | "mrna" => Some(RawKind::Transcript),
        "exon" => Some(RawKind::Exon),
        "cds" => Some(RawKind::Cds),
        "five_prime_utr" | "5utr" | "5'utr" | "5' utr" | "five_prime_utr_region" => {
            Some(RawKind::FivePrimeUtr)
        }
        "three_prime_utr" | "3utr" | "3'utr" | "3' utr" | "three_prime_utr_region" => {
            Some(RawKind::ThreePrimeUtr)
        }
        "utr" => Some(RawKind::Utr),
        _ => None,
    }
}

fn parse_attributes(text: &str) -> Attributes {
    let mut attributes = Attributes::default();
    for field in text.split(';') {
        let field = field.trim();
        if field.is_empty() {
            continue;
        }
        let (key, value) = if let Some((key, value)) = field.split_once('=') {
            (key.trim(), value.trim())
        } else if let Some((key, value)) = field.split_once(char::is_whitespace) {
            (key.trim(), value.trim())
        } else {
            continue;
        };
        let value = decode_percent(value.trim_matches('"'));
        if value.is_empty() {
            continue;
        }
        match key.to_ascii_lowercase().as_str() {
            "gene_id" | "geneid" | "gene" => attributes.gene_id = Some(value),
            "gene_name" => attributes.gene_name = Some(value),
            "transcript_id" | "transcriptid" | "transcript" => {
                attributes.transcript_id = Some(value)
            }
            "transcript_name" => attributes.transcript_name = Some(value),
            "id" => attributes.id = Some(value),
            "parent" => attributes.parent = Some(value),
            "name" => attributes.name = Some(value),
            "tag" => attributes.tags.push(value),
            _ => {}
        }
    }
    attributes
}

fn decode_percent(value: &str) -> String {
    if !value.contains('%') {
        return value.to_string();
    }
    let bytes = value.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut index = 0;
    while index < bytes.len() {
        if bytes[index] == b'%' && index + 2 < bytes.len() {
            if let (Some(high), Some(low)) =
                (hex_value(bytes[index + 1]), hex_value(bytes[index + 2]))
            {
                out.push(high << 4 | low);
                index += 3;
                continue;
            }
        }
        out.push(bytes[index]);
        index += 1;
    }
    String::from_utf8_lossy(&out).into_owned()
}

fn hex_value(byte: u8) -> Option<u8> {
    match byte {
        b'0'..=b'9' => Some(byte - b'0'),
        b'a'..=b'f' => Some(byte - b'a' + 10),
        b'A'..=b'F' => Some(byte - b'A' + 10),
        _ => None,
    }
}

#[allow(clippy::too_many_arguments)]
fn ensure_gene(
    genes: &mut Vec<GeneBuilder>,
    index: &mut HashMap<(String, String), usize>,
    contig: &str,
    id: &str,
    name: Option<&str>,
    strand: Strand,
    start: u64,
    end: u64,
) -> usize {
    let key = (contig.to_string(), id.to_string());
    let index_value = if let Some(&index_value) = index.get(&key) {
        index_value
    } else {
        let index_value = genes.len();
        genes.push(GeneBuilder {
            contig: contig.to_string(),
            name: name.unwrap_or(id).to_string(),
            strand: Some(strand),
            start,
            end,
            seen: true,
        });
        index.insert(key, index_value);
        index_value
    };
    let gene = &mut genes[index_value];
    if !gene.seen {
        gene.contig = contig.to_string();
        gene.name = name.unwrap_or(id).to_string();
        gene.start = start;
        gene.end = end;
        gene.seen = true;
    } else {
        if let Some(name) = name.filter(|name| !name.is_empty()) {
            gene.name = name.to_string();
        }
        gene.start = gene.start.min(start);
        gene.end = gene.end.max(end);
    }
    if gene.strand.is_none() {
        gene.strand = Some(strand);
    }
    index_value
}

#[allow(clippy::too_many_arguments)]
fn ensure_transcript(
    transcripts: &mut Vec<TranscriptBuilder>,
    index: &mut HashMap<(String, String), usize>,
    contig: &str,
    id: &str,
    gene_id: &str,
    name: Option<&str>,
    strand: Strand,
    start: u64,
    end: u64,
) -> usize {
    let key = (contig.to_string(), id.to_string());
    let index_value = if let Some(&index_value) = index.get(&key) {
        index_value
    } else {
        let index_value = transcripts.len();
        transcripts.push(TranscriptBuilder {
            contig: contig.to_string(),
            id: id.to_string(),
            gene_id: gene_id.to_string(),
            name: name.map(str::to_string),
            strand: Some(strand),
            start,
            end,
            seen: true,
            ..TranscriptBuilder::default()
        });
        index.insert(key, index_value);
        index_value
    };
    let transcript = &mut transcripts[index_value];
    if !transcript.seen {
        transcript.contig = contig.to_string();
        transcript.gene_id = gene_id.to_string();
        if transcript.name.is_none() {
            transcript.name = name.map(str::to_string);
        }
        transcript.start = start;
        transcript.end = end;
        transcript.seen = true;
    } else {
        if gene_id != transcript.gene_id && !gene_id.is_empty() {
            transcript.gene_id = gene_id.to_string();
        }
        if transcript.name.is_none() {
            transcript.name = name.map(str::to_string);
        }
        transcript.start = transcript.start.min(start);
        transcript.end = transcript.end.max(end);
    }
    if transcript.strand.is_none() {
        transcript.strand = Some(strand);
    }
    index_value
}

fn build_annotation(raw_features: Vec<RawFeature>) -> Result<Annotation, AnnotationError> {
    if raw_features.is_empty() {
        return Err(AnnotationError(
            "no gene, transcript, exon, CDS or UTR records found".into(),
        ));
    }
    let mut genes: Vec<GeneBuilder> = Vec::new();
    let mut gene_index: HashMap<(String, String), usize> = HashMap::new();
    let mut transcripts: Vec<TranscriptBuilder> = Vec::new();
    let mut transcript_index: HashMap<(String, String), usize> = HashMap::new();
    let mut tx_gene: HashMap<(String, String), String> = HashMap::new();

    for raw in &raw_features {
        if let (Some(transcript_id), Some(gene_id)) = (&raw.transcript_id, &raw.gene_id) {
            tx_gene
                .entry((raw.contig.clone(), transcript_id.clone()))
                .or_insert_with(|| gene_id.clone());
        }
    }

    for raw in &raw_features {
        if raw.kind == RawKind::Gene {
            if let Some(gene_id) = raw.gene_id.as_deref().or(raw.gene_name.as_deref()) {
                ensure_gene(
                    &mut genes,
                    &mut gene_index,
                    &raw.contig,
                    gene_id,
                    raw.gene_name.as_deref(),
                    raw.strand,
                    raw.start,
                    raw.end,
                );
            }
        }
        if let Some(gene_id) = raw.gene_id.as_deref() {
            ensure_gene(
                &mut genes,
                &mut gene_index,
                &raw.contig,
                gene_id,
                raw.gene_name.as_deref(),
                raw.strand,
                raw.start,
                raw.end,
            );
        }
    }

    for raw in &raw_features {
        if raw.kind == RawKind::Gene {
            continue;
        }
        let transcript_id = match &raw.transcript_id {
            Some(id) if !id.is_empty() => Some(id.clone()),
            _ if raw.kind == RawKind::Transcript => continue,
            _ => None,
        };
        let gene_id = raw
            .gene_id
            .clone()
            .or_else(|| {
                transcript_id
                    .as_ref()
                    .and_then(|id| tx_gene.get(&(raw.contig.clone(), id.clone())).cloned())
            })
            .or_else(|| raw.gene_name.clone())
            .or_else(|| transcript_id.clone());
        let Some(gene_id) = gene_id else {
            continue;
        };
        let transcript_key = transcript_id.clone().unwrap_or_else(|| gene_id.clone());
        ensure_gene(
            &mut genes,
            &mut gene_index,
            &raw.contig,
            &gene_id,
            raw.gene_name.as_deref(),
            raw.strand,
            raw.start,
            raw.end,
        );
        let transcript_index_value = ensure_transcript(
            &mut transcripts,
            &mut transcript_index,
            &raw.contig,
            &transcript_key,
            &gene_id,
            raw.transcript_name.as_deref(),
            raw.strand,
            raw.start,
            raw.end,
        );
        let transcript = &mut transcripts[transcript_index_value];
        match raw.kind {
            RawKind::Exon => transcript.exons.push((raw.start, raw.end)),
            RawKind::Cds => transcript.cds.push((raw.start, raw.end)),
            RawKind::FivePrimeUtr => transcript.utr5.push((raw.start, raw.end)),
            RawKind::ThreePrimeUtr => transcript.utr3.push((raw.start, raw.end)),
            RawKind::Utr => {
                let five = raw.tags.iter().any(|tag| {
                    matches!(
                        tag.to_ascii_lowercase().as_str(),
                        "5utr" | "five_prime_utr" | "5'utr"
                    )
                });
                let three = raw.tags.iter().any(|tag| {
                    matches!(
                        tag.to_ascii_lowercase().as_str(),
                        "3utr" | "three_prime_utr" | "3'utr"
                    )
                });
                if five {
                    transcript.utr5.push((raw.start, raw.end));
                } else if three {
                    transcript.utr3.push((raw.start, raw.end));
                } else {
                    transcript.utrs.push((raw.start, raw.end));
                }
            }
            RawKind::Gene | RawKind::Transcript => {}
        }
    }

    if genes.is_empty() {
        return Err(AnnotationError(
            "annotation has no genes to classify sites".into(),
        ));
    }

    let mut contig_features: HashMap<String, Vec<Feature>> = HashMap::new();
    let mut genes_with_promoter: HashSet<usize> = HashSet::new();

    for (transcript_index_value, transcript) in transcripts.iter().enumerate() {
        let gene_key = (transcript.contig.clone(), transcript.gene_id.clone());
        let Some(&gene_index_value) = gene_index.get(&gene_key) else {
            continue;
        };
        genes_with_promoter.insert(gene_index_value);
        let gene = &genes[gene_index_value];
        let strand = transcript.strand.or(gene.strand).unwrap_or(Strand::Plus);
        let tss = tss_of(transcript.start, transcript.end, strand);
        let mut exons = merge_intervals(transcript.exons.clone());
        if exons.is_empty() {
            exons = merge_intervals(transcript.cds.clone());
        }
        if exons.is_empty() && transcript.start < transcript.end {
            exons.push((transcript.start, transcript.end));
        }
        let cds = merge_intervals(transcript.cds.clone());
        let mut utr5 = merge_intervals(transcript.utr5.clone());
        let mut utr3 = merge_intervals(transcript.utr3.clone());
        let cds_bounds = cds
            .first()
            .map(|interval| interval.0)
            .zip(cds.last().map(|interval| interval.1));
        if let Some((cds_start, cds_end)) = cds_bounds {
            for (start, end) in &transcript.utrs {
                let is_five = match strand {
                    Strand::Plus => *end <= cds_start,
                    Strand::Minus => *start >= cds_end,
                };
                let is_three = match strand {
                    Strand::Plus => *start >= cds_end,
                    Strand::Minus => *end <= cds_start,
                };
                if is_five {
                    utr5.push((*start, *end));
                } else if is_three {
                    utr3.push((*start, *end));
                }
            }
            let (five_range, three_range) = match strand {
                Strand::Plus => ((transcript.start, cds_start), (cds_end, transcript.end)),
                Strand::Minus => ((cds_end, transcript.end), (transcript.start, cds_start)),
            };
            utr5.extend(intersect_intervals(&exons, five_range.0, five_range.1));
            utr3.extend(intersect_intervals(&exons, three_range.0, three_range.1));
        }
        utr5 = merge_intervals(utr5);
        utr3 = merge_intervals(utr3);

        let mut covered = cds.clone();
        covered.extend(utr5.iter().copied());
        covered.extend(utr3.iter().copied());
        let covered = merge_intervals(covered);
        let exon_only = subtract_intervals(&exons, &covered);

        for (start, end) in cds {
            push_feature(
                &mut contig_features,
                &transcript.contig,
                Feature {
                    start,
                    end,
                    region: Region::CodingExon,
                    gene: gene_index_value,
                    transcript: Some(transcript_index_value),
                    tss: Some(tss),
                },
            );
        }
        for (start, end) in utr5 {
            push_feature(
                &mut contig_features,
                &transcript.contig,
                Feature {
                    start,
                    end,
                    region: Region::FivePrimeUtr,
                    gene: gene_index_value,
                    transcript: Some(transcript_index_value),
                    tss: Some(tss),
                },
            );
        }
        for (start, end) in utr3 {
            push_feature(
                &mut contig_features,
                &transcript.contig,
                Feature {
                    start,
                    end,
                    region: Region::ThreePrimeUtr,
                    gene: gene_index_value,
                    transcript: Some(transcript_index_value),
                    tss: Some(tss),
                },
            );
        }
        for (start, end) in exon_only {
            push_feature(
                &mut contig_features,
                &transcript.contig,
                Feature {
                    start,
                    end,
                    region: Region::Exon,
                    gene: gene_index_value,
                    transcript: Some(transcript_index_value),
                    tss: Some(tss),
                },
            );
        }
        for pair in exons.windows(2) {
            if pair[0].1 < pair[1].0 {
                push_feature(
                    &mut contig_features,
                    &transcript.contig,
                    Feature {
                        start: pair[0].1,
                        end: pair[1].0,
                        region: Region::Intron,
                        gene: gene_index_value,
                        transcript: Some(transcript_index_value),
                        tss: Some(tss),
                    },
                );
            }
        }
        if let Some((start, end)) = promoter_interval(tss, strand) {
            push_feature(
                &mut contig_features,
                &transcript.contig,
                Feature {
                    start,
                    end,
                    region: Region::Promoter,
                    gene: gene_index_value,
                    transcript: Some(transcript_index_value),
                    tss: Some(tss),
                },
            );
        }
    }

    for (gene_index_value, gene) in genes.iter().enumerate() {
        let strand = gene.strand.unwrap_or(Strand::Plus);
        let tss = tss_of(gene.start, gene.end, strand);
        if gene.start < gene.end {
            push_feature(
                &mut contig_features,
                &gene.contig,
                Feature {
                    start: gene.start,
                    end: gene.end,
                    region: Region::GeneBody,
                    gene: gene_index_value,
                    transcript: None,
                    tss: Some(tss),
                },
            );
        }
        if !genes_with_promoter.contains(&gene_index_value) {
            if let Some((start, end)) = promoter_interval(tss, strand) {
                push_feature(
                    &mut contig_features,
                    &gene.contig,
                    Feature {
                        start,
                        end,
                        region: Region::Promoter,
                        gene: gene_index_value,
                        transcript: None,
                        tss: Some(tss),
                    },
                );
            }
        }
    }

    let mut feature_count = 0;
    let contigs = contig_features
        .into_iter()
        .filter_map(|(contig, mut features)| {
            features.sort_unstable_by(|a, b| {
                a.start
                    .cmp(&b.start)
                    .then_with(|| a.end.cmp(&b.end))
                    .then_with(|| (a.region as u8).cmp(&(b.region as u8)))
                    .then_with(|| a.gene.cmp(&b.gene))
                    .then_with(|| a.transcript.cmp(&b.transcript))
            });
            features.dedup();
            if features.is_empty() {
                return None;
            }
            let mut prefix_max_end = Vec::with_capacity(features.len());
            let mut max_end = 0;
            for feature in &features {
                max_end = max_end.max(feature.end);
                prefix_max_end.push(max_end);
            }
            feature_count += features.len();
            Some((
                contig,
                ContigAnnotation {
                    features,
                    prefix_max_end,
                },
            ))
        })
        .collect();

    Ok(Annotation {
        contigs,
        gene_names: genes.into_iter().map(|gene| gene.name).collect(),
        transcript_names: transcripts
            .into_iter()
            .map(|transcript| transcript.name.unwrap_or(transcript.id))
            .collect(),
        feature_count,
    })
}

fn tss_of(start: u64, end: u64, strand: Strand) -> u64 {
    match strand {
        Strand::Plus => start,
        Strand::Minus => end.saturating_sub(1),
    }
}

fn promoter_interval(tss: u64, strand: Strand) -> Option<(u64, u64)> {
    match strand {
        Strand::Plus => {
            if tss <= PROMOTER_UPSTREAM_END {
                return None;
            }
            let start = tss.saturating_sub(PROMOTER_UPSTREAM_START);
            let end = tss - PROMOTER_UPSTREAM_END;
            (start < end).then_some((start, end))
        }
        Strand::Minus => {
            let start = tss.checked_add(PROMOTER_UPSTREAM_END)?;
            let end = tss.checked_add(PROMOTER_UPSTREAM_START)?;
            (start < end).then_some((start, end))
        }
    }
}

fn merge_intervals(mut intervals: Vec<(u64, u64)>) -> Vec<(u64, u64)> {
    intervals.retain(|(start, end)| start < end);
    intervals.sort_unstable();
    let mut merged: Vec<(u64, u64)> = Vec::with_capacity(intervals.len());
    for (start, end) in intervals {
        if let Some(last) = merged.last_mut() {
            if start <= last.1 {
                last.1 = last.1.max(end);
                continue;
            }
        }
        merged.push((start, end));
    }
    merged
}

fn intersect_intervals(
    intervals: &[(u64, u64)],
    range_start: u64,
    range_end: u64,
) -> Vec<(u64, u64)> {
    if range_start >= range_end {
        return Vec::new();
    }
    intervals
        .iter()
        .filter_map(|&(start, end)| {
            let clipped_start = start.max(range_start);
            let clipped_end = end.min(range_end);
            (clipped_start < clipped_end).then_some((clipped_start, clipped_end))
        })
        .collect()
}

fn subtract_intervals(intervals: &[(u64, u64)], covered: &[(u64, u64)]) -> Vec<(u64, u64)> {
    let mut result = Vec::new();
    for &(start, end) in intervals {
        let mut cursor = start;
        for &(covered_start, covered_end) in covered {
            if covered_end <= cursor {
                continue;
            }
            if covered_start >= end {
                break;
            }
            if covered_start > cursor {
                result.push((cursor, covered_start.min(end)));
            }
            cursor = cursor.max(covered_end);
            if cursor >= end {
                break;
            }
        }
        if cursor < end {
            result.push((cursor, end));
        }
    }
    result
}

fn push_feature(features: &mut HashMap<String, Vec<Feature>>, contig: &str, feature: Feature) {
    if feature.start < feature.end {
        features
            .entry(contig.to_string())
            .or_default()
            .push(feature);
    }
}

/// Paths searched for a local hg38 gene annotation, in priority order.
///
/// The first existing path wins. Set [`ANNOTATION_ENV`] to point at a specific
/// file, or [`ANNOTATION_DIR_ENV`] to search an additional directory.
/// `scripts/fetch_annotation.sh` installs a GENCODE annotation at the standard
/// XDG data path, where this search finds it.
pub fn annotation_candidates() -> Vec<PathBuf> {
    let mut candidates = Vec::new();
    if let Some(path) = std::env::var_os(ANNOTATION_ENV) {
        candidates.push(PathBuf::from(path));
    }
    if let Some(dir) = std::env::var_os(ANNOTATION_DIR_ENV) {
        candidates.extend(annotation_names_in(Path::new(&dir)));
    }
    if let Some(path) = std::env::var_os("ARACHNIA_HG38") {
        let path = PathBuf::from(path);
        if path.is_dir() {
            candidates.extend(annotation_names_in(&path));
        } else if let Some(parent) = path.parent() {
            candidates.extend(annotation_names_in(parent));
        }
    }
    if let Some(dir) = std::env::var_os("ARACHNIA_HG38_DIR") {
        candidates.extend(annotation_names_in(Path::new(&dir)));
    }
    if let Some(data_home) = std::env::var_os("XDG_DATA_HOME") {
        candidates.extend(annotation_names_in(
            &PathBuf::from(data_home).join("arachnia"),
        ));
    }
    if let Some(home) = std::env::var_os("HOME") {
        let home = PathBuf::from(home);
        candidates.extend(annotation_names_in(&home.join(".local/share/arachnia")));
        candidates.extend(annotation_names_in(&home.join("genomes")));
        candidates.push(home.join("genomes/hg38.gtf"));
    }
    candidates.push(PathBuf::from("/data/genomes/hg38.gtf"));
    candidates.push(PathBuf::from("/srv/genomes/hg38.gtf"));
    candidates
}

fn annotation_names_in(dir: &Path) -> Vec<PathBuf> {
    [
        "hg38.gtf",
        "hg38.gff3",
        "annotation.gtf",
        "annotation.gff3",
        "gencode.annotation.gtf",
        "gencode.annotation.gff3",
    ]
    .into_iter()
    .map(|name| dir.join(name))
    .collect()
}

/// Locate a local hg38 GTF/GFF3 annotation file, if one is installed.
pub fn find_hg38_annotation() -> Option<PathBuf> {
    annotation_candidates()
        .into_iter()
        .find(|path| path.is_file())
}

#[cfg(test)]
mod tests {
    use super::*;

    const PLUS_GTF: &str = "\
chr1\ttest\tgene\t10001\t14000\t.\t+\t.\tgene_id \"G1\"; gene_name \"GeneA\";
chr1\ttest\ttranscript\t10001\t14000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; gene_name \"GeneA\";
chr1\ttest\texon\t10001\t11000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";
chr1\ttest\tCDS\t10501\t11000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";
chr1\ttest\texon\t12501\t14000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";
chr1\ttest\tCDS\t12501\t13000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";
";

    const MINUS_GTF: &str = "\
chr2\ttest\tgene\t10001\t14000\t.\t-\t.\tgene_id \"G2\"; gene_name \"GeneB\";
chr2\ttest\ttranscript\t10001\t14000\t.\t-\t.\tgene_id \"G2\"; transcript_id \"T2\"; gene_name \"GeneB\";
chr2\ttest\texon\t10001\t11000\t.\t-\t.\tgene_id \"G2\"; transcript_id \"T2\";
chr2\ttest\texon\t12501\t14000\t.\t-\t.\tgene_id \"G2\"; transcript_id \"T2\";
chr2\ttest\tCDS\t10501\t13500\t.\t-\t.\tgene_id \"G2\"; transcript_id \"T2\";
";

    fn annotation() -> Annotation {
        Annotation::from_text(&format!("{PLUS_GTF}{MINUS_GTF}")).unwrap()
    }

    fn region(annotation: &Annotation, start: u64, end: u64) -> Region {
        annotation.annotate_interval("chr1", start, end).region
    }

    fn minus_region(annotation: &Annotation, start: u64, end: u64) -> Region {
        annotation.annotate_interval("chr2", start, end).region
    }

    #[test]
    fn classifies_plus_and_minus_gene_features() {
        let annotation = annotation();
        assert_eq!(annotation.gene_count(), 2);
        assert_eq!(annotation.transcript_count(), 2);
        assert_eq!(region(&annotation, 6_000, 6_010), Region::Promoter);
        assert_eq!(region(&annotation, 9_499, 9_500), Region::Promoter);
        assert_eq!(region(&annotation, 9_500, 9_510), Region::Intergenic);
        assert_eq!(region(&annotation, 10_000, 10_010), Region::FivePrimeUtr);
        assert_eq!(region(&annotation, 10_600, 10_610), Region::CodingExon);
        assert_eq!(region(&annotation, 11_200, 11_210), Region::Intron);
        assert_eq!(region(&annotation, 13_200, 13_210), Region::ThreePrimeUtr);
        assert_eq!(region(&annotation, 14_000, 14_010), Region::Intergenic);
        assert_eq!(minus_region(&annotation, 15_000, 15_010), Region::Promoter);
        assert_eq!(minus_region(&annotation, 18_998, 18_999), Region::Promoter);
        assert_eq!(
            minus_region(&annotation, 18_999, 19_000),
            Region::Intergenic
        );
        assert_eq!(
            minus_region(&annotation, 20_400, 20_410),
            Region::Intergenic
        );
    }

    #[test]
    fn minus_strand_utrs_follow_transcription_direction() {
        let annotation = annotation();
        assert_eq!(
            minus_region(&annotation, 10_200, 10_210),
            Region::ThreePrimeUtr
        );
        assert_eq!(
            minus_region(&annotation, 10_700, 10_710),
            Region::CodingExon
        );
        assert_eq!(
            minus_region(&annotation, 13_700, 13_710),
            Region::FivePrimeUtr
        );
        assert_eq!(
            annotation
                .annotate_interval("chr1", 6_000, 6_010)
                .gene
                .as_deref(),
            Some("GeneA")
        );
        assert_eq!(
            annotation
                .annotate_interval("chr2", 15_000, 15_010)
                .gene
                .as_deref(),
            Some("GeneB")
        );
    }

    #[test]
    fn parses_gff3_parent_chains() {
        let gff = "\
##gff-version 3
chr2\ttest\tgene\t10001\t14000\t.\t+\t.\tID=gene:G1;Name=GeneA
chr2\ttest\tmRNA\t10001\t14000\t.\t+\t.\tID=transcript:T1;Parent=gene:G1;Name=Tx1
chr2\ttest\texon\t10001\t11000\t.\t+\t.\tParent=transcript:T1
chr2\ttest\tCDS\t10501\t11000\t.\t+\t.\tParent=transcript:T1
chr2\ttest\texon\t12501\t14000\t.\t+\t.\tParent=transcript:T1
chr2\ttest\tCDS\t12501\t13000\t.\t+\t.\tParent=transcript:T1
";
        let annotation = Annotation::from_text(gff).unwrap();
        let hit = annotation.annotate_interval("chr2", 10_600, 10_610);
        assert_eq!(hit.region, Region::CodingExon);
        assert_eq!(hit.gene.as_deref(), Some("GeneA"));
        assert_eq!(hit.transcript.as_deref(), Some("Tx1"));
        assert_eq!(
            annotation.annotate_interval("chr2", 6_000, 6_010).region,
            Region::Promoter
        );
    }

    #[test]
    fn gencode_style_utr_tags_are_classified() {
        let gtf = "\
chr3\ttest\tgene\t10001\t14000\t.\t+\t.\tgene_id \"G3\"; gene_name \"GeneC\";
chr3\ttest\ttranscript\t10001\t14000\t.\t+\t.\tgene_id \"G3\"; transcript_id \"T3\";
chr3\ttest\texon\t10001\t11000\t.\t+\t.\tgene_id \"G3\"; transcript_id \"T3\";
chr3\ttest\tUTR\t10001\t10500\t.\t+\t.\tgene_id \"G3\"; transcript_id \"T3\"; tag \"5UTR\";
chr3\ttest\tCDS\t10501\t11000\t.\t+\t.\tgene_id \"G3\"; transcript_id \"T3\";
chr3\ttest\texon\t12501\t14000\t.\t+\t.\tgene_id \"G3\"; transcript_id \"T3\";
chr3\ttest\tCDS\t12501\t13000\t.\t+\t.\tgene_id \"G3\"; transcript_id \"T3\";
chr3\ttest\tUTR\t13001\t14000\t.\t+\t.\tgene_id \"G3\"; transcript_id \"T3\"; tag \"3UTR\";
";
        let annotation = Annotation::from_text(gtf).unwrap();
        assert_eq!(
            annotation.annotate_interval("chr3", 10_200, 10_210).region,
            Region::FivePrimeUtr
        );
        assert_eq!(
            annotation.annotate_interval("chr3", 13_200, 13_210).region,
            Region::ThreePrimeUtr
        );
    }

    #[test]
    fn gene_only_annotation_labels_body_and_promoter() {
        let gtf = "\
chr4\ttest\tgene\t10001\t14000\t.\t+\t.\tgene_id \"G4\"; gene_name \"GeneD\";
";
        let annotation = Annotation::from_text(gtf).unwrap();
        assert_eq!(
            annotation.annotate_interval("chr4", 6_000, 6_010).region,
            Region::Promoter
        );
        assert_eq!(
            annotation.annotate_interval("chr4", 10_500, 10_510).region,
            Region::GeneBody
        );
    }

    #[test]
    fn empty_annotation_is_rejected() {
        assert!(Annotation::from_text("# only comments\n").is_err());
    }
}
