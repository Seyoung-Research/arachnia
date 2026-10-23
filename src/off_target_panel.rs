//! Off-target analysis window: reference loading, search controls and the
//! Cas-OFFinder-style mismatch table.
//!
//! The panel is opened from a selected on-target guide. It keeps the guide,
//! PAM orientation and default PAM pattern in sync with the active enzyme, and
//! supports two reference modes:
//!
//! * **Current target sequence** — search only the sequence in the main
//!   editor. This is useful for a quick plasmid or amplicon check.
//! * **Reference genome** — load a FASTA file or directory in a background
//!   thread, optionally using a locally discovered hg38.
//!
//! Searches run on a worker thread and communicate progress through an
//! [`arachnia::offtarget::SearchControl`]. The UI polls the control, renders a
//! progress bar, and can cancel an in-flight search. Results are held in a
//! sortable table with lower-case mismatch highlighting and can be copied or
//! exported as Cas-OFFinder-compatible TSV.
//!
//! An optional GTF/GFF3 annotation can be loaded alongside the reference. Each
//! site is then labelled with its chromosome region and overlapping gene. A
//! promoter is defined as 500 bp to 5 kb upstream of a transcription start
//! site. Install a GENCODE annotation with `scripts/fetch_annotation.sh`; the
//! field is pre-filled when one is discovered.
//!
//! The panel deliberately never loads hg38 automatically: a discovered
//! reference is preselected, but the user must click **Load** or **Use hg38**
//! because the assembly occupies about 3.1 GiB of RAM.

use std::fs;
use std::path::Path;
use std::sync::atomic::Ordering;
use std::sync::mpsc::{self, Receiver, TryRecvError};
use std::sync::Arc;
use std::time::{Duration, Instant};

use arachnia::annotation::{find_hg38_annotation, Annotation, Region};
use arachnia::offtarget::{
    expand_tilde, find_hg38, resolve_reference_path, results_tsv, search_offtargets,
    validate_search, Genome, OffTargetParams, OffTargetResult, PamSide, SearchControl, SearchError,
    HG38_ENV,
};
use arachnia::{prepare_sequence, Enzyme, GuideCandidate, Strand};
use eframe::egui::{self, Align, FontId, Layout, Margin, RichText, Stroke};
use egui_extras::{Column, TableBuilder};

use crate::theme::*;

struct SearchState {
    rx: Receiver<Result<OffTargetResult, SearchError>>,
    control: Arc<SearchControl>,
    started: Instant,
}

#[derive(Default)]
struct Actions {
    load_genome: bool,
    load_annotation: bool,
    clear_annotation: bool,
    search: bool,
    cancel: bool,
}

pub(crate) struct OffTargetPanel {
    open: bool,
    candidate: Option<GuideCandidate>,
    guide_id: String,
    max_mismatches: usize,
    pam_pattern: String,
    both_strands: bool,
    use_genome: bool,
    genome_path: String,
    genome: Option<Arc<Genome>>,
    genome_source: String,
    genome_summary: String,
    genome_load: Option<Receiver<(String, Result<Genome, String>)>>,
    genome_error: Option<String>,
    annotation_path: String,
    annotation: Option<Arc<Annotation>>,
    annotation_source: String,
    annotation_summary: String,
    annotation_load: Option<Receiver<(String, Result<Annotation, String>)>>,
    annotation_error: Option<String>,
    search: Option<SearchState>,
    result: Option<OffTargetResult>,
    result_pam_len: usize,
    result_pam_side: PamSide,
    reference_label: String,
    error: Option<String>,
    status: String,
    export_open: bool,
    export_path: String,
    export_error: Option<String>,
}

impl Default for OffTargetPanel {
    fn default() -> Self {
        let hg38 = find_hg38();
        let annotation = find_hg38_annotation();
        Self {
            open: false,
            candidate: None,
            guide_id: String::new(),
            max_mismatches: 3,
            pam_pattern: "NGG".to_string(),
            both_strands: true,
            use_genome: hg38.is_some(),
            genome_path: hg38
                .as_ref()
                .map(|path| path.display().to_string())
                .unwrap_or_default(),
            genome: None,
            genome_source: String::new(),
            genome_summary: String::new(),
            genome_load: None,
            genome_error: None,
            annotation_path: annotation
                .as_ref()
                .map(|path| path.display().to_string())
                .unwrap_or_default(),
            annotation: None,
            annotation_source: String::new(),
            annotation_summary: String::new(),
            annotation_load: None,
            annotation_error: None,
            search: None,
            result: None,
            result_pam_len: 3,
            result_pam_side: PamSide::ThreePrime,
            reference_label: String::new(),
            error: None,
            status: String::new(),
            export_open: false,
            export_path: "arachnia_offtargets.tsv".to_string(),
            export_error: None,
        }
    }
}

impl OffTargetPanel {
    pub(crate) fn open(&mut self, candidate: &GuideCandidate) {
        let same_guide = self.candidate.as_ref().is_some_and(|current| {
            current.guide == candidate.guide
                && current.start == candidate.start
                && current.strand == candidate.strand
                && current.enzyme == candidate.enzyme
        });
        self.cancel_search();
        self.search = None;
        if !same_guide {
            self.result = None;
        }
        self.error = None;
        self.export_error = None;
        self.status.clear();
        if self
            .candidate
            .as_ref()
            .is_none_or(|current| current.enzyme != candidate.enzyme)
        {
            self.pam_pattern = candidate.enzyme.pam_label().to_string();
        }
        self.candidate = Some(candidate.clone());
        self.guide_id = format!("guide_{}{}", candidate.start, candidate.strand.symbol());
        self.export_path = format!(
            "arachnia_offtargets_{}{}.tsv",
            candidate.start,
            candidate.strand.symbol()
        );
        self.open = true;
    }

    pub(crate) fn show(&mut self, ctx: &egui::Context, target_sequence: &str) {
        self.poll(ctx);
        if !self.open {
            return;
        }

        let mut open = self.open;
        let mut actions = Actions::default();
        egui::Window::new("Off-target analysis")
            .id(egui::Id::unique("offtarget_window"))
            .open(&mut open)
            .resizable(true)
            .collapsible(false)
            .default_size([1_220.0, 780.0])
            .min_width(1_080.0)
            .min_height(540.0)
            .show(ctx, |ui| {
                self.contents(ui, target_sequence, &mut actions);
            });
        self.open = open;
        if !open {
            self.cancel_search();
            self.search = None;
            self.export_open = false;
        }

        if actions.cancel {
            self.cancel_search();
        }
        if actions.load_genome {
            self.start_genome_load(ctx);
        }
        if actions.load_annotation {
            self.start_annotation_load(ctx);
        }
        if actions.clear_annotation {
            self.clear_annotation();
        }
        if actions.search {
            self.start_search(target_sequence, ctx);
        }
        if self.export_open {
            self.export_window(ctx);
        }
        if self.search.is_some() && self.open {
            ctx.request_repaint_after(Duration::from_millis(80));
        }
    }

    fn poll(&mut self, ctx: &egui::Context) {
        let loaded = self
            .genome_load
            .as_ref()
            .map(|receiver| receiver.try_recv());
        match loaded {
            Some(Ok((path, Ok(genome)))) => {
                self.genome_summary = format!(
                    "{} contig{} · {}",
                    genome.contig_count(),
                    if genome.contig_count() == 1 { "" } else { "s" },
                    format_bases(genome.total_bases() as u64)
                );
                self.genome = Some(Arc::new(genome));
                self.genome_source = path;
                self.genome_error = None;
                self.genome_load = None;
                self.status = "Reference genome loaded".to_string();
            }
            Some(Ok((path, Err(err)))) => {
                self.genome_error = Some(format!("Could not load {path}: {err}"));
                self.genome_load = None;
            }
            Some(Err(TryRecvError::Disconnected)) => {
                self.genome_error = Some("Genome loader stopped unexpectedly".to_string());
                self.genome_load = None;
            }
            Some(Err(TryRecvError::Empty)) => {
                ctx.request_repaint_after(Duration::from_millis(150));
            }
            None => {}
        }

        let annotated = self
            .annotation_load
            .as_ref()
            .map(|receiver| receiver.try_recv());
        match annotated {
            Some(Ok((path, Ok(annotation)))) => {
                self.annotation_summary = format!(
                    "{} genes · {} transcripts · {} regions",
                    annotation.gene_count(),
                    annotation.transcript_count(),
                    annotation.feature_count()
                );
                let annotation = Arc::new(annotation);
                self.annotation = Some(Arc::clone(&annotation));
                self.annotation_source = path;
                self.annotation_error = None;
                self.annotation_load = None;
                if let Some(result) = &mut self.result {
                    result.annotate(&annotation);
                }
                self.status = "Gene annotation loaded".to_string();
            }
            Some(Ok((path, Err(err)))) => {
                self.annotation_error = Some(format!("Could not load {path}: {err}"));
                self.annotation_load = None;
            }
            Some(Err(TryRecvError::Disconnected)) => {
                self.annotation_error = Some("Annotation loader stopped unexpectedly".to_string());
                self.annotation_load = None;
            }
            Some(Err(TryRecvError::Empty)) => {
                ctx.request_repaint_after(Duration::from_millis(150));
            }
            None => {}
        }

        let received = self
            .search
            .as_ref()
            .map(|state| (state.started, state.rx.try_recv()));
        match received {
            Some((started, Ok(Ok(mut result)))) => {
                if let Some(annotation) = &self.annotation {
                    result.annotate(annotation);
                }
                self.status = format!(
                    "{} site{} found in {:.2} s",
                    result.sites.len(),
                    if result.sites.len() == 1 { "" } else { "s" },
                    started.elapsed().as_secs_f64()
                );
                if result.truncated {
                    self.status.push_str(" (result limit reached)");
                }
                self.result = Some(result);
                self.search = None;
                self.error = None;
            }
            Some((_, Ok(Err(SearchError::Cancelled)))) => {
                self.status = "Search cancelled".to_string();
                self.search = None;
            }
            Some((_, Ok(Err(err)))) => {
                self.error = Some(err.to_string());
                self.search = None;
            }
            Some((_, Err(TryRecvError::Disconnected))) => {
                self.error = Some("Off-target search stopped unexpectedly".to_string());
                self.search = None;
            }
            Some((_, Err(TryRecvError::Empty))) => {
                ctx.request_repaint_after(Duration::from_millis(80));
            }
            None => {}
        }
    }

    fn contents(&mut self, ui: &mut egui::Ui, target_sequence: &str, actions: &mut Actions) {
        let Some(candidate) = self.candidate.clone() else {
            return;
        };

        ui.horizontal(|ui| {
            ui.vertical(|ui| {
                ui.horizontal(|ui| {
                    let guide_text = RichText::new(&candidate.guide)
                        .monospace()
                        .size(18.0)
                        .strong()
                        .color(TEXT);
                    let pam_text = RichText::new(&candidate.pam)
                        .monospace()
                        .size(18.0)
                        .strong()
                        .color(ACCENT_TEXT);
                    // Cpf1/Cas12a PAMs sit 5' of the guide.
                    if candidate.enzyme == Enzyme::Cpf1 {
                        ui.label(pam_text);
                        ui.label(guide_text);
                    } else {
                        ui.label(guide_text);
                        ui.label(pam_text);
                    }
                });
                ui.label(
                    RichText::new(format!(
                        "{} · {} strand · positions {}–{} · on-target score {}",
                        candidate.enzyme.label(),
                        candidate.strand.symbol(),
                        candidate.start,
                        candidate.end,
                        match candidate.enzyme {
                            Enzyme::Cas9 => format!("{:.3}", candidate.score),
                            Enzyme::Cpf1 => format!("{:.1}", candidate.score),
                        }
                    ))
                    .size(10.5)
                    .color(MUTED),
                );
            });
            ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                ui.label(
                    RichText::new("Cas-OFFinder-style exact-PAM search")
                        .size(10.5)
                        .color(MUTED),
                );
            });
        });
        ui.add_space(2.0);
        ui.separator();

        egui::Frame::new()
            .fill(CARD)
            .corner_radius(8)
            .stroke(Stroke::new(1.0, BORDER))
            .inner_margin(Margin::symmetric(12, 10))
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    ui.label(RichText::new("PAM pattern").size(12.0).color(MUTED));
                    ui.add(
                        egui::TextEdit::singleline(&mut self.pam_pattern)
                            .desired_width(64.0)
                            .font(egui::TextStyle::Monospace)
                            .char_limit(8)
                            .hint_text("NGG"),
                    );
                    ui.label(
                        RichText::new("IUPAC: N any, R A/G, Y C/T, V A/C/G…")
                            .size(9.5)
                            .color(MUTED),
                    );
                    ui.label(
                        RichText::new(match candidate.enzyme {
                            Enzyme::Cas9 => "3′ of guide",
                            Enzyme::Cpf1 => "5′ of guide",
                        })
                        .size(9.5)
                        .color(ACCENT_TEXT),
                    );
                    ui.add_space(12.0);
                    ui.label(RichText::new("Max mismatches").size(12.0).color(MUTED));
                    ui.add(egui::Slider::new(&mut self.max_mismatches, 0..=6).step_by(1.0));
                    ui.add_space(8.0);
                    ui.checkbox(&mut self.both_strands, "Both strands");
                });

                ui.add_space(4.0);
                ui.horizontal(|ui| {
                    ui.radio_value(
                        &mut self.use_genome,
                        false,
                        "Current target sequence",
                    );
                    ui.radio_value(&mut self.use_genome, true, "Reference genome");
                    if self.genome_load.is_none()
                        && ui
                            .button("Use hg38")
                            .on_hover_text(format!(
                                "Find and load hg38 from {HG38_ENV}, \
                                 ~/.local/share/arachnia/hg38.fa, or ~/genomes/hg38.fa"
                            ))
                            .clicked()
                    {
                        match find_hg38() {
                            Some(path) => {
                                self.genome_path = path.display().to_string();
                                self.use_genome = true;
                                actions.load_genome = true;
                            }
                            None => {
                                self.use_genome = true;
                                self.genome_error = Some(format!(
                                    "hg38 not found. Run scripts/fetch_hg38.sh or set {HG38_ENV}."
                                ));
                            }
                        }
                    }
                    if self.use_genome {
                        ui.add_space(6.0);
                        ui.add(
                            egui::TextEdit::singleline(&mut self.genome_path)
                                .desired_width(400.0)
                                .hint_text("/path/to/genome.fa, /path/to/fasta_dir, or hg38"),
                        );
                        if self.genome_load.is_some() {
                            ui.spinner();
                            ui.label(RichText::new("Loading…").size(10.5).color(MUTED));
                        } else if ui.button("Load").clicked() {
                            actions.load_genome = true;
                        }
                    }
                });

                if self.use_genome {
                    if let Some(error) = &self.genome_error {
                        ui.label(RichText::new(error).size(10.5).color(RED));
                    } else if !self.genome_summary.is_empty() {
                        ui.label(
                            RichText::new(format!(
                                "{} · {}",
                                self.genome_summary, self.genome_source
                            ))
                            .size(10.0)
                            .color(MUTED),
                        );
                    } else {
                        ui.label(
                            RichText::new(
                                "Load a FASTA file or a directory of chromosome FASTA files. A local hg38 is detected automatically; install it with scripts/fetch_hg38.sh.",
                            )
                            .size(10.0)
                            .color(MUTED),
                        );
                    }
                } else {
                    ui.label(
                        RichText::new(
                            "Searches only inside the target box. Load a reference genome for genome-wide analysis.",
                        )
                        .size(10.0)
                        .color(MUTED),
                    );
                }

                ui.add_space(6.0);
                ui.separator();
                ui.add_space(4.0);
                ui.horizontal(|ui| {
                    ui.label(
                        RichText::new("Gene annotation")
                            .size(12.0)
                            .color(MUTED),
                    );
                    ui.add(
                        egui::TextEdit::singleline(&mut self.annotation_path)
                            .desired_width(520.0)
                            .hint_text("/path/to/annotation.gtf or annotation.gff3"),
                    );
                    if self.annotation_load.is_some() {
                        ui.spinner();
                        ui.label(RichText::new("Loading…").size(10.5).color(MUTED));
                    } else if ui.button("Load").clicked() {
                        actions.load_annotation = true;
                    }
                    if self.annotation.is_some() && ui.button("Clear").clicked() {
                        actions.clear_annotation = true;
                    }
                });
                if let Some(error) = &self.annotation_error {
                    ui.label(RichText::new(error).size(10.5).color(RED));
                } else if !self.annotation_summary.is_empty() {
                    ui.label(
                        RichText::new(format!(
                            "{} · {} · promoter = 500 bp–5 kb upstream of a TSS",
                            self.annotation_summary, self.annotation_source
                        ))
                        .size(10.0)
                        .color(MUTED),
                    );
                } else {
                    ui.label(
                        RichText::new(
                            "Optional GTF/GFF3 gene models. Sites are labelled promoter, exon, coding exon, 5′/3′ UTR, intron, gene body or intergenic. Install a GENCODE file with scripts/fetch_annotation.sh.",
                        )
                        .size(10.0)
                        .color(MUTED),
                    );
                }
            });

        ui.add_space(4.0);
        let running = self.search.is_some();
        let genome_ready = !self.use_genome || self.genome.is_some();
        let target_ready = self.use_genome || !target_sequence.trim().is_empty();
        let can_search = !running && genome_ready && target_ready && !candidate.guide.is_empty();
        ui.horizontal(|ui| {
            let button = egui::Button::new(
                RichText::new("Find off-targets")
                    .size(13.5)
                    .strong()
                    .color(WHITE),
            )
            .fill(ACCENT)
            .corner_radius(7)
            .stroke(Stroke::NONE);
            if ui
                .add_enabled(can_search, button)
                .on_hover_text("Scan the reference for guide + PAM matches")
                .clicked()
            {
                actions.search = true;
            }
            if running && ui.button("Cancel").clicked() {
                actions.cancel = true;
            }
            if !genome_ready {
                ui.label(
                    RichText::new("Load a reference genome first.")
                        .size(10.5)
                        .color(AMBER),
                );
            }
            if let Some(search) = &self.search {
                let elapsed = search.started.elapsed().as_secs_f64();
                let progress = search.control.progress();
                ui.add(
                    egui::ProgressBar::new(progress)
                        .desired_width(250.0)
                        .text(format!("{:.0}%", progress * 100.0)),
                );
                ui.label(
                    RichText::new(format!(
                        "{} / {} windows · {:.1}s",
                        format_number(search.control.scanned.load(Ordering::Relaxed)),
                        format_number(search.control.total.load(Ordering::Relaxed)),
                        elapsed
                    ))
                    .size(10.0)
                    .color(MUTED),
                );
            }
        });

        if let Some(error) = &self.error {
            ui.add_space(2.0);
            ui.label(RichText::new(error).size(11.0).color(RED));
        }

        if let Some(result) = self.result.take() {
            ui.add_space(3.0);
            ui.separator();
            self.result_header(ui, &result, &candidate);
            ui.add_space(4.0);
            results_table(ui, &result);
            ui.add_space(3.0);
            ui.label(
                RichText::new(
                    "Locations are 0-based (Cas-OFFinder/Bowtie convention). PAM must match the pattern exactly; DNA/RNA bulges are not considered.",
                )
                .size(9.5)
                .color(MUTED),
            );
            if result.sites.iter().any(|site| site.annotation.is_some()) {
                ui.label(
                    RichText::new(
                        "Region labels use the full guide + PAM window. Promoter = 500 bp to 5 kb upstream of a transcription start site.",
                    )
                    .size(9.5)
                    .color(MUTED),
                );
            }
            self.result = Some(result);
        } else if !running {
            ui.add_space(26.0);
            ui.vertical_centered(|ui| {
                ui.label(
                    RichText::new("Find Cas-OFFinder-style off-target sites")
                        .size(18.0)
                        .strong()
                        .color(TEXT),
                );
                ui.add_space(4.0);
                ui.label(
                    RichText::new(
                        "Every site with a matching PAM and up to the selected number of guide mismatches is reported on both strands.",
                    )
                    .size(11.5)
                    .color(MUTED),
                );
                ui.add_space(8.0);
                ui.label(
                    RichText::new(
                        "Tip: load a chromosome FASTA file or a directory of FASTA files for genome-wide analysis.",
                    )
                    .size(10.5)
                    .color(ACCENT_TEXT),
                );
            });
        }

        if !self.status.is_empty() && self.error.is_none() {
            ui.add_space(2.0);
            ui.label(RichText::new(&self.status).size(10.0).color(MUTED));
        }
    }

    fn result_header(
        &mut self,
        ui: &mut egui::Ui,
        result: &OffTargetResult,
        candidate: &GuideCandidate,
    ) {
        ui.horizontal(|ui| {
            ui.vertical(|ui| {
                ui.label(
                    RichText::new(format!(
                        "{} off-target site{}",
                        result.sites.len(),
                        if result.sites.len() == 1 { "" } else { "s" }
                    ))
                    .size(18.0)
                    .strong()
                    .color(TEXT),
                );
                let mut summary = format!(
                    "{} · {} · {:.2} s · {} windows scanned",
                    self.reference_label,
                    format_bases(result.genome_bases as u64),
                    result.elapsed.as_secs_f64(),
                    format_number(result.scanned_positions)
                );
                for (mismatches, count) in result.counts_by_mismatch.iter().enumerate() {
                    if *count > 0 {
                        summary.push_str(&format!(" · {mismatches} mm: {count}"));
                    }
                }
                ui.label(RichText::new(summary).size(10.0).color(MUTED));
            });
            ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                if ui.button("Export TSV…").clicked() {
                    self.export_open = true;
                    self.export_error = None;
                }
                if ui.button("Copy TSV").clicked() {
                    let tsv = results_tsv(
                        &result.sites,
                        &candidate.guide,
                        self.result_pam_len,
                        self.result_pam_side,
                        &self.guide_id,
                    );
                    ui.ctx().copy_text(tsv);
                    self.status = format!("Copied {} off-target sites", result.sites.len());
                }
                if result.truncated {
                    ui.label(RichText::new("result cap reached").size(10.0).color(AMBER));
                }
                if result.perfect_matches() > 0 {
                    ui.label(
                        RichText::new(format!(
                            "{} perfect match{}",
                            result.perfect_matches(),
                            if result.perfect_matches() == 1 {
                                ""
                            } else {
                                "es"
                            }
                        ))
                        .size(10.0)
                        .color(ACCENT_TEXT),
                    );
                }
            });
        });
    }

    fn start_search(&mut self, target_sequence: &str, ctx: &egui::Context) {
        let Some(candidate) = &self.candidate else {
            return;
        };
        self.pam_pattern = self.pam_pattern.trim().to_ascii_uppercase();
        let pam_side = match candidate.enzyme {
            Enzyme::Cas9 => PamSide::ThreePrime,
            Enzyme::Cpf1 => PamSide::FivePrime,
        };
        let params = OffTargetParams {
            max_mismatches: self.max_mismatches,
            pam_pattern: self.pam_pattern.clone(),
            pam_side,
            both_strands: self.both_strands,
        };
        if let Err(err) = validate_search(&candidate.guide, &params) {
            self.error = Some(err.to_string());
            return;
        }

        let (genome, reference_label) = if self.use_genome {
            match &self.genome {
                Some(genome) => (
                    Arc::clone(genome),
                    if self.genome_source.is_empty() {
                        "reference genome".to_string()
                    } else {
                        self.genome_source.clone()
                    },
                ),
                None => {
                    self.error = Some("Load a reference genome first.".to_string());
                    return;
                }
            }
        } else {
            match prepare_sequence(target_sequence) {
                Ok(prepared) => {
                    let name = prepared.name.unwrap_or_else(|| "target".to_string());
                    match Genome::from_sequence(name, &prepared.bases) {
                        Ok(genome) => (Arc::new(genome), "current target sequence".to_string()),
                        Err(err) => {
                            self.error = Some(err.to_string());
                            return;
                        }
                    }
                }
                Err(err) => {
                    self.error = Some(err.to_string());
                    return;
                }
            }
        };

        let guide = candidate.guide.clone();
        let pam_len = params
            .pam_pattern
            .chars()
            .filter(|ch| !ch.is_whitespace())
            .count();
        let control = Arc::new(SearchControl::new());
        let thread_control = Arc::clone(&control);
        let (tx, rx) = mpsc::channel();
        std::thread::spawn(move || {
            let result = search_offtargets(&guide, &genome, &params, &thread_control);
            let _ = tx.send(result);
        });

        self.result = None;
        self.error = None;
        self.status = "Searching…".to_string();
        self.reference_label = reference_label;
        self.result_pam_len = pam_len;
        self.result_pam_side = pam_side;
        self.search = Some(SearchState {
            rx,
            control,
            started: Instant::now(),
        });
        ctx.request_repaint();
    }

    fn start_genome_load(&mut self, ctx: &egui::Context) {
        let raw_path = self.genome_path.trim().to_string();
        if raw_path.is_empty() {
            self.genome_error = Some("Enter a FASTA file or directory path first.".to_string());
            return;
        }
        let path = match resolve_reference_path(&raw_path) {
            Ok(path) => expand_tilde(path),
            Err(err) => {
                self.genome_error = Some(err);
                return;
            }
        };
        let path = path.display().to_string();
        let (tx, rx) = mpsc::channel();
        let path_for_thread = path.clone();
        std::thread::spawn(move || {
            let result =
                Genome::from_path(Path::new(&path_for_thread)).map_err(|err| err.to_string());
            let _ = tx.send((path, result));
        });
        self.genome_error = None;
        self.status = "Loading reference genome…".to_string();
        self.genome_load = Some(rx);
        ctx.request_repaint();
    }

    fn start_annotation_load(&mut self, ctx: &egui::Context) {
        let raw_path = self.annotation_path.trim().to_string();
        if raw_path.is_empty() {
            self.annotation_error = Some("Enter a GTF/GFF3 path first.".to_string());
            return;
        }
        let path = expand_tilde(std::path::PathBuf::from(raw_path));
        let path = path.display().to_string();
        let (tx, rx) = mpsc::channel();
        let path_for_thread = path.clone();
        std::thread::spawn(move || {
            let result =
                Annotation::from_path(Path::new(&path_for_thread)).map_err(|err| err.to_string());
            let _ = tx.send((path, result));
        });
        self.annotation_error = None;
        self.status = "Loading gene annotation…".to_string();
        self.annotation_load = Some(rx);
        ctx.request_repaint();
    }

    fn clear_annotation(&mut self) {
        self.annotation = None;
        self.annotation_load = None;
        self.annotation_source.clear();
        self.annotation_summary.clear();
        self.annotation_error = None;
        if let Some(result) = &mut self.result {
            result.clear_annotations();
        }
        self.status = "Gene annotation cleared".to_string();
    }

    fn cancel_search(&mut self) {
        if let Some(search) = &self.search {
            search.control.cancel();
        }
    }

    fn export_window(&mut self, ctx: &egui::Context) {
        let mut open = self.export_open;
        let mut save = false;
        let mut cancel = false;
        let count = self.result.as_ref().map_or(0, |result| result.sites.len());
        egui::Window::new("Export off-target sites")
            .id(egui::Id::unique("offtarget_export"))
            .open(&mut open)
            .resizable(false)
            .collapsible(false)
            .default_width(560.0)
            .show(ctx, |ui| {
                ui.label(
                    RichText::new(format!(
                        "{count} site{} in Cas-OFFinder-compatible TSV format{}:",
                        if count == 1 { "" } else { "s" },
                        if self.annotation.is_some() {
                            " with Region and Gene columns"
                        } else {
                            ""
                        }
                    ))
                    .size(11.5)
                    .color(MUTED),
                );
                ui.add(
                    egui::TextEdit::singleline(&mut self.export_path)
                        .desired_width(f32::INFINITY)
                        .hint_text("arachnia_offtargets.tsv"),
                );
                if let Some(error) = &self.export_error {
                    ui.label(RichText::new(error).size(10.5).color(RED));
                }
                ui.horizontal(|ui| {
                    if ui.button("Save TSV").clicked() {
                        save = true;
                    }
                    if ui.button("Cancel").clicked() {
                        cancel = true;
                        self.export_error = None;
                    }
                });
            });
        self.export_open = open && !cancel;
        if save {
            self.save_export();
        }
    }

    fn save_export(&mut self) {
        let Some(result) = &self.result else {
            self.export_error = Some("There are no sites to export.".to_string());
            return;
        };
        let Some(candidate) = &self.candidate else {
            return;
        };
        let path = self.export_path.trim().to_string();
        if path.is_empty() {
            self.export_error = Some("Enter a file path first.".to_string());
            return;
        }
        let tsv = results_tsv(
            &result.sites,
            &candidate.guide,
            self.result_pam_len,
            self.result_pam_side,
            &self.guide_id,
        );
        match fs::write(&path, tsv) {
            Ok(()) => {
                self.export_open = false;
                self.export_error = None;
                self.status = format!("Saved {} sites to {path}", result.sites.len());
            }
            Err(err) => {
                self.export_error = Some(format!("Could not write {path}: {err}"));
            }
        }
    }
}

fn results_table(ui: &mut egui::Ui, result: &OffTargetResult) {
    if result.sites.is_empty() {
        ui.add_space(20.0);
        ui.vertical_centered(|ui| {
            ui.label(
                RichText::new("No off-target sites found")
                    .size(17.0)
                    .strong()
                    .color(TEXT),
            );
            ui.add_space(3.0);
            ui.label(
                RichText::new(
                    "Try allowing more mismatches, a different PAM pattern, or another reference.",
                )
                .size(11.0)
                .color(MUTED),
            );
        });
        return;
    }

    let annotated = result.sites.iter().any(|site| site.annotation.is_some());
    egui::Frame::new()
        .fill(CARD)
        .corner_radius(8)
        .stroke(Stroke::new(1.0, BORDER))
        .inner_margin(Margin::same(0))
        .show(ui, |ui| {
            let mut table = TableBuilder::new(ui)
                .striped(true)
                .resizable(false)
                .auto_shrink([false, false])
                .cell_layout(Layout::left_to_right(Align::Center))
                .column(Column::exact(38.0))
                .column(Column::exact(120.0).clip(true))
                .column(Column::exact(88.0));
            if annotated {
                table = table
                    .column(Column::exact(86.0))
                    .column(Column::exact(118.0).clip(true));
            }
            table = table
                .column(Column::exact(44.0))
                .column(Column::exact(46.0))
                .column(Column::remainder().at_least(230.0).clip(true))
                .column(Column::exact(60.0));
            table
                .header(27.0, |mut header| {
                    let labels = if annotated {
                        vec![
                            "#",
                            "Contig",
                            "Location",
                            "Region",
                            "Gene",
                            "Strand",
                            "MM",
                            "Guide (guide orientation)",
                            "PAM",
                        ]
                    } else {
                        vec![
                            "#",
                            "Contig",
                            "Location",
                            "Strand",
                            "MM",
                            "Guide (guide orientation)",
                            "PAM",
                        ]
                    };
                    for label in labels {
                        header.col(|ui| {
                            ui.label(RichText::new(label).size(10.0).strong().color(MUTED));
                        });
                    }
                })
                .body(|body| {
                    body.rows(22.0, result.sites.len(), |mut row| {
                        let index = row.index();
                        let site = &result.sites[index];
                        row.col(|ui| {
                            ui.label(
                                RichText::new((index + 1).to_string())
                                    .size(10.0)
                                    .color(MUTED),
                            );
                        });
                        row.col(|ui| {
                            ui.label(RichText::new(&site.contig).size(10.5).color(MUTED));
                        });
                        row.col(|ui| {
                            ui.label(
                                RichText::new(site.location.to_string())
                                    .monospace()
                                    .size(11.0)
                                    .color(MUTED),
                            );
                        });
                        if annotated {
                            row.col(|ui| {
                                if let Some(annotation) = &site.annotation {
                                    ui.label(
                                        RichText::new(annotation.region.label())
                                            .size(10.5)
                                            .color(region_color(annotation.region)),
                                    );
                                }
                            });
                            row.col(|ui| {
                                if let Some(annotation) = &site.annotation {
                                    ui.label(
                                        RichText::new(annotation.gene.as_deref().unwrap_or("—"))
                                            .size(10.5)
                                            .color(MUTED),
                                    );
                                }
                            });
                        }
                        row.col(|ui| {
                            ui.label(
                                RichText::new(site.strand.symbol().to_string())
                                    .monospace()
                                    .size(11.5)
                                    .color(match site.strand {
                                        Strand::Plus => GREEN,
                                        Strand::Minus => AMBER,
                                    }),
                            );
                        });
                        row.col(|ui| {
                            let color = match site.mismatches {
                                0 => ACCENT_TEXT,
                                1..=2 => AMBER,
                                _ => MUTED,
                            };
                            ui.label(
                                RichText::new(site.mismatches.to_string())
                                    .monospace()
                                    .size(11.5)
                                    .strong()
                                    .color(color),
                            );
                        });
                        row.col(|ui| {
                            ui.label(site_dna_job(ui, &site.guide_dna));
                        });
                        row.col(|ui| {
                            ui.label(
                                RichText::new(&site.pam)
                                    .monospace()
                                    .size(11.5)
                                    .color(ACCENT_TEXT),
                            );
                        });
                    });
                });
        });
}

fn region_color(region: Region) -> egui::Color32 {
    match region {
        Region::Promoter => AMBER,
        Region::CodingExon | Region::Exon | Region::FivePrimeUtr | Region::ThreePrimeUtr => GREEN,
        Region::Intron | Region::GeneBody | Region::Intergenic => MUTED,
    }
}

fn site_dna_job(ui: &egui::Ui, dna: &str) -> egui::text::LayoutJob {
    use egui::text::{LayoutJob, TextFormat};
    let mut job = LayoutJob::default();
    let font = FontId::monospace(11.5);
    let matched = ui.visuals().text_color();
    for ch in dna.chars() {
        let color = if ch.is_ascii_lowercase() {
            ACCENT_TEXT
        } else {
            matched
        };
        job.append(
            &ch.to_string(),
            0.0,
            TextFormat {
                font_id: font.clone(),
                color,
                ..Default::default()
            },
        );
    }
    job
}

fn format_number(value: u64) -> String {
    let text = value.to_string();
    let mut out = String::with_capacity(text.len() + text.len() / 3);
    for (index, ch) in text.chars().enumerate() {
        if index > 0 && (text.len() - index).is_multiple_of(3) {
            out.push(',');
        }
        out.push(ch);
    }
    out
}

fn format_bases(value: u64) -> String {
    if value >= 1_000_000_000 {
        format!("{:.2} Gb", value as f64 / 1_000_000_000.0)
    } else if value >= 1_000_000 {
        format!("{:.1} Mb", value as f64 / 1_000_000.0)
    } else if value >= 1_000 {
        format!("{:.1} kb", value as f64 / 1_000.0)
    } else {
        format!("{value} bp")
    }
}
