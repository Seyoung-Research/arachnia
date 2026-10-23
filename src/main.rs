//! # Arachnia desktop application
//!
//! This binary is the `eframe`/`egui` front end for Arachnia. It embeds the
//! scoring and off-target engines from the `arachnia` library crate and adds a
//! fully offline, obsidian-and-crimson desktop workflow.
//!
//! The GUI is intentionally covered here as well as in the library crate docs
//! because the binary target shares the package name `arachnia`; `cargo doc`
//! documents the library entry point at `target/doc/arachnia/index.html`.
//!
//! # Workflow
//!
//! 1. Paste a DNA sequence or FASTA record, load a FASTA file, or press
//!    **TP53 demo** to try the built-in example.
//! 2. Choose **SpCas9** (20 nt guide, 3′ `NGG`, Rule Set 3) or **Cpf1
//!    (Cas12a)** (23 nt guide, 5′ `TTTV`, SeqDeepCpf1).
//! 3. Choose the tracrRNA for Cas9 (`Hsu2013` or `Chen2013`), whether to scan
//!    both strands, and optional GC/score/poly-T filters.
//! 4. Press **Design**. The guide table can be sorted by any column, searched
//!    by guide/PAM text, and exported or copied as TSV.
//! 5. Click a row to open the selected-guide panel. It shows the full model
//!    context, approximate cut site, score, flags, and an **Off-targets…**
//!    action.
//! 6. In the off-target window, choose the PAM pattern, mismatch limit,
//!    both-strands toggle, and search either the current target sequence or a
//!    reference genome. Genome loading and searching run in background threads
//!    with progress and cancellation.
//! 7. Optionally load a GTF/GFF3 file to annotate each hit with its chromosome
//!    region and gene. Promoters are 500 bp to 5 kb upstream of a TSS.
//! 8. Ranked hits show mismatch counts, strand, 0-based location and the
//!    guide-oriented DNA with mismatches in lower case. Copy or export the
//!    results as Cas-OFFinder-compatible TSV.
//!
//! # On-target features
//!
//! * Plain DNA, FASTA and multi-record FASTA input. The first record is used
//!   and a warning is shown when additional records are ignored.
//! * Live GC, score and poly-T filters; both-strand scanning.
//! * Sortable columns: rank, guide, PAM, strand, start/end/cut, GC, score and
//!   poly-T run.
//! * Score heat colors, strand badges and flags for GC/poly-T issues.
//! * Selected-guide details with the exact context passed to the model.
//! * Copy selected/all rows to the clipboard, or export TSV to a file.
//!
//! # Off-target features
//!
//! * Current-target or reference-genome search.
//! * Exact IUPAC PAM matching, Hamming mismatches over the guide, both strands,
//!   and Cas-OFFinder/Bowtie 0-based location conventions.
//! * Background genome loading and searching with progress and cancel.
//! * Ranked mismatch table with lower-case mismatch highlighting.
//! * Optional GTF/GFF3 region annotation (promoter, exon, UTR, intron,
//!   intergenic) with the overlapping gene.
//! * Copy/export Cas-OFFinder-compatible TSV.
//! * hg38 is detected from `ARACHNIA_HG38`, `ARACHNIA_HG38_DIR`, XDG data
//!   directories and common genome folders. Install it with
//!   `scripts/fetch_hg38.sh`; the GUI loads it only when requested, because a
//!   whole hg38 reference uses about 3.1 GiB of RAM.
//! * A GTF/GFF3 annotation is detected from `ARACHNIA_ANNOTATION`, hg38
//!   directories and common genome folders. Install a GENCODE file with
//!   `scripts/fetch_annotation.sh`.
//!
//! # Visual design
//!
//! The palette is a custom “obsidian & crimson” theme: near-black backgrounds,
//! warm off-white text, muted crimson accents, and score heat colors. See the
//! private `theme` module for the palette constants.
//!
//! # Limitations
//!
//! * DNA/RNA bulges are not supported in off-target search.
//! * Compressed FASTA is not supported.
//! * hg38 and gene annotations are not bundled; they must be installed
//!   locally.
//! * Cut positions are approximate and scores are design aids, not measured
//!   editing efficiencies.
//!
//! For the reusable engine and command-line tool, see the library crate docs
//! and the `cas-offinder` binary.

mod off_target_panel;
mod theme;

use std::fs;
use std::path::Path;

use arachnia::{design, results_tsv, DesignResult, Enzyme, Filters, GuideCandidate, Strand, Tracr};
use eframe::egui::{self, Align, Color32, FontId, Layout, Margin, RichText, Stroke};
use egui_extras::{Column, TableBuilder};
use off_target_panel::OffTargetPanel;
use theme::*;

/// TP53 coding sequence (ENST00000269305), used as the built-in demo.
const DEMO: &str = concat!(
    ">TP53 CDS (demo, ENST00000269305)\n",
    "ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCA\n",
    "GACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAAGCAATG\n",
    "GATGATTTGATGCTGTCCCCGGACGATATTGAACAATGGTTCACTGAAGACCCAGGTCCA\n",
    "GATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGCAGCTCCT\n",
    "ACACCGGCGGCCCCTGCACCAGCCCCCTCCTGGCCCCTGTCATCTTCTGTCCCTTCCCAG\n",
    "AAAACCTACCAGGGCAGCTACGGTTTCCGTCTGGGCTTCTTGCATTCTGGGACAGCCAAG\n",
    "TCTGTGACTTGCACGTACTCCCCTGCCCTCAACAAGATGTTTTGCCAACTGGCCAAGACC\n",
    "TGCCCTGTGCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATG\n",
    "GCCATCTACAAGCAGTCACAGCACATGACGGAGGTTGTGAGGCGCTGCCCCCACCATGAG\n",
    "CGCTGCTCAGATAGCGATGGTCTGGCCCCTCCTCAGCATCTTATCCGAGTGGAAGGAAAT\n",
    "TTGCGTGTGGAGTATTTGGATGACAGAAACACTTTTCGACATAGTGTGGTGGTGCCCTAT\n",
    "GAGCCGCCTGAGGTTGGCTCTGACTGTACCACCATCCACTACAACTACATGTGTAACAGT\n",
    "TCCTGCATGGGCGGCATGAACCGGAGGCCCATCCTCACCATCATCACACTGGAAGACTCC\n",
    "AGTGGTAATCTACTGGGACGGAACAGCTTTGAGGTGCGTGTTTGTGCCTGTCCTGGGAGA\n",
    "GACCGGCGCACAGAGGAAGAGAATCTCCGCAAGAAAGGGGAGCCTCACCACGAGCTGCCC\n",
    "CCAGGGAGCACTAAGCGAGCACTGCCCAACAACACCAGCTCCTCTCCCCAGCCAAAGAAG\n",
    "AAACCACTGGATGGAGAATATTTCACCCTTCAGATCCGTGGGCGTGAGCGCTTCGAGATG\n",
    "TTCCGAGAGCTGAATGAGGCCTTGGAACTCAAGGATGCCCAGGCTGGGAAGGAGCCAGGG\n",
    "GGGAGCAGGGCTCACTCCAGCCACCTGAAGTCCAAAAAGGGTCAGTCTACCTCCCGCCAT\n",
    "AAAAAACTCATGTTCAAGACAGAAGGGCCTGACTCAGACTGA\n",
);

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SortKey {
    Score,
    Guide,
    Start,
    Gc,
    Cut,
    Strand,
}

struct ArachniaApp {
    sequence_input: String,
    last_designed_input: String,
    enzyme: Enzyme,
    tracr: Tracr,
    both_strands: bool,
    gc_min: f64,
    gc_max: f64,
    min_score: f64,
    filter_poly_t: bool,
    poly_t_limit: usize,
    search: String,
    sort_key: SortKey,
    sort_desc: bool,
    result: Option<DesignResult>,
    view: Vec<usize>,
    view_dirty: bool,
    selected: Option<usize>,
    prepared_name: Option<String>,
    prepared_len: usize,
    prepared_warning: Option<String>,
    error: Option<String>,
    status: String,
    want_find: bool,
    load_open: bool,
    load_path: String,
    save_open: bool,
    save_path: String,
    dialog_error: Option<String>,
    about_open: bool,
    style_done: bool,
    frame_count: u32,
    screenshot_done: bool,
    off_target: OffTargetPanel,
}

fn main() -> eframe::Result {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1320.0, 860.0])
            .with_min_inner_size([1180.0, 700.0])
            .with_title("Arachnia — CRISPR guide RNA designer"),
        ..Default::default()
    };
    eframe::run_native(
        "Arachnia",
        options,
        Box::new(|cc| Ok(Box::new(ArachniaApp::new(cc)))),
    )
}

impl ArachniaApp {
    fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        let mut app = Self {
            sequence_input: DEMO.to_string(),
            last_designed_input: String::new(),
            enzyme: Enzyme::Cas9,
            tracr: Tracr::Hsu2013,
            both_strands: true,
            gc_min: 0.0,
            gc_max: 100.0,
            min_score: -3.0,
            filter_poly_t: false,
            poly_t_limit: 4,
            search: String::new(),
            sort_key: SortKey::Score,
            sort_desc: true,
            result: None,
            view: Vec::new(),
            view_dirty: true,
            selected: None,
            prepared_name: None,
            prepared_len: 0,
            prepared_warning: None,
            error: None,
            status: "Ready".to_string(),
            want_find: false,
            load_open: false,
            load_path: String::new(),
            save_open: false,
            save_path: "arachnia_guides.tsv".to_string(),
            dialog_error: None,
            about_open: false,
            style_done: false,
            frame_count: 0,
            screenshot_done: false,
            off_target: OffTargetPanel::default(),
        };
        // Development hooks used to regenerate the README screenshots.
        if std::env::var("ARACHNIA_ENZYME")
            .map(|value| value.eq_ignore_ascii_case("cpf1"))
            .unwrap_or(false)
        {
            app.enzyme = Enzyme::Cpf1;
            app.min_score = 0.0;
        }
        app.find();
        app.rebuild_view();
        if let Ok(row) = std::env::var("ARACHNIA_SELECT_ROW") {
            if let Ok(row) = row.parse::<usize>() {
                if let Some(index) = app.view.get(row.saturating_sub(1)) {
                    app.selected = Some(*index);
                }
            }
        }
        app
    }

    fn apply_style(&mut self, ctx: &egui::Context) {
        if self.style_done {
            return;
        }
        self.style_done = true;
        theme::apply(ctx);
    }

    fn find(&mut self) {
        self.error = None;
        self.prepared_warning = None;
        let generation = Filters {
            gc_min: 0.0,
            gc_max: 100.0,
            min_score: f64::NEG_INFINITY,
            max_poly_t: None,
            both_strands: self.both_strands,
        };
        match design(&self.sequence_input, self.enzyme, self.tracr, &generation) {
            Ok((prepared, result)) => {
                let skipped = result.skipped_ambiguous;
                self.prepared_name = prepared.name;
                self.prepared_len = prepared.bases.len();
                self.prepared_warning = prepared.warning;
                self.status = format!(
                    "Scored {} guides from {} PAM sites{}",
                    result.candidates.len(),
                    result.pam_sites,
                    if skipped > 0 {
                        format!(" ({skipped} skipped: ambiguous bases)")
                    } else {
                        String::new()
                    }
                );
                self.result = Some(result);
                self.last_designed_input = self.sequence_input.clone();
                self.selected = None;
                self.view_dirty = true;
            }
            Err(err) => {
                self.error = Some(err.to_string());
                self.result = None;
                self.view.clear();
                self.selected = None;
                self.prepared_name = None;
                self.prepared_len = 0;
                self.prepared_warning = None;
                self.status = "Ready".to_string();
            }
        }
    }

    fn rebuild_view(&mut self) {
        let Some(result) = &self.result else {
            self.view.clear();
            self.view_dirty = false;
            return;
        };
        let needle = self.search.trim().to_ascii_uppercase();
        let gc_min = self.gc_min;
        let gc_max = self.gc_max;
        let min_score = self.min_score;
        let filter_poly_t = self.filter_poly_t;
        let poly_t_limit = self.poly_t_limit;
        let mut view: Vec<usize> = result
            .candidates
            .iter()
            .enumerate()
            .filter(|(_, c)| {
                c.gc >= gc_min
                    && c.gc <= gc_max
                    && c.score >= min_score
                    && (!filter_poly_t || c.poly_t_run < poly_t_limit)
                    && (needle.is_empty() || c.guide.contains(&needle) || c.pam.contains(&needle))
            })
            .map(|(i, _)| i)
            .collect();
        let key = self.sort_key;
        let desc = self.sort_desc;
        view.sort_by(|&ia, &ib| {
            let a = &result.candidates[ia];
            let b = &result.candidates[ib];
            let ord = match key {
                SortKey::Score => a
                    .score
                    .partial_cmp(&b.score)
                    .unwrap_or(std::cmp::Ordering::Equal),
                SortKey::Guide => a.guide.cmp(&b.guide),
                SortKey::Start => a.start.cmp(&b.start),
                SortKey::Gc => a.gc.partial_cmp(&b.gc).unwrap_or(std::cmp::Ordering::Equal),
                SortKey::Cut => a.cut.cmp(&b.cut),
                SortKey::Strand => a.strand.symbol().cmp(&b.strand.symbol()),
            };
            let ord = if desc { ord.reverse() } else { ord };
            ord.then(ia.cmp(&ib))
        });
        if let Some(selected) = self.selected {
            if !view.contains(&selected) {
                self.selected = None;
            }
        }
        self.view = view;
        self.view_dirty = false;
    }

    fn tsv(&self) -> String {
        match &self.result {
            Some(result) => results_tsv(&result.candidates, &self.view),
            None => String::new(),
        }
    }

    fn score_text(&self, c: &GuideCandidate) -> String {
        match c.enzyme {
            Enzyme::Cas9 => format!("{:.3}", c.score),
            Enzyme::Cpf1 => format!("{:.1}", c.score),
        }
    }

    fn score_color(&self, c: &GuideCandidate) -> Color32 {
        let t = match c.enzyme {
            Enzyme::Cas9 => ((c.score + 1.5) / 3.0).clamp(0.0, 1.0),
            Enzyme::Cpf1 => (c.score / 100.0).clamp(0.0, 1.0),
        };
        theme::score_color(t)
    }

    fn flags(c: &GuideCandidate) -> String {
        let mut parts = Vec::new();
        if c.poly_t_run >= 4 {
            parts.push(format!("poly-T({})", c.poly_t_run));
        }
        if c.gc < 20.0 {
            parts.push("low GC".into());
        } else if c.gc > 80.0 {
            parts.push("high GC".into());
        }
        if parts.is_empty() {
            "—".into()
        } else {
            parts.join(", ")
        }
    }

    fn input_edited(&self) -> bool {
        self.result.is_some() && self.sequence_input != self.last_designed_input
    }
}

impl eframe::App for ArachniaApp {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        let ctx = ui.ctx().clone();
        self.apply_style(&ctx);
        self.frame_count = self.frame_count.saturating_add(1);
        if let Ok(path) = std::env::var("ARACHNIA_SCREENSHOT") {
            if !self.screenshot_done && self.frame_count >= 3 {
                self.screenshot_done = true;
                ctx.request_screenshot(move |image| {
                    if let Err(err) = write_ppm(Path::new(&path), &image) {
                        eprintln!("screenshot failed: {err}");
                    }
                });
                ctx.request_repaint();
            }
        }
        if ctx
            .input(|i| i.key_pressed(egui::Key::Enter) && (i.modifiers.ctrl || i.modifiers.command))
        {
            self.want_find = true;
        }

        self.header(ui);
        self.status_bar(ui);
        self.input_panel(ui);
        self.results_panel(ui);
        self.dialogs(&ctx);
        self.off_target.show(&ctx, &self.sequence_input);

        if self.want_find {
            self.want_find = false;
            self.find();
            ctx.request_repaint();
        }
    }
}

impl ArachniaApp {
    fn header(&mut self, ui: &mut egui::Ui) {
        egui::Panel::top("app_header")
            .frame(
                egui::Frame::new()
                    .fill(HEADER_BG)
                    .stroke(Stroke::new(1.0, HEADER_BORDER))
                    .inner_margin(Margin::symmetric(20, 10)),
            )
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    let (rect, _) =
                        ui.allocate_exact_size(egui::vec2(34.0, 34.0), egui::Sense::hover());
                    draw_web_icon(ui, rect);
                    ui.add_space(4.0);
                    ui.vertical(|ui| {
                        ui.label(RichText::new("ARACHNIA").size(23.0).strong().color(WHITE));
                        ui.label(
                            RichText::new("CRISPR guide RNA designer")
                                .size(11.5)
                                .color(HEADER_MUTED),
                        );
                    });
                    ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                        if ui
                            .add(
                                egui::Button::new(RichText::new("About").size(12.5).color(WHITE))
                                    .frame(false),
                            )
                            .clicked()
                        {
                            self.about_open = true;
                        }
                        ui.add_space(12.0);
                        ui.label(
                            RichText::new("SpCas9 · RuleSet3      |      Cpf1 · SeqDeepCpf1")
                                .size(11.5)
                                .color(HEADER_MUTED),
                        );
                    });
                });
            });
    }

    fn status_bar(&mut self, ui: &mut egui::Ui) {
        egui::Panel::bottom("status_bar")
            .exact_size(28.0)
            .frame(
                egui::Frame::new()
                    .fill(CARD)
                    .inner_margin(Margin::symmetric(14, 4))
                    .stroke(Stroke::new(1.0, BORDER)),
            )
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    if let Some(err) = &self.error {
                        ui.label(RichText::new(err).size(11.0).color(RED));
                    } else {
                        ui.label(RichText::new(&self.status).size(11.0).color(MUTED));
                    }
                    ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                        if ui
                            .add(
                                egui::Button::new(
                                    RichText::new("About").size(10.5).color(MUTED),
                                )
                                .frame(false),
                            )
                            .clicked()
                        {
                            self.about_open = true;
                        }
                        ui.label(
                            RichText::new(
                                "RuleSet3 · 2,279 trees   ·   SeqDeepCpf1 · 102,440 weights   ·   pure Rust",
                            )
                            .size(10.0)
                            .color(MUTED),
                        );
                    });
                });
            });
    }

    fn input_panel(&mut self, ui: &mut egui::Ui) {
        egui::Panel::left("input_panel")
            .resizable(false)
            .exact_size(362.0)
            .frame(
                egui::Frame::new()
                    .fill(CARD)
                    .inner_margin(Margin::symmetric(16, 12))
                    .stroke(Stroke::new(1.0, BORDER)),
            )
            .show(ui, |ui| {
                egui::Panel::bottom("left_find_bar")
                    .frame(egui::Frame::new().fill(CARD).inner_margin(Margin {
                        left: 0,
                        right: 0,
                        top: 6,
                        bottom: 0,
                    }))
                    .show(ui, |ui| {
                        ui.separator();
                        ui.add_space(2.0);
                        if let Some(err) = &self.error {
                            egui::Frame::new()
                                .fill(ERROR_BG)
                                .corner_radius(6)
                                .inner_margin(8)
                                .show(ui, |ui| {
                                    ui.label(RichText::new(err).size(11.0).color(RED));
                                });
                            ui.add_space(4.0);
                        }
                        let can_find = !self.sequence_input.trim().is_empty();
                        let button = egui::Button::new(
                            RichText::new("Find gRNAs").size(15.0).strong().color(WHITE),
                        )
                        .fill(ACCENT)
                        .corner_radius(7)
                        .stroke(Stroke::NONE);
                        if ui
                            .add_enabled_ui(can_find, |ui| {
                                ui.add_sized([ui.available_width(), 38.0], button)
                            })
                            .inner
                            .on_hover_text("Score every PAM site (Ctrl+Enter)")
                            .clicked()
                        {
                            self.want_find = true;
                        }
                        ui.vertical_centered(|ui| {
                            ui.label(RichText::new("Ctrl + Enter").size(9.5).color(MUTED));
                        });
                    });
                egui::ScrollArea::vertical()
                    .auto_shrink([false, false])
                    .scroll_bar_visibility(egui::scroll_area::ScrollBarVisibility::AlwaysVisible)
                    .show(ui, |ui| {
                        ui.horizontal(|ui| {
                            section_label(ui, "TARGET SEQUENCE");
                            ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                                if self.prepared_len > 0 && !self.input_edited() {
                                    ui.label(
                                        RichText::new(format!("{} bp", self.prepared_len))
                                            .size(10.5)
                                            .color(MUTED),
                                    );
                                }
                            });
                        });
                        egui::ScrollArea::vertical()
                            .id_salt("sequence_editor_scroll")
                            .max_height(160.0)
                            .auto_shrink([false, false])
                            .scroll_bar_visibility(
                                egui::scroll_area::ScrollBarVisibility::AlwaysVisible,
                            )
                            .show(ui, |ui| {
                                ui.add(
                                    egui::TextEdit::multiline(&mut self.sequence_input)
                                        .code_editor()
                                        .desired_rows(9)
                                        .desired_width(f32::INFINITY)
                                        .hint_text("Paste DNA or FASTA…"),
                                );
                            });
                        ui.horizontal(|ui| {
                            if ui.button("Load file…").clicked() {
                                self.load_open = true;
                                self.dialog_error = None;
                            }
                            if ui.button("TP53 demo").clicked() {
                                self.sequence_input = DEMO.to_string();
                                self.want_find = true;
                            }
                            if ui.button("Clear").clicked() {
                                self.sequence_input.clear();
                                self.result = None;
                                self.view.clear();
                                self.selected = None;
                                self.prepared_name = None;
                                self.prepared_len = 0;
                                self.prepared_warning = None;
                                self.last_designed_input.clear();
                                self.error = None;
                                self.status = "Ready".to_string();
                            }
                        });
                        if self.input_edited() {
                            ui.label(
                                RichText::new("Sequence edited — press Find gRNAs to rescore.")
                                    .size(10.5)
                                    .color(AMBER),
                            );
                        }
                        if let Some(name) = &self.prepared_name {
                            if !self.input_edited() {
                                ui.label(
                                    RichText::new(format!("{} · {} bp", name, self.prepared_len))
                                        .size(10.0)
                                        .color(MUTED),
                                );
                            }
                        }
                        if let Some(warning) = &self.prepared_warning {
                            ui.label(RichText::new(warning).size(10.5).color(AMBER));
                        }

                        ui.add_space(3.0);
                        ui.separator();
                        ui.add_space(3.0);
                        section_label(ui, "ENZYME");
                        let width = ui.available_width();
                        let button_w = (width - 8.0) / 2.0;
                        ui.horizontal(|ui| {
                            if enzyme_button(
                                ui,
                                self.enzyme == Enzyme::Cas9,
                                "SpCas9",
                                "RuleSet 3",
                                button_w,
                            )
                            .clicked()
                                && self.enzyme != Enzyme::Cas9
                            {
                                self.enzyme = Enzyme::Cas9;
                                self.min_score = -3.0;
                                if self.result.is_some() {
                                    self.want_find = true;
                                }
                            }
                            if enzyme_button(
                                ui,
                                self.enzyme == Enzyme::Cpf1,
                                "Cpf1 (Cas12a)",
                                "SeqDeepCpf1",
                                button_w,
                            )
                            .clicked()
                                && self.enzyme != Enzyme::Cpf1
                            {
                                self.enzyme = Enzyme::Cpf1;
                                self.min_score = 0.0;
                                if self.result.is_some() {
                                    self.want_find = true;
                                }
                            }
                        });

                        let mut tracr_changed = false;
                        if self.enzyme == Enzyme::Cas9 {
                            ui.horizontal(|ui| {
                                ui.label(RichText::new("tracrRNA").size(11.5).color(MUTED));
                                egui::ComboBox::from_id_salt("tracr_combo")
                                    .selected_text(self.tracr.label())
                                    .width(150.0)
                                    .show_ui(ui, |ui| {
                                        if ui
                                            .selectable_value(
                                                &mut self.tracr,
                                                Tracr::Hsu2013,
                                                "Hsu2013 (default)",
                                            )
                                            .changed()
                                        {
                                            tracr_changed = true;
                                        }
                                        if ui
                                            .selectable_value(
                                                &mut self.tracr,
                                                Tracr::Chen2013,
                                                "Chen2013",
                                            )
                                            .changed()
                                        {
                                            tracr_changed = true;
                                        }
                                    });
                            });
                        } else {
                            ui.label(
                                RichText::new("TTTV PAM · 23 nt guide · 34 nt context")
                                    .size(10.5)
                                    .color(MUTED),
                            );
                        }
                        if tracr_changed && self.result.is_some() {
                            self.want_find = true;
                        }

                        ui.add_space(3.0);
                        ui.separator();
                        ui.add_space(3.0);
                        section_label(ui, "FILTERS");
                        if ui
                            .checkbox(&mut self.both_strands, "Search both strands")
                            .changed()
                            && self.result.is_some()
                        {
                            self.want_find = true;
                        }
                        ui.horizontal(|ui| {
                            ui.label(RichText::new("GC").size(11.5).color(MUTED));
                            let a = ui.add(
                                egui::DragValue::new(&mut self.gc_min)
                                    .range(0.0..=100.0)
                                    .speed(1.0)
                                    .suffix("%"),
                            );
                            ui.label("–");
                            let b = ui.add(
                                egui::DragValue::new(&mut self.gc_max)
                                    .range(0.0..=100.0)
                                    .speed(1.0)
                                    .suffix("%"),
                            );
                            if a.changed() || b.changed() {
                                self.view_dirty = true;
                            }
                        });
                        if self.gc_min > self.gc_max {
                            std::mem::swap(&mut self.gc_min, &mut self.gc_max);
                            self.view_dirty = true;
                        }
                        let (lo, hi, step) = match self.enzyme {
                            Enzyme::Cas9 => (-3.0, 2.5, 0.05),
                            Enzyme::Cpf1 => (0.0, 100.0, 1.0),
                        };
                        if ui
                            .add(
                                egui::Slider::new(&mut self.min_score, lo..=hi)
                                    .text("Minimum score")
                                    .step_by(step),
                            )
                            .changed()
                        {
                            self.view_dirty = true;
                        }
                        ui.horizontal(|ui| {
                            let cb = ui.checkbox(&mut self.filter_poly_t, "Filter poly-T ≥");
                            let dv = ui.add_enabled(
                                self.filter_poly_t,
                                egui::DragValue::new(&mut self.poly_t_limit).range(2..=12),
                            );
                            if cb.changed() || dv.changed() {
                                self.view_dirty = true;
                            }
                        });
                    });
            });
    }

    fn results_panel(&mut self, ui: &mut egui::Ui) {
        egui::CentralPanel::default()
            .frame(
                egui::Frame::new()
                    .fill(BG)
                    .inner_margin(Margin::symmetric(16, 12)),
            )
            .show(ui, |ui| {
                if self.result.is_none() {
                    self.empty_state(ui);
                    return;
                }
                if self.view_dirty {
                    self.rebuild_view();
                }
                self.details_panel(ui);

                let total = self
                    .result
                    .as_ref()
                    .map(|r| r.candidates.len())
                    .unwrap_or(0);
                let pam_sites = self
                    .result
                    .as_ref()
                    .map(|r| r.pam_sites)
                    .unwrap_or(0);
                let skipped = self
                    .result
                    .as_ref()
                    .map(|r| r.skipped_ambiguous)
                    .unwrap_or(0);
                let filtered = self.view.len();

                ui.horizontal(|ui| {
                    ui.vertical(|ui| {
                        ui.label(
                            RichText::new(format!("{filtered} guides"))
                                .size(21.0)
                                .strong()
                                .color(TEXT),
                        );
                        ui.label(
                            RichText::new(format!(
                                "{filtered} of {total} scored · {pam_sites} PAM sites · {skipped} skipped"
                            ))
                            .size(10.5)
                            .color(MUTED),
                        );
                    });
                    ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                        if ui.button("Export TSV…").clicked() {
                            self.save_open = true;
                            self.dialog_error = None;
                        }
                        if ui.button("Copy TSV").clicked() {
                            if self.view.is_empty() {
                                self.status = "No guides to copy".to_string();
                            } else {
                                let text = self.tsv();
                                ui.ctx().copy_text(text);
                                self.status = format!(
                                    "Copied {} rows to the clipboard",
                                    self.view.len()
                                );
                            }
                        }
                    });
                });
                ui.add_space(5.0);
                ui.horizontal(|ui| {
                    let resp = ui.add(
                        egui::TextEdit::singleline(&mut self.search)
                            .hint_text("Filter by guide or PAM sequence…")
                            .desired_width(235.0),
                    );
                    if resp.changed() {
                        self.view_dirty = true;
                    }
                    if !self.search.is_empty() && ui.small_button("Clear").clicked() {
                        self.search.clear();
                        self.view_dirty = true;
                    }
                    let protocol = match self.enzyme {
                        Enzyme::Cas9 => format!(
                            "{} nt guide · {} PAM · {}",
                            self.enzyme.guide_len(),
                            self.enzyme.pam_label(),
                            self.tracr.label()
                        ),
                        Enzyme::Cpf1 => format!(
                            "{} nt guide · {} PAM",
                            self.enzyme.guide_len(),
                            self.enzyme.pam_label()
                        ),
                    };
                    ui.label(
                        RichText::new(protocol)
                            .size(10.5)
                            .color(MUTED),
                    );
                });
                ui.add_space(6.0);

                if self.view.is_empty() {
                    ui.vertical_centered(|ui| {
                        ui.add_space(50.0);
                        ui.label(
                            RichText::new("No guides match the current filters")
                                .size(16.0)
                                .strong()
                                .color(TEXT),
                        );
                        ui.add_space(4.0);
                        ui.label(
                            RichText::new(
                                "Relax the GC range, lower the score cutoff, or clear the search.",
                            )
                            .size(11.5)
                            .color(MUTED),
                        );
                        ui.add_space(10.0);
                        if ui.button("Reset filters").clicked() {
                            self.gc_min = 0.0;
                            self.gc_max = 100.0;
                            self.min_score = match self.enzyme {
                                Enzyme::Cas9 => -3.0,
                                Enzyme::Cpf1 => 0.0,
                            };
                            self.filter_poly_t = false;
                            self.search.clear();
                            self.view_dirty = true;
                        }
                    });
                    return;
                }

                let current_key = self.sort_key;
                let current_desc = self.sort_desc;
                let selected = self.selected;
                let mut clicked_row = None;
                let mut clicked_sort: Option<SortKey> = None;
                let enzyme = self.enzyme;
                let candidates = &self.result.as_ref().unwrap().candidates;
                let view = &self.view;
                let compact = ui.available_width() < 860.0;
                let context_header = if compact {
                    "Context"
                } else {
                    "Context (guide orientation)"
                };

                egui::Frame::new()
                    .fill(CARD)
                    .corner_radius(8)
                    .stroke(Stroke::new(1.0, BORDER))
                    .inner_margin(Margin::same(0))
                    .show(ui, |ui| {
                        TableBuilder::new(ui)
                            .striped(true)
                            .resizable(false)
                            .sense(egui::Sense::click())
                            .auto_shrink([false, false])
                            .cell_layout(Layout::left_to_right(Align::Center))
                            .column(Column::exact(if compact { 34.0 } else { 40.0 }))
                            .column(Column::exact(if compact { 182.0 } else { 188.0 }))
                            .column(Column::exact(if compact { 44.0 } else { 50.0 }))
                            .column(Column::exact(if compact { 44.0 } else { 42.0 }))
                            .column(Column::exact(if compact { 82.0 } else { 96.0 }))
                            .column(Column::exact(if compact { 40.0 } else { 46.0 }))
                            .column(Column::exact(if compact { 50.0 } else { 56.0 }))
                            .column(Column::exact(if compact { 72.0 } else { 80.0 }))
                            .column(Column::exact(if compact { 84.0 } else { 88.0 }))
                            .column(
                                Column::remainder()
                                    .at_least(if compact { 70.0 } else { 80.0 })
                                    .clip(true),
                            )
                            .header(28.0, |mut header| {
                                header.col(|ui| {
                                    ui.label(
                                        RichText::new("#")
                                            .size(10.5)
                                            .strong()
                                            .color(MUTED),
                                    );
                                });
                                header.col(|ui| {
                                    if sort_header(
                                        ui,
                                        "Guide",
                                        SortKey::Guide,
                                        current_key,
                                        current_desc,
                                    ) {
                                        clicked_sort = Some(SortKey::Guide);
                                    }
                                });
                                header.col(|ui| {
                                    ui.label(
                                        RichText::new("PAM")
                                            .size(10.5)
                                            .strong()
                                            .color(MUTED),
                                    );
                                });
                                header.col(|ui| {
                                    if sort_header(
                                        ui,
                                        "Strand",
                                        SortKey::Strand,
                                        current_key,
                                        current_desc,
                                    ) {
                                        clicked_sort = Some(SortKey::Strand);
                                    }
                                });
                                header.col(|ui| {
                                    if sort_header(
                                        ui,
                                        "Position",
                                        SortKey::Start,
                                        current_key,
                                        current_desc,
                                    ) {
                                        clicked_sort = Some(SortKey::Start);
                                    }
                                });
                                header.col(|ui| {
                                    if sort_header(
                                        ui,
                                        "Cut",
                                        SortKey::Cut,
                                        current_key,
                                        current_desc,
                                    ) {
                                        clicked_sort = Some(SortKey::Cut);
                                    }
                                });
                                header.col(|ui| {
                                    if sort_header(
                                        ui,
                                        "GC",
                                        SortKey::Gc,
                                        current_key,
                                        current_desc,
                                    ) {
                                        clicked_sort = Some(SortKey::Gc);
                                    }
                                });
                                header.col(|ui| {
                                    if sort_header(
                                        ui,
                                        score_header(enzyme),
                                        SortKey::Score,
                                        current_key,
                                        current_desc,
                                    ) {
                                        clicked_sort = Some(SortKey::Score);
                                    }
                                });
                                header.col(|ui| {
                                    ui.label(
                                        RichText::new("Flags")
                                            .size(10.5)
                                            .strong()
                                            .color(MUTED),
                                    );
                                });
                                header.col(|ui| {
                                    ui.label(
                                        RichText::new(context_header)
                                            .size(10.5)
                                            .strong()
                                            .color(MUTED),
                                    );
                                });
                            })
                            .body(|body| {
                                body.rows(23.0, view.len(), |mut row| {
                                    let row_index = row.index();
                                    let candidate_index = view[row_index];
                                    let c = &candidates[candidate_index];
                                    row.set_selected(selected == Some(candidate_index));
                                    row.col(|ui| {
                                        ui.label(
                                            RichText::new((row_index + 1).to_string())
                                                .size(10.5)
                                                .color(MUTED),
                                        );
                                    });
                                    row.col(|ui| {
                                        ui.label(
                                            RichText::new(&c.guide)
                                                .monospace()
                                                .size(12.5),
                                        );
                                    });
                                    row.col(|ui| {
                                        ui.label(
                                            RichText::new(&c.pam)
                                                .monospace()
                                                .size(12.0)
                                                .color(ACCENT_TEXT),
                                        );
                                    });
                                    row.col(|ui| {
                                        ui.label(
                                            RichText::new(c.strand.symbol().to_string())
                                                .monospace()
                                                .size(12.0)
                                                .color(match c.strand {
                                                    Strand::Plus => GREEN,
                                                    Strand::Minus => AMBER,
                                                }),
                                        );
                                    });
                                    row.col(|ui| {
                                        ui.label(
                                            RichText::new(format!(
                                                "{}–{}",
                                                c.start, c.end
                                            ))
                                            .monospace()
                                            .size(11.0)
                                            .color(MUTED),
                                        );
                                    });
                                    row.col(|ui| {
                                        ui.label(
                                            RichText::new(c.cut.to_string())
                                                .monospace()
                                                .size(11.0)
                                                .color(MUTED),
                                        );
                                    });
                                    row.col(|ui| {
                                        ui.label(
                                            RichText::new(format!("{:.1}%", c.gc))
                                                .size(11.0),
                                        );
                                    });
                                    row.col(|ui| {
                                        let score = score_text_for(enzyme, c);
                                        ui.label(
                                            RichText::new(score)
                                                .monospace()
                                                .size(12.0)
                                                .strong()
                                                .color(score_color_for(
                                                    enzyme, c.score,
                                                )),
                                        );
                                    });
                                    row.col(|ui| {
                                        let flags = Self::flags(c);
                                        let color = if flags == "—" {
                                            MUTED
                                        } else {
                                            AMBER
                                        };
                                        ui.label(
                                            RichText::new(flags)
                                                .size(10.0)
                                                .color(color),
                                        );
                                    });
                                    row.col(|ui| {
                                        ui.label(
                                            RichText::new(&c.context)
                                                .monospace()
                                                .size(9.5)
                                                .color(MUTED),
                                        );
                                    });
                                    let response = row
                                        .response()
                                        .on_hover_cursor(egui::CursorIcon::PointingHand);
                                    if response.clicked() {
                                        clicked_row = Some(candidate_index);
                                    }
                                });
                            });
                    });

                if let Some(key) = clicked_sort {
                    if self.sort_key == key {
                        self.sort_desc = !self.sort_desc;
                    } else {
                        self.sort_key = key;
                        self.sort_desc = default_sort_desc(key);
                    }
                    self.view_dirty = true;
                }
                if let Some(index) = clicked_row {
                    self.selected = Some(index);
                }
            });
    }

    fn details_panel(&mut self, ui: &mut egui::Ui) {
        if self.selected.is_none() {
            return;
        }
        let Some(candidate) = self
            .result
            .as_ref()
            .and_then(|r| self.selected.and_then(|i| r.candidates.get(i)))
            .cloned()
        else {
            self.selected = None;
            return;
        };
        let mut close = false;
        let mut copy_guide: Option<String> = None;
        let mut copy_context: Option<String> = None;
        let mut open_offtarget = false;

        egui::Panel::bottom("selected_guide")
            .resizable(true)
            .default_size(190.0)
            .min_size(120.0)
            .frame(
                egui::Frame::new()
                    .fill(CARD)
                    .inner_margin(Margin::symmetric(16, 8))
                    .stroke(Stroke::new(1.0, BORDER)),
            )
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    section_label(ui, "SELECTED GUIDE");
                    ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                        if ui
                            .add(
                                egui::Button::new(RichText::new("✕").size(13.0).color(MUTED))
                                    .frame(false),
                            )
                            .clicked()
                        {
                            close = true;
                        }
                    });
                });
                let score_text = self.score_text(&candidate);
                let score_color = self.score_color(&candidate);
                let wide = ui.available_width() >= 860.0;
                ui.horizontal_top(|ui| {
                    guide_details_header(
                        ui,
                        &candidate,
                        &mut copy_guide,
                        &mut copy_context,
                        &mut open_offtarget,
                    );
                    ui.add_space(10.0);
                    guide_details_meta(ui, &candidate, &score_text, score_color);
                    if wide {
                        ui.add_space(10.0);
                        guide_details_context(ui, &candidate, 260.0);
                    }
                });
                if !wide {
                    ui.separator();
                    guide_details_context(ui, &candidate, 0.0);
                }
            });

        if close {
            self.selected = None;
        }
        if let Some(text) = copy_guide {
            ui.ctx().copy_text(text);
            self.status = "Copied guide + PAM to the clipboard".to_string();
        }
        if let Some(text) = copy_context {
            ui.ctx().copy_text(text);
            self.status = "Copied model context to the clipboard".to_string();
        }
        if open_offtarget {
            self.off_target.open(&candidate);
        }
    }

    fn empty_state(&mut self, ui: &mut egui::Ui) {
        ui.vertical_centered(|ui| {
            ui.add_space(80.0);
            ui.label(
                RichText::new("Design CRISPR guides")
                    .size(30.0)
                    .strong()
                    .color(TEXT),
            );
            ui.add_space(6.0);
            ui.label(
                RichText::new(
                    "Paste a target sequence or load a FASTA file, choose an enzyme, and press Find gRNAs.",
                )
                .size(13.0)
                .color(MUTED),
            );
            ui.add_space(18.0);
            if ui
                .add_sized(
                    [210.0, 38.0],
                    egui::Button::new(
                        RichText::new("Load TP53 demo").color(WHITE),
                    )
                    .fill(ACCENT)
                    .corner_radius(7),
                )
                .clicked()
            {
                self.sequence_input = DEMO.to_string();
                self.want_find = true;
            }
            ui.add_space(28.0);
            ui.label(
                RichText::new(
                    "On-target scoring · Cas-OFFinder-style off-target analysis",
                )
                .size(11.0)
                .color(MUTED),
            );
        });
    }

    fn dialogs(&mut self, ctx: &egui::Context) {
        if self.about_open {
            let mut open = self.about_open;
            egui::Window::new("About Arachnia")
                .open(&mut open)
                .resizable(true)
                .collapsible(false)
                .default_width(580.0)
                .show(ctx, |ui| {
                    ui.label(
                        RichText::new("Arachnia")
                            .size(22.0)
                            .strong()
                            .color(TEXT),
                    );
                    ui.label(
                        RichText::new(
                            "A CRISPick-style guide designer with Cas-OFFinder-style off-target analysis, written entirely in Rust.",
                        )
                        .size(12.5)
                        .color(MUTED),
                    );
                    ui.separator();
                    ui.label(
                        RichText::new(
                            "Cas9 guides are scored with Rule Set 3, the LightGBM model distributed with gpp-rnd/rs3, using a faithful Rust port of the sglearn feature extraction. Cpf1 guides are scored with SeqDeepCpf1, a 1-D CNN from Kim et al. (2018), evaluated by a pure-Rust forward pass.",
                        )
                        .size(12.0),
                    );
                    ui.add_space(6.0);
                    ui.label(
                        RichText::new(
                            "Select a guide and press Off-targets… to scan the current sequence or a loaded FASTA genome for every exact-PAM site within a chosen Hamming distance. A local hg38 reference is detected and preselected when available (install it with scripts/fetch_hg38.sh). Load a GTF/GFF3 file in the same window to add Region and Gene labels; install a GENCODE annotation with scripts/fetch_annotation.sh. The search follows Cas-OFFinder's mismatch model and reports 0-based locations on both strands.",
                        )
                        .size(12.0),
                    );
                    ui.add_space(8.0);
                    ui.label(RichText::new("Limitations").strong().size(12.5));
                    ui.label(
                        RichText::new(
                            "• Off-target search is mismatch-only: the PAM must match exactly and DNA/RNA bulges are not considered.\n• Rule Set 3 is used in sequence mode (30 nt context).\n• Cpf1 scores are predictions from a model trained on the original paired-library data.\n• Cut positions are approximate.",
                        )
                        .size(11.5)
                        .color(MUTED),
                    );
                    ui.add_space(8.0);
                    ui.label(RichText::new("Sources & licenses").strong().size(12.5));
                    ui.hyperlink_to(
                        "gpp-rnd/rs3 — Rule Set 3 (Apache-2.0)",
                        "https://github.com/gpp-rnd/rs3",
                    );
                    ui.hyperlink_to(
                        "alisakumbara/SeqDeepCpf1-Py3 — SeqDeepCpf1 (MIT)",
                        "https://github.com/alisakumbara/SeqDeepCpf1-Py3",
                    );
                    ui.hyperlink_to(
                        "Lattice-Automation/seqfold — RNA folding (MIT)",
                        "https://github.com/Lattice-Automation/seqfold",
                    );
                    ui.hyperlink_to(
                        "snugel/cas-offinder — Cas-OFFinder (BSD-3-Clause)",
                        "https://github.com/snugel/cas-offinder",
                    );
                    ui.add_space(6.0);
                    ui.label(
                        RichText::new(
                            "Kim, H. et al. Deep learning improves prediction of CRISPR–Cpf1 guide RNA activity. Nat. Biotechnol. 36, 239–241 (2018). doi:10.1038/nbt.4061",
                        )
                        .size(10.5)
                        .color(MUTED),
                    );
                    ui.label(
                        RichText::new(
                            "Bae, S., Park, J. & Kim, J.-S. Cas-OFFinder: a fast and versatile algorithm that searches for potential off-target sites of Cas9 RNA-guided endonucleases. Bioinformatics 30, 1473–1475 (2014). doi:10.1093/bioinformatics/btu048",
                        )
                        .size(10.5)
                        .color(MUTED),
                    );
                });
            self.about_open = open;
        }

        let mut load_requested = false;
        if self.load_open {
            let mut open = self.load_open;
            egui::Window::new("Load sequence")
                .open(&mut open)
                .resizable(false)
                .collapsible(false)
                .default_width(520.0)
                .show(ctx, |ui| {
                    ui.label(
                        RichText::new("Path to a FASTA or plain-text sequence file:")
                            .size(12.0)
                            .color(MUTED),
                    );
                    let resp = ui.add(
                        egui::TextEdit::singleline(&mut self.load_path)
                            .desired_width(f32::INFINITY)
                            .hint_text("/path/to/sequence.fa"),
                    );
                    if resp.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)) {
                        load_requested = true;
                    }
                    if let Some(err) = &self.dialog_error {
                        ui.label(RichText::new(err).size(11.0).color(RED));
                    }
                    ui.horizontal(|ui| {
                        if ui.button("Load").clicked() {
                            load_requested = true;
                        }
                        if ui.button("Cancel").clicked() {
                            self.load_open = false;
                            self.dialog_error = None;
                        }
                    });
                });
            self.load_open = open;
            if load_requested {
                self.load_sequence_file();
            }
        }

        let mut save_requested = false;
        if self.save_open {
            let mut open = self.save_open;
            egui::Window::new("Export guides as TSV")
                .open(&mut open)
                .resizable(false)
                .collapsible(false)
                .default_width(520.0)
                .show(ctx, |ui| {
                    ui.label(
                        RichText::new(format!("Exporting {} visible guides:", self.view.len()))
                            .size(12.0)
                            .color(MUTED),
                    );
                    let resp = ui.add(
                        egui::TextEdit::singleline(&mut self.save_path)
                            .desired_width(f32::INFINITY)
                            .hint_text("arachnia_guides.tsv"),
                    );
                    if resp.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)) {
                        save_requested = true;
                    }
                    if let Some(err) = &self.dialog_error {
                        ui.label(RichText::new(err).size(11.0).color(RED));
                    }
                    ui.horizontal(|ui| {
                        if ui.button("Save TSV").clicked() {
                            save_requested = true;
                        }
                        if ui.button("Cancel").clicked() {
                            self.save_open = false;
                            self.dialog_error = None;
                        }
                    });
                });
            self.save_open = open;
            if save_requested {
                self.save_results();
            }
        }
    }

    fn load_sequence_file(&mut self) {
        let path = self.load_path.trim().to_string();
        if path.is_empty() {
            self.dialog_error = Some("Enter a file path first.".to_string());
            return;
        }
        match fs::read_to_string(&path) {
            Ok(text) => {
                self.sequence_input = text;
                self.load_open = false;
                self.dialog_error = None;
                self.status = format!("Loaded {path}");
                self.want_find = true;
            }
            Err(err) => {
                self.dialog_error = Some(format!("Could not read {path}: {err}"));
            }
        }
    }

    fn save_results(&mut self) {
        if self.view_dirty {
            self.rebuild_view();
        }
        let path = self.save_path.trim().to_string();
        if path.is_empty() {
            self.dialog_error = Some("Enter a file path first.".to_string());
            return;
        }
        if self.view.is_empty() {
            self.dialog_error = Some("There are no guides to export.".to_string());
            return;
        }
        let tsv = self.tsv();
        match fs::write(&path, tsv) {
            Ok(()) => {
                self.save_open = false;
                self.dialog_error = None;
                self.status = format!("Saved {} guides to {path}", self.view.len());
            }
            Err(err) => {
                self.dialog_error = Some(format!("Could not write {path}: {err}"));
            }
        }
    }
}

fn section_label(ui: &mut egui::Ui, text: &str) {
    ui.label(RichText::new(text).size(10.0).strong().color(MUTED));
}

fn enzyme_button(
    ui: &mut egui::Ui,
    selected: bool,
    title: &str,
    subtitle: &str,
    width: f32,
) -> egui::Response {
    let fill = if selected { ACCENT } else { WEAK_BG };
    let color = if selected { WHITE } else { TEXT };
    let stroke = if selected {
        Stroke::new(1.0, ACCENT)
    } else {
        Stroke::new(1.0, BORDER)
    };
    ui.add_sized(
        [width, 48.0],
        egui::Button::new(
            RichText::new(format!("{title}\n{subtitle}"))
                .size(11.5)
                .color(color),
        )
        .fill(fill)
        .stroke(stroke)
        .corner_radius(7),
    )
}

fn sort_header(ui: &mut egui::Ui, label: &str, key: SortKey, current: SortKey, desc: bool) -> bool {
    let arrow = if current == key {
        if desc {
            " ▼"
        } else {
            " ▲"
        }
    } else {
        ""
    };
    ui.add(
        egui::Button::new(
            RichText::new(format!("{label}{arrow}"))
                .size(10.5)
                .strong()
                .color(MUTED),
        )
        .frame(false),
    )
    .clicked()
}

fn default_sort_desc(key: SortKey) -> bool {
    matches!(key, SortKey::Score | SortKey::Gc)
}

fn score_header(enzyme: Enzyme) -> &'static str {
    match enzyme {
        Enzyme::Cas9 => "RS3 score",
        Enzyme::Cpf1 => "Cpf1 score",
    }
}

fn score_text_for(enzyme: Enzyme, c: &GuideCandidate) -> String {
    match enzyme {
        Enzyme::Cas9 => format!("{:.3}", c.score),
        Enzyme::Cpf1 => format!("{:.1}", c.score),
    }
}

fn score_color_for(enzyme: Enzyme, score: f64) -> Color32 {
    let t = match enzyme {
        Enzyme::Cas9 => ((score + 1.5) / 3.0).clamp(0.0, 1.0),
        Enzyme::Cpf1 => (score / 100.0).clamp(0.0, 1.0),
    };
    theme::score_color(t)
}

fn meta_cell(ui: &mut egui::Ui, label: &str, value: &str, color: Color32) {
    ui.vertical(|ui| {
        ui.set_min_width(48.0);
        ui.label(RichText::new(label).size(9.0).strong().color(MUTED));
        ui.label(RichText::new(value).size(11.5).color(color));
    });
}

fn guide_details_header(
    ui: &mut egui::Ui,
    candidate: &GuideCandidate,
    copy_guide: &mut Option<String>,
    copy_context: &mut Option<String>,
    open_offtarget: &mut bool,
) {
    ui.vertical(|ui| {
        ui.set_min_width(250.0);
        ui.horizontal(|ui| {
            ui.label(
                RichText::new(&candidate.guide)
                    .monospace()
                    .size(17.0)
                    .strong()
                    .color(TEXT),
            );
            ui.label(
                RichText::new(&candidate.pam)
                    .monospace()
                    .size(17.0)
                    .strong()
                    .color(ACCENT_TEXT),
            );
        });
        ui.label(
            RichText::new(format!(
                "{} · {} strand · {} nt guide",
                candidate.enzyme.label(),
                candidate.strand.symbol(),
                candidate.guide.len()
            ))
            .size(10.5)
            .color(MUTED),
        );
        ui.horizontal(|ui| {
            if ui.small_button("Copy guide").clicked() {
                *copy_guide = Some(format!("{}{}", candidate.guide, candidate.pam));
            }
            if ui.small_button("Copy context").clicked() {
                *copy_context = Some(candidate.context.clone());
            }
            if ui
                .add(
                    egui::Button::new(RichText::new("Off-targets…").size(11.0).color(WHITE))
                        .fill(ACCENT)
                        .corner_radius(5),
                )
                .on_hover_text("Cas-OFFinder-style mismatch search")
                .clicked()
            {
                *open_offtarget = true;
            }
        });
    });
}

fn guide_details_meta(
    ui: &mut egui::Ui,
    candidate: &GuideCandidate,
    score_text: &str,
    score_color: Color32,
) {
    ui.vertical(|ui| {
        egui::Grid::new("guide_props")
            .num_columns(4)
            .spacing([12.0, 5.0])
            .show(ui, |ui| {
                meta_cell(ui, "Score", score_text, score_color);
                meta_cell(ui, "GC", &format!("{:.1}%", candidate.gc), TEXT);
                meta_cell(
                    ui,
                    "Position",
                    &format!("{}–{}", candidate.start, candidate.end),
                    TEXT,
                );
                meta_cell(ui, "Cut", &candidate.cut.to_string(), TEXT);
                ui.end_row();
                meta_cell(
                    ui,
                    "Strand",
                    &candidate.strand.symbol().to_string(),
                    match candidate.strand {
                        Strand::Plus => GREEN,
                        Strand::Minus => AMBER,
                    },
                );
                meta_cell(ui, "PAM", &candidate.pam, ACCENT);
                meta_cell(
                    ui,
                    "Poly-T",
                    &candidate.poly_t_run.to_string(),
                    if candidate.poly_t_run >= 4 {
                        AMBER
                    } else {
                        TEXT
                    },
                );
                meta_cell(ui, "Enzyme", candidate.enzyme.label(), TEXT);
                ui.end_row();
            });
    });
}

fn guide_details_context(ui: &mut egui::Ui, candidate: &GuideCandidate, min_width: f32) {
    ui.vertical(|ui| {
        ui.set_min_width(min_width);
        ui.label(
            RichText::new("CONTEXT (guide orientation)")
                .size(9.5)
                .strong()
                .color(MUTED),
        );
        ui.add(egui::Label::new(context_job(ui, candidate)).wrap());
        if candidate.poly_t_run >= 4 {
            ui.label(
                RichText::new(format!(
                    "poly-T run of {} may terminate U6 transcription",
                    candidate.poly_t_run
                ))
                .size(10.0)
                .color(AMBER),
            );
        } else if candidate.gc < 20.0 || candidate.gc > 80.0 {
            ui.label(
                RichText::new("GC content is outside 20–80%")
                    .size(10.0)
                    .color(AMBER),
            );
        } else {
            ui.label(
                RichText::new("No basic sequence warnings.")
                    .size(10.0)
                    .color(GREEN),
            );
        }
        ui.label(
            RichText::new(
                "Cut site is approximate (Cas9: 3 bp upstream of PAM; Cpf1: staggered 18/23).",
            )
            .size(9.5)
            .color(MUTED),
        );
    });
}

fn context_job(ui: &egui::Ui, c: &GuideCandidate) -> egui::text::LayoutJob {
    use egui::text::{LayoutJob, TextFormat};
    let mut job = LayoutJob::default();
    let font = FontId::monospace(13.0);
    let flank = ui.visuals().weak_text_color();
    let guide = ui.visuals().text_color();
    let pam = ACCENT_TEXT;
    let (a, b, d) = match c.enzyme {
        Enzyme::Cas9 => (4, 24, 27),
        Enzyme::Cpf1 => (8, 31, 4),
    };
    let fmt = |color| TextFormat {
        font_id: font.clone(),
        color,
        ..Default::default()
    };
    let ctx = &c.context;
    if b == 24 {
        job.append(&ctx[..a], 0.0, fmt(flank));
        job.append(&ctx[a..b], 0.0, fmt(guide));
        job.append(&ctx[b..d], 0.0, fmt(pam));
        job.append(&ctx[d..], 0.0, fmt(flank));
    } else {
        job.append(&ctx[..4], 0.0, fmt(flank));
        job.append(&ctx[4..8], 0.0, fmt(pam));
        job.append(&ctx[8..31], 0.0, fmt(guide));
        job.append(&ctx[31..], 0.0, fmt(flank));
    }
    job
}

fn draw_web_icon(ui: &mut egui::Ui, rect: egui::Rect) {
    let painter = ui.painter();
    let center = rect.center();
    let radius = rect.width() * 0.45;
    let stroke = Stroke::new(1.2, ACCENT_TEXT);
    painter.circle_stroke(center, radius, stroke);
    painter.circle_stroke(center, radius * 0.55, stroke);
    for i in 0..6 {
        let angle = i as f32 * std::f32::consts::TAU / 6.0;
        let p = center + egui::vec2(angle.cos(), angle.sin()) * radius;
        painter.line_segment([center, p], stroke);
    }
    painter.circle_filled(center, 2.5, ACCENT_TEXT);
}

fn write_ppm(path: &Path, image: &egui::ColorImage) -> std::io::Result<()> {
    let mut bytes = format!("P6\n{} {}\n255\n", image.width(), image.height()).into_bytes();
    bytes.reserve(image.pixels.len() * 3);
    for px in &image.pixels {
        bytes.extend_from_slice(&[px.r(), px.g(), px.b()]);
    }
    std::fs::write(path, bytes)
}
