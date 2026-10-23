//! In-app file explorer for sequence, genome, annotation and TSV paths.
//!
//! Arachnia is intentionally self-contained and does not depend on a native
//! file-dialog toolkit. This module provides a small `egui` explorer that is
//! shared by every path field in the application:
//!
//! * **Open sequence** — the target-sequence dialog ([`FileKind::Sequence`]).
//! * **Open reference genome** — the off-target window
//!   ([`FileKind::Genome`]); a FASTA file or a folder of chromosome FASTA
//!   files can be selected.
//! * **Open gene annotation** — the off-target window
//!   ([`FileKind::Annotation`]); GTF and GFF3 files are listed.
//! * **Save TSV** — the on-target and off-target export dialogs
//!   ([`FileKind::Tsv`]).
//!
//! The explorer is deliberately conservative:
//!
//! * It only reads directory metadata. It never opens or modifies a file
//!   until the caller performs the load or save.
//! * It never creates directories. A save path must already have an existing
//!   parent folder.
//! * It does not follow compressed files specially; `.gz` inputs are rejected
//!   by the loaders themselves because Arachnia does not support compressed
//!   FASTA or GTF/GFF3.
//! * Directory listings are capped at [`MAX_ENTRIES`] entries so that opening
//!   an enormous folder cannot stall the UI for an unbounded amount of time.
//!
//! # Navigation
//!
//! The path bar accepts an absolute path, a relative path, or a leading `~`.
//! Press `Enter` or click **Go** to navigate. Files can be selected with a
//! single click and opened with a double click. Folders are entered with a
//! single click. **Up**, **Home** and **Data** provide quick navigation; the
//! **Data** button points at the XDG Arachnia folder used by
//! `scripts/fetch_hg38.sh` and `scripts/fetch_annotation.sh`.
//!
//! When a local hg38 reference or GENCODE annotation is detected, a quick-pick
//! button selects it immediately. Hidden entries and non-matching file types
//! are hidden by default; both can be revealed with checkboxes. The selected
//! file's full path and size are shown above the action row.
//!
//! # Save mode
//!
//! Save mode shows the current folder and a file-name field. A relative file
//! name is resolved against the current folder; an absolute path is used as
//! entered. For [`FileKind::Tsv`], a missing extension is completed to
//! `.tsv`. If the destination already exists, the action button changes to
//! **Overwrite** and requires a second click, which prevents accidental data
//! loss without a separate modal dialog.
//!
//! # Integration
//!
//! Call [`FileExplorer::open_for_open`] or [`FileExplorer::open_for_save`] to
//! start a session, then call [`FileExplorer::show`] once per frame. It returns
//! `Some(path)` when the user confirms a selection, and `None` otherwise. The
//! caller is responsible for applying the returned path to its own text field;
//! this keeps the explorer independent of the main application state.

use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use arachnia::annotation::find_hg38_annotation;
use arachnia::offtarget::{expand_tilde, find_hg38, is_hg38_alias};
use eframe::egui::{self, Align, Layout, Margin, RichText, Stroke};

use crate::theme::*;

/// Maximum number of directory entries read for one listing.
///
/// The cap keeps the explorer responsive on directories containing hundreds
/// of thousands of files. When it is reached, a notice is shown and the
/// remaining entries are ignored until the user narrows the folder.
pub(crate) const MAX_ENTRIES: usize = 4_000;

/// What the explorer is being used for.
///
/// The kind controls the window title, the accepted file extensions and
/// whether folders can be selected as the final result. It is also the value
/// the caller uses to route the returned path.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum FileKind {
    /// A FASTA or plain-text DNA sequence for on-target design.
    Sequence,
    /// A reference FASTA file or a folder of chromosome FASTA files.
    Genome,
    /// A GTF/GFF3 gene annotation.
    Annotation,
    /// A TSV export destination.
    Tsv,
}

impl FileKind {
    /// Window title for the given explorer mode.
    fn title(self, mode: ExplorerMode) -> &'static str {
        match (self, mode) {
            (FileKind::Sequence, ExplorerMode::Open) => "Open sequence file",
            (FileKind::Genome, ExplorerMode::Open) => "Select reference genome",
            (FileKind::Annotation, ExplorerMode::Open) => "Select gene annotation",
            (FileKind::Tsv, ExplorerMode::Save) => "Save TSV",
            (FileKind::Tsv, ExplorerMode::Open) => "Open TSV",
            (_, ExplorerMode::Save) => "Save file",
        }
    }

    /// One-line explanation shown at the top of the explorer.
    fn description(self, mode: ExplorerMode) -> &'static str {
        match (self, mode) {
            (FileKind::Sequence, _) => {
                "Choose a FASTA or plain-text DNA file. The first FASTA record is used for design."
            }
            (FileKind::Genome, _) => {
                "Choose a reference FASTA file or a folder of chromosome FASTA files. The whole reference is loaded into memory."
            }
            (FileKind::Annotation, _) => {
                "Choose a GTF or GFF3 gene annotation. Compressed files must be decompressed first."
            }
            (FileKind::Tsv, ExplorerMode::Save) => {
                "Choose a folder and file name for the TSV export."
            }
            (FileKind::Tsv, ExplorerMode::Open) => "Choose a TSV file.",
        }
    }

    /// Extensions highlighted by the filter, all lower-case.
    fn extensions(self) -> &'static [&'static str] {
        match self {
            FileKind::Sequence => &["fa", "fasta", "fna", "fas", "ffn", "seq", "txt"],
            FileKind::Genome => &["fa", "fasta", "fna", "fas", "ffn", "seq"],
            FileKind::Annotation => &["gtf", "gff3", "gff", "txt"],
            FileKind::Tsv => &["tsv"],
        }
    }

    /// Default file name used when save mode starts without a path.
    fn default_file_name(self) -> &'static str {
        match self {
            FileKind::Tsv => "arachnia_export.tsv",
            FileKind::Sequence => "sequence.fa",
            FileKind::Genome => "genome.fa",
            FileKind::Annotation => "annotation.gtf",
        }
    }

    /// Whether the current folder itself can be selected as the result.
    fn allows_folders(self) -> bool {
        matches!(self, FileKind::Genome)
    }

    /// Whether a regular file has an extension appropriate for this kind.
    ///
    /// The comparison is case-insensitive. Files without an extension do not
    /// match; the **All files** checkbox can still reveal them.
    fn matches_file(self, path: &Path) -> bool {
        let Some(ext) = path.extension().and_then(|ext| ext.to_str()) else {
            return false;
        };
        let ext = ext.to_ascii_lowercase();
        self.extensions().iter().any(|candidate| *candidate == ext)
    }
}

/// Whether the explorer is choosing an existing path or a new one.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ExplorerMode {
    Open,
    Save,
}

/// A detected file or folder offered as a one-click shortcut.
#[derive(Clone, Debug)]
struct QuickPick {
    label: &'static str,
    path: PathBuf,
}

/// One row in the directory listing.
#[derive(Clone, Debug)]
struct Entry {
    name: String,
    path: PathBuf,
    is_dir: bool,
    size: u64,
}

/// A reusable in-app file explorer window.
///
/// See the [module documentation](self) for behavior and integration notes.
pub(crate) struct FileExplorer {
    /// Unique `egui` window id for this instance. The application has one
    /// explorer for the main window and one for the off-target panel; both
    /// may exist at the same time.
    instance: u64,
    open: bool,
    mode: ExplorerMode,
    kind: FileKind,
    /// Current directory shown in the list.
    dir: PathBuf,
    /// Text shown in the editable path bar.
    path_input: String,
    /// Highlighted file, if any.
    selected: Option<PathBuf>,
    /// File name used in save mode.
    file_name: String,
    entries: Vec<Entry>,
    quick_picks: Vec<QuickPick>,
    show_hidden: bool,
    show_all_files: bool,
    /// Error from navigation or validation, shown above the list.
    error: Option<String>,
    /// Error from reading the current directory.
    list_error: Option<String>,
    /// Whether the entry cap truncated the current listing.
    truncated: bool,
    /// Set after the first save click when the target already exists.
    confirm_overwrite: Option<PathBuf>,
}

static NEXT_INSTANCE: AtomicU64 = AtomicU64::new(1);

impl Default for FileExplorer {
    fn default() -> Self {
        let dir = std::env::current_dir().unwrap_or_else(|_| PathBuf::from("."));
        Self {
            instance: NEXT_INSTANCE.fetch_add(1, Ordering::Relaxed),
            open: false,
            mode: ExplorerMode::Open,
            kind: FileKind::Sequence,
            path_input: dir.display().to_string(),
            dir,
            selected: None,
            file_name: String::new(),
            entries: Vec::new(),
            quick_picks: Vec::new(),
            show_hidden: false,
            show_all_files: false,
            error: None,
            list_error: None,
            truncated: false,
            confirm_overwrite: None,
        }
    }
}

impl FileExplorer {
    /// The kind of path this explorer is currently choosing.
    pub(crate) fn kind(&self) -> FileKind {
        self.kind
    }

    /// Close the explorer without returning a path.
    pub(crate) fn close(&mut self) {
        self.open = false;
        self.confirm_overwrite = None;
    }

    /// Start an open-file session.
    ///
    /// `current` is the path already present in the caller's text field. It
    /// may be empty, an existing file or folder, a relative path, a leading
    /// `~` path, or the `hg38`/`GRCh38` shorthand for a genome. If it points at
    /// a file, the explorer starts in its parent folder with the file
    /// highlighted. If no usable path is supplied, a sensible default folder
    /// is chosen (the detected hg38/annotation folder for those kinds, and the
    /// working directory otherwise).
    pub(crate) fn open_for_open(&mut self, kind: FileKind, current: &str) {
        self.begin(kind, ExplorerMode::Open, current);
    }

    /// Start a save-file session.
    ///
    /// `current` is treated like [`Self::open_for_open`], except that a
    /// non-existent path is split into its parent folder and file name. For
    /// [`FileKind::Tsv`], the file name defaults to `arachnia_export.tsv` when
    /// no name is supplied.
    pub(crate) fn open_for_save(&mut self, kind: FileKind, current: &str) {
        self.begin(kind, ExplorerMode::Save, current);
    }

    /// Render the explorer and return a confirmed path.
    ///
    /// Returns `Some(path)` exactly once when the user clicks **Select**,
    /// **Use this folder**, or **Save** (including the overwrite
    /// confirmation). Returns `None` while the window is open or after it is
    /// cancelled.
    pub(crate) fn show(&mut self, ctx: &egui::Context) -> Option<PathBuf> {
        if !self.open {
            return None;
        }

        let mut open = self.open;
        let mut chosen = None;
        let mut cancel = false;
        egui::Window::new(self.kind.title(self.mode))
            .id(egui::Id::unique(("arachnia_file_explorer", self.instance)))
            .open(&mut open)
            .resizable(true)
            .collapsible(false)
            .default_size([760.0, 560.0])
            .min_width(560.0)
            .min_height(380.0)
            .show(ctx, |ui| {
                if ui.input(|input| input.key_pressed(egui::Key::Escape)) {
                    cancel = true;
                }
                if let Some(path) = self.body(ui) {
                    chosen = Some(path);
                }
            });

        if cancel || !open {
            self.close();
        }
        if chosen.is_some() {
            self.open = false;
            self.confirm_overwrite = None;
        }
        chosen
    }

    /// Start a session and resolve the initial folder and selection.
    fn begin(&mut self, kind: FileKind, mode: ExplorerMode, current: &str) {
        self.open = true;
        self.mode = mode;
        self.kind = kind;
        self.error = None;
        self.list_error = None;
        self.truncated = false;
        self.confirm_overwrite = None;
        self.show_hidden = false;
        self.show_all_files = false;
        self.selected = None;
        self.file_name.clear();

        let trimmed = current.trim();
        let supplied = if trimmed.is_empty() {
            match kind {
                FileKind::Genome => find_hg38(),
                FileKind::Annotation => find_hg38_annotation(),
                FileKind::Sequence | FileKind::Tsv => None,
            }
        } else if kind == FileKind::Genome && is_hg38_alias(trimmed) {
            find_hg38()
        } else {
            Some(expand_tilde(PathBuf::from(trimmed)))
        };

        match supplied {
            Some(path) if path.is_file() => {
                self.dir = parent_dir(&path).unwrap_or_else(|| fallback_dir(kind));
                self.selected = Some(path.clone());
                self.file_name = file_name_of(&path);
            }
            Some(path) if path.is_dir() => {
                self.dir = path;
                if mode == ExplorerMode::Save {
                    self.file_name = kind.default_file_name().to_string();
                }
            }
            Some(path) => {
                // A non-existent save target (or a stale open path): start in
                // its parent when possible and keep the file name for saving.
                let parent = parent_dir(&path);
                if let Some(parent) = parent.filter(|parent| parent.is_dir()) {
                    self.dir = parent;
                    if mode == ExplorerMode::Save {
                        self.file_name = file_name_of(&path);
                    }
                } else {
                    self.dir = fallback_dir(kind);
                    if mode == ExplorerMode::Save && !trimmed.is_empty() {
                        self.file_name = file_name_of(&path);
                    }
                }
            }
            None => {
                self.dir = fallback_dir(kind);
                if mode == ExplorerMode::Save {
                    self.file_name = kind.default_file_name().to_string();
                }
            }
        }

        if self.dir.as_os_str().is_empty() {
            self.dir = fallback_dir(kind);
        }
        if mode == ExplorerMode::Save && self.file_name.is_empty() {
            self.file_name = kind.default_file_name().to_string();
        }

        self.path_input = self.dir.display().to_string();
        self.quick_picks = quick_picks(kind);
        self.refresh();
    }

    /// Draw the window contents. Returns a path when the user confirms one.
    fn body(&mut self, ui: &mut egui::Ui) -> Option<PathBuf> {
        ui.label(
            RichText::new(self.kind.description(self.mode))
                .size(11.5)
                .color(MUTED),
        );

        if !self.quick_picks.is_empty() {
            ui.add_space(2.0);
            let mut picked = None;
            ui.horizontal_wrapped(|ui| {
                for pick in &self.quick_picks {
                    let button = ui
                        .button(RichText::new(pick.label).size(10.5).color(ACCENT_TEXT))
                        .on_hover_text(pick.path.display().to_string());
                    if button.clicked() {
                        picked = Some(pick.path.clone());
                    }
                }
            });
            if let Some(path) = picked {
                self.select_path(path);
            }
        }

        ui.add_space(3.0);
        ui.separator();
        ui.add_space(3.0);

        let mut navigate_to = None;
        let mut refresh = false;
        ui.horizontal(|ui| {
            let up = self
                .dir
                .parent()
                .filter(|parent| !parent.as_os_str().is_empty())
                .map(Path::to_path_buf);
            if ui
                .add_enabled(up.is_some(), egui::Button::new("Up"))
                .on_hover_text("Go to the parent folder")
                .clicked()
            {
                navigate_to = up;
            }
            if ui
                .button("Home")
                .on_hover_text("Go to your home folder")
                .clicked()
            {
                navigate_to = home_dir();
            }
            if ui
                .button("Data")
                .on_hover_text(
                    "Open the Arachnia data folder used by fetch_hg38.sh and fetch_annotation.sh",
                )
                .clicked()
            {
                navigate_to = data_dir();
            }
            if ui
                .button("Refresh")
                .on_hover_text("Re-read this folder")
                .clicked()
            {
                refresh = true;
            }

            let width = (ui.available_width() - 46.0).max(120.0);
            let response = ui.add(
                egui::TextEdit::singleline(&mut self.path_input)
                    .desired_width(width)
                    .hint_text("/path/to/folder"),
            );
            let enter = response.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter));
            if enter || ui.button("Go").clicked() {
                navigate_to = Some(expand_tilde(PathBuf::from(self.path_input.trim())));
            }
        });
        if let Some(path) = navigate_to {
            self.navigate(path);
        } else if refresh {
            self.refresh();
        }

        ui.horizontal(|ui| {
            ui.checkbox(&mut self.show_hidden, "Hidden")
                .on_hover_text("Show names beginning with a dot");
            if self.mode == ExplorerMode::Open {
                ui.checkbox(&mut self.show_all_files, "All files")
                    .on_hover_text("Show files with any extension");
            }
            ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                if self.truncated {
                    ui.label(
                        RichText::new(format!("showing first {MAX_ENTRIES} entries"))
                            .size(10.0)
                            .color(AMBER),
                    );
                }
            });
        });

        if let Some(error) = &self.error {
            ui.label(RichText::new(error).size(10.5).color(RED));
        }

        ui.add_space(2.0);
        let mut select = None;
        let mut enter_dir = None;
        egui::Frame::new()
            .fill(CARD)
            .corner_radius(8)
            .stroke(Stroke::new(1.0, BORDER))
            .inner_margin(Margin::symmetric(8, 6))
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    ui.label(RichText::new("Name").size(10.0).strong().color(MUTED));
                    ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                        ui.label(RichText::new("Type").size(10.0).strong().color(MUTED));
                        ui.add_space(54.0);
                        ui.label(RichText::new("Size").size(10.0).strong().color(MUTED));
                    });
                });
                ui.separator();

                if let Some(error) = &self.list_error {
                    ui.add_space(8.0);
                    ui.label(RichText::new(error).size(10.5).color(RED));
                    return;
                }

                let mut visible = 0usize;
                egui::ScrollArea::vertical()
                    .id_salt("file_explorer_list")
                    .max_height(310.0)
                    .min_scrolled_height(180.0)
                    .auto_shrink([false, false])
                    .show(ui, |ui| {
                        let total = self.entries.len();
                        for index in 0..total {
                            let (path, name, is_dir, size) = {
                                let entry = &self.entries[index];
                                (
                                    entry.path.clone(),
                                    entry.name.clone(),
                                    entry.is_dir,
                                    entry.size,
                                )
                            };
                            if !self.show_hidden && name.starts_with('.') {
                                continue;
                            }
                            if !is_dir
                                && self.mode == ExplorerMode::Open
                                && !self.show_all_files
                                && !self.kind.matches_file(&path)
                            {
                                continue;
                            }
                            visible += 1;
                            let selected = self.selected.as_deref() == Some(path.as_path());
                            let label = if is_dir {
                                format!("{name}/")
                            } else {
                                name.clone()
                            };
                            let color = if is_dir { ACCENT_TEXT } else { TEXT };
                            let type_label = entry_type(&path, is_dir);
                            let name_width = (ui.available_width() - 170.0).max(180.0);
                            let response = ui
                                .horizontal(|ui| {
                                    let response = ui.add_sized(
                                        [name_width, 20.0],
                                        egui::Button::selectable(
                                            selected,
                                            RichText::new(label).size(12.0).color(color),
                                        ),
                                    );
                                    ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                                        ui.label(
                                            RichText::new(if is_dir {
                                                String::new()
                                            } else {
                                                format_size(size)
                                            })
                                            .size(10.0)
                                            .color(MUTED),
                                        );
                                        ui.add_space(12.0);
                                        ui.label(
                                            RichText::new(type_label).size(10.0).color(MUTED),
                                        );
                                    });
                                    response
                                })
                                .inner
                                .on_hover_text(path.display().to_string());

                            if response.clicked() {
                                if is_dir {
                                    enter_dir = Some(path.clone());
                                } else {
                                    self.selected = Some(path.clone());
                                    self.confirm_overwrite = None;
                                    if self.mode == ExplorerMode::Save {
                                        self.file_name = name.clone();
                                    }
                                }
                            }
                            if response.double_clicked() && !is_dir {
                                if self.mode == ExplorerMode::Open {
                                    select = Some(path.clone());
                                } else {
                                    self.file_name = name.clone();
                                }
                            }
                        }
                        if visible == 0 {
                            ui.add_space(24.0);
                            ui.vertical_centered(|ui| {
                                ui.label(
                                    RichText::new(if self.mode == ExplorerMode::Open {
                                        "No matching files in this folder"
                                    } else {
                                        "This folder is empty"
                                    })
                                    .size(12.0)
                                    .color(MUTED),
                                );
                                if self.mode == ExplorerMode::Open
                                    && !self.show_all_files
                                    && !self.kind.extensions().is_empty()
                                {
                                    ui.label(
                                        RichText::new(
                                            "Tick “All files” to show every file, or navigate to another folder.",
                                        )
                                        .size(10.5)
                                        .color(MUTED),
                                    );
                                }
                            });
                        }
                    });
            });

        if let Some(path) = enter_dir {
            self.navigate(path);
        }

        ui.add_space(4.0);

        match self.mode {
            ExplorerMode::Open => self.open_footer(ui, select),
            ExplorerMode::Save => self.save_footer(ui),
        }
    }

    /// Footer for open mode: selected path, folder selection and actions.
    fn open_footer(&mut self, ui: &mut egui::Ui, select: Option<PathBuf>) -> Option<PathBuf> {
        let mut chosen = select;
        let mut cancel = false;
        if let Some(selected) = &self.selected {
            let size = fs::metadata(selected).map(|meta| meta.len()).unwrap_or(0);
            ui.add(
                egui::Label::new(
                    RichText::new(selected.display().to_string())
                        .size(11.0)
                        .monospace()
                        .color(TEXT),
                )
                .truncate(),
            );
            ui.label(
                RichText::new(format!(
                    "{} · {}",
                    entry_type(selected, selected.is_dir()),
                    if selected.is_dir() {
                        "folder".to_string()
                    } else {
                        format_size(size)
                    }
                ))
                .size(10.0)
                .color(MUTED),
            );
        } else {
            ui.label(
                RichText::new("Select a file from the list, or use the current folder.")
                    .size(11.0)
                    .color(MUTED),
            );
        }
        ui.horizontal(|ui| {
            ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                if ui.button("Cancel").clicked() {
                    cancel = true;
                }
                let can_select = self.selected.as_ref().is_some_and(|path| {
                    path.is_file() || (self.kind.allows_folders() && path.is_dir())
                });
                if ui
                    .add_enabled(can_select, egui::Button::new("Select"))
                    .clicked()
                {
                    chosen = self.selected.clone();
                }
                if self.kind.allows_folders() && ui.button("Use this folder").clicked() {
                    chosen = Some(self.dir.clone());
                }
            });
        });
        if cancel {
            self.close();
            return None;
        }
        if let Some(path) = chosen {
            self.open = false;
            self.confirm_overwrite = None;
            return Some(path);
        }
        None
    }

    /// Footer for save mode: file name, resolved path and overwrite check.
    fn save_footer(&mut self, ui: &mut egui::Ui) -> Option<PathBuf> {
        ui.separator();
        ui.add_space(2.0);
        let mut save = false;
        let mut cancel = false;
        let mut chosen = None;
        ui.horizontal(|ui| {
            ui.label(RichText::new("File name").size(11.5).color(MUTED));
            let response = ui.add(
                egui::TextEdit::singleline(&mut self.file_name)
                    .desired_width((ui.available_width() - 210.0).max(120.0))
                    .hint_text(self.kind.default_file_name()),
            );
            let enter = response.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter));
            if ui.button("Cancel").clicked() {
                cancel = true;
            }

            let target = resolve_save_path(&self.dir, &self.file_name, self.kind);
            let exists = target.as_ref().ok().is_some_and(|path| path.exists());
            let confirmed = target
                .as_ref()
                .ok()
                .is_some_and(|path| self.confirm_overwrite.as_ref() == Some(path));
            let label = if exists && confirmed {
                "Overwrite"
            } else {
                "Save"
            };
            let button = egui::Button::new(label)
                .fill(ACCENT)
                .corner_radius(6)
                .stroke(Stroke::NONE);
            if ui.add_enabled(target.is_ok(), button).clicked() {
                save = true;
            }
            if enter {
                save = true;
            }

            if save {
                match target {
                    Ok(path) => {
                        if path.exists() && self.confirm_overwrite.as_ref() != Some(&path) {
                            self.confirm_overwrite = Some(path);
                            self.error = Some(
                                "That file already exists. Click Overwrite to replace it."
                                    .to_string(),
                            );
                        } else {
                            chosen = Some(path);
                        }
                    }
                    Err(err) => {
                        self.error = Some(err);
                    }
                }
            }
        });

        if let Ok(path) = resolve_save_path(&self.dir, &self.file_name, self.kind) {
            ui.label(
                RichText::new(format!("Save to {}", path.display()))
                    .size(10.0)
                    .monospace()
                    .color(MUTED),
            );
        }
        if let Some(error) = &self.error {
            ui.label(RichText::new(error).size(10.5).color(RED));
        }
        if cancel {
            self.close();
            return None;
        }
        if let Some(path) = chosen {
            self.open = false;
            self.confirm_overwrite = None;
            return Some(path);
        }
        None
    }

    /// Re-read the current folder and rebuild the visible entry list.
    fn refresh(&mut self) {
        self.entries.clear();
        self.list_error = None;
        self.truncated = false;

        let read_dir = match fs::read_dir(&self.dir) {
            Ok(read_dir) => read_dir,
            Err(err) => {
                self.list_error = Some(format!("Could not read {}: {err}", self.dir.display()));
                return;
            }
        };

        for entry in read_dir {
            if self.entries.len() >= MAX_ENTRIES {
                self.truncated = true;
                break;
            }
            let Ok(entry) = entry else {
                continue;
            };
            let path = entry.path();
            let is_dir = path.is_dir();
            let name = entry.file_name().to_string_lossy().into_owned();
            let size = if is_dir {
                0
            } else {
                entry.metadata().map(|meta| meta.len()).unwrap_or(0)
            };
            self.entries.push(Entry {
                name,
                path,
                is_dir,
                size,
            });
        }

        self.entries.sort_by(|a, b| {
            b.is_dir
                .cmp(&a.is_dir)
                .then_with(|| a.name.to_lowercase().cmp(&b.name.to_lowercase()))
                .then_with(|| a.name.cmp(&b.name))
        });
    }

    /// Change folders, clearing the previous selection.
    fn navigate(&mut self, path: PathBuf) {
        let path = expand_tilde(path);
        if path.is_dir() {
            self.dir = path;
            self.path_input = self.dir.display().to_string();
            self.selected = None;
            self.error = None;
            self.confirm_overwrite = None;
            self.refresh();
        } else if path.is_file() {
            self.select_path(path);
        } else {
            self.error = Some(format!("Folder not found: {}", path.display()));
        }
    }

    /// Select a file (or navigate into a folder) without confirming it.
    fn select_path(&mut self, path: PathBuf) {
        if path.is_dir() {
            self.navigate(path);
        } else if path.is_file() {
            if let Some(parent) = parent_dir(&path) {
                self.dir = parent;
                self.path_input = self.dir.display().to_string();
            }
            self.file_name = file_name_of(&path);
            self.selected = Some(path);
            self.error = None;
            self.confirm_overwrite = None;
            self.refresh();
        } else {
            self.error = Some(format!("File not found: {}", path.display()));
        }
    }
}

/// Resolve a save destination from the current folder and a file-name field.
///
/// An absolute file name is used as entered; a relative one is joined to
/// `dir`. A leading `~` is expanded. For [`FileKind::Tsv`], a missing
/// extension is completed to `.tsv`. Empty names and missing parent folders
/// are rejected with a user-facing message.
fn resolve_save_path(dir: &Path, file_name: &str, kind: FileKind) -> Result<PathBuf, String> {
    let trimmed = file_name.trim();
    if trimmed.is_empty() {
        return Err("Enter a file name first.".to_string());
    }
    let mut path = expand_tilde(PathBuf::from(trimmed));
    if !path.is_absolute() {
        path = dir.join(path);
    }
    if path.file_name().is_none() {
        return Err("Enter a file name, not a folder.".to_string());
    }
    if kind == FileKind::Tsv && path.extension().is_none() {
        path.set_extension("tsv");
    }
    if path.is_dir() {
        return Err("That path is a folder; enter a file name.".to_string());
    }
    let Some(parent) = path.parent() else {
        return Err("The destination folder does not exist.".to_string());
    };
    if !parent.is_dir() {
        return Err(format!(
            "The destination folder does not exist: {}",
            parent.display()
        ));
    }
    Ok(path)
}

/// A sensible starting folder for a kind when the caller supplies no path.
fn fallback_dir(kind: FileKind) -> PathBuf {
    let detected = match kind {
        FileKind::Genome => find_hg38(),
        FileKind::Annotation => find_hg38_annotation(),
        FileKind::Sequence | FileKind::Tsv => None,
    };
    if let Some(path) = detected {
        if path.is_dir() {
            return path;
        }
        if let Some(parent) = path.parent() {
            return parent.to_path_buf();
        }
    }
    if matches!(kind, FileKind::Sequence | FileKind::Tsv) {
        if let Ok(cwd) = std::env::current_dir() {
            return cwd;
        }
    }
    if let Some(data) = data_dir().filter(|path| path.is_dir()) {
        return data;
    }
    home_dir()
        .or_else(|| std::env::current_dir().ok())
        .unwrap_or_else(|| PathBuf::from("."))
}

/// Quick-pick shortcuts shown for the current kind.
fn quick_picks(kind: FileKind) -> Vec<QuickPick> {
    let mut picks = Vec::new();
    match kind {
        FileKind::Genome => {
            if let Some(path) = find_hg38() {
                picks.push(QuickPick {
                    label: "Use detected hg38",
                    path,
                });
            }
        }
        FileKind::Annotation => {
            if let Some(path) = find_hg38_annotation() {
                picks.push(QuickPick {
                    label: "Use detected annotation",
                    path,
                });
            }
        }
        FileKind::Sequence | FileKind::Tsv => {}
    }
    if let Some(data) = data_dir().filter(|path| path.is_dir()) {
        picks.push(QuickPick {
            label: "Arachnia data folder",
            path: data,
        });
    }
    picks
}

/// `$HOME`, if set.
fn home_dir() -> Option<PathBuf> {
    std::env::var_os("HOME").map(PathBuf::from)
}

/// `${XDG_DATA_HOME:-$HOME/.local/share}/arachnia`.
///
/// The folder is returned even when it does not exist; callers check
/// [`Path::is_dir`] before navigating there.
fn data_dir() -> Option<PathBuf> {
    data_dir_from(
        std::env::var_os("XDG_DATA_HOME").map(PathBuf::from),
        home_dir(),
    )
}

/// Pure helper behind [`data_dir`], kept separate for unit testing.
fn data_dir_from(xdg_data_home: Option<PathBuf>, home: Option<PathBuf>) -> Option<PathBuf> {
    xdg_data_home
        .or_else(|| home.map(|home| home.join(".local/share")))
        .map(|base| base.join("arachnia"))
}

/// File name component as lossy UTF-8, or an empty string.
fn file_name_of(path: &Path) -> String {
    path.file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_default()
}

/// Parent folder of a path, resolving an empty relative parent to the current
/// working directory.
fn parent_dir(path: &Path) -> Option<PathBuf> {
    let parent = path.parent()?;
    if parent.as_os_str().is_empty() {
        std::env::current_dir().ok()
    } else {
        Some(parent.to_path_buf())
    }
}

/// Short type label for the listing.
fn entry_type(path: &Path, is_dir: bool) -> String {
    if is_dir {
        return "folder".to_string();
    }
    match path
        .extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| ext.to_ascii_lowercase())
        .as_deref()
    {
        Some("fa" | "fasta" | "fna" | "fas" | "ffn" | "seq") => "FASTA".to_string(),
        Some("gtf") => "GTF".to_string(),
        Some("gff" | "gff3") => "GFF3".to_string(),
        Some("tsv") => "TSV".to_string(),
        Some("txt") => "text".to_string(),
        Some(ext) if !ext.is_empty() => ext.to_ascii_uppercase(),
        _ => "file".to_string(),
    }
}

/// Human-readable binary file size (B, KiB, MiB, GiB, TiB).
fn format_size(bytes: u64) -> String {
    const UNITS: [&str; 5] = ["B", "KiB", "MiB", "GiB", "TiB"];
    if bytes < 1024 {
        return format!("{bytes} B");
    }
    let mut value = bytes as f64;
    let mut unit = 0;
    while value >= 1024.0 && unit + 1 < UNITS.len() {
        value /= 1024.0;
        unit += 1;
    }
    format!("{value:.1} {}", UNITS[unit])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn file_kind_matching_is_case_insensitive() {
        assert!(FileKind::Genome.matches_file(Path::new("hg38.fa")));
        assert!(FileKind::Genome.matches_file(Path::new("HG38.FASTA")));
        assert!(!FileKind::Genome.matches_file(Path::new("annotation.gtf")));
        assert!(FileKind::Annotation.matches_file(Path::new("genes.GFF3")));
        assert!(FileKind::Sequence.matches_file(Path::new("target.txt")));
        assert!(!FileKind::Genome.matches_file(Path::new("README")));
        assert!(FileKind::Tsv.matches_file(Path::new("out.tsv")));
    }

    #[test]
    fn save_path_adds_tsv_extension() {
        let dir = std::env::temp_dir();
        assert_eq!(
            resolve_save_path(&dir, "guides", FileKind::Tsv).unwrap(),
            dir.join("guides.tsv")
        );
        assert_eq!(
            resolve_save_path(&dir, "guides.tsv", FileKind::Tsv).unwrap(),
            dir.join("guides.tsv")
        );
        let absolute = dir.join("absolute.tsv");
        assert_eq!(
            resolve_save_path(&dir, absolute.to_str().unwrap(), FileKind::Tsv).unwrap(),
            absolute
        );
    }

    #[test]
    fn save_path_rejects_empty_names_and_missing_parents() {
        let dir = std::env::temp_dir();
        assert!(resolve_save_path(&dir, "  ", FileKind::Tsv).is_err());
        assert!(resolve_save_path(&dir, "missing/sub/out.tsv", FileKind::Tsv).is_err());
    }

    #[test]
    fn size_formatting_uses_binary_units() {
        assert_eq!(format_size(0), "0 B");
        assert_eq!(format_size(1023), "1023 B");
        assert_eq!(format_size(1024), "1.0 KiB");
        assert_eq!(format_size(3 * 1024 * 1024), "3.0 MiB");
        assert_eq!(format_size(3 * 1024 * 1024 * 1024), "3.0 GiB");
    }

    #[test]
    fn data_dir_prefers_xdg_data_home() {
        assert_eq!(
            data_dir_from(
                Some(PathBuf::from("/tmp/arachnia-test-data")),
                Some(PathBuf::from("/home/tester"))
            )
            .unwrap(),
            PathBuf::from("/tmp/arachnia-test-data/arachnia")
        );
        assert_eq!(
            data_dir_from(None, Some(PathBuf::from("/home/tester"))).unwrap(),
            PathBuf::from("/home/tester/.local/share/arachnia")
        );
        assert_eq!(data_dir_from(None, None), None);
    }

    #[test]
    fn explorer_renders_in_a_headless_context() {
        let ctx = egui::Context::default();
        let mut explorer = FileExplorer::default();
        explorer.open_for_open(FileKind::Genome, "");
        let mut output = ctx.run_ui(egui::RawInput::default(), |ui| {
            let ctx = ui.ctx().clone();
            let _ = explorer.show(&ctx);
        });
        output.textures_delta.clear();
        explorer.open_for_save(FileKind::Tsv, "arachnia-test-output.tsv");
        let mut output = ctx.run_ui(egui::RawInput::default(), |ui| {
            let ctx = ui.ctx().clone();
            let _ = explorer.show(&ctx);
        });
        output.textures_delta.clear();
    }
}
