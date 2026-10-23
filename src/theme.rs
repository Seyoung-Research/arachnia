//! Obsidian & crimson palette for the Arachnia UI.
//!
//! The palette uses near-black obsidian backgrounds, warm off-white text,
//! muted crimson accents, and restrained green/amber/red status colors. The
//! constants here are shared by the main window and the off-target panel;
//! [`apply`] installs the matching `egui` visuals and text styles once per
//! frame context.

use eframe::egui::{self, Color32, FontId, Stroke, TextStyle};

pub(crate) const HEADER_BG: Color32 = Color32::from_rgb(26, 10, 14);
pub(crate) const HEADER_BORDER: Color32 = Color32::from_rgb(74, 24, 33);
pub(crate) const HEADER_MUTED: Color32 = Color32::from_rgb(199, 154, 162);
pub(crate) const ACCENT: Color32 = Color32::from_rgb(198, 42, 65);
pub(crate) const ACCENT_TEXT: Color32 = Color32::from_rgb(255, 138, 155);
pub(crate) const ACCENT_SOFT: Color32 = Color32::from_rgb(60, 29, 36);
pub(crate) const BG: Color32 = Color32::from_rgb(13, 10, 11);
pub(crate) const CARD: Color32 = Color32::from_rgb(22, 17, 20);
pub(crate) const CARD_ALT: Color32 = Color32::from_rgb(28, 21, 24);
pub(crate) const INPUT_BG: Color32 = Color32::from_rgb(16, 12, 13);
pub(crate) const WEAK_BG: Color32 = Color32::from_rgb(33, 24, 27);
pub(crate) const TEXT: Color32 = Color32::from_rgb(243, 235, 237);
pub(crate) const MUTED: Color32 = Color32::from_rgb(166, 143, 148);
pub(crate) const BORDER: Color32 = Color32::from_rgb(50, 36, 40);
pub(crate) const GREEN: Color32 = Color32::from_rgb(125, 199, 154);
pub(crate) const AMBER: Color32 = Color32::from_rgb(226, 181, 108);
pub(crate) const RED: Color32 = Color32::from_rgb(255, 125, 139);
pub(crate) const WHITE: Color32 = Color32::from_rgb(255, 247, 248);
pub(crate) const ERROR_BG: Color32 = Color32::from_rgb(43, 16, 21);

pub(crate) fn apply(ctx: &egui::Context) {
    let mut style = (*ctx.style_of(egui::Theme::Dark)).clone();
    style.visuals = egui::Visuals::dark();
    style.visuals.override_text_color = Some(TEXT);
    style.visuals.panel_fill = BG;
    style.visuals.window_fill = CARD;
    style.visuals.extreme_bg_color = INPUT_BG;
    style.visuals.faint_bg_color = CARD_ALT;
    style.visuals.code_bg_color = INPUT_BG;
    style.visuals.hyperlink_color = ACCENT_TEXT;
    style.visuals.warn_fg_color = AMBER;
    style.visuals.error_fg_color = RED;
    style.visuals.selection.bg_fill = ACCENT_SOFT;
    style.visuals.selection.stroke = Stroke::new(1.0, ACCENT_TEXT);
    style.visuals.window_stroke = Stroke::new(1.0, BORDER);
    style.visuals.widgets.noninteractive.bg_fill = CARD;
    style.visuals.widgets.noninteractive.weak_bg_fill = CARD_ALT;
    style.visuals.widgets.noninteractive.fg_stroke = Stroke::new(1.0, TEXT);
    style.visuals.widgets.noninteractive.bg_stroke = Stroke::new(1.0, BORDER);
    style.visuals.widgets.noninteractive.corner_radius = egui::CornerRadius::same(6);
    style.visuals.widgets.inactive.bg_fill = WEAK_BG;
    style.visuals.widgets.inactive.weak_bg_fill = WEAK_BG;
    style.visuals.widgets.inactive.fg_stroke = Stroke::new(1.0, TEXT);
    style.visuals.widgets.inactive.bg_stroke = Stroke::new(1.0, BORDER);
    style.visuals.widgets.inactive.corner_radius = egui::CornerRadius::same(6);
    style.visuals.widgets.hovered.bg_fill = ACCENT_SOFT;
    style.visuals.widgets.hovered.weak_bg_fill = ACCENT_SOFT;
    style.visuals.widgets.hovered.fg_stroke = Stroke::new(1.0, WHITE);
    style.visuals.widgets.hovered.bg_stroke = Stroke::new(1.0, ACCENT);
    style.visuals.widgets.hovered.corner_radius = egui::CornerRadius::same(6);
    style.visuals.widgets.active.bg_fill = ACCENT_SOFT;
    style.visuals.widgets.active.weak_bg_fill = ACCENT_SOFT;
    style.visuals.widgets.active.fg_stroke = Stroke::new(1.0, WHITE);
    style.visuals.widgets.active.bg_stroke = Stroke::new(1.0, ACCENT);
    style.visuals.widgets.active.corner_radius = egui::CornerRadius::same(6);
    style.visuals.widgets.open.weak_bg_fill = ACCENT_SOFT;
    style.visuals.window_corner_radius = egui::CornerRadius::same(10);
    style.spacing.item_spacing = egui::vec2(8.0, 7.0);
    style.spacing.button_padding = egui::vec2(10.0, 5.0);
    style.spacing.interact_size.y = 24.0;
    style
        .text_styles
        .insert(TextStyle::Body, FontId::proportional(14.0));
    style
        .text_styles
        .insert(TextStyle::Button, FontId::proportional(14.0));
    style
        .text_styles
        .insert(TextStyle::Small, FontId::proportional(11.5));
    style
        .text_styles
        .insert(TextStyle::Heading, FontId::proportional(20.0));
    style
        .text_styles
        .insert(TextStyle::Monospace, FontId::monospace(13.0));
    ctx.set_theme(egui::ThemePreference::Dark);
    ctx.set_style_of(egui::Theme::Dark, style);
}

/// Red heat scale for on-target model scores: muted mauve (low) to bright
/// rose (high).
pub(crate) fn score_color(t: f64) -> Color32 {
    let t = t.clamp(0.0, 1.0);
    let start = (155.0, 124.0, 130.0);
    let end = (255.0, 138.0, 155.0);
    Color32::from_rgb(
        (start.0 + (end.0 - start.0) * t) as u8,
        (start.1 + (end.1 - start.1) * t) as u8,
        (start.2 + (end.2 - start.2) * t) as u8,
    )
}
