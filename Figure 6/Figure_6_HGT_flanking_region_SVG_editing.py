#!/usr/bin/env python3
"""
Edit final anti-stress HGT flanking-region SVGs into an Illustrator-friendly
style requested by the user.

Input:
  /Volumes/Haleo_DK01/WD_Takakia_F/metagenome-data/Figure5/
    HGT_flanking_plots_final_antistress/*.SVG

Output:
  HGT_flanking_plots_final_antistress/edited_ai_style/*_edited.SVG

Edits:
  - keep original gene arrows and homology ribbons;
  - remove original text labels in the transfer plot region;
  - add a grey title bar with italic gene name, Identity score, and Bin transfer;
  - add R/E colour bars at the left of the two flanking rows
    (Bin 1-118 = Rhizomatousphere/R; >118 = Endophytic/E).
  - restyle ordinary genes as grey and transferred genes as green;
  - thicken gene arrows and move the two flanking rows closer together.
"""

from __future__ import annotations

import csv
import html
import re
from pathlib import Path


ROOT = Path("/Volumes/Haleo_DK01/WD_Takakia_F")
FIG5 = ROOT / "metagenome-data" / "Figure5"
SVG_DIR = FIG5 / "HGT_flanking_plots_final_antistress"
OUT_DIR = SVG_DIR / "edited_ai_style"

SUMMARY_TSV = FIG5 / "Figure5_antistress_final_HGT_event_overlap_summary.tsv"
RECOMMENDED_TSV = FIG5 / "Figure5_antistress_final_HGT_main_text_recommended_12.tsv"
HGT_TSV = ROOT / "all_g_detected_HGTs.txt"

R_COLOUR = "#82B3D7"
E_COLOUR = "#CF6C72"
TITLE_FILL = "#D3D3D3"
BLACK = "#000000"
GENE_GREY = "#666666"
TRANSFER_GREEN = "#00E000"
HOMOLOGY_RED = "#FF5A5A"

ROW_MID_Y = 354.0
LOWER_ROW_CENTER_Y = 210.8268
UPPER_ROW_CENTER_Y = 497.8346
LOWER_ROW_MIN_Y = 169.8256
LOWER_ROW_MAX_Y = 251.8279
UPPER_ROW_MIN_Y = 456.8335
UPPER_ROW_MAX_Y = 538.8358
ROW_SHIFT = 35.0
TRACK_THICKNESS_FACTOR = 1.25
AXIS_STROKE_WIDTH = "3"

TITLE_X_PAD = 38.0
TITLE_Y = 38.0
TITLE_H = 74.0
CROP_X_PAD = 30.0
CROP_TOP_PAD = 18.0
CROP_BOTTOM_PAD = 30.0


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def bin_number(gene_or_bin: str) -> int | None:
    m = re.search(r"Bin[_ ](\d+)", gene_or_bin)
    if not m:
        return None
    return int(m.group(1))


def bin_label(bin_or_gene: str) -> str:
    n = bin_number(bin_or_gene)
    return f"Bin {n}" if n is not None else bin_or_gene


def habitat_colour(bin_or_gene: str) -> str:
    n = bin_number(bin_or_gene)
    if n is None:
        return "#BDBDBD"
    return R_COLOUR if 1 <= n <= 118 else E_COLOUR


def find_identity(hgt_rows: list[dict[str, str]], gene1: str, gene2: str) -> tuple[float | None, str]:
    for row in hgt_rows:
        g1 = row.get("Gene_1", "")
        g2 = row.get("Gene_2", "")
        if {g1, g2} == {gene1, gene2}:
            identity = row.get("Identity", "")
            direction = row.get("direction", "")
            try:
                return float(identity), direction
            except ValueError:
                return None, direction
    return None, ""


def direction_title(direction: str, fallback_gene1: str, fallback_gene2: str) -> str:
    if "-->" in direction:
        left, right = direction.split("-->", 1)
        return f"{bin_label(left)} -&gt; {bin_label(right)}"
    return f"{bin_label(fallback_gene1)} -&gt; {bin_label(fallback_gene2)}"


def compact_gene_label(gene: str) -> str:
    gene = gene.strip()
    gene = gene.replace("; ", "/").replace(";", "/")
    gene = gene.replace(" (", "(")
    return gene


def svg_size(svg_text: str) -> tuple[float, float]:
    m = re.search(r'<svg[^>]*\bviewBox="([^"]+)"', svg_text)
    if m:
        vals = [float(x) for x in re.split(r"[,\s]+", m.group(1).strip()) if x]
        if len(vals) == 4:
            return vals[2], vals[3]
    mw = re.search(r'<svg[^>]*\bwidth="([0-9.]+)"', svg_text)
    mh = re.search(r'<svg[^>]*\bheight="([0-9.]+)"', svg_text)
    return (float(mw.group(1)) if mw else 1417.0, float(mh.group(1)) if mh else 708.0)


def remove_original_text(svg_text: str) -> str:
    # Original SVGs place all transfer-plot labels in <text> elements.
    # Keep <title>/<desc> metadata but remove visible <text> labels.
    return re.sub(r"\s*<text\b[^>]*>.*?</text>", "", svg_text, flags=re.S)


def format_style(style: str, fill: str | None = None, stroke: str | None = None,
                 stroke_width: str | None = None) -> str:
    parts = []
    seen = set()
    for raw in style.split(";"):
        if ":" not in raw:
            continue
        key, value = raw.split(":", 1)
        key = key.strip()
        value = value.strip()
        if not key:
            continue
        if key == "fill" and fill is not None:
            value = fill
        elif key == "stroke" and stroke is not None:
            value = stroke
        elif key == "stroke-width" and stroke_width is not None:
            value = stroke_width
        parts.append(f"{key}: {value}")
        seen.add(key)
    if fill is not None and "fill" not in seen:
        parts.append(f"fill: {fill}")
    if stroke is not None and "stroke" not in seen:
        parts.append(f"stroke: {stroke}")
    if stroke_width is not None and "stroke-width" not in seen:
        parts.append(f"stroke-width: {stroke_width}")
    return "; ".join(parts) + ";"


def parse_points(points: str) -> list[tuple[float, float]]:
    nums = [float(x) for x in re.split(r"[,\s]+", points.strip()) if x]
    return list(zip(nums[0::2], nums[1::2]))


def serialize_points(points: list[tuple[float, float]]) -> str:
    return ", ".join(f"{x:.4f} {y:.4f}" for x, y in points)


def transform_row_y(y: float) -> float:
    """Thicken each row around its own centre and move both rows inward.

    The MetaCHIP drawing group is vertically flipped by SVG transform. These
    coordinates are still the original group coordinates; visual y positions for
    overlay elements are computed separately below.
    """
    if y < ROW_MID_Y:
        return LOWER_ROW_CENTER_Y + (y - LOWER_ROW_CENTER_Y) * TRACK_THICKNESS_FACTOR + ROW_SHIFT
    if y > ROW_MID_Y:
        return UPPER_ROW_CENTER_Y + (y - UPPER_ROW_CENTER_Y) * TRACK_THICKNESS_FACTOR - ROW_SHIFT
    return y


def display_y_from_group_y(group_y: float, svg_height: float) -> float:
    return svg_height - group_y


def point_bbox(points: list[tuple[float, float]]) -> tuple[float, float, float, float]:
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    return min(xs), min(ys), max(xs), max(ys)


def visible_text_gene_groups(svg_text: str) -> list[list[int]]:
    genes = []
    for m in re.finditer(r"<text\b[^>]*>.*?</text>", svg_text, flags=re.S):
        content = re.sub(r"<[^>]+>", "", m.group(0))
        gm = re.search(r"Gene_(\d+)", content)
        if gm:
            genes.append(int(gm.group(1)))
    groups: list[list[int]] = []
    for g in genes:
        if not groups or g != groups[-1][-1] + 1:
            groups.append([g])
        else:
            groups[-1].append(g)
    return groups[:2]


def highlighted_transfer_gene_numbers(svg_text: str) -> set[int]:
    """Return all Gene numbers that MetaCHIP labelled with an identity score."""
    out = set()
    for m in re.finditer(r"<text\b[^>]*>.*?</text>", svg_text, flags=re.S):
        content = re.sub(r"<[^>]+>", "", m.group(0))
        # Transferred genes are labelled as e.g. "Gene_02576 (81.4)" in blue.
        gm = re.search(r"Gene_(\d+)\s*\([0-9.]+\)", content)
        if gm:
            out.add(int(gm.group(1)))
    return out


def local_gene_number(gene_id: str) -> int | None:
    m = re.search(r"_(\d+)$", gene_id)
    if not m:
        return None
    return int(m.group(1))


def classify_gene_polygons(svg_text: str, gene1: str, gene2: str) -> tuple[set[int], str, str]:
    groups = visible_text_gene_groups(svg_text)
    highlighted_numbers = highlighted_transfer_gene_numbers(svg_text)
    fallback_numbers = {n for n in (local_gene_number(gene1), local_gene_number(gene2)) if n is not None}
    target_numbers = highlighted_numbers if highlighted_numbers else fallback_numbers

    polygon_matches = list(re.finditer(r"<polygon\b[^>]*>", svg_text))
    gene_poly_positions = []
    for poly_idx, m in enumerate(polygon_matches):
        pts_m = re.search(r'points="([^"]+)"', m.group(0))
        if not pts_m:
            continue
        points = parse_points(pts_m.group(1))
        if len(points) not in (5, 7):
            continue
        xmin, ymin, xmax, ymax = point_bbox(points)
        height = ymax - ymin
        # Gene arrows in MetaCHIP plots are 5/7-point polygons, roughly half
        # the background row height. Background and homology polygons have 4
        # points and are excluded above.
        if 20 <= height <= 60:
            gene_poly_positions.append(poly_idx)

    target_poly_indices: set[int] = set()
    row_to_gene_number: list[list[int]] = groups
    offset = 0
    row_target_gene: list[str | None] = []
    for group in row_to_gene_number:
        current_target = None
        for i, gene_num in enumerate(group):
            if offset + i >= len(gene_poly_positions):
                continue
            if gene_num in target_numbers:
                target_poly_indices.add(gene_poly_positions[offset + i])
                if gene_num == local_gene_number(gene1):
                    current_target = gene1
                elif gene_num == local_gene_number(gene2):
                    current_target = gene2
        row_target_gene.append(current_target)
        offset += len(group)

    # In the original transformed coordinate system, the first gene-label row
    # corresponds to the lower displayed row, while the second corresponds to
    # the upper displayed row. Fall back to gene1/gene2 if parsing fails.
    bottom_gene = row_target_gene[0] if len(row_target_gene) > 0 and row_target_gene[0] else gene2
    top_gene = row_target_gene[1] if len(row_target_gene) > 1 and row_target_gene[1] else gene1
    return target_poly_indices, top_gene, bottom_gene


def restyle_and_relayout_original_svg(svg_text: str, transfer_poly_indices: set[int]) -> str:
    poly_counter = -1

    def repl_polygon(match: re.Match[str]) -> str:
        nonlocal poly_counter
        poly_counter += 1
        tag = match.group(0)
        pts_m = re.search(r'points="([^"]+)"', tag)
        style_m = re.search(r'style="([^"]+)"', tag)
        if not pts_m:
            return tag
        points = parse_points(pts_m.group(1))
        is_gene = len(points) in (5, 7)
        # Move all original flanking plot geometry toward the middle and thicken
        # the full tracks, not only the arrows. This includes grey backgrounds,
        # red homology ribbons, and gene arrows.
        points = [(x, transform_row_y(y)) for x, y in points]
        tag = tag[:pts_m.start(1)] + serialize_points(points) + tag[pts_m.end(1):]
        style_m = re.search(r'style="([^"]+)"', tag)
        if style_m:
            style = style_m.group(1)
            if is_gene:
                colour = TRANSFER_GREEN if poly_counter in transfer_poly_indices else GENE_GREY
                style = format_style(style, fill=colour, stroke=colour, stroke_width="1")
            elif "100%" in style and "fill" in style:
                # Homology ribbons: keep them red but make the colour match the
                # compact reference more closely.
                style = format_style(style, fill=HOMOLOGY_RED, stroke=HOMOLOGY_RED, stroke_width="1")
            tag = tag[:style_m.start(1)] + style + tag[style_m.end(1):]
        return tag

    def repl_path(match: re.Match[str]) -> str:
        tag = match.group(0)
        d_m = re.search(r'd="([^"]+)"', tag)
        if not d_m:
            return tag
        d = d_m.group(1)
        # The original paths are simple M/L line segments with comma-separated
        # coordinates. Adjust every y coordinate.
        def repl_coord(cm: re.Match[str]) -> str:
            cmd = cm.group(1)
            x = float(cm.group(2))
            y = float(cm.group(3))
            return f"{cmd} {x:.6f},{transform_row_y(y):.6f}"

        d2 = re.sub(r"([ML])\s*([0-9.\-]+),([0-9.\-]+)", repl_coord, d)
        tag = tag[:d_m.start(1)] + d2 + tag[d_m.end(1):]
        style_m = re.search(r'style="([^"]+)"', tag)
        if style_m:
            style = format_style(style_m.group(1), stroke_width=AXIS_STROKE_WIDTH)
            tag = tag[:style_m.start(1)] + style + tag[style_m.end(1):]
        return tag

    svg_text = re.sub(r"<polygon\b[^>]*>", repl_polygon, svg_text)
    svg_text = re.sub(r"<path\b[^>]*>", repl_path, svg_text)
    return svg_text


def make_overlay(width: float, height: float, title_gene: str, score: float | None,
                 direction: str, top_gene: str, bottom_gene: str) -> str:
    # Coordinates tuned to the MetaCHIP flanking plots: two flanking rows are
    # around y=170-252 and y=456-538 in normal SVG coordinates.
    header_x = TITLE_X_PAD
    header_y = TITLE_Y
    header_w = width - TITLE_X_PAD * 2
    header_h = TITLE_H
    text_y = header_y + header_h * 0.66
    font_size = 36 if len(title_gene) < 10 else 32
    if score is None:
        score_text = "NA"
    else:
        score_text = f"{score:.1f}"
    escaped_gene = html.escape(title_gene)
    escaped_direction = direction
    plain_title = f"{title_gene} ({score_text}), {direction.replace('&gt;', '>')}"
    text_length = min(header_w - 160, max(520, font_size * 0.58 * len(plain_title)))
    text_x = width / 2 - text_length / 2
    title_text = (
        f'<tspan font-style="italic">{escaped_gene}</tspan>'
        f'<tspan> ({score_text}), {escaped_direction}</tspan>'
    )

    top_colour = habitat_colour(top_gene)
    bottom_colour = habitat_colour(bottom_gene)
    bar_x = 42
    bar_w = 38
    # Overlay coordinates are outside the flipped MetaCHIP group, so convert
    # transformed group coordinates back to visible SVG coordinates.
    top_group_min = transform_row_y(UPPER_ROW_MIN_Y)
    top_group_max = transform_row_y(UPPER_ROW_MAX_Y)
    bottom_group_min = transform_row_y(LOWER_ROW_MIN_Y)
    bottom_group_max = transform_row_y(LOWER_ROW_MAX_Y)
    top_y = display_y_from_group_y(top_group_max, height)
    bottom_y = display_y_from_group_y(bottom_group_max, height)
    bar_h = top_group_max - top_group_min

    return f"""
  <g id="codex_hgt_annotation_layer">
    <rect x="{header_x:.1f}" y="{header_y:.1f}" width="{header_w:.1f}" height="{header_h:.1f}"
          fill="{TITLE_FILL}" stroke="none"/>
    <text x="{text_x:.1f}" y="{text_y:.1f}" text-anchor="start"
          textLength="{text_length:.1f}" lengthAdjust="spacingAndGlyphs"
          font-family="Helvetica, Arial, sans-serif" font-size="{font_size}"
          font-weight="700" fill="{BLACK}">{title_text}</text>
    <rect x="{bar_x:.1f}" y="{top_y:.1f}" width="{bar_w:.1f}" height="{bar_h:.1f}"
          fill="{top_colour}" stroke="none"/>
    <rect x="{bar_x:.1f}" y="{bottom_y:.1f}" width="{bar_w:.1f}" height="{bar_h:.1f}"
          fill="{bottom_colour}" stroke="none"/>
  </g>
"""


def cropped_viewbox(width: float, height: float) -> tuple[float, float, float, float]:
    """Crop away the unused white canvas while keeping the SVG coordinate system."""
    bottom_group_min = transform_row_y(LOWER_ROW_MIN_Y)
    bottom_visible_y = display_y_from_group_y(bottom_group_min, height)
    crop_x = CROP_X_PAD
    crop_y = CROP_TOP_PAD
    crop_w = width - CROP_X_PAD * 2
    crop_h = bottom_visible_y + CROP_BOTTOM_PAD - crop_y
    return crop_x, crop_y, crop_w, crop_h


def crop_svg_canvas(svg_text: str, crop_x: float, crop_y: float, crop_w: float, crop_h: float) -> str:
    """Set the exported artboard to the actual figure area instead of the full MetaCHIP canvas."""

    def repl_svg(match: re.Match[str]) -> str:
        tag = match.group(0)
        attrs = {
            "width": f"{crop_w:.1f}",
            "height": f"{crop_h:.1f}",
            "viewBox": f"{crop_x:.1f} {crop_y:.1f} {crop_w:.1f} {crop_h:.1f}",
            "preserveAspectRatio": "xMinYMin meet",
        }
        for key, value in attrs.items():
            if re.search(rf'\b{key}="[^"]*"', tag):
                tag = re.sub(rf'\b{key}="[^"]*"', f'{key}="{value}"', tag, count=1)
            else:
                tag = tag[:-1] + f' {key}="{value}">'
        return tag

    return re.sub(r"<svg\b[^>]*>", repl_svg, svg_text, count=1)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary_rows = {row["flanking_plot"]: row for row in read_tsv(SUMMARY_TSV)}
    hgt_rows = read_tsv(HGT_TSV)

    edited_manifest: list[dict[str, str]] = []
    for flanking_plot, row in sorted(summary_rows.items()):
        src = SVG_DIR / flanking_plot
        if not src.exists():
            continue
        gene1 = row["gene1"]
        gene2 = row["gene2"]
        score, direction_raw = find_identity(hgt_rows, gene1, gene2)
        direction = direction_title(direction_raw, gene1, gene2)
        gene_label = compact_gene_label(row["final_genes"])

        svg = src.read_text(encoding="utf-8", errors="ignore")
        width, height = svg_size(svg)
        transfer_poly_indices, top_gene_for_colour, bottom_gene_for_colour = classify_gene_polygons(svg, gene1, gene2)
        edited = restyle_and_relayout_original_svg(svg, transfer_poly_indices)
        edited = remove_original_text(edited)
        overlay = make_overlay(
            width, height, gene_label, score, direction,
            top_gene_for_colour, bottom_gene_for_colour
        )
        if "</svg>" not in edited:
            raise RuntimeError(f"Could not find closing </svg> in {src}")
        edited = edited.replace("</svg>", overlay + "\n</svg>")
        edited = crop_svg_canvas(edited, *cropped_viewbox(width, height))

        out_name = src.stem + "_edited.SVG"
        out = OUT_DIR / out_name
        out.write_text(edited, encoding="utf-8")
        edited_manifest.append({
            "original_svg": str(src),
            "edited_svg": str(out),
            "flanking_plot": flanking_plot,
            "gene_label": gene_label,
            "identity_score": "" if score is None else f"{score:.3f}",
            "direction_raw": direction_raw,
            "direction_display": direction.replace("&gt;", ">"),
            "gene1": gene1,
            "gene2": gene2,
            "n_transfer_arrows_detected": str(len(transfer_poly_indices)),
            "top_row_colour_gene": top_gene_for_colour,
            "bottom_row_colour_gene": bottom_gene_for_colour,
            "top_bar_colour": habitat_colour(top_gene_for_colour),
            "bottom_bar_colour": habitat_colour(bottom_gene_for_colour),
            "gene1_habitat": "R" if habitat_colour(gene1) == R_COLOUR else "E",
            "gene2_habitat": "R" if habitat_colour(gene2) == R_COLOUR else "E",
        })

    manifest_path = OUT_DIR / "edited_hgt_flanking_svg_manifest.tsv"
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = list(edited_manifest[0].keys()) if edited_manifest else []
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(edited_manifest)
    print(f"Wrote {len(edited_manifest)} edited SVG files to {OUT_DIR}")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
