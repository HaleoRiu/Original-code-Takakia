#!/usr/bin/env python3
"""Build final Figure 5 anti-stress KO tables and final HGT subset.

Final user decisions implemented here:
  - 12 figure panels: total + 11 categories.
  - Keep the settled base categories from Figure5_AS_base_v4:
    Amino acids, Cold stress, Flavonoids, Iron carrier, Terpenoids,
    Trehalose, VB12.
  - Carotenoids: KEGG modules M00097, M01038, M01044, M01045.
  - DNA repair: 14 UV/ROS-relevant bacterial KOs only.
  - Oxidative stress: 23 conservative ROS detox / oxidative damage repair /
    iron-Fenton protection KOs.
  - Colonization: 24 Biofilm/QS-mediated colonization KOs.
  - HGT: keep only events whose annotation KOs are present in this final
    anti-stress KO table, copy all matched flanking plots, and recommend up to
    12 representative events for the main text.
"""

from __future__ import annotations

import re
import shutil
from pathlib import Path

import pandas as pd


BASE_DIR = Path("/Volumes/Haleo_DK01/WD_Takakia_F/metagenome-data")
FIG5_DIR = BASE_DIR / "Figure5"
TPM_FILE = BASE_DIR / "tpm_renormalized.txt"
BASE_V4_CATALOG = FIG5_DIR / "Figure5_AS_base_v4_KO_catalog.tsv"
LIST_FILE = BASE_DIR / "list.csv"
DONOR_FILE = BASE_DIR / "donor.tsv"
RECIPIENT_FILE = BASE_DIR / "recipient.tsv"
FLANKING_DIR = BASE_DIR / "flanking-region" / "all_g115_Flanking_region_plots"
KEGG_KO_LIST = Path("work/kegg_ko_list.tsv")

PREFIX = "Figure5_antistress_final"
HGT_OUT_DIR = FIG5_DIR / "HGT_flanking_plots_final_antistress"

CATEGORY_ORDER = [
    "All anti-stress genes",
    "Amino acids",
    "Carotenoids",
    "Cold stress",
    "Colonization",
    "DNA repair",
    "Flavonoids",
    "Iron carrier",
    "Oxidative stress",
    "Terpenoids",
    "Trehalose",
    "VB12",
]

BASE_CATEGORIES = [
    "Amino acids",
    "Cold stress",
    "Flavonoids",
    "Iron carrier",
    "Terpenoids",
    "Trehalose",
    "VB12",
]

CAROTENOID_MODULES = {
    "M00097": {
        "name": "Lycopene biosynthesis, geranylgeranyl-PP => lycopene",
        "kos": ["K02291", "K17841", "K02293", "K15744", "K00514", "K09835", "K10027", "K15745"],
    },
    "M01038": {
        "name": "Staphyloxanthin biosynthesis, bacteria, farnesyl-PP => staphyloxanthin",
        "kos": ["K10208", "K10209", "K10210", "K25072", "K10211", "K10212"],
    },
    "M01044": {
        "name": "Astaxanthin biosynthesis, beta-carotene => astaxanthin",
        "kos": ["K09836", "K15746", "K23037"],
    },
    "M01045": {
        "name": "Zeaxanthin biosynthesis, lycopene => beta-carotene => zeaxanthin",
        "kos": ["K06443", "K22502", "K14605", "K14606", "K02294", "K15747", "K15746", "K22492"],
    },
}

DNA_REPAIR_14 = {
    "K01669": ("Photoreactivation/direct UV repair", "photolyase direct repair of UV-induced CPD lesions"),
    "K06876": ("Photoreactivation/direct UV repair", "6-4 photolyase direct repair of UV photoproducts"),
    "K03701": ("NER/UvrABC", "UvrA nucleotide-excision repair damage recognition"),
    "K03702": ("NER/UvrABC", "UvrB nucleotide-excision repair helicase/damage verification"),
    "K03703": ("NER/UvrABC", "UvrC nucleotide-excision repair incision"),
    "K03657": ("NER/UvrABC", "UvrD/PcrA repair helicase for NER/BER/MMR turnover"),
    "K03723": ("Transcription-coupled NER", "Mfd transcription-repair coupling factor"),
    "K10563": ("Oxidative BER", "MutM/Fpg oxidative-base DNA glycosylase"),
    "K03575": ("Oxidative BER", "MutY A/G-specific adenine glycosylase"),
    "K03574": ("Oxidative BER", "MutT 8-oxo-dGTP sanitizing hydrolase"),
    "K05522": ("Oxidative BER", "Nei/endonuclease VIII oxidized-base repair"),
    "K10773": ("Oxidative BER", "Nth/endonuclease III oxidized pyrimidine repair"),
    "K01142": ("BER AP-endonuclease", "XthA/exodeoxyribonuclease III AP-endonuclease"),
    "K01151": ("BER AP-endonuclease", "Nfo/deoxyribonuclease IV AP-endonuclease"),
}

OXIDATIVE_STRESS_23 = {
    "K00428": ("ROS/peroxide detox", "cytochrome c peroxidase"),
    "K00432": ("ROS/peroxide detox", "glutathione peroxidase"),
    "K00518": ("ROS detox", "nickel superoxide dismutase"),
    "K03386": ("ROS/peroxide detox", "AhpC/peroxiredoxin alkyl-hydroperoxide reductase"),
    "K03387": ("ROS/peroxide detox", "AhpF reductase partner for AhpC"),
    "K03564": ("ROS/peroxide detox", "thioredoxin-dependent peroxiredoxin"),
    "K03594": ("Iron/Fenton protection", "bacterioferritin limits free-iron driven Fenton chemistry"),
    "K03781": ("ROS/peroxide detox", "catalase"),
    "K03782": ("ROS/peroxide detox", "catalase-peroxidase"),
    "K04063": ("ROS/peroxide detox", "OsmC/Ohr organic-hydroperoxide peroxiredoxin"),
    "K04564": ("ROS detox", "Fe/Mn superoxide dismutase"),
    "K04565": ("ROS detox", "Cu/Zn superoxide dismutase"),
    "K04756": ("ROS/peroxide detox", "AhpD lipoyl-dependent peroxiredoxin subunit"),
    "K05919": ("ROS detox", "superoxide reductase"),
    "K07147": ("Oxidative protein-damage repair", "MsrP methionine sulfoxide reductase catalytic subunit"),
    "K07217": ("ROS/peroxide detox", "manganese catalase"),
    "K07304": ("Oxidative protein-damage repair", "MsrA peptide-methionine S-oxide reductase"),
    "K07305": ("Oxidative protein-damage repair", "MsrB peptide-methionine R-oxide reductase"),
    "K08968": ("Oxidative protein-damage repair", "MsrC L-methionine R-oxide reductase"),
    "K11065": ("ROS/peroxide detox", "thiol peroxidase / atypical peroxiredoxin"),
    "K12267": ("Oxidative protein-damage repair", "bifunctional MsrAB"),
    "K17247": ("Oxidative protein-damage repair", "MsrQ heme-binding partner of MsrP/Q system"),
    "K22336": ("Iron/Fenton protection", "bacterioferritin B limits free-iron driven Fenton chemistry"),
}

COLONIZATION_24 = {
    "K07173": ("AI-2 quorum sensing", "LuxS autoinducer-2 synthesis"),
    "K11216": ("AI-2 quorum sensing", "LsrK autoinducer-2 kinase"),
    "K10909": ("AI-2 quorum sensing", "LuxQ autoinducer sensor kinase/phosphatase"),
    "K10912": ("AI-2 quorum sensing", "LuxO quorum-sensing response regulator"),
    "K13060": ("AHL quorum sensing", "LasI acyl-homoserine lactone synthase"),
    "K13061": ("AHL quorum sensing", "RhlI/PhzI/SolI/CepI/TofI acyl-homoserine lactone synthase"),
    "K18096": ("AHL quorum sensing", "BjaI/RpaI/BraI/RhiI AHL synthase"),
    "K20250": ("AHL quorum sensing", "CciI AHL synthase"),
    "K07782": ("LuxR-family QS regulator", "SdiA LuxR-family quorum-sensing regulator"),
    "K18098": ("LuxR-family QS regulator", "BjaR1/RpaR/RhiR LuxR-family quorum-sensing regulator"),
    "K18099": ("LuxR-family QS regulator", "RhlR/PhzR LuxR-family quorum-sensing regulator"),
    "K18304": ("LuxR-family QS regulator", "LasR LuxR-family quorum-sensing regulator"),
    "K19666": ("LuxR-family QS regulator", "SolR/CepR/TofR LuxR-family quorum-sensing regulator"),
    "K19731": ("LuxR-family QS regulator", "CciR quorum-sensing regulator"),
    "K19734": ("LuxR-family QS regulator", "ExpR quorum-sensing regulator"),
    "K11935": ("Biofilm/EPS", "PgaA biofilm PGA synthesis protein"),
    "K11937": ("Biofilm/EPS", "PgaD biofilm PGA synthesis protein"),
    "K16566": ("Biofilm/EPS", "ExoY exopolysaccharide production protein"),
    "K16568": ("Biofilm/EPS", "ExoZ exopolysaccharide production protein"),
    "K19421": ("Biofilm/EPS", "EpsC polysaccharide biosynthesis protein"),
    "K19422": ("Biofilm/EPS", "EpsD exopolysaccharide glycosyltransferase"),
    "K20918": ("Biofilm/EPS", "VpsT positive regulator of biofilm formation"),
    "K20952": ("Biofilm/EPS", "RbmC biofilm structure modulator"),
    "K21154": ("Biofilm/EPS", "PssM exopolysaccharide transferase"),
}

HGT_CATEGORY_PRIORITY = {
    "Cold stress": 0,
    "Oxidative stress": 1,
    "DNA repair": 2,
    "Trehalose": 3,
    "Carotenoids": 4,
    "Iron carrier": 5,
    "Terpenoids": 6,
    "VB12": 7,
    "Amino acids": 8,
    "Colonization": 9,
    "Flavonoids": 10,
}

REPRESENTATIVE_GENE_PRIORITY = [
    "cspA",
    "katG", "katE", "SOD2", "SOD1", "sodN", "gpx", "ahpC", "msrA", "msrB", "bfr",
    "phr", "phrB", "uvrA", "uvrB", "uvrC", "uvrD", "mutM", "mutY", "mutT", "nei", "nth",
    "otsA", "otsB",
    "crtB", "crtI", "crtY", "crtZ",
    "entA", "entB", "entC", "entF", "fepA",
    "luxS", "lasI", "rhlI", "lasR", "rhlR", "pgaA", "epsC",
]


def parse_kegg_ko_list(path: Path) -> pd.DataFrame:
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip() or "\t" not in line:
            continue
        kid, raw = line.split("\t", 1)
        k = kid.replace("ko:", "").strip()
        if "; " in raw:
            aliases, definition = raw.split("; ", 1)
        elif ";" in raw:
            aliases, definition = raw.split(";", 1)
            definition = definition.strip()
        else:
            aliases, definition = raw, ""
        rows.append({"K": k, "KEGG_aliases": aliases, "KEGG_definition": definition, "KEGG_raw": raw})
    return pd.DataFrame(rows)


def compact_gene_name(k: str, aliases: str, definition: str = "") -> str:
    parts = [x.strip() for x in str(aliases).split(",") if x.strip()]
    useful = []
    for x in parts:
        if x == k or re.match(r"^K\d{5}$", x) or re.match(r"^E\d", x):
            continue
        useful.append(x)
    if useful:
        return useful[0]
    m = re.match(r"([A-Za-z0-9_./-]+);", str(definition))
    return m.group(1) if m else k


def sample_metadata(sample_cols: list[str]) -> pd.DataFrame:
    rows = []
    for sample in sample_cols:
        m = re.match(r"^([SP])(\d{2})_(\d+)$", sample)
        if not m:
            continue
        prefix, elev2, rep = m.groups()
        rows.append(
            {
                "sample": sample,
                "Group": "R" if prefix == "S" else "E",
                "Habitat": "Rhizomatousphere" if prefix == "S" else "Endophytic",
                "Elevation": int(elev2) * 100,
                "Replicate": int(rep),
                "PairID": f"{elev2}_{rep}",
            }
        )
    return pd.DataFrame(rows)


def split_kos(val) -> list[str]:
    if pd.isna(val):
        return []
    return [x.strip().replace("ko:", "") for x in re.split(r"[,;]\s*", str(val)) if re.match(r"^ko:K\d{5}$|^K\d{5}$", x.strip())]


def parse_hgt_events() -> pd.DataFrame:
    events = []
    for line in LIST_FILE.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if not re.search(r"\.SVG$", line, flags=re.I):
            line = f"{line}.SVG"
        stem = re.sub(r"\.SVG$", "", line, flags=re.I)
        parts = stem.split("___")
        if len(parts) == 2:
            events.append({"flanking_plot": line, "gene1": parts[0], "gene2": parts[1]})
    return pd.DataFrame(events)


def gene_priority(gene: str) -> int:
    gene_lower = str(gene).split(",")[0].strip()
    for i, candidate in enumerate(REPRESENTATIVE_GENE_PRIORITY):
        if gene_lower.lower() == candidate.lower():
            return i
    return 1000


def main() -> None:
    FIG5_DIR.mkdir(parents=True, exist_ok=True)

    ko_info = parse_kegg_ko_list(KEGG_KO_LIST).set_index("K")
    tpm = pd.read_csv(TPM_FILE, sep="\t").rename(columns={"feature-id": "K"})
    tpm["K"] = tpm["K"].astype(str)
    sample_cols = [c for c in tpm.columns if re.match(r"^[SP]\d{2}_\d+$", c)]
    tpm_present = set(tpm["K"])

    base_catalog = pd.read_csv(BASE_V4_CATALOG, sep="\t")
    base_catalog = base_catalog[base_catalog["category"].isin(BASE_CATEGORIES)].copy()
    records: dict[tuple[str, str], dict[str, str]] = {}

    def add_record(category: str, k: str, evidence: str, reason: str, gene_override: str | None = None) -> None:
        if k not in tpm_present:
            return
        if k in ko_info.index:
            aliases = str(ko_info.loc[k, "KEGG_aliases"])
            definition = str(ko_info.loc[k, "KEGG_definition"])
        else:
            aliases = ""
            definition = ""
        gene = gene_override or compact_gene_name(k, aliases, definition)
        records[(category, k)] = {
            "category": category,
            "gene": gene,
            "K": k,
            "KEGG_aliases": aliases,
            "KEGG_definition": definition,
            "KEGG_evidence": evidence,
            "inclusion_reason": reason,
            "final_selection_note": "retained in final Figure5 anti-stress KO set",
        }

    for _, row in base_catalog.iterrows():
        # K07173 is LuxS: biochemically it appears in sulfur/amino-acid
        # metabolism, but in this final Figure 5 curation it is interpreted as
        # AI-2 quorum-sensing/colonization. Keep categories mutually cleaner by
        # assigning it only to Colonization.
        if str(row["category"]) == "Amino acids" and str(row["K"]) == "K07173":
            continue
        add_record(
            str(row["category"]),
            str(row["K"]),
            str(row.get("KEGG_evidence", "")),
            str(row.get("inclusion_reason", "")),
            gene_override=str(row.get("gene", "")) or None,
        )

    for mid, meta in CAROTENOID_MODULES.items():
        for k in meta["kos"]:
            add_record(
                "Carotenoids",
                k,
                f"{mid}: {meta['name']}",
                "user-selected bacterial/prokaryote-compatible carotenoid module from ko00906",
            )

    for k, (branch, reason) in DNA_REPAIR_14.items():
        add_record("DNA repair", k, branch, reason)

    for k, (branch, reason) in OXIDATIVE_STRESS_23.items():
        add_record("Oxidative stress", k, branch, reason)

    for k, (branch, reason) in COLONIZATION_24.items():
        add_record("Colonization", k, branch, reason)

    catalog = pd.DataFrame(records.values())
    category_rank = {c: i for i, c in enumerate(CATEGORY_ORDER)}
    catalog["_category_order"] = catalog["category"].map(category_rank)
    catalog = catalog.sort_values(["_category_order", "K"]).drop(columns="_category_order").reset_index(drop=True)
    duplicated = catalog.duplicated(["category", "gene"], keep=False)
    catalog.loc[duplicated, "gene"] = catalog.loc[duplicated].apply(lambda r: f"{r['gene']} ({r['K']})", axis=1)

    ko_tpm = catalog.merge(tpm, on="K", how="left")
    parts = []
    for category in CATEGORY_ORDER[1:]:
        sub = ko_tpm[ko_tpm["category"].eq(category)].drop_duplicates("K")
        summed = sub.set_index("K")[sample_cols].sum(axis=0) if len(sub) else pd.Series(0, index=sample_cols)
        parts.append(pd.DataFrame({"sample": sample_cols, "category": category, "TPM": [float(summed[s]) for s in sample_cols], "n_KO": sub["K"].nunique()}))

    total = ko_tpm.drop_duplicates("K").set_index("K")[sample_cols].sum(axis=0)
    parts.insert(0, pd.DataFrame({"sample": sample_cols, "category": "All anti-stress genes", "TPM": [float(total[s]) for s in sample_cols], "n_KO": catalog["K"].nunique()}))

    abundance = pd.concat(parts, ignore_index=True).merge(sample_metadata(sample_cols), on="sample", how="left")
    abundance["category"] = pd.Categorical(abundance["category"], CATEGORY_ORDER, ordered=True)
    abundance = abundance.sort_values(["category", "Group", "Elevation", "Replicate"])

    summary = (
        catalog.groupby("category", as_index=False)
        .agg(
            n_KO=("K", "nunique"),
            K_numbers=("K", lambda x: "; ".join(sorted(set(x)))),
            genes=("gene", lambda x: "; ".join(sorted(set(x)))),
            evidence=("KEGG_evidence", lambda x: " | ".join(sorted(set(x)))),
        )
    )
    total_summary = pd.DataFrame(
        [{
            "category": "All anti-stress genes",
            "n_KO": catalog["K"].nunique(),
            "K_numbers": "; ".join(sorted(set(catalog["K"]))),
            "genes": "; ".join(sorted(set(catalog["gene"]))),
            "evidence": "sum of all unique final anti-stress KOs",
        }]
    )
    summary = pd.concat([total_summary, pd.DataFrame({"category": CATEGORY_ORDER[1:]}).merge(summary, on="category", how="left")], ignore_index=True)

    catalog.to_csv(FIG5_DIR / f"{PREFIX}_KO_catalog.tsv", sep="\t", index=False)
    ko_tpm.to_csv(FIG5_DIR / f"{PREFIX}_KO_TPM_table.tsv", sep="\t", index=False)
    abundance.to_csv(FIG5_DIR / f"{PREFIX}_category_TPM_by_sample.tsv", sep="\t", index=False)
    summary.to_csv(FIG5_DIR / f"{PREFIX}_category_summary.tsv", sep="\t", index=False)

    # HGT screening against final KO set.
    final_kos = set(catalog["K"])
    final_kos_by_k = (
        catalog.groupby("K", as_index=False)
        .agg(
            category=("category", lambda x: " | ".join(sorted(set(map(str, x))))),
            gene=("gene", lambda x: "; ".join(sorted(set(map(str, x))))),
            KEGG_definition=("KEGG_definition", lambda x: " | ".join(sorted(set(map(str, x))))),
            KEGG_evidence=("KEGG_evidence", lambda x: " | ".join(sorted(set(map(str, x))))),
            inclusion_reason=("inclusion_reason", lambda x: " | ".join(sorted(set(map(str, x))))),
        )
        .set_index("K")
        .to_dict("index")
    )
    events = parse_hgt_events()
    ann = []
    for role, fn in [("donor", DONOR_FILE), ("recipient", RECIPIENT_FILE)]:
        df = pd.read_csv(fn, sep="\t")
        df = df.rename(columns={"#query": "query"})
        df["annotation_file"] = role
        ann.append(df)
    ann = pd.concat(ann, ignore_index=True)

    long = []
    for side in ["gene1", "gene2"]:
        tmp = events[["flanking_plot", "gene1", "gene2"]].copy()
        tmp["query"] = tmp[side]
        tmp = tmp.merge(ann, on="query", how="left")
        tmp["event_side"] = side
        long.append(tmp)
    long = pd.concat(long, ignore_index=True)
    long["K"] = long.get("KEGG_ko", pd.Series(index=long.index, dtype=object)).apply(split_kos)
    long = long.explode("K")
    hgt_hits = long[long["K"].isin(final_kos)].copy()
    for col in ["category", "gene", "KEGG_definition", "KEGG_evidence", "inclusion_reason"]:
        hgt_hits[f"final_{col}"] = hgt_hits["K"].map(lambda x: final_kos_by_k.get(x, {}).get(col, ""))

    hgt_cols = [
        "flanking_plot", "gene1", "gene2", "event_side", "query", "annotation_file",
        "K", "final_category", "final_gene", "final_KEGG_definition",
        "Preferred_name", "Description", "KEGG_ko", "KEGG_Pathway", "PFAMs",
        "final_KEGG_evidence", "final_inclusion_reason",
    ]
    hgt_hits = hgt_hits[[c for c in hgt_cols if c in hgt_hits.columns]].drop_duplicates()

    event_summary = (
        hgt_hits.groupby(["flanking_plot", "gene1", "gene2"], dropna=False)
        .agg(
            final_categories=("final_category", lambda x: " | ".join(sorted(set(map(str, x))))),
            final_K_numbers=("K", lambda x: "; ".join(sorted(set(map(str, x))))),
            final_genes=("final_gene", lambda x: "; ".join(sorted(set(map(str, x))))),
            KEGG_definitions=("final_KEGG_definition", lambda x: " | ".join(sorted(set(map(str, x))))),
            annotated_genes=("query", lambda x: "; ".join(sorted(set(str(v) for v in x if pd.notna(v))))),
            annotation_files=("annotation_file", lambda x: "; ".join(sorted(set(str(v) for v in x if pd.notna(v))))),
        )
        .reset_index()
    )

    # Recommend up to 12 events: first preserve category diversity, then use
    # mechanistic specificity and non-redundant genes.
    event_summary["primary_category"] = event_summary["final_categories"].str.split(" \\| ").str[0]
    event_summary["category_rank"] = event_summary["primary_category"].map(HGT_CATEGORY_PRIORITY).fillna(99).astype(int)
    event_summary["gene_rank"] = event_summary["final_genes"].apply(lambda s: min(gene_priority(g.strip()) for gene in str(s).split(";") for g in gene.split(",")))
    event_summary["n_categories"] = event_summary["final_categories"].apply(lambda s: len(set(str(s).split(" | "))))
    event_summary["n_KO"] = event_summary["final_K_numbers"].apply(lambda s: len(set(str(s).split("; "))))

    sorted_events = event_summary.sort_values(["category_rank", "gene_rank", "flanking_plot"]).copy()
    selected_indices = []
    used_categories = set()
    used_genes = set()

    for idx, row in sorted_events.iterrows():
        cats = set(str(row["final_categories"]).split(" | "))
        genes = set(g.strip() for g in str(row["final_genes"]).split(";") if g.strip())
        if cats - used_categories:
            selected_indices.append(idx)
            used_categories.update(cats)
            used_genes.update(genes)
        if len(selected_indices) >= 12:
            break

    for idx, row in sorted_events.iterrows():
        if len(selected_indices) >= 12:
            break
        genes = set(g.strip() for g in str(row["final_genes"]).split(";") if g.strip())
        if idx not in selected_indices and genes - used_genes:
            selected_indices.append(idx)
            used_genes.update(genes)

    for idx, row in sorted_events.iterrows():
        if len(selected_indices) >= 12:
            break
        if idx not in selected_indices:
            selected_indices.append(idx)

    recommended = event_summary.loc[selected_indices].copy()
    recommended["recommendation_rank"] = range(1, len(recommended) + 1)
    recommended["recommendation_reason"] = recommended.apply(
        lambda r: (
            "prioritized for main-text visualization: representative final anti-stress category, "
            "explicit KO/gene annotation, and reduced redundancy across genes/categories"
        ),
        axis=1,
    )

    hgt_hits.to_csv(FIG5_DIR / f"{PREFIX}_HGT_gene_hit_overlap.tsv", sep="\t", index=False)
    event_summary.drop(columns=["category_rank", "gene_rank"], errors="ignore").to_csv(
        FIG5_DIR / f"{PREFIX}_HGT_event_overlap_summary.tsv", sep="\t", index=False
    )
    recommended.drop(columns=["category_rank", "gene_rank"], errors="ignore").to_csv(
        FIG5_DIR / f"{PREFIX}_HGT_main_text_recommended_12.tsv", sep="\t", index=False
    )

    if HGT_OUT_DIR.exists():
        for old in HGT_OUT_DIR.iterdir():
            if old.is_file() and (old.suffix.lower() in {".svg", ".tsv", ".csv"}):
                old.unlink()
    HGT_OUT_DIR.mkdir(parents=True, exist_ok=True)
    copied = []
    missing = []
    for plot in sorted(event_summary["flanking_plot"].unique()):
        src = FLANKING_DIR / plot
        if not src.exists():
            src = FLANKING_DIR / re.sub(r"\.SVG$", ".svg", plot, flags=re.I)
        if src.exists():
            dst = HGT_OUT_DIR / src.name
            shutil.copy2(src, dst)
            copied.append(dst.name)
        else:
            missing.append(plot)
    pd.DataFrame({"flanking_plot": copied}).to_csv(HGT_OUT_DIR / "manifest_all_final_antistress_hgt_flanking_plots.tsv", sep="\t", index=False)
    pd.DataFrame({"missing_flanking_plot": missing}).to_csv(HGT_OUT_DIR / "missing_final_antistress_hgt_flanking_plots.tsv", sep="\t", index=False)
    recommended[["flanking_plot", "recommendation_rank", "final_categories", "final_genes", "final_K_numbers"]].to_csv(
        HGT_OUT_DIR / "manifest_main_text_recommended_12.tsv",
        sep="\t",
        index=False,
    )

    print("Final KO catalog:", len(catalog), "KOs")
    print(summary[["category", "n_KO"]].to_string(index=False))
    print("Final HGT events:", len(event_summary))
    print("Copied flanking plots:", len(copied), "missing:", len(missing))
    print("Recommended main-text HGT events:", len(recommended))
    print(recommended[["recommendation_rank", "flanking_plot", "final_categories", "final_genes", "final_K_numbers"]].to_string(index=False))


if __name__ == "__main__":
    main()
