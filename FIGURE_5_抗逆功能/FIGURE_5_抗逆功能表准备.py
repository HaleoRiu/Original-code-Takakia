#!/usr/bin/env python3
"""Rebuild Figure 5 anti-stress KO tables after GSH/amino-acid and trehalose edits.

This script deliberately starts from the previously settled final KO catalog and
changes only the parts requested in the latest revision:

1. "Amino acids" is replaced by a glutathione-related amino-acid category:
   - remove Pro/proline related KOs from the old final amino-acid set;
   - keep Glu/Gly/Cys and sulfur-amino-acid precursor KOs from the old set;
   - add glutathione synthesis/recycling/GST detoxification KOs requested by the
     user, including GST/fosA-family GST-like KOs.
2. Trehalose is rebuilt in three versions:
   - main/original: the previously used strict M00565 trehalose module;
   - full_reference: central trehalose synthesis, interconversion, degradation,
     and bacterial transport KOs from KEGG ko00500.
   - merged: strict M00565 plus the full_reference KO set; this is the
     preferred revised "Trehalose metabolism" category because it keeps the
     original module and adds the new synthesis-metabolism KOs without deleting
     any previous KO.
   - bacterial_conservative: merged set minus very-low-abundance trehalose KOs
     and the eukaryote-biased/non-canonical bacterial bifunctional TPS KO.

The output tables are KO-level TPM summaries compatible with the existing
Figure 5 plotting style.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd


BASE_DIR = Path("/Volumes/Haleo_DK01/WD_Takakia_F/metagenome-data")
FIG5_DIR = BASE_DIR / "Figure5"
TPM_FILE = BASE_DIR / "tpm_renormalized.txt"
FINAL_CATALOG = FIG5_DIR / "Figure5_antistress_final_KO_catalog.tsv"
KEGG_KO_LIST = Path("/Users/haleoriu/Documents/Codex/2026-06-21/zai/work/kegg_ko_list.tsv")

AMINO_NEW_CATEGORY = "Glutathione and related amino acids"

PROLINE_REMOVE = {"K00147", "K00286", "K00931"}
GSH_ADDITION_KOS = [
    "K01919",  # gshA
    "K01920",  # gshB
    "K06048",  # gshA/ybdK
    "K28983",  # gshAB; absent from current TPM table, kept in audit
    "K00383",  # gor/GSR
    "K00799",  # GST
    "K11208",  # GST-like
    "K21253",  # fosA
    "K21264",  # fosA2
    "K21265",  # fosA5
]

# Central trehalose synthesis/interconversion/degradation/transport KOs from
# ko00500. Mycobacterial trehalose-mycolate/sulfation/acyltransferase branches
# are intentionally excluded because they describe cell-wall lipid decoration
# rather than the general trehalose synthesis-metabolism route requested here.
TREHALOSE_FULL_KOS = [
    "K00697",  # otsA
    "K01087",  # otsB
    "K01194",  # treA/treF
    "K01226",  # treC
    "K01236",  # treZ/glgZ
    "K02819",  # treB/treP PTS
    "K03731",  # trehalose-6P phosphorylase; absent from current TPM
    "K05342",  # treP alpha,alpha-trehalose phosphorylase
    "K05343",  # treS
    "K06044",  # treY/glgY
    "K10236",  # thuE/lpqY
    "K10237",  # thuF/sugA
    "K10238",  # thuG/sugB
    "K13057",  # treT
    "K16055",  # bifunctional TPS
    "K17311",  # trehalose ABC substrate-binding; absent from current TPM
    "K17312",  # trehalose ABC permease; absent
    "K17313",  # trehalose ABC permease; absent
    "K17314",  # trehalose ABC ATP-binding; absent
    "K22248",  # retaining trehalose phosphorylase; absent
    "K22934",  # generic alpha,alpha-trehalase; absent
]

TREHALOSE_BACTERIAL_CONSERVATIVE_EXCLUDE = {
    "K02819",  # trehalose PTS; present in only 1/60 samples in this dataset
    "K01226",  # trehalose-6-phosphate hydrolase; very low abundance here
    "K16055",  # bifunctional TPS/TPP; eukaryote-heavy/non-canonical for bacterial MAG interpretation
}

CATEGORY_ORDER = [
    "All anti-stress genes",
    "DNA repair",
    AMINO_NEW_CATEGORY,
    "Cold stress",
    "Flavonoids",
    "Terpenoids",
    "VB12",
    "Carotenoids",
    "Iron carrier",
    "Colonization",
    "Oxidative stress",
    "Trehalose",
]


def parse_kegg_ko_list(path: Path) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
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
        rows.append(
            {
                "K": k,
                "KEGG_aliases": aliases,
                "KEGG_definition": definition,
                "KEGG_raw": raw,
            }
        )
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


def make_catalog_row(
    k: str,
    category: str,
    kegg: pd.DataFrame,
    evidence: str,
    reason: str,
    note: str,
) -> dict[str, str]:
    hit = kegg[kegg["K"] == k]
    if hit.empty:
        aliases = ""
        definition = ""
        gene = k
    else:
        aliases = hit.iloc[0]["KEGG_aliases"]
        definition = hit.iloc[0]["KEGG_definition"]
        gene = compact_gene_name(k, aliases, definition)
    return {
        "category": category,
        "gene": gene,
        "K": k,
        "KEGG_aliases": aliases,
        "KEGG_definition": definition,
        "KEGG_evidence": evidence,
        "inclusion_reason": reason,
        "final_selection_note": note,
    }


def build_catalog(trehalose_mode: str, present_kos: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    original = pd.read_csv(FINAL_CATALOG, sep="\t")
    kegg = parse_kegg_ko_list(KEGG_KO_LIST)

    kept = original[~original["category"].isin(["Amino acids", "Trehalose"])].copy()

    amino_old = original[original["category"] == "Amino acids"].copy()
    amino_old = amino_old[~amino_old["K"].isin(PROLINE_REMOVE)].copy()
    amino_old["category"] = AMINO_NEW_CATEGORY
    amino_old["inclusion_reason"] = amino_old["inclusion_reason"].map(
        lambda x: f"{x}; retained after removing Pro-related KOs"
    )
    amino_old["final_selection_note"] = (
        "Glu/Gly/Cys and sulfur-amino-acid precursor KO retained; Pro-related KOs removed"
    )

    added_rows = [
        make_catalog_row(
            k,
            AMINO_NEW_CATEGORY,
            kegg,
            "KEGG ko00480 / user-requested GSH-GST supplement",
            "glutathione synthesis, recycling, or glutathione S-transferase detoxification",
            "added in GSH/amino-acid revision; absent KOs are recorded in audit but not summed",
        )
        for k in GSH_ADDITION_KOS
    ]
    amino_add = pd.DataFrame(added_rows)
    amino_all = pd.concat([amino_old, amino_add], ignore_index=True)
    amino_all = amino_all[amino_all["K"].isin(present_kos)].drop_duplicates("K", keep="first")

    if trehalose_mode == "original_m00565":
        tre = original[original["category"] == "Trehalose"].copy()
        tre["final_selection_note"] = "strict M00565 trehalose module retained for main comparison"
    elif trehalose_mode == "full_reference":
        tre_rows = [
            make_catalog_row(
                k,
                "Trehalose",
                kegg,
                "KEGG ko00500 central trehalose synthesis/metabolism/transport",
                "central trehalose synthesis, interconversion, degradation, or bacterial transport",
                "full trehalose synthesis-metabolism reference route",
            )
            for k in TREHALOSE_FULL_KOS
        ]
        tre = pd.DataFrame(tre_rows)
        tre = tre[tre["K"].isin(present_kos)].drop_duplicates("K", keep="first")
    elif trehalose_mode == "merged_m00565_plus_full":
        tre_original = original[original["category"] == "Trehalose"].copy()
        tre_original["final_selection_note"] = (
            "strict M00565 module member retained in merged trehalose metabolism category"
        )
        tre_rows = [
            make_catalog_row(
                k,
                "Trehalose",
                kegg,
                "KEGG ko00500 central trehalose synthesis/metabolism/transport",
                "central trehalose synthesis, interconversion, degradation, or bacterial transport",
                "added to M00565 in merged trehalose metabolism category",
            )
            for k in TREHALOSE_FULL_KOS
        ]
        tre_full = pd.DataFrame(tre_rows)
        tre = pd.concat([tre_original, tre_full], ignore_index=True)
        tre = tre[tre["K"].isin(present_kos)].drop_duplicates("K", keep="first")
    elif trehalose_mode == "bacterial_conservative":
        tre_original = original[original["category"] == "Trehalose"].copy()
        tre_original["final_selection_note"] = (
            "strict M00565 module member retained in bacteria-focused trehalose metabolism category"
        )
        tre_rows = [
            make_catalog_row(
                k,
                "Trehalose",
                kegg,
                "KEGG ko00500 central trehalose synthesis/metabolism/transport",
                "central bacterial-compatible trehalose synthesis, interconversion, degradation, or transport",
                "added to M00565 in bacteria-focused trehalose metabolism category",
            )
            for k in TREHALOSE_FULL_KOS
        ]
        tre_full = pd.DataFrame(tre_rows)
        tre = pd.concat([tre_original, tre_full], ignore_index=True)
        tre = tre[~tre["K"].isin(TREHALOSE_BACTERIAL_CONSERVATIVE_EXCLUDE)]
        tre = tre[tre["K"].isin(present_kos)].drop_duplicates("K", keep="first")
    else:
        raise ValueError(f"unknown trehalose_mode: {trehalose_mode}")

    catalog = pd.concat([kept, amino_all, tre], ignore_index=True)
    catalog = catalog[catalog["K"].isin(present_kos)].drop_duplicates(["category", "K"], keep="first")
    catalog["category"] = pd.Categorical(catalog["category"], categories=CATEGORY_ORDER[1:], ordered=True)
    catalog = catalog.sort_values(["category", "K"]).reset_index(drop=True)

    audit_rows = []
    for kind, category, kos in [
        ("GSH/amino-acid candidate", AMINO_NEW_CATEGORY, GSH_ADDITION_KOS),
        ("Trehalose full candidate", "Trehalose", TREHALOSE_FULL_KOS),
        (
            "Trehalose merged candidate",
            "Trehalose",
            sorted(set(original[original["category"] == "Trehalose"]["K"].tolist()) | set(TREHALOSE_FULL_KOS)),
        ),
        (
            "Trehalose bacterial-conservative excluded",
            "Trehalose",
            sorted(TREHALOSE_BACTERIAL_CONSERVATIVE_EXCLUDE),
        ),
        ("Removed Pro-related KO", "removed_from_amino_acids", sorted(PROLINE_REMOVE)),
    ]:
        for k in kos:
            hit = kegg[kegg["K"] == k]
            audit_rows.append(
                {
                    "candidate_set": kind,
                    "category": category,
                    "K": k,
                    "present_in_TPM": k in present_kos,
                    "KEGG_aliases": "" if hit.empty else hit.iloc[0]["KEGG_aliases"],
                    "KEGG_definition": "" if hit.empty else hit.iloc[0]["KEGG_definition"],
                }
            )
    audit = pd.DataFrame(audit_rows)
    return catalog, audit


def write_abundance_tables(prefix: str, catalog: pd.DataFrame) -> None:
    tpm = pd.read_csv(TPM_FILE, sep="\t")
    sample_cols = [c for c in tpm.columns if re.match(r"^[SP]\d{2}_\d+$", c)]
    meta = sample_metadata(sample_cols)

    ko_tpm = (
        catalog.merge(tpm, left_on="K", right_on="feature-id", how="left")
        .drop(columns=["feature-id"])
        .sort_values(["category", "K"])
    )
    ko_tpm.to_csv(FIG5_DIR / f"{prefix}_KO_TPM_table.tsv", sep="\t", index=False)
    catalog.to_csv(FIG5_DIR / f"{prefix}_KO_catalog.tsv", sep="\t", index=False)

    long_rows = []
    for category, sub in ko_tpm.groupby("category", observed=False):
        if pd.isna(category) or len(sub) == 0:
            continue
        sums = sub[sample_cols].sum(axis=0)
        for sample, value in sums.items():
            long_rows.append(
                {
                    "sample": sample,
                    "category": str(category),
                    "TPM": float(value),
                    "n_KO": int(sub["K"].nunique()),
                }
            )

    unique_kos = catalog["K"].drop_duplicates().tolist()
    all_sums = tpm[tpm["feature-id"].isin(unique_kos)][sample_cols].sum(axis=0)
    for sample, value in all_sums.items():
        long_rows.append(
            {
                "sample": sample,
                "category": "All anti-stress genes",
                "TPM": float(value),
                "n_KO": int(len(unique_kos)),
            }
        )

    cat_long = pd.DataFrame(long_rows).merge(meta, on="sample", how="left")
    cat_long["category"] = pd.Categorical(cat_long["category"], categories=CATEGORY_ORDER, ordered=True)
    cat_long = cat_long.sort_values(["category", "Group", "Elevation", "Replicate"])
    cat_long.to_csv(FIG5_DIR / f"{prefix}_category_TPM_by_sample.tsv", sep="\t", index=False)

    summary_rows = [
        {
            "category": "All anti-stress genes",
            "n_KO": int(len(unique_kos)),
            "KOs": "; ".join(sorted(unique_kos)),
            "genes": "; ".join(sorted(catalog["gene"].dropna().astype(str).unique())),
        }
    ]
    for category, sub in catalog.groupby("category", observed=False):
        if pd.isna(category) or len(sub) == 0:
            continue
        summary_rows.append(
            {
                "category": str(category),
                "n_KO": int(sub["K"].nunique()),
                "KOs": "; ".join(sorted(sub["K"].unique())),
                "genes": "; ".join(sorted(sub["gene"].dropna().astype(str).unique())),
            }
        )
    pd.DataFrame(summary_rows).to_csv(FIG5_DIR / f"{prefix}_category_summary.tsv", sep="\t", index=False)


def main() -> None:
    FIG5_DIR.mkdir(parents=True, exist_ok=True)
    tpm = pd.read_csv(TPM_FILE, sep="\t", usecols=["feature-id"])
    present_kos = set(tpm["feature-id"].astype(str))

    outputs = {
        "Figure5_antistress_GSHaa_main": "original_m00565",
        "Figure5_antistress_GSHaa_trehalose_full_reference": "full_reference",
        "Figure5_antistress_GSHaa_trehalose_merged": "merged_m00565_plus_full",
        "Figure5_antistress_GSHaa_trehalose_bacterial_conservative": "bacterial_conservative",
    }
    combined_audits = []
    for prefix, mode in outputs.items():
        catalog, audit = build_catalog(mode, present_kos)
        audit.insert(0, "output_prefix", prefix)
        combined_audits.append(audit)
        write_abundance_tables(prefix, catalog)
        print(f"{prefix}: wrote {catalog['K'].nunique()} unique KOs")
        print(catalog["category"].value_counts().sort_index().to_string())

    pd.concat(combined_audits, ignore_index=True).to_csv(
        FIG5_DIR / "Figure5_antistress_GSHaa_candidate_audit.tsv",
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    main()
