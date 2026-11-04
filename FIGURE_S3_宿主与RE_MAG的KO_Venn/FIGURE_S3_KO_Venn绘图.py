import json
from pathlib import Path
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor, black, white

output_dir = Path("output/takakia_ko_comparison")
counts = json.loads((output_dir / "venn_counts.json").read_text())
regions = counts["regions"]
total = counts["union"]

def n(key):
    return regions[key]

pdf_path = output_dir / "Takakia_R_E_MAGs_KO_Venn.pdf"
c = canvas.Canvas(str(pdf_path), pagesize=(800, 760))
c.setTitle("Takakia and R/E MAG KEGG ortholog overlap")
c.setFillColor(white)
c.rect(0, 0, 800, 760, fill=1, stroke=0)
c.setFillColor(black)
c.setFont("Helvetica-Bold", 16)
c.drawCentredString(400, 730, "KEGG ortholog overlap: Takakia vs. R/E MAGs")
c.setFont("Helvetica", 9)
c.drawCentredString(400, 711, "Presence/absence at the KO level; each MAG group is the union across its MAGs")

circles = [
    (255, 440, HexColor("#29a879")),  # Takakia
    (545, 440, HexColor("#cf6570")),  # Endophyte
    (400, 245, HexColor("#6da8d6")),  # Rhizoidsphere
]
for x, y, color in circles:
    c.setFillColor(color)
    c.circle(x, y, 185, fill=1, stroke=0)

c.setFillColor(black)
c.setFont("Helvetica-Bold", 14)
c.drawCentredString(185, 543, "Takakia")
c.drawCentredString(615, 543, "Endophytic")
c.drawCentredString(615, 525, "MAGs")
c.drawCentredString(400, 110, "Rhizomatousphere MAGs")

def label(x, y, value):
    c.setFont("Helvetica-Bold", 12)
    c.drawCentredString(x, y, str(value))
    c.setFont("Helvetica", 10)
    c.drawCentredString(x, y - 15, f"({value / total * 100:.1f}%)")

label(160, 410, n("Takakia"))
label(400, 535, n("Takakia & Endophytic_MAGs"))
label(281, 315, n("Takakia & Rhizomatousphere_MAGs"))
label(400, 375, n("Takakia & Rhizomatousphere_MAGs & Endophytic_MAGs"))
label(640, 410, n("Endophytic_MAGs"))
label(520, 315, n("Rhizomatousphere_MAGs & Endophytic_MAGs"))
label(400, 205, n("Rhizomatousphere_MAGs"))

c.setFont("Helvetica", 8.5)
c.drawCentredString(400, 52, f"Takakia: {counts['set_sizes']['Takakia']:,} KOs; R MAGs: {counts['set_sizes']['Rhizomatousphere_MAGs']:,}; E MAGs: {counts['set_sizes']['Endophytic_MAGs']:,}; union: {total:,}.")
c.drawCentredString(400, 38, f"MAGs included: R = {counts['R_MAGs']}, E = {counts['E_MAGs']}. A KO is present when detected in at least one member of the set.")
c.save()
print(pdf_path)
