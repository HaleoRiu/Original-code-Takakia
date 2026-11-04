import fs from "node:fs/promises";

const takakiaPath = "/Volumes/Haleo_DK01/WD_Takakia/pathview/takakia_ko_list.txt";
const magPath = "/Volumes/Haleo_DK01/WD_Takakia/next_steps/MAGs_KO_presence_absence_matrix.txt";
const groupPath = "/Volumes/Haleo_DK01/WD_Takakia/metabolic_70/group_MAGs.txt";

const takakia = new Set((await fs.readFile(takakiaPath, "utf8")).match(/K\d{5}/g) || []);
const groupRows = (await fs.readFile(groupPath, "utf8")).trim().split(/\r?\n/).slice(1).map(x => x.split("\t"));
const groupMap = new Map(groupRows.map(([id, group]) => [id.replaceAll(" ", "_"), group]));
const [header, ...rows] = (await fs.readFile(magPath, "utf8")).trim().split(/\r?\n/).map(x => x.split("\t"));
const mags = header.slice(1);
const groups = { Rhizoidsphere: new Set(), Endophyte: new Set() };
for (const row of rows) {
  const ko = row[0];
  row.slice(1).forEach((v, i) => {
    const group = groupMap.get(mags[i]);
    if (!group) throw new Error(`Unclassified MAG: ${mags[i]}`);
    if (Number(v) > 0) groups[group].add(ko);
  });
}
const sets = { Takakia: takakia, Rhizomatousphere_MAGs: groups.Rhizoidsphere, Endophytic_MAGs: groups.Endophyte };
const names = Object.keys(sets);
const regions = {};
const union = new Set(names.flatMap(n => [...sets[n]]));
for (const ko of union) {
  const key = names.filter(n => sets[n].has(ko)).join(" & ");
  regions[key] = (regions[key] || 0) + 1;
}
console.log(JSON.stringify({
  total_MAGs: mags.length,
  R_MAGs: mags.filter(m => groupMap.get(m) === "Rhizoidsphere").length,
  E_MAGs: mags.filter(m => groupMap.get(m) === "Endophyte").length,
  set_sizes: Object.fromEntries(names.map(n => [n, sets[n].size])),
  regions,
  union: union.size,
}, null, 2));
