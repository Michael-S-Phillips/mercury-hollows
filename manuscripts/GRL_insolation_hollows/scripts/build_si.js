// Build the Supporting Information (.docx) from si.json and figs/figS*.png.
const fs = require("fs");
const { Document, Packer, Paragraph, TextRun, ImageRun, AlignmentType, Table, TableRow, TableCell, WidthType, BorderStyle,
        PageOrientation } = require("docx");
const S = JSON.parse(fs.readFileSync(__dirname + "/si.json", "utf8"));
const FONT = "Times New Roman", SZ = 24, LINE = 276;
const run = (text, o = {}) => new TextRun({ text, font: FONT, size: o.size || SZ, ...o });
function runs(text, o = {}) {
  return text.split(/(\[[^\]]*\])/g).filter(s => s.length).map(p => /^\[(VERIFY|To be|DOI|co-authors)/.test(p) ? run(p, { ...o, highlight: "yellow" }) : run(p, o));
}
const P = (text, o = {}) => new Paragraph({ spacing: { line: LINE, after: 160 }, alignment: o.align, children: runs(text, o.run || {}) });
const H = (text) => new Paragraph({ spacing: { line: LINE, before: 240, after: 120 }, keepNext: true, children: [run(text, { bold: true })] });
const b = [];
b.push(new Paragraph({ spacing: { after: 120 }, children: [run("Geophysical Research Letters", { italics: true })] }));
b.push(new Paragraph({ spacing: { after: 120 }, children: [run(S.title, { bold: true })] }));
b.push(new Paragraph({ spacing: { after: 0 }, children: runs(S.authors) }));
b.push(new Paragraph({ spacing: { after: 240 }, children: [run(S.affil)] }));
b.push(H("Contents of this file"));
S.contents.forEach(c => b.push(new Paragraph({ spacing: { after: 0 }, indent: { left: 360 }, children: [run(c)] })));
b.push(H("Introduction"));
b.push(P(S.intro));
for (const [h, ps] of S.texts) { b.push(H(h)); ps.forEach(p => b.push(P(p))); }
for (const [file, cap] of S.figures) {
  const png = fs.readFileSync(__dirname + "/figs/" + file);
  const w = png.readUInt32BE(16), h = png.readUInt32BE(20); const wpx = 6.5 * 96, hpx = Math.round(wpx * h / w);
  b.push(new Paragraph({ spacing: { before: 240, after: 120 }, alignment: AlignmentType.CENTER, keepNext: true,
    children: [new ImageRun({ type: "png", data: png, transformation: { width: wpx, height: hpx } })] }));
  b.push(P(cap));
}
// Table S1
b.push(H(S.table_title));
const widths = [1350, 2100, 2350, 2150, 1410];   // twips, sum 9360 = 6.5 in
const border = { style: BorderStyle.SINGLE, size: 4, color: "999999" };
const cell = (t, i, head) => new TableCell({ width: { size: widths[i], type: WidthType.DXA },
  borders: { top: border, bottom: border, left: border, right: border }, margins: { top: 40, bottom: 40, left: 80, right: 80 },
  shading: head ? { fill: "EDEBE6", type: "clear", color: "auto" } : undefined,
  children: [new Paragraph({ spacing: { after: 0 }, children: [run(t, { size: 18, bold: !!head })] })] });
const rows = [new TableRow({ tableHeader: true, children: S.table_header.map((t, i) => cell(t, i, true)) })]
  .concat(S.table.map(r => new TableRow({ children: r.map((t, i) => cell(t, i, false)) })));
b.push(new Table({ width: { size: 9360, type: WidthType.DXA }, columnWidths: widths, rows }));
b.push(new Paragraph({ spacing: { before: 240 }, children: runs(S.references_note, { size: 20 }) }));
const doc = new Document({ creator: "Michael S. Phillips", title: "Supporting Information",
  styles: { default: { document: { run: { font: FONT, size: SZ } } } },
  sections: [{ properties: { page: { size: { width: 12240, height: 15840 }, margin: { top: 1440, bottom: 1440, left: 1440, right: 1440 } } }, children: b }] });
const out = process.argv[2] || (__dirname + "/GRL_insolation_hollows_SI_draft_v2.docx");
Packer.toBuffer(doc).then(buf => { fs.writeFileSync(out, buf); console.log("written", out); });
