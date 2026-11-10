// Build the GRL draft (.docx) from paper_text.json and figs/*.png.
// AGU initial-submission layout: US Letter, 1" margins, Times New Roman 12 pt, double spacing, continuous line numbers,
// figures placed after the paragraph that first cites them.
const fs = require("fs");
const { Document, Packer, Paragraph, TextRun, ImageRun, AlignmentType, LineNumberRestartFormat, PageNumber, Footer,
        LevelFormat } = require("docx");
const T = JSON.parse(fs.readFileSync(__dirname + "/paper_text.json", "utf8"));
const FONT = "Times New Roman", SZ = 24, LINE = 480;

const run = (text, o = {}) => new TextRun({ text, font: FONT, size: SZ, ...o });
// Highlight [VERIFY ...] and [To be completed ...] / [DOI ...] / [co-authors ...] markers in yellow so they are easy to find.
function runs(text, o = {}) {
  const parts = text.split(/(\[[^\]]*\])/g).filter(s => s.length);
  return parts.map(p => /^\[(VERIFY|To be|DOI|co-authors)/.test(p) ? run(p, { ...o, highlight: "yellow" }) : run(p, o));
}
const P = (text, o = {}) => new Paragraph({ spacing: { line: LINE, after: 0 }, indent: o.noIndent ? undefined : { firstLine: 720 },
  alignment: o.align, children: runs(text, o.run || {}) });
const H = (text) => new Paragraph({ spacing: { line: LINE, before: 240, after: 0 }, keepNext: true, children: [run(text, { bold: true })] });
const Label = (text) => new Paragraph({ spacing: { line: LINE, before: 240, after: 0 }, keepNext: true, children: [run(text, { bold: true })] });

const b = [];
b.push(new Paragraph({ spacing: { line: LINE, after: 0 }, children: [run("Geophysical Research Letters", { italics: true })] }));
b.push(new Paragraph({ spacing: { line: LINE, before: 240, after: 0 }, children: [run(T.title, { bold: true })] }));
b.push(new Paragraph({ spacing: { line: LINE, after: 0 }, children: runs(T.authors) }));
T.affil.forEach(a => b.push(new Paragraph({ spacing: { line: LINE, after: 0 }, children: [run(a)] })));
b.push(new Paragraph({ spacing: { line: LINE, before: 240, after: 0 }, children: [run(T.corresponding)] }));
b.push(Label("Key Points:"));
T.key_points.forEach(k => b.push(new Paragraph({ spacing: { line: LINE, after: 0 }, numbering: { reference: "kp", level: 0 }, children: [run(k)] })));
b.push(Label("Abstract"));
b.push(P(T.abstract, { noIndent: true }));
b.push(Label("Plain Language Summary"));
b.push(P(T.pls, { noIndent: true }));

const figPara = (n) => {
  const [file, cap] = T.figures[String(n)];
  const png = fs.readFileSync(__dirname + "/figs/" + file);
  // image width 6.5 in; height from the PNG header
  const w = png.readUInt32BE(16), h = png.readUInt32BE(20);
  const wpx = 6.5 * 96, hpx = Math.round(wpx * h / w);
  return [
    new Paragraph({ spacing: { before: 240, after: 120 }, alignment: AlignmentType.CENTER, keepNext: true,
      children: [new ImageRun({ type: "png", data: png, transformation: { width: wpx, height: hpx } })] }),
    new Paragraph({ spacing: { line: LINE, after: 240 }, children: runs(cap) }),
  ];
};

for (const [h, ps] of T.sections) {
  if (h === "FIG") { figPara(ps).forEach(x => b.push(x)); continue; }
  b.push(H(h));
  ps.forEach((p, i) => b.push(P(p, { noIndent: i === 0 })));
}
b.push(Label("Acknowledgments"));
b.push(P(T.ack, { noIndent: true }));
b.push(Label("Open Research"));
b.push(P(T.open_research, { noIndent: true }));
b.push(Label("References"));
T.references.forEach(r => b.push(new Paragraph({ spacing: { line: LINE, after: 0 }, indent: { left: 720, hanging: 720 }, children: runs(r) })));

const doc = new Document({
  creator: "Michael S. Phillips", title: T.title,
  styles: { default: { document: { run: { font: FONT, size: SZ } } } },
  numbering: { config: [{ reference: "kp", levels: [{ level: 0, format: LevelFormat.BULLET, text: "•", alignment: AlignmentType.LEFT,
    style: { paragraph: { indent: { left: 720, hanging: 360 } } } }] }] },
  sections: [{
    properties: {
      page: { size: { width: 12240, height: 15840 }, margin: { top: 1440, bottom: 1440, left: 1440, right: 1440 } },
      lineNumbers: { countBy: 1, restart: LineNumberRestartFormat.CONTINUOUS, distance: 360 },
    },
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.CENTER,
      children: [new TextRun({ children: [PageNumber.CURRENT], font: FONT, size: 20 })] })] }) },
    children: b,
  }],
});
const out = process.argv[2] || (__dirname + "/GRL_insolation_hollows_draft_v1.docx");
Packer.toBuffer(doc).then(buf => { fs.writeFileSync(out, buf); console.log("written", out); });
