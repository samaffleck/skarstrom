#!/usr/bin/env python3
"""Convert the paper PDFs in papers/<paper>/assets/ into Markdown.

Uses Docling (https://github.com/docling-project/docling):

  - body text and section headings as Markdown
  - tables as Markdown tables
  - figures exported as PNGs into ``<stem>_artifacts/``, referenced relatively
  - equations as ``<!-- formula-not-decoded -->`` placeholders, or as LaTeX
    (``$$...$$``) with ``--formulas``

This is a LOCAL reading aid only. PDFs, the Markdown extractions and the
extracted figures are copyrighted paper material and are git-ignored (see
papers/README.md): the repo keeps only the transcribed numbers
(reference_values.csv) and the digitised curves (adsim/data/).

Usage (docling installed via pipx / uv / venv):

    python3 -m pip install docling            # or: pipx install docling
    python3 adsim/test/papers/extract_pdfs.py            # only missing papers
    python3 adsim/test/papers/extract_pdfs.py --force    # re-extract everything
    python3 adsim/test/papers/extract_pdfs.py perez_2019 # a single paper folder

First run downloads the Docling layout / table models (~500 MB).

``--formulas`` turns on Docling's formula enrichment, which recovers the model
equations as LaTeX. It is worth it for a supplementary-information PDF full of
model equations, but it is *very* slow on CPU (minutes per equation-dense page,
and it has been seen to stall for hours on a whole paper) and it mangles
equations that sit next to two-column body text. Run it per paper, not over the
whole folder.
"""

from __future__ import annotations

import argparse
import re
import sys
import time
from pathlib import Path

from docling.datamodel.base_models import InputFormat
from docling.datamodel.pipeline_options import PdfPipelineOptions, TableFormerMode
from docling.document_converter import DocumentConverter, PdfFormatOption
from docling_core.types.doc import ImageRefMode

PAPERS = Path(__file__).resolve().parent

# PDF stem -> markdown stem. Anything not listed keeps the PDF stem.
NAME_OVERRIDES = {
    "arora_2021_paper": "2021_Arora",
    "arora_2021_si": "2021_Arora_SI",
}

MARKER = "<!-- extracted from {pdf} by adsim/test/papers/extract_pdfs.py (docling) -->"
MARKER_PREFIX = "<!-- extracted from "
PLACEHOLDER = "<!-- formula-not-decoded -->"
FORMULA_RE = re.compile(r"\$\$(.*?)\$\$", re.DOTALL)


def make_converter(formulas: bool, force_ocr: bool) -> DocumentConverter:
    opts = PdfPipelineOptions()
    opts.images_scale = 2.0
    opts.generate_picture_images = True
    opts.generate_page_images = False
    opts.do_table_structure = True
    opts.table_structure_options.mode = TableFormerMode.ACCURATE
    opts.do_ocr = True
    opts.ocr_options.force_full_page_ocr = force_ocr
    opts.do_formula_enrichment = formulas
    return DocumentConverter(
        format_options={InputFormat.PDF: PdfFormatOption(pipeline_options=opts)}
    )


def sanitize_formulas(text: str) -> tuple[str, int]:
    """Revert mangled LaTeX blocks to the not-decoded placeholder."""
    dropped = 0

    def repl(m: re.Match) -> str:
        nonlocal dropped
        latex = m.group(1)
        mangled = (
            len(latex) > 800
            or latex.count(r"\quad") > 15
            or latex.count(r"{ \quad }") > 4
        )
        if mangled:
            dropped += 1
            return PLACEHOLDER
        return m.group(0)

    return FORMULA_RE.sub(repl, text), dropped


def md_target(pdf: Path) -> Path:
    return pdf.parent / f"{NAME_OVERRIDES.get(pdf.stem, pdf.stem)}.md"


def already_extracted(md_path: Path) -> bool:
    return md_path.is_file() and MARKER_PREFIX in md_path.read_text(encoding="utf-8")


def convert(pdf: Path, converter: DocumentConverter) -> tuple[Path, int, int]:
    md_path = md_target(pdf)
    artifacts = pdf.parent / f"{md_path.stem}_artifacts"

    doc = converter.convert(pdf).document
    doc.save_as_markdown(md_path, image_mode=ImageRefMode.REFERENCED, artifacts_dir=artifacts)

    text = md_path.read_text(encoding="utf-8")
    # Docling writes absolute artifact paths; make them relative to the .md.
    text = text.replace(f"{md_path.parent}/", "")
    text, dropped = sanitize_formulas(text)
    text = f"{text.rstrip()}\n\n{MARKER.format(pdf=pdf.name)}\n"
    md_path.write_text(text, encoding="utf-8")

    kept = text.count("$$") // 2
    return md_path, kept, dropped


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("papers", nargs="*", help="paper folder names (default: all)")
    ap.add_argument("--force", action="store_true", help="re-extract already-converted PDFs")
    ap.add_argument(
        "--formulas",
        action="store_true",
        help="convert equations to LaTeX (accurate but very slow — see module docstring)",
    )
    ap.add_argument(
        "--ocr",
        action="store_true",
        help="ignore the PDF text layer and OCR every page — use when a PDF's fonts "
        "are broken (mangled superscripts/units); costs some accuracy elsewhere",
    )
    args = ap.parse_args(argv)

    folders = args.papers or ["*"]
    pdfs = sorted({p for f in folders for p in PAPERS.glob(f"{f}/assets/*.pdf")})
    if not pdfs:
        print("no PDFs found under papers/*/assets/", file=sys.stderr)
        return 1

    todo = [p for p in pdfs if args.force or not already_extracted(md_target(p))]
    for pdf in pdfs:
        if pdf not in todo:
            print(f"--- {pdf.parent.parent.name}/{pdf.name}: up to date, skipping", flush=True)
    if not todo:
        return 0

    converter = make_converter(args.formulas, args.ocr)
    summary = []
    for pdf in todo:
        t0 = time.time()
        print(f"==> {pdf.parent.parent.name}/{pdf.name}", flush=True)
        md_path, kept, dropped = convert(pdf, converter)
        text = md_path.read_text(encoding="utf-8")
        images = len(list((md_path.parent / f"{md_path.stem}_artifacts").glob("*.png")))
        line = (
            f"{md_path.parent.parent.name}/{md_path.name}: images={images} "
            f"table_rows={text.count(chr(10) + '|')} formulas={kept} "
            f"dropped={dropped} {time.time() - t0:.0f}s"
        )
        summary.append(line)
        print(f"    {line}", flush=True)

    print("\n=== summary ===")
    for line in summary:
        print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
