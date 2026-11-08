"""MkDocs hook: write llms.txt and llms-full.txt into the site root.

Follows https://llmstxt.org/. Uses only the nav, each page's Markdown source and
``config.site_url``, so it works on any release tree (the Pages root is built from the
latest release tag with this hook added by scripts/docs/build_site.sh).

* ``llms.txt``: ``# <site_name>``, the site description, one ``##`` section per top-level
  nav section with ``- [Title](absolute URL): one line`` entries, and a final section with
  the full-text file, the repository, PyPI and AGENTS.md.
* ``llms-full.txt``: every nav page's Markdown source in nav order, each preceded by
  ``# <title>`` and its canonical URL, with HTML comments (``<!-- doc-setup ... -->``, run
  markers) removed.

The one-line text of an entry is the page's ``description`` front matter if present, else
its first paragraph (admonitions, tables, lists, code and generated notices skipped),
with Markdown links and emphasis flattened and trimmed to about 200 characters.
"""

import os
import re

_COMMENT = re.compile(r"<!--.*?-->[ \t]*\n?", re.S)
_FENCE = re.compile(r"^\s*(```|~~~)")
_LINK = re.compile(r"\[([^\]]*)\]\([^)]*\)")
_MAX = 200

_pages = {}  # src_uri -> (title, url, markdown, meta)


def _strip_comments(text):
    return _COMMENT.sub("", text)


def _first_paragraph(md):
    para, in_code = [], False
    for line in md.splitlines():
        if _FENCE.match(line):
            in_code = not in_code
            continue
        s = line.strip()
        if in_code:
            continue
        if not s:
            if para:
                break
            continue
        # skip headings, admonitions and their indented bodies, tables, lists, raw HTML, math
        if s.startswith(("#", "!!!", "???", "![", "|", "- ", "* ", "<", "$$", ">")) or re.match(r"\d+\. ", s):
            if para:
                break
            continue
        if line.startswith(("    ", "\t")) and not para:
            continue
        para.append(s)
    text = " ".join(para)
    text = _LINK.sub(r"\1", text)
    text = re.sub(r"\{[^}]*\}", "", text) if "{:" in text or '{ ' in text else text
    text = re.sub(r"[*`]+", "", text)
    text = re.sub(r"(?<!\w)_+|_+(?!\w)", "", text)
    text = re.sub(r"\s+", " ", text).strip()
    if len(text) > _MAX:
        cut = text[:_MAX]
        sentence = max(cut.rfind(". "), -1)
        text = cut[: sentence + 1] if sentence > 80 else cut.rsplit(" ", 1)[0] + "…"
    return text


def on_page_markdown(markdown, page, config, files):
    # page.title is resolved (nav title, else first heading); the markdown is the source
    # without front matter.
    _pages[page.file.src_uri] = (page.title, page.canonical_url, markdown, dict(page.meta or {}))
    return markdown


def _nav_sections(nav):
    """[(section title or None, [Page, ...])] in nav order, top-level sections only."""
    out, loose = [], []

    def pages(item):
        if getattr(item, "is_page", False):
            return [item]
        found = []
        for child in getattr(item, "children", None) or []:
            found += pages(child)
        return found

    for item in nav.items:
        if getattr(item, "is_page", False):
            loose.append(item)
        else:
            if loose:
                out.append((None, loose))
                loose = []
            out.append((item.title, pages(item)))
    if loose:
        out.append((None, loose))
    return out


def _url(config, page):
    return (config.get("site_url") or "").rstrip("/") + "/" + (page.url or "")


_CURRENT_NAV = {}


def on_nav(nav, config, files):
    _CURRENT_NAV["nav"] = nav
    return nav


def on_post_build(config):
    _write(config, config["site_dir"], _CURRENT_NAV.get("nav"))


def _write(config, site_dir, nav):
    if nav is None:
        return
    base = (config.get("site_url") or "").rstrip("/") + "/"
    name = config.get("site_name") or "Skarstrom"
    desc = config.get("site_description") or ""

    index = [f"# {name}", "", f"> {desc}", ""]
    full = [f"# {name}", "", f"> {desc}", "", f"Site: {base}", ""]
    for section, pages in _nav_sections(nav):
        entries = []
        for page in pages:
            src = page.file.src_uri
            if src not in _pages:
                continue
            title, _, md, meta = _pages[src]
            url = _url(config, page)
            line = meta.get("description") or _first_paragraph(_strip_comments(md))
            entries.append(f"- [{title}]({url})" + (f": {line}" if line else ""))
            body = _strip_comments(md).strip()
            full += ["", "---", "", f"# {title}", "", f"Canonical URL: {url}", "", body, ""]
        if entries:
            index += [f"## {section or 'Pages'}", ""] + entries + [""]

    repo = config.get("repo_url") or "https://github.com/samaffleck/skarstrom"
    index += [
        "## Optional",
        "",
        f"- [Full documentation in one file]({base}llms-full.txt): every page's Markdown in nav order",
        f"- [GitHub repository]({repo}): source, issues and releases",
        "- [PyPI](https://pypi.org/project/skarstrom/): `pip install skarstrom`",
        f"- [AGENTS.md]({repo}/blob/main/AGENTS.md): guide for coding agents using and working on Skarstrom",
        "",
    ]
    os.makedirs(site_dir, exist_ok=True)
    with open(os.path.join(site_dir, "llms.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(index))
    with open(os.path.join(site_dir, "llms-full.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(full).rstrip() + "\n")
