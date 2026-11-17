"""Single-file HTML source of truth. No project copies or sidecar exports."""
from __future__ import annotations

import json
from html.parser import HTMLParser
from pathlib import Path

from kbiosuite.version import __version__

from .dependencies import resource_root
from .model import validate

MARKER = "kstack-project"


def safe_json(value: dict) -> str:
    return (json.dumps(value, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
            .replace("&", "\\u0026").replace("<", "\\u003c").replace(">", "\\u003e")
            .replace("\u2028", "\\u2028").replace("\u2029", "\\u2029"))


def to_html(project: dict) -> str:
    validate(project)
    root = resource_root()
    css = (root / "renderer/styles.css").read_text()
    scripts = [(root / path).read_text() for path in
               ("renderer/model.js", "export/svg.js", "export/fasta.js",
                "renderer/gene-menu.js", "renderer/export-preview.js", "renderer/app.js")]
    template = (root / "renderer/template.html").read_text()
    return (template.replace("<!--KSTACK_VERSION-->", __version__)
            .replace("<!--KSTACK_STYLE-->", f"<style>{css}</style>")
            .replace("<!--KSTACK_DATA-->",
                     f'<script id="{MARKER}" type="application/json" data-schema="1">{safe_json(project)}</script>')
            .replace("<!--KSTACK_SCRIPTS-->", "\n".join(f"<script>{s}</script>" for s in scripts)))


class ProjectReader(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=False)
        self.depth = False
        self.chunks = []
        self.count = 0

    def handle_starttag(self, tag, attrs):
        attributes = dict(attrs)
        if tag == "script" and attributes.get("id") == MARKER:
            if attributes.get("type") != "application/json":
                raise ValueError("Invalid Kstack project marker.")
            self.depth = True
            self.count += 1

    def handle_endtag(self, tag):
        if tag == "script":
            self.depth = False

    def handle_data(self, data):
        if self.depth:
            self.chunks.append(data)


def from_html(html: str) -> dict:
    parser = ProjectReader()
    parser.feed(html)
    if parser.count != 1:
        raise ValueError("Expected one Kstack project marker in this HTML file.")
    return validate(json.loads("".join(parser.chunks)))


def load(path: str | Path) -> dict:
    return from_html(Path(path).read_text())
