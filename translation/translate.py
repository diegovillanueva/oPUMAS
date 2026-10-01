#!/usr/bin/env python3
"""🔤🔁 Translate identifiers in free-form Fortran between the languages of dictionary.csv.

🧭📘 Usage: translate.py --from PUMAS --to CF-CCPP [--no-units] [--wrap] [--root DIR] [files...]
📂🗂️ Without files, every file listed in files.csv is translated in place.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
META_COLUMNS = {"id", "kind", "scope", "units", "meaning", "cf_ccpp_source", "standard_name", "notes"}
TOKEN = re.compile(
    r"""(?P<str>'(?:[^']|'')*'|"(?:[^"]|"")*")
       |(?P<num>(?:\d+\.?\d*|\.\d+)(?:[eEdD][+-]?\d+)?(?:_\w+)?)
       |(?P<dotop>\.[A-Za-z]+\.)
       |(?P<id>[A-Za-z]\w*)
       |(?P<other>.)""",
    re.X | re.S,
)
DIRECTIVE = re.compile(r"^(\s*!\$(?:acc|omp)\b)(.*)$", re.I | re.S)
TYPE_BEGIN = re.compile(r"^\s*type\b(?!\s*\()(?:\s*,[^:]*)?(?:\s*::)?\s*(\w+)", re.I)
TYPE_END = re.compile(r"^\s*end\s*type\b", re.I)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def language_columns(rows: list[dict[str, str]]) -> list[str]:
    return [c for c in rows[0] if c not in META_COLUMNS]


def unit_suffix(units: str) -> str:
    """🧮 'kg kg-1 s-1' -> '_kg_per_kg_per_s'; '1' or '' -> ''."""
    parts = []
    for tok in units.split():
        m = re.fullmatch(r"([A-Za-z]+)(-?\d+)?", tok)
        if not m:
            return ""
        unit, power = m.group(1), int(m.group(2) or 1)
        text = unit + (str(abs(power)) if abs(power) > 1 else "")
        parts.append(("per_" if power < 0 else "") + text)
    return "_" + "_".join(parts) if parts else ""


def spelled(row: dict[str, str], lang: str, suffix_langs: set[str]) -> str:
    name = row[lang].strip()
    if name and lang in suffix_langs:
        name += unit_suffix(row["units"])
    return name


def suffix_languages(units: bool) -> set[str]:
    if not units:
        return set()
    return {r["language"] for r in read_csv(HERE / "languages.csv") if r["unit_suffix"] == "yes"}


class Maps:
    """🗺️ Lookup tables for one translation direction."""

    def __init__(self, rows: list[dict[str, str]], src: str, dst: str, units: bool = True) -> None:
        suf = suffix_languages(units)
        self.glob: dict[str, tuple[str, str]] = {}
        self.per_file: dict[str, dict[str, tuple[str, str]]] = {}
        self.comp: dict[str, dict[str, tuple[str, str]]] = {}
        for row in rows:
            a, b = spelled(row, src, suf), spelled(row, dst, suf)
            if not a or not b:
                continue
            scope = row["scope"].strip()
            if scope == "global":
                table = self.glob
            elif scope.startswith("type:"):
                table = self.comp.setdefault(scope[5:].lower(), {})
            elif scope.startswith("file:"):
                table = self.per_file.setdefault(scope[5:], {})
            else:
                raise ValueError(f"unknown scope {scope!r} in row {row['id']}")
            if a.lower() in table:
                raise ValueError(f"duplicate {src} name {a!r} in scope {scope}")
            table[a.lower()] = (a, b)
        self.any_comp: dict[str, tuple[str, str]] = {}
        for table in self.comp.values():
            self.any_comp.update(table)


def swap(token: str, table: dict[str, tuple[str, str]], problems: list[str] | None) -> str:
    hit = table.get(token.lower())
    if hit is None:
        return token
    a, b = hit
    if token == a:
        return b
    if token == a.upper() and b != b.upper():
        return b.upper()
    if problems is not None:
        problems.append(token)
    return b


def translate_code(code: str, table: dict[str, tuple[str, str]], comp: dict[str, tuple[str, str]],
                   problems: list[str] | None) -> str:
    out, prev = [], ""
    pos = 0
    for m in TOKEN.finditer(code):
        text = m.group(0)
        if m.lastgroup == "other" and text == "!":
            pos = m.start()
            out.append(code[pos:])
            return "".join(out)
        if m.lastgroup == "id":
            text = swap(text, comp if prev == "%" else table, problems)
        if not text.isspace():
            prev = text
        out.append(text)
    return "".join(out)


def translate_text(text: str, maps: Maps, filename: str, mode: str = "full",
                   problems: list[str] | None = None) -> str:
    """🔁 Translate one file's text. mode 'components' only touches derived-type components."""
    table: dict[str, tuple[str, str]] = {}
    if mode == "full":
        table = dict(maps.glob)
        table.update(maps.per_file.get(filename, {}))
    out, in_type = [], None
    for line in text.splitlines(keepends=True):
        stripped = line.lstrip()
        if stripped.startswith("#"):
            out.append(line)
            continue
        m = DIRECTIVE.match(line)
        if m:
            out.append(m.group(1) + translate_code(m.group(2), table, maps.any_comp, problems))
            continue
        if in_type is None:
            t = TYPE_BEGIN.match(line)
            if t and t.group(1).lower() != "is":
                in_type = t.group(1).lower()
            out.append(translate_code(line, table, maps.any_comp, problems))
            continue
        if TYPE_END.match(line):
            in_type = None
            out.append(line)
            continue
        # 🧱 inside a type definition only that type's own components are translated
        out.append(translate_code(line, maps.comp.get(in_type, {}), maps.any_comp, problems))
    return "".join(out)


def wrap(text: str, width: int = 132) -> str:
    """✂️ One-way export step: break over-long code lines at a comma or blank with '&'."""
    out = []
    for line in text.splitlines(keepends=True):
        body = line.rstrip("\n")
        if len(body) <= width or body.lstrip().startswith(("!", "#")):
            out.append(line)
            continue
        sentinel = ""
        while len(body) > width:
            cut = max(body.rfind(",", 0, width - 2), body.rfind(" ", 20, width - 2))
            if cut < 20 or body.count("'", 0, cut + 1) % 2 or body.count('"', 0, cut + 1) % 2 or "!" in body[:cut]:
                break
            out.append(body[: cut + 1] + " &\n")
            indent = " " * (len(body) - len(body.lstrip()) + 5)
            body = indent + sentinel + body[cut + 1:].lstrip()
        out.append(body + "\n")
    return "".join(out)


def file_modes(root: Path) -> dict[str, str]:
    return {r["path"]: r["mode"] for r in read_csv(root / "translation" / "files.csv")}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--from", dest="src", required=True)
    ap.add_argument("--to", dest="dst", required=True)
    ap.add_argument("--no-units", action="store_true", help="do not append unit suffixes")
    ap.add_argument("--wrap", action="store_true", help="break lines longer than 132 characters (one-way export)")
    ap.add_argument("--root", type=Path, default=HERE.parent, help="repository root holding translation/")
    ap.add_argument("files", nargs="*", help="paths relative to the root; default: all of files.csv")
    args = ap.parse_args()
    root = args.root.resolve()
    rows = read_csv(root / "translation" / "dictionary.csv")
    langs = language_columns(rows)
    for lang in (args.src, args.dst):
        if lang not in langs:
            sys.exit(f"unknown language {lang!r}; available: {', '.join(langs)}")
    maps = Maps(rows, args.src, args.dst, units=not args.no_units)
    modes = file_modes(root)
    written = []
    for rel in args.files or list(modes):
        path = root / rel
        if not path.is_file():
            continue
        old = path.read_text()
        new = translate_text(old, maps, Path(rel).name, modes.get(rel, "full"))
        if args.wrap:
            new = wrap(new)
        if new != old:
            path.write_text(new)
            written.append(path)
    for path in written:
        print(path)


if __name__ == "__main__":
    main()
