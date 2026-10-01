#!/usr/bin/env python3
"""✅🔎 Static checks of dictionary.csv against the Fortran sources.

🧭📘 Usage: check.py [--root DIR] [--branches]
🧪📋 Checks: legal and unique names, protected public names, collisions, case variants, byte-exact round trip.
🌿🔁 With --branches: translate(tree of branch A) must equal the tree of branch B for every language pair.
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

from translate import DIRECTIVE, TOKEN, Maps, file_modes, language_columns, read_csv, spelled, suffix_languages, translate_text

RESERVED = set("""abs acos aimag aint all allocatable allocate allocated anint any asin associated atan atan2 call case
ceiling character class close cmplx complex contains cos cosh count cycle data dble deallocate dim dimension do
dot_product double elemental else elsewhere end enddo endif epsilon exit exp floor function gamma huge if implicit in
index inout int integer intent interface kind len log log10 logical matmul max maxloc maxval merge min minloc minval
mod module modulo nint none null open optional out pack parameter pointer present print private procedure product
program public pure read real recursive reshape result return save select shape sign sin sinh size spread sqrt stop
subroutine sum tan tanh target then tiny transfer trim type ubound lbound unpack use where while write erf erfc""".split())
USE_ONLY = re.compile(r"^\s*use\s+(\w+)\s*(?:,\s*only\s*:\s*(.*))?$", re.I)
PROC = re.compile(r"^\s*(?:(?:pure|elemental|impure|recursive)\s+)*(?:subroutine|function)\s+(\w+)\s*\(([^)]*)\)", re.I)
MODULE = re.compile(r"^\s*module\s+(?!procedure\b)(\w+)", re.I)


def statements(text: str) -> list[str]:
    """🧵 Join continuation lines, drop comments (strings in the sources never contain '!')."""
    out, cur = [], ""
    for line in text.splitlines():
        if line.lstrip().startswith(("!", "#")):
            continue
        code = line.split("!")[0].strip()
        if not code:
            continue
        cur += code.lstrip("&").rstrip("&")
        if not code.endswith("&"):
            out.append(cur)
            cur = ""
    return out


def identifiers(text: str) -> tuple[set[str], set[str]]:
    """🔤 Lower-case identifiers of a file: (plain, used as %component)."""
    plain, comp = set(), set()
    for line in text.splitlines():
        if line.lstrip().startswith("#"):
            continue
        m = DIRECTIVE.match(line)
        code = m.group(2) if m else line
        prev = ""
        for t in TOKEN.finditer(code):
            if t.lastgroup == "other" and t.group(0) == "!":
                break
            if t.lastgroup == "id":
                (comp if prev == "%" else plain).add(t.group(0).lower())
            if not t.group(0).isspace():
                prev = t.group(0)
    return plain, comp


def protected_names(root: Path, modes: dict[str, str]) -> tuple[set[str], list[str]]:
    """🛡️ Names of the fully translated modules that untranslated code can see."""
    errors: list[str] = []
    full = [p for p, m in modes.items() if m == "full"]
    full_stmts = {p: statements((root / p).read_text()) for p in full}
    modules, dummies = set(), {}
    for stmts in full_stmts.values():
        for s in stmts:
            if m := MODULE.match(s):
                modules.add(m.group(1).lower())
            if m := PROC.match(s):
                dummies[m.group(1).lower()] = {a.strip().lower() for a in m.group(2).split(",") if a.strip()}
    others = [p for p in list(root.glob("*.F90")) + list(root.glob("interface/*/*.F90"))
              if not p.is_symlink() and str(p.relative_to(root)) not in full]
    imported = set()
    for p in others:
        for s in statements(p.read_text(errors="replace")):
            m = USE_ONLY.match(s)
            if not m or m.group(1).lower() not in modules:
                continue
            if m.group(2) is None:
                errors.append(f"{p.name}: 'use {m.group(1)}' without only-list, cannot determine protected names")
                continue
            for item in m.group(2).split(","):
                imported.add(item.split("=>")[-1].strip().lower())
    prot = set(modules) | imported
    for stmts in full_stmts.values():
        for s in stmts:
            m = USE_ONLY.match(s)
            if m and m.group(1).lower() not in modules:
                prot.add(m.group(1).lower())
                for item in (m.group(2) or "").split(","):
                    prot |= {x.strip().lower() for x in item.split("=>") if x.strip()}
    for name in imported:
        prot |= dummies.get(name, set())
    for stmts in full_stmts.values():
        for s in stmts:
            if m := re.match(r"^\s*namelist\s*/\s*\w+\s*/(.*)$", s, re.I):
                prot |= {a.strip().lower() for a in m.group(1).split(",")}
    return prot, errors


def run_checks(root: Path) -> list[str]:
    errors: list[str] = []
    rows = read_csv(root / "translation" / "dictionary.csv")
    langs = language_columns(rows)
    modes = file_modes(root)
    texts = {p: (root / p).read_text() for p in modes}
    prot, errs = protected_names(root, modes)
    errors += errs
    suf = suffix_languages(True)
    ids = {p: identifiers(t) for p, t in texts.items()}
    plain_all = set().union(*(ids[p][0] for p, m in modes.items() if m == "full"))
    comp_all = set().union(*(c for _, c in ids.values()))
    ref = langs[0]
    comp_names = {r[ref].lower() for r in rows if r["scope"].startswith("type:")}
    for p in list(root.glob("*.F90")) + list(root.glob("interface/*/*.F90")):
        rel = str(p.relative_to(root))
        if p.is_symlink() or rel in modes:
            continue
        hit = identifiers(p.read_text(errors="replace"))[1] & comp_names
        if hit:
            errors.append(f"{rel}: uses translated components {sorted(hit)[:5]} but is not in files.csv")
    seen: set[str] = set()
    for r in rows:
        if r["id"] in seen:
            errors.append(f"duplicate id {r['id']}")
        seen.add(r["id"])
        src = r[ref].lower()
        is_comp = r["scope"].startswith("type:")
        if not is_comp and src in prot:
            errors.append(f"{r['id']}: {r[ref]} is a protected public name")
        if src not in (comp_all if is_comp else plain_all):
            errors.append(f"{r['id']}: {r[ref]} does not occur in the translated files")
        for lang in langs[1:]:
            name = spelled(r, lang, suf)
            if not name:
                continue
            low = name.lower()
            if not re.fullmatch(r"[A-Za-z]\w*", name) or len(name) > 63:
                errors.append(f"{r['id']}: {lang} name {name!r} is illegal or longer than 63 ({len(name)})")
            if low in RESERVED or low in prot:
                errors.append(f"{r['id']}: {lang} name {name!r} is a keyword, intrinsic or protected name")
            if low in (comp_all if is_comp else plain_all):
                errors.append(f"{r['id']}: {lang} name {name!r} already exists in the translated files")
    for lang in langs[1:]:
        try:
            fwd, back = Maps(rows, ref, lang), Maps(rows, lang, ref)
        except ValueError as exc:
            errors.append(str(exc))
            continue
        if set(fwd.glob) & set().union(*fwd.per_file.values(), set()):
            errors.append(f"{lang}: a name is both global and file-scoped")
        for p, text in texts.items():
            problems: list[str] = []
            name = Path(p).name
            new = translate_text(text, fwd, name, modes[p], problems)
            if problems:
                errors.append(f"{p}: unsupported case variants {sorted(set(problems))[:8]}")
            if translate_text(new, back, name, modes[p]) != text:
                errors.append(f"{p}: round trip {ref} -> {lang} -> {ref} is not byte-identical")
    return errors


def check_branches(root: Path) -> list[str]:
    """🌿 The tree of every language branch must equal the translation of the reference branch."""
    errors: list[str] = []
    rows = read_csv(root / "translation" / "dictionary.csv")
    langs = read_csv(root / "translation" / "languages.csv")
    modes = file_modes(root)

    def git(*a: str) -> str:
        return subprocess.run(["git", "-C", str(root), *a], check=True, capture_output=True, text=True).stdout

    ref = langs[0]
    for lang in langs[1:]:
        changed = set(git("diff", "--name-only", ref["branch"], lang["branch"]).split())
        extra = changed - set(modes)
        if extra:
            errors.append(f"{lang['branch']}: differs from {ref['branch']} in untranslated files {sorted(extra)[:5]}")
        rows_b = rows
        maps = Maps(rows_b, ref["language"], lang["language"])
        for p in modes:
            a = git("show", f"{ref['branch']}:{p}")
            b = git("show", f"{lang['branch']}:{p}")
            if translate_text(a, maps, Path(p).name, modes[p]) != b:
                errors.append(f"{p}: {lang['branch']} is not the translation of {ref['branch']}")
    return errors


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent)
    ap.add_argument("--branches", action="store_true")
    args = ap.parse_args()
    root = args.root.resolve()
    errors = check_branches(root) if args.branches else run_checks(root)
    for e in errors:
        print("ERROR:", e)
    print(f"{len(errors)} error(s)")
    sys.exit(1 if errors else 0)


if __name__ == "__main__":
    main()
