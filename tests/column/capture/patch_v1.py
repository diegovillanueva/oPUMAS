#!/usr/bin/env python3
"""🔬🧷 Copy micro_pumas_v1.F90 and micro_pumas_utils.F90 with capture calls.

🧪 Every process call in micro_pumas_tend is recorded before and after the call, and the
   module-level constants of both modules are recorded once (site 'module_state').
🔁 Regenerated on every capture build, so it follows upstream merges automatically.
Usage: patch_v1.py <src_dir> <out_dir>
"""

import re
import sys
from pathlib import Path

ROUTINES = {
    "size_dist_param_liq", "size_dist_param_basic", "kk2000_liq_autoconversion",
    "sb2001v2_liq_autoconversion", "sb2001v2_accre_cld_water_rain", "ice_autoconversion",
    "immersion_freezing", "contact_freezing", "snow_self_aggregation",
    "accrete_cloud_water_snow", "secondary_ice_production", "accrete_rain_snow",
    "heterogeneous_rain_freezing", "accrete_cloud_water_rain", "self_collection_rain",
    "accrete_cloud_ice_snow", "bergeron_process_snow", "vapor_deposition_onto_snow",
    "ice_deposition_sublimation", "graupel_collecting_snow", "graupel_collecting_cld_water",
    "graupel_riming_liquid_snow", "graupel_collecting_rain", "graupel_rain_riming_snow",
    "graupel_rime_splintering", "evaporate_sublimate_precip_graupel",
    "evaporate_sublimate_precip", "qsat_water", "qsat_ice", "calc_rercld",
    "avg_diameter_vec", "Sedimentation", "Sedimentation_implicit",
}
CALL = re.compile(r"^(\s*)call\s+(\w+)\s*\(", re.IGNORECASE)


def strip_comment(line: str) -> str:
    return line.split("!", 1)[0].rstrip()


def split_args(text: str) -> list[str]:
    """✂️🧮 Split an argument list at top-level commas."""
    args, depth, cur = [], 0, ""
    for ch in text:
        if ch in "([":
            depth += 1
        elif ch in ")]":
            depth -= 1
        if ch == "," and depth == 0:
            args.append(cur.strip())
            cur = ""
        else:
            cur += ch
    if cur.strip():
        args.append(cur.strip())
    return args


PROPS_FIELDS = ["rho", "eff_dim", "shape_coef", "lambda_bounds", "min_mean_mass"]


def capture_lines(indent: str, key: str, phase: str, expr: str) -> list[str]:
    """🧊📏 One capture call, or one per field for a hydrometeor property struct."""
    if re.fullmatch(r"mg_\w+_props", expr):
        return [f"{indent}call capture('{key}_{f}', '{phase}', '{expr}%{f}', {expr}%{f})\n"
                for f in PROPS_FIELDS]
    return [f"{indent}call capture('{key}', '{phase}', '{expr}', {expr})\n"]


def captures(indent: str, site: str, args: list[str], phase: str) -> list[str]:
    out = [f"{indent}call capture_site('{site}', mgncol, nlev, '{phase}')\n"]
    for i, a in enumerate(args, 1):
        m = re.match(r"^(\w+)\s*=(?!=)\s*(.+)$", a)
        key, expr = (m.group(1), m.group(2)) if m else (f"a{i:02d}", a)
        out += capture_lines(indent, key, phase, expr)
    return out


def module_scalars(lines: list[str]) -> list[str]:
    """🗂️🔢 Names of module-level real(r8)/logical variables and parameters of rank 0 or 1."""
    names, in_type, i = [], False, 0
    while i < len(lines):
        text = strip_comment(lines[i])
        while text.endswith("&"):
            i += 1
            text = text[:-1] + " " + strip_comment(lines[i]).lstrip().lstrip("&")
        i += 1
        low = text.lower().strip()
        if low == "contains":
            break
        if re.match(r"^type\s*(::|\w)", low) and not low.startswith("type("):
            in_type = True
        if re.match(r"^end\s*type", low):
            in_type = False
            continue
        m = re.match(r"^\s*(real\(r8\)|logical)(?=[\s,:])([^:]*)::(.*)$", text, re.IGNORECASE)
        if in_type or not m:
            continue
        for item in split_args(m.group(3)):
            name = item.split("=")[0].strip()
            dims = re.match(r"^(\w+)\s*(\((.*)\))?$", name)
            if dims and (dims.group(3) is None or "," not in dims.group(3)):
                names.append(dims.group(1))
    return names


def state_subroutine(module: str, names: list[str], props: list[str]) -> str:
    body = [f"subroutine capture_{module}_state()\n",
            "  use pumas_column_capture, only: capture\n"]
    for n in names:
        body.append(f"  call capture('{module}_{n}', 'in', '{n}', {n})\n")
    for n in props:
        body += capture_lines("  ", f"{module}_{n}", "in", n)
    body.append(f"end subroutine capture_{module}_state\n")
    return "".join(body)


def patch_utils(src: Path, dst: Path) -> None:
    """🧰🔬 Add capture_utils_state() to micro_pumas_utils so v1 can record its constants."""
    lines = src.read_text().splitlines(keepends=True)
    names = module_scalars(lines)
    props = [f"mg_{h}_props" for h in ("liq", "ice", "rain", "snow", "graupel", "hail")]
    text = "".join(lines)
    i = re.search(r"^contains\s*$", text, re.MULTILINE).start()
    text = text[:i] + "public :: capture_utils_state\n\n" + text[i:]
    j = re.search(r"^end module micro_pumas_utils", text, re.MULTILINE).start()
    text = text[:j] + state_subroutine("utils", names, props) + "\n" + text[j:]
    dst.write_text(text)


def main() -> None:
    srcdir, outdir = Path(sys.argv[1]), Path(sys.argv[2])
    patch_utils(srcdir / "micro_pumas_utils.F90", outdir / "micro_pumas_utils.F90")
    src, dst = srcdir / "micro_pumas_v1.F90", outdir / "micro_pumas_v1.F90"
    lines = src.read_text().splitlines(keepends=True)
    v1_names = module_scalars(lines)
    out, count, in_tend, i, used, first = [], {}, False, 0, False, True
    while i < len(lines):
        line = lines[i]
        low = line.lower()
        if re.match(r"^\s*subroutine\s+micro_pumas_tend\b", low):
            in_tend = True
        if re.match(r"^\s*end\s+subroutine\s+micro_pumas_tend\b", low):
            in_tend = False
        if not used and re.match(r"^\s*use\s+micro_pumas_utils\b", low):
            out.append("use pumas_column_capture, only: capture, capture_site\n")
            out.append("use micro_pumas_utils, only: capture_utils_state\n")
            used = True
        m = CALL.match(line)
        if in_tend and m and m.group(2) in ROUTINES:
            stmt = [line]
            while strip_comment(lines[i]).endswith("&"):
                i += 1
                stmt.append(lines[i])
            text = " ".join(strip_comment(s).rstrip("&").strip() for s in stmt)
            text = re.sub(r"^\s*&", "", text)
            body = text[text.index("(") + 1: text.rindex(")")].replace("& ", " ").replace(" &", " ")
            args = split_args(body)
            name = m.group(2)
            count[name] = count.get(name, 0) + 1
            site = f"{name}_{count[name]}"
            if first:
                ind = m.group(1)
                out += [f"{ind}call capture_site('module_state', mgncol, nlev, 'in')\n",
                        f"{ind}call capture_v1_state()\n", f"{ind}call capture_utils_state()\n"]
                first = False
            out += captures(m.group(1), site, args, "in")
            out += stmt
            out += captures(m.group(1), site, args, "out")
        else:
            out.append(line)
        i += 1
    if not used:
        raise SystemExit("no 'use micro_pumas_utils' line found")
    text = "".join(out)
    j = re.search(r"^end module micro_pumas_v1", text, re.MULTILINE).start()
    text = text[:j] + state_subroutine("v1", v1_names, []) + "\n" + text[j:]
    dst.write_text(text)
    print(f"{sum(count.values())} call sites in {len(count)} routines -> {dst.resolve()}")


if __name__ == "__main__":
    main()
