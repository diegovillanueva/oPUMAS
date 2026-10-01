"""🧪🔤 Tests of translate.py on small Fortran snippets and of the real dictionary."""
from pathlib import Path

import check
from translate import Maps, translate_text, unit_suffix, wrap

ROWS = [
    {"id": "pra", "scope": "global", "units": "kg kg-1 s-1", "PUMAS": "pra", "CF-CCPP": "local_cldliq_mass_tend_accretion_by_rain"},
    {"id": "rho", "scope": "global", "units": "kg m-3", "PUMAS": "rho", "CF-CCPP": "air_density"},
    {"id": "sc", "scope": "global", "units": "1", "PUMAS": "sc", "CF-CCPP": "schmidt_number"},
    {"id": "t%pratot", "scope": "type:proc_rates_type", "units": "kg kg-1 s-1", "PUMAS": "pratot",
     "CF-CCPP": "cldliq_mass_tend_accretion_by_rain"},
    {"id": "t%umr", "scope": "type:proc_rates_type", "units": "m s-1", "PUMAS": "umr", "CF-CCPP": "rain_mass_weighted_fall_speed"},
    {"id": "f", "scope": "file:a.F90", "units": "", "PUMAS": "tmp", "CF-CCPP": "scratch"},
]
FWD, BACK = Maps(ROWS, "PUMAS", "CF-CCPP"), Maps(ROWS, "CF-CCPP", "PUMAS")
PRA = "local_cldliq_mass_tend_accretion_by_rain_kg_per_kg_per_s"


def fwd(text: str, name: str = "a.F90", mode: str = "full") -> str:
    return translate_text(text, FWD, name, mode)


def test_unit_suffix() -> None:
    assert unit_suffix("kg kg-1 s-1") == "_kg_per_kg_per_s"
    assert unit_suffix("kg-1 s-1") == "_per_kg_per_s"
    assert unit_suffix("m-4") == "_per_m4"
    assert unit_suffix("1") == "" and unit_suffix("") == ""


def test_code_is_translated_comments_and_strings_are_not() -> None:
    src = "  pra(i) = rho(i)*2._r8 ! pra uses rho\n  call f('pra rho', \"it's pra\")\n"
    assert fwd(src) == f"  {PRA}(i) = air_density_kg_per_m3(i)*2._r8 ! pra uses rho\n  call f('pra rho', \"it's pra\")\n"


def test_whole_identifier_only_and_kind_suffix() -> None:
    src = "  prax = pra_x + rhow + 1.0_rho + sc\n"
    assert fwd(src) == "  prax = pra_x + rhow + 1.0_rho + schmidt_number\n"


def test_directives_and_preprocessor() -> None:
    assert fwd("!$acc parallel copyin(pra,rho)\n") == f"!$acc parallel copyin({PRA},air_density_kg_per_m3)\n"
    assert fwd("#ifdef pra\n") == "#ifdef pra\n"
    assert fwd("! pra\n") == "! pra\n"


def test_components_and_type_block() -> None:
    src = ("type, public :: proc_rates_type\n  real(r8), allocatable :: pratot(:,:), umr(:,:)\nend type proc_rates_type\n"
           "x = proc_rates%pratot(i,k) + umr + pratot\n")
    out = fwd(src)
    assert "allocatable :: cldliq_mass_tend_accretion_by_rain_kg_per_kg_per_s(:,:), rain_mass_weighted_fall_speed_m_per_s(:,:)" in out
    assert out.endswith("x = proc_rates%cldliq_mass_tend_accretion_by_rain_kg_per_kg_per_s(i,k) + umr + pratot\n")


def test_components_mode_leaves_plain_names() -> None:
    assert fwd("pra = this % pratot\n", mode="components") == "pra = this % cldliq_mass_tend_accretion_by_rain_kg_per_kg_per_s\n"


def test_file_scope_and_upper_case() -> None:
    assert fwd("tmp = RHO\n") == "scratch = AIR_DENSITY_KG_PER_M3\n"
    assert fwd("tmp = RHO\n", name="b.F90") == "tmp = AIR_DENSITY_KG_PER_M3\n"


def test_round_trip_and_no_units() -> None:
    src = "  pra(i) = rho(i)*sc ! c\n!$acc loop private(pra)\nx = p%umr\n"
    assert translate_text(fwd(src), BACK, "a.F90") == src
    plain = Maps(ROWS, "PUMAS", "CF-CCPP", units=False)
    assert translate_text("x = rho\n", plain, "a.F90") == "x = air_density\n"


def test_wrap_keeps_lines_short() -> None:
    line = "  x = " + " + ".join(f"very_long_name_number_{i}(i,k)" for i in range(12)) + "\n"
    out = wrap(line)
    assert all(len(part) <= 132 for part in out.splitlines())
    assert out.replace(" &\n", "").replace(" ", "") == line.replace(" ", "")


def test_real_dictionary_passes_static_checks() -> None:
    assert check.run_checks(Path(__file__).resolve().parent.parent) == []


def test_components_of_other_types_keep_their_names() -> None:
    src = "type :: MGHydrometeorProps\n  real(r8) :: rho\nend type MGHydrometeorProps\nx = props%rho*rho\nselect type (a)\ntype is (real)\n  y = rho\nend select\n"
    assert fwd(src) == src.replace("%rho*rho", "%rho*air_density_kg_per_m3").replace("y = rho", "y = air_density_kg_per_m3")
