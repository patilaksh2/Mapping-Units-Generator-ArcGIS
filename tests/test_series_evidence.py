"""Tests for the SLUSI two-stage driver core in mu_generator_arcgis.py:
build_series_evidence(), rank_series_candidates(),
primary_driver_confidence() and the Stage 2 completion helpers.

Ported from the QGIS repo's tests/test_series_evidence.py. The pure
functions are extracted from the script source with ast and exec'd in
isolation, so these tests run in plain Python - no ArcGIS Pro install
needed. Run directly (`python3 tests/test_series_evidence.py`) or under
pytest.
"""
import ast
from collections import Counter, defaultdict
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "mu_generator_arcgis.py"

PURE_FUNCTIONS = {
    "_modal_value",
    "build_series_evidence",
    "rank_series_candidates",
    "primary_driver_confidence",
    "complete_mu_components",
    "assemble_mu_code",
    "apply_stage2",
}


def _load_pure_functions():
    source = SCRIPT.read_text(encoding="utf-8")
    tree = ast.parse(source)
    namespace = {"Counter": Counter, "defaultdict": defaultdict}
    found = set()
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name in PURE_FUNCTIONS:
            exec(ast.get_source_segment(source, node), namespace)
            found.add(node.name)
    missing = PURE_FUNCTIONS - found
    assert not missing, f"functions not found in mu_generator_arcgis.py: {missing}"
    return namespace


_ns = _load_pure_functions()
build_series_evidence = _ns["build_series_evidence"]
rank_series_candidates = _ns["rank_series_candidates"]
primary_driver_confidence = _ns["primary_driver_confidence"]
apply_stage2 = _ns["apply_stage2"]


def _row(series, geology, physio, depth, surtex="r", slope="B", erosion="2"):
    return {"series": series, "geology": geology, "physiography": physio,
            "depth": depth, "surtex": surtex, "slope": slope,
            "erosion": erosion}


def _evidence():
    rows = []
    # Series AD - Basalt; Plateau (depths 3,3,3,4) + Piedmont (depth 3).
    rows += [_row("AD", "Basalt", "Plateau", "3"),
             _row("AD", "Basalt", "Plateau", "3"),
             _row("AD", "Basalt", "Plateau", "3"),
             _row("AD", "Basalt", "Plateau", "4", slope="C"),
             _row("AD", "Basalt", "Piedmont", "3", erosion="1")]
    # Series CH - Basalt; Plateau only (depths 5,5,5,5,3).
    rows += [_row("CH", "Basalt", "Plateau", "5", surtex="s", slope="A",
                  erosion="1"),
             _row("CH", "Basalt", "Plateau", "5", surtex="s", slope="A",
                  erosion="1"),
             _row("CH", "Basalt", "Plateau", "5", surtex="s", slope="A",
                  erosion="1"),
             _row("CH", "Basalt", "Plateau", "5", surtex="s", slope="A",
                  erosion="1"),
             _row("CH", "Basalt", "Plateau", "3", surtex="s", slope="B",
                  erosion="1")]
    return build_series_evidence(rows)


def test_single_candidate_match():
    # Only AD occurs on Basalt + Piedmont, so it is the sole candidate.
    candidates = rank_series_candidates(_evidence(), "Basalt", "Piedmont")
    assert [c["series"] for c in candidates] == ["AD"]


def test_depth_agreement_breaks_tie():
    evidence = _evidence()
    # No known depth: AD and CH tie on modal-depth agreement (4/5 each)
    # and on support (5 each) - the deterministic order puts AD first.
    ranked = rank_series_candidates(evidence, "Basalt", "Plateau")
    assert [c["series"] for c in ranked] == ["AD", "CH"]
    assert ranked[0]["depth_agreement"] == 0.8
    assert ranked[1]["depth_agreement"] == 0.8
    # With a known depth class of 5, CH (4 of its 5 polygons) beats AD
    # (none of its polygons) - depth agreement breaks the tie.
    ranked = rank_series_candidates(evidence, "Basalt", "Plateau",
                                    known_depth="5")
    assert ranked[0]["series"] == "CH"
    assert ranked[0]["depth_agreement"] == 0.8
    assert ranked[1]["series"] == "AD"
    assert ranked[1]["depth_agreement"] == 0.0


def test_empty_candidates():
    evidence = _evidence()
    # Unseen geology, and unseen physiography: no candidates at all.
    assert rank_series_candidates(evidence, "Granite", "Plateau") == []
    assert rank_series_candidates(evidence, "Basalt", "Valley floor") == []
    assert rank_series_candidates({}, "Basalt", "Plateau") == []


def test_confidence_formula():
    assert primary_driver_confidence(1.0, 5, 10) == 90.0   # 50 + 30 + 10
    assert primary_driver_confidence(0.0, 0, 10) == 50.0    # base only
    assert primary_driver_confidence(1.0, 10, 10) == 100.0  # clamped


def test_stage2_completion():
    evidence = _evidence()
    # Modal components on the Plateau: depth 3, tex r, slope B, erosion 2.
    comps, mu_new, note = apply_stage2(evidence, "AD", "Plateau")
    assert mu_new == "AD3rB2"
    assert note == ""
    # Same-physiography preference: AD's only Piedmont polygon has
    # erosion 1, so the Piedmont MU takes erosion 1, not the series-wide
    # modal erosion 2.
    comps, mu_new, note = apply_stage2(evidence, "AD", "Piedmont")
    assert comps["ErosionNum"] == "1"
    assert mu_new == "AD3rB1"
    # A series with no mapped polygons: nothing is invented, everything
    # is reported missing instead.
    comps, mu_new, note = apply_stage2(evidence, "ZZ", "Plateau")
    assert mu_new == "ZZ"
    assert note.startswith("MU incomplete:")


if __name__ == "__main__":
    test_single_candidate_match()
    test_depth_agreement_breaks_tie()
    test_empty_candidates()
    test_confidence_formula()
    test_stage2_completion()
    print("test_series_evidence: all tests passed")
