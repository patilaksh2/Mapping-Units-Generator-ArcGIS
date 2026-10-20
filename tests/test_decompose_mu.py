"""Tests for decompose_mu() in mu_generator_arcgis.py.

Ported from the QGIS repo's tests/test_decompose_mu.py. The function
(and the module-level _MU_PATTERN it uses) is extracted from the script
source with ast and exec'd in isolation, so these tests run in plain
Python - no ArcGIS Pro install needed. Run directly
(`python3 tests/test_decompose_mu.py`) or under pytest.
"""
import ast
import re
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "mu_generator_arcgis.py"


def _load_decompose_mu():
    source = SCRIPT.read_text(encoding="utf-8")
    tree = ast.parse(source)
    namespace = {"re": re}
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
                isinstance(target, ast.Name) and target.id == "_MU_PATTERN"
                for target in node.targets):
            exec(ast.get_source_segment(source, node), namespace)
        elif isinstance(node, ast.FunctionDef) and node.name == "decompose_mu":
            exec(ast.get_source_segment(source, node), namespace)
    if "decompose_mu" not in namespace:
        raise AssertionError("decompose_mu not found in mu_generator_arcgis.py")
    return namespace["decompose_mu"]


decompose_mu = _load_decompose_mu()


def test_full_codes():
    assert decompose_mu("AD3rB2") == {
        "SeriesCode": "AD", "DepthNum": "3", "SurTexCode": "r",
        "SlopeCode": "B", "ErosionNum": "2",
    }
    assert decompose_mu("CH5sA1") == {
        "SeriesCode": "CH", "DepthNum": "5", "SurTexCode": "s",
        "SlopeCode": "A", "ErosionNum": "1",
    }
    assert decompose_mu("MR4tC3") == {
        "SeriesCode": "MR", "DepthNum": "4", "SurTexCode": "t",
        "SlopeCode": "C", "ErosionNum": "3",
    }


def test_empty_input():
    assert decompose_mu(None) is None
    assert decompose_mu("") is None
    assert decompose_mu("   ") is None


def test_non_matching_codes():
    assert decompose_mu("UNK") is None        # no grammar structure at all
    assert decompose_mu("AD3rB") is None      # too short (missing erosion)
    assert decompose_mu("AD33rB2") is None    # multi-digit depth not in grammar
    assert decompose_mu("A13rB2") is None     # series must be exactly 2 letters
    assert decompose_mu("AD3rB22") is None    # too long


if __name__ == "__main__":
    test_full_codes()
    test_empty_input()
    test_non_matching_codes()
    print("test_decompose_mu: all tests passed")
