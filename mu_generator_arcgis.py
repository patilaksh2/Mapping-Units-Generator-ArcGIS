r"""
Soil Mapping Unit (MU) Generator - ArcGIS Pro (arcpy) port
============================================================
ArcGIS Pro Python script that fills in blank Mapping Unit (MU) codes on a
soil survey polygon layer, following the SLUSI two-stage driver method:

    Tier 0 : IDW interpolation from ground-truth points   (highest authority)
    Stage 1: confirm the soil SERIES from the PRIMARY drivers (geology,
             physiography, depth) - a series-evidence table built from
             the already-mapped polygons; the trained classifier (Random
             Forest, or k-NN fallback) works within that evidence
    Stage 2: complete the MAPPING UNIT from the SECONDARY drivers
             (surface texture, slope, erosion) using the series' modal
             components

    Fallbacks when Stage 1 finds no primary-driver evidence, in order:
    classifier alone, rule-based keyword matching, nearest mapped
    polygon, flagged default.

This is the ArcGIS port of the QGIS "Mapping-Units-Generator"
(mu_generator.py). Same tier design, same CONFIGURATION concept and the
same output fields (SeriesCode, DepthNum, SurTexCode, SlopeCode,
ErosionNum, MU_New, MU_Source, MU_Conf, MU_Method, MU_Flag); the QGIS /
Processing calls are replaced by their arcpy equivalents (Project,
Repair Geometry, Spatial Join with LARGEST_OVERLAP, arcpy.da cursors).

The port also carries FIXED logic rather than the QGIS original's known
bugs:
  * decompose_mu() is a strict regex parser (Series + Depth + Surface
    texture + Slope + Erosion), not positional letter/digit slicing
    (which corrupted the slope code);
  * EVERY input layer (soil, physiography, geology and the ground-truth
    points) is reprojected to the target CRS, and CRSs are asserted
    equal before any join or distance calculation;
  * classifier coordinates are z-score normalised, so the k-NN fallback
    is not dominated by raw map coordinates;
  * the adjacency step only FLAGS conflicting polygons in MU_Flag - it
    never invents or rewrites MU codes;
  * Tier 0 (IDW) results are decomposed into components exactly like
    existing codes, instead of being stored inconsistently.

Requirements
------------
- ArcGIS Pro 3.x (arcpy). No Spatial Analyst licence is needed - the
  IDW is implemented manually in pure Python.
- numpy + scikit-learn are OPTIONAL: ArcGIS Pro's Python ships numpy;
  if scikit-learn is importable it is used for the Random Forest tier,
  otherwise the script falls back to a pure-Python weighted k-NN (the
  same graceful-degradation strategy as the QGIS original).

How to run
----------
Option A - ArcGIS Pro Python window:
    exec(open(r"C:\path\to\mu_generator_arcgis.py").read())

Option B - register the file as a script tool in a custom toolbox and
    run it from the Geoprocessing pane (it takes no parameters - edit
    the CONFIGURATION block below first).

Option C - standalone with ArcGIS Pro's Python:
    "C:\Program Files\ArcGIS\Pro\bin\Python\Scripts\propy.bat" mu_generator_arcgis.py

Everything you need to change lives in the CONFIGURATION block below.
The logic beneath it should not need editing for normal use.
"""

import math
import os
import random
import re
import statistics
from collections import Counter, defaultdict

import arcpy

# Optional dependencies - the script works without them.
try:
    import numpy as np
    HAVE_NUMPY = True
except ImportError:
    HAVE_NUMPY = False

try:
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import StratifiedKFold, train_test_split
    from sklearn.metrics import confusion_matrix, accuracy_score
    HAVE_SKLEARN = True
except ImportError:
    HAVE_SKLEARN = False


# =============================================================================
# CONFIGURATION - edit everything in this block for your own project.
# Nothing below this block should need to change for normal use.
# =============================================================================

# ---- Paths ------------------------------------------------------------------
# Inputs may be shapefiles or feature classes inside a file geodatabase.
# Prefer absolute paths: the ArcGIS Pro Python window's working folder is
# not guaranteed to be this script's folder.
SOIL_LAYER    = r"./data/prefieldpoly.shp"
PHYSIO_LAYER  = r"./data/physiography.shp"
GEOLOGY_LAYER = r"./data/geology.shp"
OUTPUT_FC     = r"./output/MappingUnit_Output.shp"

# Ground-truth point layers used for IDW interpolation (Tier 0 evidence).
# Leave this list empty if you don't have field point data.
GROUND_TRUTH_POINT_LAYERS = [
    # r"./data/field_points_2026.shp",
]

# ---- Field names --------------------------------------------------------
PHYSIO_NAME_FLD  = "legend_sho"   # name/legend field on the physiography layer
GEOLOGY_NAME_FLD = "lithologic"   # lithology field on the geology layer
EXISTING_MU_FLD  = "MU"           # field on the soil layer holding known MU codes
POINT_MU_FLD     = "MU"           # field on the ground-truth points holding MU codes
DEPTH_FLD        = ""             # OPTIONAL: field on the SOIL layer holding a known
                                  # depth class per polygon (a primary driver - feeds
                                  # Stage 1 depth agreement and Stage 2). "" = disabled.
SLOPE_FLD        = ""             # OPTIONAL: field on the SOIL layer holding a slope
                                  # code per polygon (used by Stage 2 before the
                                  # series' modal slope). "" = disabled.

# ---- CRS ------------------------------------------------------------------
# EPSG well-known ID (WKID) of the projected CRS for your study area
# (32643 = WGS 84 / UTM zone 43N). Using a projected CRS matters here
# because distances (IDW search radius, nearest-neighbour fallback) are
# computed in map units.
TARGET_CRS_WKID = 32643

# ---- Non-soil polygon codes --------------------------------------------
# MU codes (or existing values in EXISTING_MU_FLD) that represent non-soil
# land uses and should be excluded from classification entirely.
NON_SOIL_CODES = {"WB", "HAB", "ROAD", "QUARRY"}  # water body, habitation, etc.

# ---- Series lookup table -------------------------------------------------
# This is the single most project-specific part of the script. It maps a
# short series code to attributes used for validation/reporting. Rebuild
# this table for your own survey area - do not assume the example values
# below apply to your district.
#
#   code : (typical_depth_class, description)
SERIES_LOOKUP = {
    # "AD": (3, "example series - replace with your own"),
    # "CH": (5, "example series - replace with your own"),
}

# ---- Tuning parameters ---------------------------------------------------
IDW_SEARCH_RADIUS_M   = 500.0   # metres
IDW_POWER             = 2.0     # inverse-distance weighting power
IDW_CONFIDENCE_MIN    = 80.0    # % - auto-accept IDW result above this
CLASSIFIER_CONF_MIN   = 60.0    # % - accept classifier result above this
KNN_NEIGHBOURS        = 5
RANDOM_SEED           = 42

random.seed(RANDOM_SEED)
if HAVE_NUMPY:
    np.random.seed(RANDOM_SEED)


# =============================================================================
# PIPELINE LOGIC - shouldn't normally need edits below this line.
# =============================================================================

def log(step, total, msg):
    print(f"[{step}/{total}] {msg}")


def safe_str(value, max_len):
    """Defensive truncation before writing a string attribute (fixes the
    field-width save failure some shapefile drivers hit on long strings)."""
    if value is None:
        return ""
    s = str(value)
    return s[:max_len] if len(s) > max_len else s


def _clean_code_value(value):
    """Normalise an optional-field value (DEPTH_FLD / SLOPE_FLD) to a
    clean code string; "" for None/empty."""
    if value is None:
        return ""
    if isinstance(value, float) and value.is_integer():
        value = int(value)
    return str(value).strip()


def check_exists(path, label):
    if not arcpy.Exists(path):
        raise RuntimeError(f"{label} not found: {path}")


# ---- Step 0: geometry fix + CRS alignment ----------------------------------

def project_to_target(path, out_name, target_sr):
    """Copy/project a dataset into the in_memory workspace in the target
    CRS. Returns the in_memory path."""
    out_fc = os.path.join("in_memory", out_name)
    if arcpy.Exists(out_fc):
        arcpy.management.Delete(out_fc)
    in_sr = arcpy.Describe(path).spatialReference
    if in_sr is not None and in_sr.factoryCode == TARGET_CRS_WKID:
        arcpy.management.CopyFeatures(path, out_fc)
    else:
        arcpy.management.Project(path, out_fc, target_sr)
    return out_fc


def assert_crs(fc, label):
    sr = arcpy.Describe(fc).spatialReference
    code = sr.factoryCode if sr is not None else None
    if code != TARGET_CRS_WKID:
        raise RuntimeError(
            f"CRS mismatch for {label}: expected EPSG:{TARGET_CRS_WKID}, "
            f"got {sr.name if sr is not None else 'unknown'} (code {code})"
        )


def prepare_polygon_layer(path, out_name, target_sr):
    """Project a polygon layer to the target CRS (if needed), repair its
    geometry, and verify the CRS. FIX vs the QGIS original: this is
    applied to ALL polygon layers, not just the soil layer - joining or
    measuring across mismatched CRSs otherwise fails silently."""
    fc = project_to_target(path, out_name, target_sr)
    arcpy.management.RepairGeometry(fc, "DELETE_NULL")
    assert_crs(fc, out_name)
    return fc


# ---- Steps 3-4: spatial joins (largest overlap) ----------------------------

def spatial_join_largest_overlap(target_fc, join_fc, join_field, out_fc):
    """One-to-one spatial join that keeps, for every target feature, the
    attributes of the join feature with the LARGEST OVERLAP - the arcpy
    equivalent of QGIS joinattributesbylocation with METHOD=2. (The QGIS
    original shipped METHOD=0, a one-to-many join that duplicated
    polygons; do not 'fix' this back.)

    Only `join_field` is carried over from the join features, under its
    original name - so the target must not already have a field with
    that name (true for the physiography/geology name fields here)."""
    if arcpy.Exists(out_fc):
        arcpy.management.Delete(out_fc)
    fms = arcpy.FieldMappings()
    fms.addTable(target_fc)
    fm = arcpy.FieldMap()
    fm.addInputField(join_fc, join_field)
    fms.addFieldMap(fm)
    arcpy.analysis.SpatialJoin(
        target_fc, join_fc, out_fc,
        join_operation="JOIN_ONE_TO_ONE",
        join_type="KEEP_ALL",
        field_mapping=fms,
        match_option="LARGEST_OVERLAP",
    )
    return out_fc


def drop_join_scaffolding(fc):
    """Remove the bookkeeping fields Spatial Join adds (Join_Count,
    TARGET_FID) so the output schema stays clean."""
    names = {f.name for f in arcpy.ListFields(fc)}
    junk = [n for n in ("Join_Count", "TARGET_FID",
                        "Join_Count_1", "TARGET_FID_1") if n in names]
    if junk:
        try:
            arcpy.management.DeleteField(fc, junk)
        except Exception:
            pass  # cosmetic only - never fail the run over this


# ---- Step 5: decompose existing MU codes / build training set -------------

# MU code grammar: SeriesCode (2 letters) + DepthNum (1 digit) +
# SurfaceTexCode (1 letter) + SlopeCode (1 letter) + ErosionNum (1 digit),
# e.g. "AD3rB2".
_MU_PATTERN = re.compile(r"^([A-Za-z]{2})(\d)([A-Za-z])([A-Za-z])(\d)$")


def decompose_mu(mu_code):
    """
    Split a composite MU code into its components, or return None if the
    code does not match the grammar above.

    FIX vs the QGIS original: the original sliced letters/digits
    positionally and ended up writing the surface-texture character into
    SlopeCode as well (AD3rB2 -> SlopeCode 'r' instead of 'B'). This
    regex parser cannot do that, and it rejects malformed codes instead
    of silently mis-parsing them.
    """
    if not mu_code:
        return None
    m = _MU_PATTERN.match(str(mu_code).strip())
    if m is None:
        return None
    series, depth, surtex, slope, erosion = m.groups()
    return {
        "SeriesCode": series,
        "DepthNum": depth,
        "SurTexCode": surtex,
        "SlopeCode": slope,
        "ErosionNum": erosion,
    }


def read_polygons(fc):
    """Read the (joined) soil layer into a list of record dicts with a
    SearchCursor: OID, existing MU value, physiography/geology values,
    centroid coordinates and the geometry (kept for the adjacency
    check in Step 8)."""
    oid_fld = arcpy.Describe(fc).OIDFieldName
    field_names = {f.name for f in arcpy.ListFields(fc)}
    cursor_fields = [oid_fld, EXISTING_MU_FLD, "SHAPE@"]
    if PHYSIO_NAME_FLD in field_names:
        cursor_fields.append(PHYSIO_NAME_FLD)
    if GEOLOGY_NAME_FLD in field_names:
        cursor_fields.append(GEOLOGY_NAME_FLD)
    if DEPTH_FLD and DEPTH_FLD in field_names:
        cursor_fields.append(DEPTH_FLD)
    if SLOPE_FLD and SLOPE_FLD in field_names:
        cursor_fields.append(SLOPE_FLD)

    records = []
    with arcpy.da.SearchCursor(fc, cursor_fields) as cur:
        for row in cur:
            data = dict(zip(cursor_fields, row))
            geom = data["SHAPE@"]
            centroid = geom.centroid if geom is not None else None
            records.append({
                "oid": data[oid_fld],
                "mu": data.get(EXISTING_MU_FLD),
                "physio": data.get(PHYSIO_NAME_FLD) or "",
                "geology": data.get(GEOLOGY_NAME_FLD) or "",
                "depth": _clean_code_value(data.get(DEPTH_FLD)) if DEPTH_FLD else "",
                "slope": _clean_code_value(data.get(SLOPE_FLD)) if SLOPE_FLD else "",
                "x": centroid.X if centroid is not None else 0.0,
                "y": centroid.Y if centroid is not None else 0.0,
                "geom": geom,
            })
    return records


def build_training_set(records, non_soil_codes):
    rows, labels = [], []
    for rec in records:
        mu = rec["mu"]
        if not mu or str(mu).strip().upper() in non_soil_codes:
            continue
        comps = decompose_mu(mu)
        if not comps or not comps["SeriesCode"]:
            continue
        rows.append({
            "physio": rec["physio"],
            "geology": rec["geology"],
            "x": rec["x"],
            "y": rec["y"],
            "fid": rec["oid"],
        })
        labels.append(comps["SeriesCode"])
    return rows, labels


def build_evidence_rows(records, non_soil_codes):
    """One evidence dict per already-mapped record (see
    build_series_evidence for the dict shape)."""
    rows = []
    for rec in records:
        mu = rec["mu"]
        if not mu or str(mu).strip().upper() in non_soil_codes:
            continue
        comps = decompose_mu(mu)
        if not comps:
            continue
        rows.append({
            "series": comps["SeriesCode"],
            "geology": rec["geology"],
            "physiography": rec["physio"],
            "depth": comps["DepthNum"],
            "surtex": comps["SurTexCode"],
            "slope": comps["SlopeCode"],
            "erosion": comps["ErosionNum"],
        })
    return rows


# ---- Series-evidence table: the SLUSI two-stage core -----------------------
# Stage 1 confirms the soil SERIES from the primary drivers (geology,
# physiography, depth); Stage 2 completes the MAPPING UNIT from the
# secondary drivers (surface texture, slope, erosion). These functions
# are pure (plain dicts/lists in and out) and identical in behaviour to
# the QGIS original, where they are unit-tested outside QGIS
# (tests/test_series_evidence.py in that repo).

def _modal_value(counter):
    """Most frequent key of a Counter; ties break on the value itself so
    results are deterministic. None for an empty counter."""
    if not counter:
        return None
    return sorted(counter.items(), key=lambda kv: (-kv[1], str(kv[0])))[0][0]


def build_series_evidence(evidence_rows):
    """
    Build the series-evidence table used by Stage 1 and Stage 2.

    Pure function: `evidence_rows` is a list of dicts with keys
    series / geology / physiography / depth / surtex / slope / erosion
    (one dict per already-mapped polygon, components from decompose_mu).
    Returns {series: {...}} holding the geology value set, the
    physiography value set, the depth-class Counter, the support
    (polygon count), the modal depth class, and the modal MU components
    overall and per physiography.
    """
    acc = {}
    for row in evidence_rows:
        entry = acc.setdefault(row["series"], {
            "geologies": set(),
            "physiographies": set(),
            "support": 0,
            "depth_counts": Counter(),
            "_comp_all": defaultdict(Counter),
            "_comp_physio": defaultdict(lambda: defaultdict(Counter)),
        })
        entry["geologies"].add(row["geology"])
        entry["physiographies"].add(row["physiography"])
        entry["support"] += 1
        entry["depth_counts"][row["depth"]] += 1
        for key, comp in (("depth", "DepthNum"), ("surtex", "SurTexCode"),
                          ("slope", "SlopeCode"), ("erosion", "ErosionNum")):
            entry["_comp_all"][comp][row[key]] += 1
            entry["_comp_physio"][row["physiography"]][comp][row[key]] += 1

    evidence = {}
    for series, entry in acc.items():
        evidence[series] = {
            "geologies": entry["geologies"],
            "physiographies": entry["physiographies"],
            "support": entry["support"],
            "depth_counts": entry["depth_counts"],
            "depth_modal": _modal_value(entry["depth_counts"]),
            "components": {comp: _modal_value(counter)
                           for comp, counter in entry["_comp_all"].items()},
            "components_by_physio": {
                physio: {comp: _modal_value(counter)
                         for comp, counter in comps.items()}
                for physio, comps in entry["_comp_physio"].items()
            },
        }
    return evidence


def rank_series_candidates(evidence, geology, physiography,
                           known_depth=None, classifier_probs=None):
    """
    Stage 1 candidate ranking (pure function).

    Candidates are the series whose observed geology set contains
    `geology` AND whose observed physiography set contains
    `physiography`. Ranked by, in order:
      1. depth agreement - if `known_depth` is given, the share of the
         candidate's mapped polygons carrying that depth class;
         otherwise the share carrying the candidate's modal depth class;
      2. classifier probability for the series (`classifier_probs` maps
         series -> probability %; a series missing from it counts as -1);
      3. support (number of mapped polygons);
      4. series code, ascending (deterministic tie-break).
    Returns a list of candidate dicts, best first; empty when no series
    matches both primary drivers.
    """
    probs = classifier_probs or {}
    candidates = []
    for series, entry in evidence.items():
        if geology not in entry["geologies"]:
            continue
        if physiography not in entry["physiographies"]:
            continue
        support = entry["support"]
        if known_depth:
            agreed_depth = known_depth
        else:
            agreed_depth = entry["depth_modal"]
        depth_agreement = (entry["depth_counts"].get(agreed_depth, 0) / support
                           if support else 0.0)
        candidates.append({
            "series": series,
            "depth_agreement": depth_agreement,
            "agreed_depth": agreed_depth,
            "modal_depth": entry["depth_modal"],
            "classifier_prob": probs.get(series),
            "support": support,
        })
    candidates.sort(key=lambda c: (
        -c["depth_agreement"],
        -(c["classifier_prob"] if c["classifier_prob"] is not None else -1.0),
        -c["support"],
        c["series"],
    ))
    return candidates


def primary_driver_confidence(depth_agreement, support, total_support):
    """
    Confidence (%) for a primary_driver_match assignment:

        50 + 30 * depth-agreement share + 20 * relative support

    where relative support = the candidate's polygon count / the total
    mapped-polygon count across all series. Clamped to 0-100. The base
    50 reflects that both primary drivers (geology + physiography)
    already matched exactly; depth agreement and prevalence refine it.
    """
    relative_support = (support / total_support) if total_support else 0.0
    score = 50.0 + 30.0 * depth_agreement + 20.0 * relative_support
    return max(0.0, min(100.0, score))


def complete_mu_components(entry, physiography, depth_num=None, slope_code=None):
    """
    Stage 2 component completion (pure function).

    `entry` is one series' evidence-table entry (or None when the series
    has no mapped polygons). Modal component values are taken from the
    series' polygons in the SAME physiography first, falling back to the
    series-wide modal values. `depth_num` (the polygon's agreed depth,
    e.g. from DEPTH_FLD) and `slope_code` (from SLOPE_FLD) override the
    modal values when provided. Returns (components, missing):
    components maps DepthNum / SurTexCode / SlopeCode / ErosionNum to a
    value ("" when unknown) and `missing` lists the components left
    blank - values are never invented.
    """
    overall = entry["components"] if entry else {}
    by_physio = (entry["components_by_physio"].get(physiography, {})
                 if entry else {})
    comps = {}
    for comp in ("DepthNum", "SurTexCode", "SlopeCode", "ErosionNum"):
        comps[comp] = by_physio.get(comp) or overall.get(comp) or ""
    if depth_num:
        comps["DepthNum"] = depth_num
    if slope_code:
        comps["SlopeCode"] = slope_code
    missing = [comp for comp, val in comps.items() if not val]
    return comps, missing


def assemble_mu_code(series, comps):
    """Join series + components into the SLUSI-format MU code. Blank
    components contribute nothing (they are reported via MU_Flag)."""
    return (series + comps.get("DepthNum", "") + comps.get("SurTexCode", "")
            + comps.get("SlopeCode", "") + comps.get("ErosionNum", ""))


def apply_stage2(evidence, series, physiography, known_depth="", slope_code=""):
    """
    Stage 2: complete the mapping unit for a series assigned in Stage 1.
    Returns (components dict including SeriesCode, assembled MU code,
    incomplete-components note, "" when the MU is complete).
    """
    comp_vals, missing = complete_mu_components(
        evidence.get(series), physiography,
        depth_num=known_depth or None, slope_code=slope_code or None)
    comps = {"SeriesCode": series}
    comps.update(comp_vals)
    note = ""
    if missing:
        note = "MU incomplete: " + ", ".join(missing) + " unknown"
    return comps, assemble_mu_code(series, comp_vals), note


# ---- Feature encoding for the classifier -----------------------------------

def coordinate_stats(rows):
    """Mean/std of the training centroids, for z-score normalisation."""
    def _stats(vals):
        mean = statistics.mean(vals)
        std = statistics.pstdev(vals) if len(vals) > 1 else 0.0
        return mean, (std if std > 0 else 1.0)
    return _stats([r["x"] for r in rows]), _stats([r["y"] for r in rows])


def encode_features(rows, physio_vocab, geology_vocab, stats):
    """One-hot encode physiography/geology + z-score normalised
    coordinates. FIX vs the QGIS original, whose docstring claimed
    normalised coordinates but appended raw ones: with raw UTM values
    (~1e5-1e6) beside 0/1 one-hot features, Euclidean distance in the
    k-NN fallback was dominated by position and physiography/geology
    barely mattered."""
    (mx, sx), (my, sy) = stats
    encoded = []
    for r in rows:
        vec = [1.0 if r["physio"] == p else 0.0 for p in physio_vocab]
        vec += [1.0 if r["geology"] == g else 0.0 for g in geology_vocab]
        vec += [(r["x"] - mx) / sx, (r["y"] - my) / sy]
        encoded.append(vec)
    return encoded


# ---- Step 6: train & validate classifier -----------------------------------

def train_classifier(rows, labels):
    physio_vocab = sorted({r["physio"] for r in rows})
    geology_vocab = sorted({r["geology"] for r in rows})
    stats = coordinate_stats(rows)
    X = encode_features(rows, physio_vocab, geology_vocab, stats)
    y = labels

    class_counts = Counter(y)
    print(f"  Training set: {len(y)} samples, {len(class_counts)} series classes")
    print(f"  Class distribution: {dict(class_counts)}")

    if HAVE_SKLEARN and HAVE_NUMPY and len(y) >= 10:
        X_arr = np.array(X)
        y_arr = np.array(y)

        # Stratified holdout
        try:
            X_train, X_test, y_train, y_test = train_test_split(
                X_arr, y_arr, test_size=0.2, stratify=y_arr, random_state=RANDOM_SEED
            )
        except ValueError:
            # some class has too few members to stratify - fall back to plain split
            X_train, X_test, y_train, y_test = train_test_split(
                X_arr, y_arr, test_size=0.2, random_state=RANDOM_SEED
            )

        # NOTE: trains on the true empirical class prior (no artificial
        # rebalancing), which avoids over-predicting rare classes relative
        # to their real prevalence.
        clf = RandomForestClassifier(
            n_estimators=300, random_state=RANDOM_SEED, class_weight=None
        )
        clf.fit(X_train, y_train)
        holdout_acc = accuracy_score(y_test, clf.predict(X_test))
        print(f"  Holdout accuracy: {holdout_acc * 100:.1f}%")

        # k-fold CV
        try:
            skf = StratifiedKFold(n_splits=min(5, min(class_counts.values())),
                                   shuffle=True, random_state=RANDOM_SEED)
            cv_scores = []
            for train_idx, test_idx in skf.split(X_arr, y_arr):
                fold_clf = RandomForestClassifier(
                    n_estimators=300, random_state=RANDOM_SEED
                )
                fold_clf.fit(X_arr[train_idx], y_arr[train_idx])
                cv_scores.append(fold_clf.score(X_arr[test_idx], y_arr[test_idx]))
            print(f"  k-fold CV accuracy: {statistics.mean(cv_scores) * 100:.1f}% "
                  f"(+/- {statistics.pstdev(cv_scores) * 100:.1f})")
        except ValueError:
            print("  k-fold CV skipped (too few samples in a class)")

        # Confusion matrix + feature importance
        cm = confusion_matrix(y_test, clf.predict(X_test), labels=sorted(class_counts))
        print(f"  Confusion matrix (labels={sorted(class_counts)}):\n{cm}")
        feature_names = ([f"physio:{p}" for p in physio_vocab] +
                         [f"geology:{g}" for g in geology_vocab] +
                         ["x_norm", "y_norm"])
        top_features = np.argsort(clf.feature_importances_)[::-1][:5]
        print("  Top features by importance: "
              + ", ".join(feature_names[i] for i in top_features))

        # Retrain on all data for final model
        clf.fit(X_arr, y_arr)
        return ("random_forest", clf, physio_vocab, geology_vocab, stats)

    else:
        print("  numpy/scikit-learn not available (or too few samples) - "
              "using pure-Python weighted k-NN instead")
        return ("knn", (X, y), physio_vocab, geology_vocab, stats)


def classify_knn_probabilities(point_vec, training_X, training_y, k):
    """Distance-weighted k-NN class distribution: {label: probability %}."""
    dists = []
    for i, xv in enumerate(training_X):
        d = math.sqrt(sum((a - b) ** 2 for a, b in zip(point_vec, xv)))
        dists.append((d, training_y[i]))
    dists.sort(key=lambda t: t[0])
    weights = defaultdict(float)
    for d, label in dists[:k]:
        weights[label] += 1.0 / (d + 1e-6)
    total = sum(weights.values())
    if not total:
        return {}
    return {label: (w / total) * 100 for label, w in weights.items()}


def classify_knn(point_vec, training_X, training_y, k):
    probs = classify_knn_probabilities(point_vec, training_X, training_y, k)
    if not probs:
        return None, 0.0
    best_label = max(probs, key=probs.get)
    return best_label, probs[best_label]


def predict_series(model_info, physio, geology, x, y):
    kind, model, physio_vocab, geology_vocab, stats = model_info
    (mx, sx), (my, sy) = stats
    vec = [1.0 if physio == p else 0.0 for p in physio_vocab]
    vec += [1.0 if geology == g else 0.0 for g in geology_vocab]
    vec += [(x - mx) / sx, (y - my) / sy]

    if kind == "random_forest":
        proba = model.predict_proba([vec])[0]
        classes = model.classes_
        best_idx = int(proba.argmax())
        return classes[best_idx], float(proba[best_idx] * 100)
    else:
        train_X, train_y = model
        return classify_knn(vec, train_X, train_y, KNN_NEIGHBOURS)


def predict_series_probabilities(model_info, physio, geology, x, y):
    """Full {series: probability %} distribution for one polygon, using
    the same encoding/model as predict_series. Stage 1 uses it to rank
    primary-driver candidates; the classifier's top series is its max."""
    kind, model, physio_vocab, geology_vocab, stats = model_info
    (mx, sx), (my, sy) = stats
    vec = [1.0 if physio == p else 0.0 for p in physio_vocab]
    vec += [1.0 if geology == g else 0.0 for g in geology_vocab]
    vec += [(x - mx) / sx, (y - my) / sy]

    if kind == "random_forest":
        proba = model.predict_proba([vec])[0]
        return {str(c): float(p * 100) for c, p in zip(model.classes_, proba)}
    train_X, train_y = model
    return classify_knn_probabilities(vec, train_X, train_y, KNN_NEIGHBOURS)


# ---- Tier 0: IDW from ground-truth points ----------------------------------

def load_ground_truth_points(paths, mu_field, target_sr):
    """Load ground-truth points as (x, y, mu) tuples. FIX vs the QGIS
    original: point layers are reprojected to the target CRS first -
    the original measured the search radius in whatever units the
    points happened to be stored in (degrees vs metres)."""
    points = []
    for i, path in enumerate(paths):
        check_exists(path, "Ground-truth point layer")
        fc = project_to_target(path, f"gt_points_{i}", target_sr)
        assert_crs(fc, path)
        with arcpy.da.SearchCursor(fc, [mu_field, "SHAPE@XY"]) as cur:
            for mu, (x, y) in cur:
                if not mu:
                    continue
                points.append((x, y, str(mu)))
    return points


def idw_predict(x, y, points, radius, power):
    candidates = [(px, py, mu) for px, py, mu in points
                  if math.hypot(px - x, py - y) <= radius]
    if not candidates:
        return None, 0.0

    weights = defaultdict(float)
    total_weight = 0.0
    for px, py, mu in candidates:
        d = math.hypot(px - x, py - y)
        w = 1.0 / (d ** power + 1e-6)
        weights[mu] += w
        total_weight += w

    best_mu = max(weights, key=weights.get)
    confidence = (weights[best_mu] / total_weight) * 100 if total_weight else 0.0
    return best_mu, confidence


# ---- Tier 2: rule-based keyword matching -----------------------------------

def rule_based_match(physio, geology, series_lookup):
    """
    Very simple keyword fallback: looks for a series whose description
    shares a keyword with the physiography/geology text. Replace this with
    real domain rules for your own survey area.
    """
    text = f"{physio} {geology}".lower()
    for code, (depth, description) in series_lookup.items():
        if description and any(word in text for word in description.lower().split()):
            return code, 40.0  # low-ish confidence, rule-based tier
    return None, 0.0


# ---- Tier 3: nearest mapped polygon -----------------------------------------

def build_mapped_list(records, non_soil_codes):
    """(centroid_x, centroid_y, mu) for every already-mapped polygon."""
    mapped = []
    for rec in records:
        mu = rec["mu"]
        if mu and str(mu).strip().upper() not in non_soil_codes:
            mapped.append((rec["x"], rec["y"], str(mu)))
    return mapped


def nearest_neighbor_mu(x, y, mapped):
    """MU code of the mapped polygon whose centroid is nearest to (x, y).
    Computed manually - no Near tool / licence dependency."""
    best_mu, best_d = None, None
    for mx, my, mu in mapped:
        d = math.hypot(mx - x, my - y)
        if best_d is None or d < best_d:
            best_d, best_mu = d, mu
    return best_mu


# ---- Step 7: assignment ------------------------------------------------------

# (field name, arcpy type, text length) - names/lengths mirror the QGIS
# original; all names fit the 10-character shapefile limit.
OUTPUT_FIELDS = [
    ("SeriesCode", "TEXT", 20),
    ("DepthNum", "TEXT", 5),
    ("SurTexCode", "TEXT", 10),
    ("SlopeCode", "TEXT", 10),
    ("ErosionNum", "TEXT", 5),
    ("MU_New", "TEXT", 30),
    ("MU_Source", "TEXT", 20),
    ("MU_Conf", "DOUBLE", None),
    ("MU_Method", "TEXT", 40),
    ("MU_Flag", "TEXT", 100),
]

_EMPTY_COMPS = {"SeriesCode": "", "DepthNum": "", "SurTexCode": "",
                "SlopeCode": "", "ErosionNum": ""}


def ensure_fields(fc):
    """Add the output fields to the joined layer if they don't exist."""
    existing = {f.name for f in arcpy.ListFields(fc)}
    for name, ftype, length in OUTPUT_FIELDS:
        if name not in existing:
            if ftype == "TEXT":
                arcpy.management.AddField(fc, name, "TEXT", field_length=length)
            else:
                arcpy.management.AddField(fc, name, "DOUBLE")


def assign_mapping_units(records, gt_points, model_info, mapped,
                         evidence, total_support):
    """Run the two-stage assignment over the in-memory records, filling
    in comps / mu_new / source / conf / method / flag on each processed
    record. Non-soil polygons are left untouched. Returns
    (method_counts, confidences)."""
    method_counts = Counter()
    confidences = []

    for rec in records:
        existing_mu = rec["mu"]

        if existing_mu and str(existing_mu).strip().upper() not in NON_SOIL_CODES:
            # Surveyor-verified - never altered, only decomposed.
            comps = decompose_mu(existing_mu)
            method, conf, flag = "existing", 100.0, ""
            mu_new = str(existing_mu)
            if comps is None:
                comps = dict(_EMPTY_COMPS)
                flag = ("existing MU code does not match the expected "
                        "pattern - kept as-is")
        elif existing_mu and str(existing_mu).strip().upper() in NON_SOIL_CODES:
            continue  # leave non-soil polygons untouched
        else:
            series, conf, method, flag = None, 0.0, "fallback", ""
            full_code = None

            # Tier 0: IDW from ground-truth points
            if gt_points:
                idw_mu, idw_conf = idw_predict(
                    rec["x"], rec["y"], gt_points,
                    IDW_SEARCH_RADIUS_M, IDW_POWER
                )
                if idw_mu and idw_conf >= IDW_CONFIDENCE_MIN:
                    series, conf, method = idw_mu, idw_conf, "idw_ground_truth"
                    full_code = idw_mu

            # Stage 1: confirm the series from the PRIMARY drivers
            # (geology + physiography + depth). The classifier works
            # WITHIN the primary-driver evidence: its top series is
            # accepted only when it is one of the evidence candidates;
            # otherwise the top-ranked candidate wins.
            if series is None:
                probs = {}
                pred_series, pred_conf = None, 0.0
                if model_info:
                    probs = predict_series_probabilities(
                        model_info, rec["physio"], rec["geology"],
                        rec["x"], rec["y"]
                    )
                    if probs:
                        pred_series = max(probs, key=probs.get)
                        pred_conf = probs[pred_series]
                candidates = rank_series_candidates(
                    evidence, rec["geology"], rec["physio"],
                    known_depth=rec["depth"] or None,
                    classifier_probs=probs or None,
                )
                candidate_series = {c["series"] for c in candidates}
                if candidates and pred_series in candidate_series:
                    series, conf, method = pred_series, pred_conf, "ml_classifier"
                elif candidates:
                    top = candidates[0]
                    series, method = top["series"], "primary_driver_match"
                    conf = primary_driver_confidence(
                        top["depth_agreement"], top["support"], total_support
                    )
                    if pred_series is not None:
                        flag = ("series set by primary-driver evidence "
                                "over classifier")
                elif (pred_series is not None
                        and pred_conf >= CLASSIFIER_CONF_MIN):
                    # No primary-driver evidence (unseen geology /
                    # physiography combination): classifier alone.
                    series, conf, method = pred_series, pred_conf, "ml_classifier"

            # Fallback: rule-based keyword matching
            if series is None:
                rule_series, rule_conf = rule_based_match(
                    rec["physio"], rec["geology"], SERIES_LOOKUP
                )
                if rule_series:
                    series, conf, method = rule_series, rule_conf, "rule_based"

            # Fallback: nearest already-mapped polygon (series only -
            # its components are completed in Stage 2 like every other
            # Stage 1 path)
            if series is None:
                nn_mu = nearest_neighbor_mu(rec["x"], rec["y"], mapped)
                if nn_mu:
                    nn_comps = decompose_mu(nn_mu)
                    series = nn_comps["SeriesCode"] if nn_comps else nn_mu
                    conf, method = 30.0, "nearest_neighbor"

            # Last resort: flagged default
            if series is None:
                series, conf, method = "UNK", 0.0, "fallback"
                flag = "needs manual review - no evidence matched"

            if full_code is not None:
                # Tier 0 votes for a COMPLETE MU code. Decompose it into
                # components exactly like an existing code (FIX: the QGIS
                # original stored the whole code in SeriesCode too).
                comps = decompose_mu(full_code) or dict(_EMPTY_COMPS,
                                                        SeriesCode=series)
                mu_new = full_code
            elif method in ("primary_driver_match", "ml_classifier",
                            "rule_based", "nearest_neighbor"):
                # Stage 2: complete the mapping unit from the SECONDARY
                # drivers - the confirmed series' modal components
                # (same-physiography polygons preferred).
                comps, mu_new, incomplete_note = apply_stage2(
                    evidence, series, rec["physio"],
                    known_depth=rec["depth"], slope_code=rec["slope"],
                )
                if incomplete_note:
                    flag = (f"{flag}; {incomplete_note}" if flag
                            else incomplete_note)
            else:
                comps = dict(_EMPTY_COMPS, SeriesCode=series)
                mu_new = series

        rec["comps"] = comps
        rec["mu_new"] = mu_new
        rec["source"] = "existing" if method == "existing" else "predicted"
        rec["conf"] = conf
        rec["method"] = method
        rec["flag"] = flag
        method_counts[method] += 1
        if method != "existing":
            confidences.append(conf)

    return method_counts, confidences


def write_results(fc, records):
    """Write the computed fields back to the joined layer with an
    UpdateCursor, keyed by OID."""
    oid_fld = arcpy.Describe(fc).OIDFieldName
    by_oid = {r["oid"]: r for r in records if "mu_new" in r}
    fields = [oid_fld] + [name for name, _, _ in OUTPUT_FIELDS]
    with arcpy.da.UpdateCursor(fc, fields) as cur:
        for row in cur:
            rec = by_oid.get(row[0])
            if rec is None:
                continue  # non-soil polygon - untouched
            comps = rec["comps"]
            cur.updateRow([row[0],
                           safe_str(comps.get("SeriesCode"), 20),
                           safe_str(comps.get("DepthNum"), 5),
                           safe_str(comps.get("SurTexCode"), 10),
                           safe_str(comps.get("SlopeCode"), 10),
                           safe_str(comps.get("ErosionNum"), 5),
                           safe_str(rec["mu_new"], 30),
                           rec["source"],
                           float(rec["conf"]),
                           safe_str(rec["method"], 40),
                           safe_str(rec["flag"], 100)])


# ---- Step 8: adjacency conflict check (flag-only) ---------------------------

def flag_adjacency_conflicts(records):
    """
    Flag predicted polygons that TOUCH another predicted polygon with an
    identical MU_New, so a human can verify them.

    FIX vs the QGIS original: the original 'resolved' such conflicts by
    silently incrementing the MU code's last digit - i.e. it invented
    soil data (AD -> AD1). This version never alters a code; it writes
    the conflict to MU_Flag instead. Surveyor-verified ('existing')
    polygons are never flagged, and a predicted polygon agreeing with a
    verified neighbour is not a conflict.
    """
    predicted = [r for r in records
                 if r.get("source") == "predicted" and r.get("geom") is not None]
    flagged = 0
    for rec in predicted:
        geom = rec["geom"]
        extent = geom.extent
        for other in predicted:
            if other["oid"] == rec["oid"] or other["mu_new"] != rec["mu_new"]:
                continue
            o_ext = other["geom"].extent
            if (extent.XMax < o_ext.XMin or o_ext.XMax < extent.XMin or
                    extent.YMax < o_ext.YMin or o_ext.YMax < extent.YMin):
                continue  # bounding boxes don't even overlap
            if geom.touches(other["geom"]):
                note = "adjacent predicted polygon has identical MU_New - verify"
                rec["flag"] = f"{rec['flag']}; {note}" if rec["flag"] else note
                flagged += 1
                break
    return flagged


# ---- Step 9: save output -----------------------------------------------------

def save_output(joined_fc, output_path):
    parent = os.path.dirname(os.path.abspath(output_path))
    if ".gdb" not in output_path.lower() and parent and not os.path.exists(parent):
        os.makedirs(parent, exist_ok=True)
    if arcpy.Exists(output_path):
        arcpy.management.Delete(output_path)
    arcpy.management.CopyFeatures(joined_fc, output_path)
    return output_path


def add_to_map(output_path):
    """Add the output to the current ArcGIS Pro project's map, if a
    project is open. Silently skipped for standalone runs."""
    try:
        aprx = arcpy.mp.ArcGISProject("CURRENT")
        maps = aprx.listMaps()
        if not maps:
            return
        (aprx.activeMap or maps[0]).addDataFromPath(output_path)
        print("  Output layer added to the current map.")
    except Exception as exc:
        print(f"  (Output not added to a map: {exc})")


# ---- Main pipeline -----------------------------------------------------------

def run_pipeline():
    TOTAL_STEPS = 10
    target_sr = arcpy.SpatialReference(TARGET_CRS_WKID)
    arcpy.env.overwriteOutput = True

    log(0, TOTAL_STEPS, "Fixing geometry and aligning CRS (all layers)...")
    for path, label in ((SOIL_LAYER, "Soil layer"),
                        (PHYSIO_LAYER, "Physiography layer"),
                        (GEOLOGY_LAYER, "Geology layer")):
        check_exists(path, label)
    soil = prepare_polygon_layer(SOIL_LAYER, "soil_prep", target_sr)
    physio = prepare_polygon_layer(PHYSIO_LAYER, "physio_prep", target_sr)
    geology = prepare_polygon_layer(GEOLOGY_LAYER, "geology_prep", target_sr)

    log(1, TOTAL_STEPS, "Layers loaded (all reprojected to target CRS).")

    log(2, TOTAL_STEPS, "Loading ground-truth points for IDW...")
    gt_points = load_ground_truth_points(
        GROUND_TRUTH_POINT_LAYERS, POINT_MU_FLD, target_sr)
    print(f"  Loaded {len(gt_points)} ground-truth points")

    log(3, TOTAL_STEPS, "Spatial join: soil <- physiography (largest overlap)...")
    joined = spatial_join_largest_overlap(
        soil, physio, PHYSIO_NAME_FLD, os.path.join("in_memory", "joined_physio"))

    log(4, TOTAL_STEPS, "Spatial join: result <- geology (largest overlap)...")
    joined = spatial_join_largest_overlap(
        joined, geology, GEOLOGY_NAME_FLD, os.path.join("in_memory", "joined_final"))
    drop_join_scaffolding(joined)

    log(5, TOTAL_STEPS, "Decomposing existing MU codes / building training set...")
    records = read_polygons(joined)
    rows, labels = build_training_set(records, NON_SOIL_CODES)
    print(f"  {len(rows)} verified polygons available for training")
    evidence = build_series_evidence(
        build_evidence_rows(records, NON_SOIL_CODES))
    total_support = sum(e["support"] for e in evidence.values())
    print(f"  Series-evidence table: {len(evidence)} series from "
          f"{total_support} mapped polygons")

    log(6, TOTAL_STEPS, "Training & validating classifier...")
    model_info = None
    if rows:
        model_info = train_classifier(rows, labels)
    else:
        print("  No training data available - classifier tier will be skipped")

    log(7, TOTAL_STEPS, "Assigning series to unmapped polygons...")
    ensure_fields(joined)
    mapped = build_mapped_list(records, NON_SOIL_CODES)
    method_counts, confidences = assign_mapping_units(
        records, gt_points, model_info, mapped, evidence, total_support)

    log(8, TOTAL_STEPS, "Checking adjacency conflicts (flag-only)...")
    adj_flagged = flag_adjacency_conflicts(records)
    print(f"  {adj_flagged} polygons flagged for adjacency review")
    write_results(joined, records)

    log(9, TOTAL_STEPS, "Saving output...")
    save_output(joined, OUTPUT_FC)
    print(f"  Saved: {OUTPUT_FC}")

    log(10, TOTAL_STEPS, "Summary report")
    print("  MU assignment method breakdown:")
    for method, count in method_counts.most_common():
        print(f"    {method}: {count}")
    if confidences:
        print(f"  Mean confidence (predicted polygons): "
              f"{statistics.mean(confidences):.1f}%")
    flagged_count = sum(1 for r in records if r.get("flag"))
    print(f"  Polygons flagged for manual review: {flagged_count}")

    add_to_map(OUTPUT_FC)
    print("\nDone. Output saved to the path in OUTPUT_FC.")


if __name__ == "__main__":
    run_pipeline()
