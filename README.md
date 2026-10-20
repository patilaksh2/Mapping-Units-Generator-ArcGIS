# Mapping Units Generator — ArcGIS Pro (arcpy) port

Fills in blank **Soil Mapping Unit (MU) codes** on a semi-detailed soil
survey polygon layer, using the evidence already on the map instead of
guessing: nearby field observations, landscape position (physiography),
geology, and the pattern of already-mapped neighbouring polygons.

This is the **ArcGIS Pro port** of the
[QGIS Mapping-Units-Generator](https://github.com/patilaksh2/Mapping-Units-Generator).
Same two-stage method, same configuration concept, same output fields —
the QGIS/Processing calls are replaced by their `arcpy` equivalents, and
the port carries the *fixed* logic (see "Differences from the QGIS
version").

## How it works

Every unmapped polygon is assigned by the first step that produces an
answer:

| Step | Method | When it wins | Confidence |
|------|--------|--------------|------------|
| Tier 0 | IDW interpolation from ground-truth field points | Field points within `IDW_SEARCH_RADIUS_M` agree strongly | ≥ `IDW_CONFIDENCE_MIN` (default 80%) |
| Stage 1 | Primary-driver evidence + classifier (see *Method*) | The polygon's geology + physiography match mapped series evidence; the classifier's top series is accepted if it is among the evidence candidates, otherwise the top-ranked candidate wins (`primary_driver_match`) | Classifier probability, or the heuristic score 50 + 30 × depth term + 20 × relative support |
| Fallback | Classifier alone | No evidence candidates, but the predicted series probability is high enough | ≥ `CLASSIFIER_CONF_MIN` (default 60%) |
| Fallback | Rule-based keyword matching of physiography/geology text against `SERIES_LOOKUP` | A series description shares a keyword with the polygon's landscape/geology text | flat 40% |
| Fallback | Nearest already-mapped polygon (by centroid distance) | Nothing else matched | flat 30% |
| Last resort | Flagged default `UNK` | No evidence at all | 0%, flagged for manual review |

Before assignment, the soil layer is spatially joined with the
physiography and geology layers (**largest-overlap** join), so every
polygon knows its landscape unit and lithology. Polygons that already
have a surveyor-verified MU code are **never altered** — their codes are
only decomposed into components. Non-soil polygons (water bodies,
habitation, roads, quarries — see `NON_SOIL_CODES`) are left completely
untouched.

After assignment, a final **adjacency check** flags (never rewrites)
predicted polygons that touch another predicted polygon with an
identical new MU, so a surveyor can verify them.

## Method: SLUSI two-stage driver hierarchy

Assignment follows the two-stage driver hierarchy used in SLUSI-style soil survey — the **primary drivers confirm the soil series first**; only then do the **secondary drivers complete the mapping unit**.

**Stage 1 — confirm the series (primary drivers: Geology, Physiography, Depth).** From the already-mapped polygons the script builds a *series-evidence table*: for every series, the geology values and physiography values it occurs with, its depth-class distribution, its polygon count (support), and its modal MU components. For an unmapped polygon, candidate series are those whose observed geology **and** physiography both match the polygon's joined values. Candidates are ranked by (1) the depth term — when a depth field is configured (`DEPTH_FLD`), true *depth agreement*: the share of the candidate's polygons carrying the polygon's own depth class; when `DEPTH_FLD` is not set, a *depth-concentration weight*: the share of the candidate series' mapped polygons that sit at its modal depth class, which measures how depth-consistent the series itself is rather than agreement with this polygon; (2) the classifier's probability for that series; (3) support. The trained classifier (Random Forest / k-NN) still predicts a series, but it operates *within* this evidence: if its top series is among the candidates it is accepted (`ml_classifier`); if candidates exist but the classifier's pick is not one of them, the top-ranked candidate wins (`primary_driver_match`, score = 50 + 30 × depth-term share + 20 × relative support) and the polygon is flagged "series set by primary-driver evidence over classifier". Only when the candidate set is empty (an unseen geology/physiography combination) does the pipeline fall through to the classifier alone (if it clears `CLASSIFIER_CONF_MIN`), then keyword rules, then the nearest mapped polygon, then the flagged `UNK` default. Tier 0 (IDW from ground-truth points, which returns a complete MU code) stays above Stage 1 as direct field evidence. One deliberate design rule: when the classifier's top series is inside the candidate set, it is accepted at its model probability **without** the `CLASSIFIER_CONF_MIN` gate — membership in the primary-driver evidence is itself the gate. `CLASSIFIER_CONF_MIN` applies only on the no-candidate path, where the classifier stands alone.

**Stage 2 — complete the mapping unit (secondary drivers).** Once the series is fixed, `DepthNum`, `SurTexCode`, `SlopeCode` and `ErosionNum` are filled from the evidence table's modal component values for that series — preferring polygons in the same physiography, falling back to series-wide modal values. A configured `DEPTH_FLD` value supplies the depth class, and `SLOPE_FLD` the slope code, when populated. `MU_New` is then assembled as the full code `Series + Depth + SurfaceTex + Slope + Erosion`. A component with no modal value anywhere is left blank and flagged (`MU incomplete: <component> unknown`) — values are never invented.

## Configuration

Everything project-specific lives in the `CONFIGURATION` block at the
top of `mu_generator_arcgis.py`:

- **Paths** — `SOIL_LAYER`, `PHYSIO_LAYER`, `GEOLOGY_LAYER`,
  `OUTPUT_FC` (shapefiles or feature classes in a file geodatabase;
  absolute paths recommended), plus `GROUND_TRUTH_POINT_LAYERS`
  (field point layers for Tier 0; may be empty).
- **Field names** — `PHYSIO_NAME_FLD` (`legend_sho`), `GEOLOGY_NAME_FLD`
  (`lithologic`), `EXISTING_MU_FLD` and `POINT_MU_FLD` (`MU`), plus the
  optional `DEPTH_FLD` / `SLOPE_FLD` (soil-layer fields holding a known
  depth class / slope code per polygon; blank = disabled).
- **CRS** — `TARGET_CRS_WKID`, the EPSG WKID of a projected CRS
  (default `32643` = WGS 84 / UTM zone 43N). **All** input layers are
  reprojected to it, and the script raises an error if any layer does
  not end up in it — distances are computed in map units, so a projected
  CRS is essential.
- **`NON_SOIL_CODES`** — MU values excluded from classification.
- **`SERIES_LOOKUP`** — your survey area's series table:
  `code: (typical_depth_class, description)`. The descriptions feed the
  Tier-2 keyword matcher, so rebuild this for your own district.
- **Tuning** — IDW radius/power and confidence threshold, classifier
  confidence threshold, k-NN neighbours, `RANDOM_SEED`.

## Inputs

| Input | Geometry | Required content |
|-------|----------|------------------|
| Soil polygons | Polygon | `MU` field: filled for surveyed polygons, blank for unmapped ones |
| Physiography | Polygon | `legend_sho` landscape/legend text |
| Geology | Polygon | `lithologic` lithology text |
| Ground-truth points (optional) | Point | `MU` field with observed codes |

MU codes must follow the grammar
`Series(2 letters) + Depth(1 digit) + SurfaceTexture(1 letter) + Slope(1 letter) + Erosion(1 digit)`,
e.g. `AD3rB2`. Codes that don't match are kept as-is and flagged.

## Output fields

The output is a copy of the soil layer (originals are never modified)
with these fields added:

| Field | Type | Meaning |
|-------|------|---------|
| `SeriesCode` | Text (20) | Soil series code — decomposed from the existing code, or predicted |
| `DepthNum` | Text (5) | Depth-class digit from the MU code |
| `SurTexCode` | Text (10) | Surface-texture code letter |
| `SlopeCode` | Text (10) | Slope-class code letter |
| `ErosionNum` | Text (5) | Erosion-class digit |
| `MU_New` | Text (30) | Final MU: the full original code for existing polygons and Tier-0 (IDW) assignments; for Stage-1 assignments, the assembled Series+Depth+Texture+Slope+Erosion code (components completed in Stage 2) |
| `MU_Source` | Text (20) | `existing` (surveyor-verified) or `predicted` |
| `MU_Conf` | Double | Heuristic confidence **score** (0–100) — *not* a calibrated probability. 100 for existing; for predictions the value mixes model probabilities, IDW vote shares and heuristic values depending on `MU_Method`, so scores are not strictly comparable across methods |
| `MU_Method` | Text (40) | `existing` / `idw_ground_truth` / `primary_driver_match` / `ml_classifier` / `rule_based` / `nearest_neighbor` / `fallback` |
| `MU_Flag` | Text (100) | Review flags: unmatched codes, Tier-4 defaults, adjacency conflicts — empty means nothing to check |

The console ends with a summary: counts per method, mean confidence of
predicted polygons, and how many polygons are flagged for manual review.

## How to run

**Requirements:** ArcGIS Pro 3.x. No Spatial Analyst licence is needed
(the IDW is implemented manually). `numpy` ships with ArcGIS Pro's
Python; `scikit-learn` is optional — without it, Tier 1 falls back to
the pure-Python weighted k-NN.

1. Edit the `CONFIGURATION` block in `mu_generator_arcgis.py`.
2. Run it any of these ways:
   - **Python window** in ArcGIS Pro:
     `exec(open(r"C:\path\to\mu_generator_arcgis.py").read())`
   - As a **script tool** in a custom toolbox (the script takes no
     parameters).
   - **Standalone**:
     `"C:\Program Files\ArcGIS\Pro\bin\Python\Scripts\propy.bat" mu_generator_arcgis.py`
3. The output is written to `OUTPUT_FC` and, if a project is open,
   added to the current map.

## Known limitations

- Validation is random, not spatial: the classifier's holdout / k-fold splits are stratified but random, so neighbouring polygons of the same series appear on both sides of every split; under spatial autocorrelation the printed accuracy is optimistic (an upper bound, not field accuracy).
- Coordinates are classifier features: the (normalised) polygon centroid is part of the model input, so part of what the classifier learns is *where* each series occurs, not only transferable landscape relationships. Performance may not carry over to a new sheet or area.
- Stage 2 modal imputation homogenises within-series variation: every predicted polygon of a series in a physiography receives the same texture/slope/erosion, erasing real variation within the series (see also the Stage 2 note below).
- Geology/physiography matching is by literal name: candidate matching compares the joined legend/lithology strings exactly, so inconsistent naming between layers (synonyms, case differences, spacing) yields no candidates and the polygon falls back to lower tiers.
- The classifier is trained on physiography + geology + (normalised)
  position only — the same evidence set as the QGIS original. It cannot
  use covariates the inputs don't carry.
- Tier 2's keyword matcher is deliberately naive (any shared word in a
  series description counts). Treat its 40% confidence accordingly, and
  replace it with real domain rules for production use.
- Stage 2 completes depth/texture/slope/erosion from the assigned
  series' modal components; a series with no mapped polygons (e.g. one
  matched only by a keyword rule) has no modal values, so its
  components stay blank and the polygon is flagged `MU incomplete`.
- The adjacency check and the Tier-3 nearest search are O(n²)-ish over
  in-memory centroids/geometries; fine for survey-sheet sized layers,
  slow on very large ones.
- If `OUTPUT_FC` points inside a file geodatabase, the `.gdb` must
  already exist.
- Shapefile outputs inherit shapefile limits (10-character field names
  — all output names fit; 254-character text max).

## Differences from the QGIS version

| Aspect | QGIS original | This port |
|--------|---------------|-----------|
| Runtime | QGIS Python console, `qgis.core` + Processing | ArcGIS Pro, `arcpy` |
| Geometry fix | `native:fixgeometries` | `RepairGeometry` |
| Reprojection | `native:reprojectlayer` — **soil layer only** | `Project` applied to **every** layer (incl. ground-truth points), with a CRS assertion before joins |
| Spatial join | `native:joinattributesbylocation` (shipped with `METHOD: 0` = one-to-many, duplicating polygons) | `SpatialJoin` with `JOIN_ONE_TO_ONE` + `LARGEST_OVERLAP` |
| MU parsing | Positional letter/digit slicing (corrupted `SlopeCode`) | Strict regex parser; non-matching codes kept as-is and flagged |
| Classifier features | Raw coordinates appended (docstring claimed normalised) | Z-score normalised coordinates — the k-NN fallback now actually uses physiography/geology |
| Tier 0 result | Full MU code stored in `SeriesCode` as well | Decomposed into components like every other tier |
| Adjacency step | Silently incremented MU digits to break ties (invented data) | **Flag-only**: writes the conflict to `MU_Flag`, never alters a code |
| Output | Shapefile via `QgsVectorFileWriter` | Shapefile or geodatabase feature class via `CopyFeatures`; added to the Pro map when a project is open |

## License

MIT — see [LICENSE](LICENSE).
