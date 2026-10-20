# Changelog

All notable changes to the ArcGIS Pro port of Mapping-Units-Generator.

## [1.1.0] - 2026-10-20 - SLUSI two-stage driver method

### Changed
- Assignment restructured into the SLUSI two-stage driver hierarchy (mirrors QGIS v10). **Stage 1 (primary drivers — geology, physiography, depth)** confirms the soil series: a series-evidence table is built from the already-mapped polygons (per series: observed geology/physiography values, depth-class distribution, support, modal MU components overall and per physiography), and an unmapped polygon's candidate series must match its geology AND physiography. Candidates are ranked by depth agreement, then classifier probability, then support.
- The classifier now works within the evidence: its top series is accepted (`ml_classifier`) only when it is among the candidates; otherwise the top-ranked candidate wins as `primary_driver_match` (confidence = 50 + 30 × depth-agreement share + 20 × relative support) and the polygon is flagged "series set by primary-driver evidence over classifier". With no candidates (unseen geology/physiography combination), the previous fall-through order is preserved: classifier alone above `CLASSIFIER_CONF_MIN` → keyword rules → nearest mapped polygon → flagged `UNK`.
- **Stage 2 (secondary drivers)** completes the mapping unit for every Stage-1 assignment: DepthNum/SurTexCode/SlopeCode/ErosionNum are filled with the series' modal component values (same-physiography polygons preferred, series-wide modal as fallback) and `MU_New` is assembled as the full Series+Depth+Texture+Slope+Erosion code. Components with no modal value anywhere stay blank and are flagged `MU incomplete: <component> unknown`; nothing is invented.
- New optional configuration: `DEPTH_FLD` (soil-layer field with a known depth class per polygon — feeds Stage 1 depth agreement and supplies the Stage 2 depth) and `SLOPE_FLD` (soil-layer slope-code field — overrides the Stage 2 modal slope). Empty (the default) disables each; the fields are read via the record cursor when configured and present.
- Tier 0 (IDW ground truth) unchanged and still the highest authority; surveyor-verified polygons unchanged.

## [1.0.0] - 2026-10-20

### Added
- Initial ArcGIS Pro (arcpy) port of the QGIS Mapping-Units-Generator
  (`mu_generator.py`), as `mu_generator_arcgis.py`.
- Same five-tier assignment pipeline: Tier 0 IDW from ground-truth
  points (manual implementation, no Spatial Analyst licence needed),
  Tier 1 Random Forest via scikit-learn with a pure-Python weighted
  k-NN fallback, Tier 2 keyword rules, Tier 3 nearest mapped polygon,
  Tier 4 flagged default.
- Same output fields as the QGIS original (`SeriesCode`, `DepthNum`,
  `SurTexCode`, `SlopeCode`, `ErosionNum`, `MU_New`, `MU_Source`,
  `MU_Conf`, `MU_Method`, `MU_Flag`) and the same end-of-run summary
  (method counts, mean confidence, flagged count).
- arcpy equivalents throughout: `Project` for reprojection,
  `RepairGeometry` for geometry fixing, `SpatialJoin`
  (`JOIN_ONE_TO_ONE` + `LARGEST_OVERLAP`) for the physiography/geology
  joins, `arcpy.da` cursors for reading/assignment/write-back, and
  automatic add-to-map via `arcpy.mp` when a project is open.

### Fixed (relative to the QGIS original)
- `decompose_mu()` is now a strict regex parser for the
  Series+Depth+Texture+Slope+Erosion grammar; the original's positional
  slicing wrote the wrong character into `SlopeCode`.
- All input layers — including ground-truth points — are reprojected to
  the target CRS and asserted equal before joining/measuring; the
  original reprojected only the soil layer, so mismatched CRSs failed
  silently.
- Classifier coordinates are z-score normalised, so the k-NN fallback
  is no longer dominated by raw map coordinates.
- The adjacency step only flags conflicts in `MU_Flag`; it no longer
  rewrites MU codes to break ties.
- Tier 0 (IDW) results are decomposed into MU components consistently,
  instead of storing the full code in `SeriesCode`.
