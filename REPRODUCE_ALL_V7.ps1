param(
    [Parameter(Mandatory = $true)]
    [string]$PythonSpatial,
    [Parameter(Mandatory = $true)]
    [string]$PythonSingleCell
)

$ErrorActionPreference = "Stop"
$project = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $project

& $PythonSpatial scripts/run_unified_primary_pipeline_v2.py
& $PythonSpatial scripts/run_gse278687_block_null_v3_999.py
1..4 | ForEach-Object { $env:SPATIAL_NULL_BATCH = "$_"; & $PythonSpatial scripts/run_spatial_nulls_v3_999.py all }
& $PythonSpatial scripts/combine_spatial_null_batches_v3.py
1..4 | ForEach-Object { $env:MASK_MORAN_BATCH = "$_"; & $PythonSpatial scripts/run_mask_moran_sensitivity_v4.py gse278687; & $PythonSpatial scripts/run_mask_moran_sensitivity_v4.py gse277116 }
& $PythonSpatial scripts/combine_mask_moran_batches_v4.py
& $PythonSpatial scripts/run_scoreability_threshold_sensitivity_v4.py
& $PythonSpatial scripts/run_gse217845_marker_gated_single_cell_validation.py
& $PythonSpatial scripts/run_conditional_scale_joint_sensitivities_v7.py
& $PythonSpatial scripts/run_spatial_holdout_validation.py
& $PythonSingleCell scripts/audit_gse202051_author_reference.py
& $PythonSingleCell scripts/run_gse202051_author_annotated_attribution.py
& $PythonSpatial scripts/generate_main_figures_v4.py
& $PythonSpatial scripts/generate_spatial_and_supplementary_figures_v6.py
& $PythonSpatial scripts/generate_conditional_scale_figures_v7.py
& $PythonSpatial scripts/generate_holdout_figure_v1.py
& $PythonSpatial scripts/assemble_submission_supplements_v7.py
& $PythonSpatial scripts/build_manuscript_docx_v7.py
& $PythonSpatial scripts/build_release_manifest_v7.py
& $PythonSpatial scripts/verify_release_v7.py
