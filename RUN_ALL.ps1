param(
    [Parameter(Mandatory = $true)]
    [string]$Python,
    [switch]$DownloadInputs
)

$ErrorActionPreference = "Stop"
$project = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $project

if ($DownloadInputs) {
    & $Python scripts/fetch_public_inputs.py --download
}
& $Python scripts/fetch_public_inputs.py

# These two preparation steps unpack only the GSE278687 files read by the
# locked pipeline. Raw archives themselves remain untouched and untracked.
& powershell -NoProfile -ExecutionPolicy Bypass -File scripts/prepare_spatial_validation_expression.ps1
& powershell -NoProfile -ExecutionPolicy Bypass -File scripts/prepare_gse278687_spatial_coordinates.ps1

& $Python scripts/qc_gse277116_full_packages.py
& $Python scripts/run_unified_primary_pipeline_v2.py
& $Python scripts/run_gse278687_block_null_v3_999.py
1..4 | ForEach-Object { $env:SPATIAL_NULL_BATCH = "$_"; & $Python scripts/run_spatial_nulls_v3_999.py all }
Remove-Item Env:SPATIAL_NULL_BATCH -ErrorAction SilentlyContinue
& $Python scripts/combine_spatial_null_batches_v3.py
1..4 | ForEach-Object {
    $env:MASK_MORAN_BATCH = "$_"
    & $Python scripts/run_mask_moran_sensitivity_v4.py gse278687
    & $Python scripts/run_mask_moran_sensitivity_v4.py gse277116
}
Remove-Item Env:MASK_MORAN_BATCH -ErrorAction SilentlyContinue
& $Python scripts/combine_mask_moran_batches_v4.py
& $Python scripts/run_scoreability_threshold_sensitivity_v4.py
& $Python scripts/run_gse217845_marker_gated_single_cell_validation.py
& $Python scripts/run_conditional_scale_joint_sensitivities_v7.py
& $Python scripts/run_spatial_holdout_validation.py
& $Python scripts/audit_gse202051_author_reference.py
& $Python scripts/run_gse202051_author_annotated_attribution.py
& $Python scripts/generate_main_figures_v4.py
& $Python scripts/generate_spatial_and_supplementary_figures_v6.py
& $Python scripts/generate_conditional_scale_figures_v7.py
& $Python scripts/generate_holdout_figure_v1.py
& $Python scripts/assemble_submission_supplements_v7.py
& $Python scripts/build_manuscript_docx_v7.py
& $Python scripts/build_release_manifest.py
& $Python scripts/verify_release.py --require-inputs

