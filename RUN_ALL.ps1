param(
    [Parameter(Mandatory = $true)]
    [string]$Python,
    [switch]$DownloadInputs
)

$ErrorActionPreference = "Stop"
$project = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $project

function Invoke-Analysis {
    param([Parameter(Mandatory = $true)][string]$Script, [string[]]$Arguments = @())
    & $Python $Script @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "Analysis failed: $Script (exit code $LASTEXITCODE)"
    }
}

function Invoke-Preparation {
    param([Parameter(Mandatory = $true)][string]$Script)
    & powershell -NoProfile -ExecutionPolicy Bypass -File $Script
    if ($LASTEXITCODE -ne 0) {
        throw "Preparation failed: $Script (exit code $LASTEXITCODE)"
    }
}

if ($DownloadInputs) {
    Invoke-Analysis scripts/fetch_public_inputs.py --download
}
Invoke-Analysis scripts/fetch_public_inputs.py

# These two preparation steps unpack only the GSE278687 files read by the
# locked pipeline. Raw archives themselves remain untouched and untracked.
Invoke-Preparation scripts/prepare_spatial_validation_expression.ps1
Invoke-Preparation scripts/prepare_gse278687_spatial_coordinates.ps1

Invoke-Analysis scripts/qc_gse277116_full_packages.py
Invoke-Analysis scripts/run_unified_primary_pipeline_v2.py
Invoke-Analysis scripts/run_gse278687_block_null_v3_999.py
Invoke-Analysis scripts/run_gse278687_block_size_sensitivity_v1.py
1..4 | ForEach-Object { $env:SPATIAL_NULL_BATCH = "$_"; Invoke-Analysis scripts/run_spatial_nulls_v3_999.py all }
Remove-Item Env:SPATIAL_NULL_BATCH -ErrorAction SilentlyContinue
Invoke-Analysis scripts/combine_spatial_null_batches_v3.py
1..4 | ForEach-Object {
    $env:MASK_MORAN_BATCH = "$_"
    Invoke-Analysis scripts/run_mask_moran_sensitivity_v4.py gse278687
    Invoke-Analysis scripts/run_mask_moran_sensitivity_v4.py gse277116
}
Remove-Item Env:MASK_MORAN_BATCH -ErrorAction SilentlyContinue
Invoke-Analysis scripts/combine_mask_moran_batches_v4.py
Invoke-Analysis scripts/run_program_scoreability_v3.py
Invoke-Analysis scripts/run_scoreability_threshold_sensitivity_v4.py
Invoke-Analysis scripts/run_gse217845_marker_gated_single_cell_validation.py
Invoke-Analysis scripts/run_conditional_scale_joint_sensitivities_v7.py
Invoke-Analysis scripts/run_spatial_holdout_validation.py
Invoke-Analysis scripts/audit_gse202051_author_reference.py
Invoke-Analysis scripts/run_gse202051_author_annotated_attribution.py
Invoke-Analysis scripts/generate_gse277116_representative_spatial_maps_v1.py
Invoke-Analysis scripts/generate_main_figures_v4.py
Invoke-Analysis scripts/generate_spatial_and_supplementary_figures_v6.py
Invoke-Analysis scripts/generate_conditional_scale_figures_v7.py
Invoke-Analysis scripts/generate_holdout_figure_v1.py
Invoke-Analysis scripts/assemble_submission_supplements_v7.py
Invoke-Analysis scripts/build_manuscript_docx_v7.py
Invoke-Analysis scripts/build_release_manifest.py
Invoke-Analysis scripts/verify_release.py --require-inputs
Invoke-Analysis scripts/build_public_release_zip.py
