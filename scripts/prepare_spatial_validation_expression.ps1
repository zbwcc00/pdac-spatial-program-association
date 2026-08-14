$projectRoot = Split-Path -Parent $PSScriptRoot
$rawSpatial = Join-Path $projectRoot "data\00_raw\spatial"
$unpackedSpatial = Join-Path $projectRoot "data\01_unpacked\spatial"

$targets = @(
    @{ archive = "GSE278687_RAW.tar"; destination = "GSE278687"; pattern = "GSM*_filtered_feature_bc_matrix.h5" }
)

foreach ($target in $targets) {
    $archive = Join-Path $rawSpatial $target.archive
    $destination = Join-Path $unpackedSpatial $target.destination
    if (-not (Test-Path -LiteralPath $archive)) {
        throw "Missing archive: $archive"
    }
    New-Item -ItemType Directory -Force -Path $destination | Out-Null
    $files = tar -tf $archive | Where-Object { $_ -like $target.pattern }
    if (@($files).Count -eq 0) {
        throw "No matching files in $archive"
    }
    & tar -xf $archive -C $destination @files
    Write-Host "$($target.destination): extracted $(@($files).Count) expression files"
}
