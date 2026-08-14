$projectRoot = Split-Path -Parent $PSScriptRoot
$raw = Join-Path $projectRoot "data\00_raw\spatial\GSE278687_RAW.tar"
$destination = Join-Path $projectRoot "data\01_unpacked\spatial\GSE278687_spatial"
New-Item -ItemType Directory -Force -Path $destination | Out-Null

$nested = tar -tf $raw | Where-Object { $_ -match '_spatial\.tar\.gz$' }
foreach ($entry in $nested) {
    $sample = $entry -replace '_spatial\.tar\.gz$', ''
    $temp = Join-Path $destination "_temp_$sample"
    New-Item -ItemType Directory -Force -Path $temp | Out-Null
    & tar -xf $raw -C $temp $entry
    $target = Join-Path $destination $sample
    New-Item -ItemType Directory -Force -Path $target | Out-Null
    & tar -xzf (Join-Path $temp $entry) -C $target
    Remove-Item -LiteralPath $temp -Recurse -Force
}
Write-Host "Extracted spatial coordinate packages: $($nested.Count)"
