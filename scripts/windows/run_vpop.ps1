param(
    [string]$RepoRoot = ".",
    [string]$PythonExe = "",
    [string]$JuliaExe = "julia",
    [string]$RunTag = "vpop_windows_run",
    [int]$NCandidates = 72,
    [int]$NSelect = 50,
    [int]$NRandomSubsets = 1500,
    [string]$SpreadMode = "primary",
    [int]$IncludeStandardRegimens = 1,
    [switch]$SkipSim,
    [string]$MetricsCsv = ""
)

$ErrorActionPreference = "Stop"
$RepoRoot = (Resolve-Path $RepoRoot).Path
$LocalPython = Join-Path $RepoRoot ".venv\Scripts\python.exe"
$PythonResolved = $PythonExe
if ([string]::IsNullOrWhiteSpace($PythonResolved)) {
    if (Test-Path $LocalPython) {
        $PythonResolved = $LocalPython
    } else {
        $PythonResolved = "python"
    }
}

$JuliaProject = Join-Path $RepoRoot "julia"

$ArgsList = @(
    "scripts/sample_and_prune_vpop.py",
    "--run-tag", $RunTag,
    "--n-candidates", "$NCandidates",
    "--n-select", "$NSelect",
    "--n-random-subsets", "$NRandomSubsets",
    "--spread-mode", $SpreadMode,
    "--include-standard-regimens", "$IncludeStandardRegimens",
    "--julia-bin", $JuliaExe,
    "--julia-project", $JuliaProject
)

if ($SkipSim) {
    if ([string]::IsNullOrWhiteSpace($MetricsCsv)) {
        throw "-SkipSim was set, but -MetricsCsv is empty."
    }
    $ArgsList += @("--skip-sim", "--metrics-csv", $MetricsCsv)
}

Push-Location $RepoRoot
try {
    Write-Host "Running VPop generation with run tag: $RunTag"
    & $PythonResolved @ArgsList
} finally {
    Pop-Location
}

$Summary = Join-Path $RepoRoot "generated\vpop_pruning\$RunTag\pruning_summary.json"
Write-Host "Done. Summary: $Summary"
