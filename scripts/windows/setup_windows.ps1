param(
    [string]$RepoRoot = ".",
    [string]$PythonExe = "python",
    [string]$JuliaExe = "julia",
    [switch]$SkipPython,
    [switch]$SkipJulia
)

$ErrorActionPreference = "Stop"
$RepoRoot = (Resolve-Path $RepoRoot).Path
$VenvPython = Join-Path $RepoRoot ".venv\Scripts\python.exe"
$Requirements = Join-Path $RepoRoot "requirements.txt"
$JuliaProject = Join-Path $RepoRoot "julia"

Write-Host "RepoRoot: $RepoRoot"

if (-not $SkipPython) {
    Write-Host "[Python] creating/updating .venv and installing requirements"
    & $PythonExe -m venv (Join-Path $RepoRoot ".venv")
    & $VenvPython -m pip install --upgrade pip
    & $VenvPython -m pip install -r $Requirements
}

if (-not $SkipJulia) {
    Write-Host "[Julia] instantiating environment at $JuliaProject"
    & $JuliaExe "--project=$JuliaProject" -e "using Pkg; Pkg.instantiate(); Pkg.precompile()"
}

Write-Host "Setup completed."
Write-Host "Use scripts/windows/run_vpop.ps1 to execute VPop generation."
