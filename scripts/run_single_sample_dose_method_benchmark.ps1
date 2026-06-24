param(
    [string]$Root = "",
    [string]$OutRoot = "generated\figures\optimization\single_sample_dose_method_benchmark_20260311",
    [string]$Solver = "rodas4p",
    [string]$SaveAtDt = "1.0",
    [string]$RsEvals = "180",
    [string]$McSteps = "260",
    [string]$NmIters = "60",
    [string]$LbfgsIters = "25",
    [string]$ThreadsPerJob = "4",
    [string]$PrintEvery = "5",
    [bool]$HideWindows = $true,
    [switch]$Foreground,
    [switch]$Resume,
    [switch]$SkipPlot,
    [string]$MethodFilter = "",
    [string]$RandomSeed = "20260311",
    [string[]]$Objectives = @("simple", "clinical", "tracking")
)

$ErrorActionPreference = "Stop"

$Root = if ([string]::IsNullOrWhiteSpace($Root)) {
    (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
} else {
    (Resolve-Path $Root).Path
}
$resolvedOutRoot = if ([System.IO.Path]::IsPathRooted($OutRoot)) { $OutRoot } else { Join-Path $Root $OutRoot }
New-Item -ItemType Directory -Force -Path $resolvedOutRoot | Out-Null

$runner = Join-Path $Root "julia\benchmark_single_sample_dose_methods_mosun.jl"
if (-not (Test-Path $runner)) {
    throw "Runner not found: $runner"
}

$jobs = foreach ($obj in $Objectives) {
    $stdout = Join-Path $resolvedOutRoot "$obj.stdout.log"
    $stderr = Join-Path $resolvedOutRoot "$obj.stderr.log"
    $cmd = @"
`$env:JULIA_NUM_THREADS = '$ThreadsPerJob'
`$env:SINGLE_BENCH_OBJECTIVE = '$obj'
`$env:SINGLE_BENCH_OUT_ROOT = '$resolvedOutRoot'
`$env:SINGLE_BENCH_SOLVER = '$Solver'
`$env:SINGLE_BENCH_SAVEAT_DT = '$SaveAtDt'
`$env:SINGLE_BENCH_RS_EVALS = '$RsEvals'
`$env:SINGLE_BENCH_MC_STEPS = '$McSteps'
`$env:SINGLE_BENCH_NM_ITERS = '$NmIters'
`$env:SINGLE_BENCH_LBFGS_ITERS = '$LbfgsIters'
`$env:SINGLE_BENCH_PRINT_EVERY = '$PrintEvery'
`$env:SINGLE_BENCH_RANDOM_SEED = '$RandomSeed'
`$env:SINGLE_BENCH_METHOD_FILTER = '$MethodFilter'
`$env:SINGLE_BENCH_RESUME = '$(if ($Resume.IsPresent) { 1 } else { 0 })'
Set-Location '$Root'
julia --project=./julia julia/benchmark_single_sample_dose_methods_mosun.jl
"@
    if ($Foreground.IsPresent) {
        Write-Host "running objective=$obj in foreground"
        & powershell -NoProfile -ExecutionPolicy Bypass -Command $cmd
        continue
    }
    $startArgs = @{
        FilePath = "powershell"
        ArgumentList = @("-NoProfile", "-ExecutionPolicy", "Bypass", "-Command", $cmd)
        WorkingDirectory = $Root
        RedirectStandardOutput = $stdout
        RedirectStandardError = $stderr
        PassThru = $true
    }
    if ($HideWindows) {
        $startArgs["WindowStyle"] = "Hidden"
    }
    Start-Process @startArgs
}

if (-not $Foreground.IsPresent) {
    $jobs | Wait-Process
}

$shouldPlot = -not $SkipPlot.IsPresent -and @($Objectives) -contains "simple" -and @($Objectives) -contains "clinical" -and @($Objectives) -contains "tracking"
if ($shouldPlot) {
    $env:SINGLE_BENCH_RESULTS_ROOT = $resolvedOutRoot
    python (Join-Path $Root "scripts\plot_single_sample_dose_method_benchmark.py")
}
