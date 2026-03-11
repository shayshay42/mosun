param(
    [string]$Root = "C:\Users\shaya\mosun",
    [string]$OutRoot = "C:\Users\shaya\mosun\generated\figures\optimization\single_sample_dose_method_benchmark_20260311",
    [string]$Solver = "rodas4p",
    [string]$SaveAtDt = "1.0",
    [string]$RsEvals = "180",
    [string]$McSteps = "260",
    [string]$NmIters = "60",
    [string]$LbfgsIters = "25",
    [string]$ThreadsPerJob = "4",
    [string]$PrintEvery = "5",
    [bool]$HideWindows = $true,
    [string[]]$Objectives = @("simple", "clinical", "tracking")
)

$ErrorActionPreference = "Stop"

New-Item -ItemType Directory -Force -Path $OutRoot | Out-Null

$runner = Join-Path $Root "julia\benchmark_single_sample_dose_methods_mosun.jl"
if (-not (Test-Path $runner)) {
    throw "Runner not found: $runner"
}

$jobs = foreach ($obj in $Objectives) {
    $stdout = Join-Path $OutRoot "$obj.stdout.log"
    $stderr = Join-Path $OutRoot "$obj.stderr.log"
    $cmd = @"
`$env:JULIA_NUM_THREADS = '$ThreadsPerJob'
`$env:SINGLE_BENCH_OBJECTIVE = '$obj'
`$env:SINGLE_BENCH_SOLVER = '$Solver'
`$env:SINGLE_BENCH_SAVEAT_DT = '$SaveAtDt'
`$env:SINGLE_BENCH_RS_EVALS = '$RsEvals'
`$env:SINGLE_BENCH_MC_STEPS = '$McSteps'
`$env:SINGLE_BENCH_NM_ITERS = '$NmIters'
`$env:SINGLE_BENCH_LBFGS_ITERS = '$LbfgsIters'
`$env:SINGLE_BENCH_PRINT_EVERY = '$PrintEvery'
Set-Location '$Root'
julia --project=./julia julia/benchmark_single_sample_dose_methods_mosun.jl
"@
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

$jobs | Wait-Process

python (Join-Path $Root "scripts\plot_single_sample_dose_method_benchmark.py")
