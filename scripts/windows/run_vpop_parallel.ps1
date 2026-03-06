param(
    [string]$RepoRoot = ".",
    [string]$PythonExe = "",
    [string]$JuliaExe = "julia",
    [string]$RunTagPrefix = "vpop_windows_parallel",
    [int]$NCandidates = 384,
    [int]$NSelect = 140,
    [int]$NRandomSubsets = 100000,
    [int]$NumWorkers = 0,
    [int]$Seed = 20260228,
    [string]$SpreadMode = "primary",
    [int]$IncludeStandardRegimens = 1,
    [switch]$PublishBestToAssets
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

if ($NumWorkers -le 0) {
    $NumWorkers = [Math]::Max(1, [Environment]::ProcessorCount - 1)
}
$NumWorkers = [Math]::Max(1, $NumWorkers)
$DrawsPerWorker = [int][Math]::Ceiling($NRandomSubsets / [double]$NumWorkers)
$BaseRunTag = "${RunTagPrefix}_base"

Write-Host "[parallel-vpop] RepoRoot = $RepoRoot"
Write-Host "[parallel-vpop] Python = $PythonResolved"
Write-Host "[parallel-vpop] Julia = $JuliaExe"
Write-Host "[parallel-vpop] Workers = $NumWorkers (draws/worker = $DrawsPerWorker)"
Write-Host "[parallel-vpop] Base run tag = $BaseRunTag"

Push-Location $RepoRoot
try {
    # 1) Build a single candidate/simulation metric table once.
    $BaseArgs = @(
        "scripts/sample_and_prune_vpop.py",
        "--run-tag", $BaseRunTag,
        "--n-candidates", "$NCandidates",
        "--n-select", "$NSelect",
        "--n-random-subsets", "1",
        "--seed", "$Seed",
        "--spread-mode", $SpreadMode,
        "--include-standard-regimens", "$IncludeStandardRegimens",
        "--julia-bin", $JuliaExe,
        "--julia-project", $JuliaProject
    )
    Write-Host "[parallel-vpop] Running base simulation pass..."
    & $PythonResolved @BaseArgs
    if ($LASTEXITCODE -ne 0) {
        throw "Base VPop run failed with exit code $LASTEXITCODE"
    }

    $MetricsCsv = Join-Path $RepoRoot "generated\vpop_pruning\$BaseRunTag\sim\phase1_metrics_julia.csv"
    if (!(Test-Path $MetricsCsv)) {
        throw "Expected metrics csv not found: $MetricsCsv"
    }

    # 2) Launch parallel subset-search workers over the same candidate/sim table.
    $jobs = @()
    $WorkerTags = @()
    for ($w = 1; $w -le $NumWorkers; $w++) {
        $tag = "{0}_w{1:D2}" -f $RunTagPrefix, $w
        $subsetSeed = $Seed + 1000 + $w
        $WorkerTags += $tag
        Write-Host "[parallel-vpop] Launch worker $w/$NumWorkers tag=$tag subset_seed=$subsetSeed"

        $job = Start-Job -Name $tag -ArgumentList @(
            $RepoRoot,
            $PythonResolved,
            $tag,
            $NCandidates,
            $NSelect,
            $DrawsPerWorker,
            $Seed,
            $subsetSeed,
            $SpreadMode,
            $IncludeStandardRegimens,
            $JuliaExe,
            $JuliaProject,
            $MetricsCsv
        ) -ScriptBlock {
            param(
                $RepoRootArg,
                $PyArg,
                $RunTagArg,
                $NCandArg,
                $NSelArg,
                $NDrawArg,
                $SeedArg,
                $SubsetSeedArg,
                $SpreadArg,
                $IncludeStdArg,
                $JuliaExeArg,
                $JuliaProjArg,
                $MetricsArg
            )
            Set-Location $RepoRootArg
            $Args = @(
                "scripts/sample_and_prune_vpop.py",
                "--run-tag", "$RunTagArg",
                "--n-candidates", "$NCandArg",
                "--n-select", "$NSelArg",
                "--n-random-subsets", "$NDrawArg",
                "--seed", "$SeedArg",
                "--subset-seed", "$SubsetSeedArg",
                "--spread-mode", "$SpreadArg",
                "--include-standard-regimens", "$IncludeStdArg",
                "--julia-bin", "$JuliaExeArg",
                "--julia-project", "$JuliaProjArg",
                "--skip-sim",
                "--metrics-csv", "$MetricsArg"
            )
            & $PyArg @Args
            if ($LASTEXITCODE -ne 0) {
                throw "Worker run failed for $RunTagArg exit code=$LASTEXITCODE"
            }
        }
        $jobs += $job
    }

    Write-Host "[parallel-vpop] Waiting for workers..."
    $null = Wait-Job -Job $jobs
    foreach ($job in $jobs) {
        Receive-Job -Job $job
    }

    # 3) Pick best worker by objective.best_score.
    $workerSummaries = @()
    foreach ($tag in $WorkerTags) {
        $summaryPath = Join-Path $RepoRoot "generated\vpop_pruning\$tag\pruning_summary.json"
        if (!(Test-Path $summaryPath)) {
            throw "Missing worker summary: $summaryPath"
        }
        $obj = Get-Content $summaryPath -Raw | ConvertFrom-Json
        $score = [double]$obj.objective.best_score
        $workerSummaries += [PSCustomObject]@{
            run_tag = $tag
            best_score = $score
            summary_path = $summaryPath
            selected_patients_csv = (Join-Path $RepoRoot "generated\vpop_pruning\$tag\selected_patients.csv")
        }
    }
    $best = $workerSummaries | Sort-Object -Property best_score | Select-Object -First 1

    $combinedSummaryPath = Join-Path $RepoRoot "generated\vpop_pruning\${RunTagPrefix}_parallel_summary.json"
    $combined = [PSCustomObject]@{
        created_utc = (Get-Date).ToUniversalTime().ToString("o")
        run_tag_prefix = $RunTagPrefix
        base_run_tag = $BaseRunTag
        n_workers = $NumWorkers
        draws_per_worker = $DrawsPerWorker
        total_requested_draws = $NRandomSubsets
        seed = $Seed
        spread_mode = $SpreadMode
        metrics_csv = $MetricsCsv
        best_worker = $best
        workers = $workerSummaries
    }
    ($combined | ConvertTo-Json -Depth 6) | Set-Content $combinedSummaryPath
    Write-Host "[parallel-vpop] Wrote $combinedSummaryPath"
    Write-Host "[parallel-vpop] Best run: $($best.run_tag) score=$($best.best_score)"

    if ($PublishBestToAssets) {
        $assetsDir = Join-Path $RepoRoot "assets\generated_vpop"
        if (!(Test-Path $assetsDir)) {
            New-Item -ItemType Directory -Path $assetsDir | Out-Null
        }
        Copy-Item -Path $best.selected_patients_csv -Destination (Join-Path $assetsDir "selected_patients.csv") -Force
        Copy-Item -Path $best.summary_path -Destination (Join-Path $assetsDir "pruning_summary.json") -Force
        Write-Host "[parallel-vpop] Published best selected cohort to assets/generated_vpop/"
    }
}
finally {
    Pop-Location
    if ($jobs) {
        $jobs | Remove-Job -Force -ErrorAction SilentlyContinue | Out-Null
    }
}
