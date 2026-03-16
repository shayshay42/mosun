param(
    [string]$RepoRoot = ".",
    [string]$JuliaExe = "julia",
    [string]$PythonExe = "",
    [ValidateSet("both", "auc_btumor", "auc_il6combo")]
    [string]$Only = "both",
    [int]$N = 4000,
    [int]$ThreadsPerJob = 0,
    [int]$M = 4,
    [double]$HorizonDays = 168.0,
    [double]$SaveDtDays = 0.25,
    [string]$DoseDays = "0,7,14,21,42,63,84,105,126,147",
    [string]$DoseMg = "1,2,60,60,30,30,30,30,30,30",
    [switch]$Sequential,
    [switch]$SkipPlots
)

$ErrorActionPreference = "Stop"
$RepoRoot = (Resolve-Path $RepoRoot).Path
$JuliaProject = Join-Path $RepoRoot "julia"
$Runner = Join-Path $RepoRoot "julia\run_efast_aucsum_dlbcl_ranges.jl"
$Plotter = Join-Path $RepoRoot "scripts\plot_efast_best_spd.py"
$LocalPython = Join-Path $RepoRoot ".venv\Scripts\python.exe"

if (!(Test-Path $Runner)) {
    throw "Missing Julia runner: $Runner"
}

$PythonResolved = $PythonExe
if ([string]::IsNullOrWhiteSpace($PythonResolved)) {
    if (Test-Path $LocalPython) {
        $PythonResolved = $LocalPython
    } else {
        $PythonResolved = "python"
    }
}

if ($ThreadsPerJob -le 0) {
    $ThreadsPerJob = [Math]::Max(1, [int][Math]::Floor(([Environment]::ProcessorCount - 2) / 2.0))
}
$ThreadsPerJob = [Math]::Max(1, $ThreadsPerJob)

$TotalEvals = 56 * $N
$RunDefs = @()
if ($Only -eq "both" -or $Only -eq "auc_btumor") {
    $RunDefs += [PSCustomObject]@{
        Name = "auc_btumor"
        OutDir = (Join-Path $RepoRoot "generated\figures\sensitivity\auc_btumor_dlbcl_ranges_n$N")
        IndicesCsv = "efast_auc_btumor_indices.csv"
        MetaJson = "efast_auc_btumor_meta.json"
        OutPng = "efast_auc_btumor_indices.png"
    }
}
if ($Only -eq "both" -or $Only -eq "auc_il6combo") {
    $RunDefs += [PSCustomObject]@{
        Name = "auc_il6combo"
        OutDir = (Join-Path $RepoRoot "generated\figures\sensitivity\auc_il6combo_dlbcl_ranges_n$N")
        IndicesCsv = "efast_auc_il6combo_indices.csv"
        MetaJson = "efast_auc_il6combo_meta.json"
        OutPng = "efast_auc_il6combo_indices.png"
    }
}

foreach ($run in $RunDefs) {
    if (!(Test-Path $run.OutDir)) {
        New-Item -ItemType Directory -Path $run.OutDir | Out-Null
    }
}

Write-Host "[efast-n4000] RepoRoot = $RepoRoot"
Write-Host "[efast-n4000] Julia = $JuliaExe"
Write-Host "[efast-n4000] Python = $PythonResolved"
Write-Host "[efast-n4000] Mode = $Only"
Write-Host "[efast-n4000] N = $N, M = $M, expected evals per output = $TotalEvals"
Write-Host "[efast-n4000] ThreadsPerJob = $ThreadsPerJob"
Write-Host "[efast-n4000] Sequential = $Sequential"
Write-Host "[efast-n4000] DoseDays = $DoseDays"
Write-Host "[efast-n4000] DoseMg = $DoseMg"

function Invoke-EfastRun {
    param(
        [string]$RepoRootArg,
        [string]$JuliaExeArg,
        [string]$JuliaProjectArg,
        [string]$RunnerArg,
        [string]$OutputKindArg,
        [string]$OutDirArg,
        [int]$NArg,
        [int]$MArg,
        [int]$ThreadsArg,
        [double]$HorizonArg,
        [double]$SaveDtArg,
        [string]$DoseDaysArg,
        [string]$DoseMgArg
    )

    $logPath = Join-Path $OutDirArg "run.log"
    Push-Location $RepoRootArg
    try {
        $env:JULIA_NUM_THREADS = "$ThreadsArg"
        $env:EFAST_N = "$NArg"
        $env:EFAST_M = "$MArg"
        $env:EFAST_OUTPUT_KIND = $OutputKindArg
        $env:EFAST_OUT_DIR = $OutDirArg
        $env:EFAST_HORIZON_DAYS = "$HorizonArg"
        $env:EFAST_SAVE_DT = "$SaveDtArg"
        $env:EFAST_DOSE_DAYS = $DoseDaysArg
        $env:EFAST_DOSE_MG = $DoseMgArg
        & $JuliaExeArg "--project=$JuliaProjectArg" $RunnerArg *> $logPath
        if ($LASTEXITCODE -ne 0) {
            throw "Julia run failed for $OutputKindArg with exit code $LASTEXITCODE. See $logPath"
        }
    }
    finally {
        Pop-Location
    }
}

$jobs = @()
try {
    if ($Sequential -or $RunDefs.Count -eq 1) {
        foreach ($run in $RunDefs) {
            Write-Host "[efast-n4000] Running $($run.Name) sequentially..."
            Invoke-EfastRun `
                -RepoRootArg $RepoRoot `
                -JuliaExeArg $JuliaExe `
                -JuliaProjectArg $JuliaProject `
                -RunnerArg $Runner `
                -OutputKindArg $run.Name `
                -OutDirArg $run.OutDir `
                -NArg $N `
                -MArg $M `
                -ThreadsArg $ThreadsPerJob `
                -HorizonArg $HorizonDays `
                -SaveDtArg $SaveDtDays `
                -DoseDaysArg $DoseDays `
                -DoseMgArg $DoseMg
        }
    } else {
        foreach ($run in $RunDefs) {
            Write-Host "[efast-n4000] Launching $($run.Name) in parallel..."
            $job = Start-Job -Name $run.Name -ArgumentList @(
                $RepoRoot,
                $JuliaExe,
                $JuliaProject,
                $Runner,
                $run.Name,
                $run.OutDir,
                $N,
                $M,
                $ThreadsPerJob,
                $HorizonDays,
                $SaveDtDays,
                $DoseDays,
                $DoseMg
            ) -ScriptBlock {
                param(
                    $RepoRootArg,
                    $JuliaExeArg,
                    $JuliaProjectArg,
                    $RunnerArg,
                    $OutputKindArg,
                    $OutDirArg,
                    $NArg,
                    $MArg,
                    $ThreadsArg,
                    $HorizonArg,
                    $SaveDtArg,
                    $DoseDaysArg,
                    $DoseMgArg
                )
                $ErrorActionPreference = "Stop"
                $logPath = Join-Path $OutDirArg "run.log"
                Set-Location $RepoRootArg
                $env:JULIA_NUM_THREADS = "$ThreadsArg"
                $env:EFAST_N = "$NArg"
                $env:EFAST_M = "$MArg"
                $env:EFAST_OUTPUT_KIND = $OutputKindArg
                $env:EFAST_OUT_DIR = $OutDirArg
                $env:EFAST_HORIZON_DAYS = "$HorizonArg"
                $env:EFAST_SAVE_DT = "$SaveDtArg"
                $env:EFAST_DOSE_DAYS = $DoseDaysArg
                $env:EFAST_DOSE_MG = $DoseMgArg
                & $JuliaExeArg "--project=$JuliaProjectArg" $RunnerArg *> $logPath
                if ($LASTEXITCODE -ne 0) {
                    throw "Julia run failed for $OutputKindArg with exit code $LASTEXITCODE. See $logPath"
                }
            }
            $jobs += $job
        }

        Write-Host "[efast-n4000] Waiting for parallel jobs..."
        $null = Wait-Job -Job $jobs
        foreach ($job in $jobs) {
            Receive-Job -Job $job
        }
    }
}
finally {
    if ($jobs.Count -gt 0) {
        $jobs | Remove-Job -Force -ErrorAction SilentlyContinue | Out-Null
    }
}

if (-not $SkipPlots) {
    if (!(Test-Path $Plotter)) {
        Write-Warning "Missing plot script: $Plotter"
    } else {
        foreach ($run in $RunDefs) {
            $indicesCsv = Join-Path $run.OutDir $run.IndicesCsv
            $metaJson = Join-Path $run.OutDir $run.MetaJson
            $outPng = Join-Path $run.OutDir $run.OutPng
            if ((Test-Path $indicesCsv) -and (Test-Path $metaJson)) {
                Write-Host "[efast-n4000] Plotting $($run.Name)..."
                & $PythonResolved $Plotter --indices-csv $indicesCsv --meta-json $metaJson --out $outPng
                if ($LASTEXITCODE -ne 0) {
                    throw "Plot generation failed for $($run.Name) with exit code $LASTEXITCODE"
                }
            } else {
                Write-Warning "Skipping plot for $($run.Name): missing $indicesCsv or $metaJson"
            }
        }
    }
}

Write-Host "[efast-n4000] Completed."
foreach ($run in $RunDefs) {
    Write-Host "  $($run.Name): $($run.OutDir)"
}
