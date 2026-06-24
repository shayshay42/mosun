param(
    [string]$RemoteAlias = "emad-cpu",
    [string]$CampaignRoot = "",
    [string]$SessionName = "susilo_efast_hosseini_fig5_vpop1000_20260516",
    [int]$JuliaThreads = 100,
    [int]$PruneWorkers = 24,
    [int]$NiceLevel = 5,
    [int]$ShardSize = 100,
    [int]$NResim = 60000,
    [int]$ResimLimit = 0,
    [switch]$Bootstrap
)

$ErrorActionPreference = "Stop"
if ([string]::IsNullOrWhiteSpace($CampaignRoot)) {
    $remoteUser = if ([string]::IsNullOrWhiteSpace($env:MOSUN_REMOTE_USER)) { $env:USERNAME } else { $env:MOSUN_REMOTE_USER }
    $CampaignRoot = "/projects/$remoteUser/mosun_susilo_fig5_efast"
}
$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$RemoteRepo = "$CampaignRoot/repo"
$RemoteScript = "$CampaignRoot/run_susilo_efast_hosseini_fig5_vpop1000_tmux.sh"
$OutDir = "$RemoteRepo/generated/figures/vpop_pruning/susilo_efast_hosseini2020_fig5_vpop1000_20260516"

$SyncPaths = @(
    "julia",
    "scripts",
    "assets/digitization_hosseini_2020_vpop",
    "README.md",
    "generated/model_tables",
    "generated/model_reconstruction_bundle",
    "generated/figures/reference/dlbcl_stack_parameter_ranges_all"
) | Where-Object { Test-Path (Join-Path $RepoRoot $_) }

$LargeSyncPaths = @(
    "generated/server_results/susilo_fig5_il6_dummy_efast_20260502",
    "generated/server_results/susilo_fig5_il6_day84_checkpointed_efast_20260504"
) | Where-Object { Test-Path (Join-Path $RepoRoot $_) }

$LargeSyncRequiredFiles = @{
    "generated/server_results/susilo_fig5_il6_dummy_efast_20260502" = "susilo_fig5_il6_dummy_efast_20260502/simulation_metrics_long.csv"
    "generated/server_results/susilo_fig5_il6_day84_checkpointed_efast_20260504" = "susilo_fig5_il6_day84_checkpointed_efast_20260504/significant_parameter_screen.csv"
}

if ($SyncPaths.Count -eq 0) {
    throw "No sync paths resolved from $RepoRoot"
}

ssh $RemoteAlias "mkdir -p '$RemoteRepo'"
if ($LASTEXITCODE -ne 0) { throw "remote mkdir failed with exit code $LASTEXITCODE" }

Push-Location $RepoRoot
try {
    $tmpFile = [System.IO.Path]::ChangeExtension([System.IO.Path]::GetTempFileName(), ".tar")
    try {
        & tar -cf $tmpFile @SyncPaths
        if ($LASTEXITCODE -ne 0) { throw "tar archive creation failed with exit code $LASTEXITCODE" }

        $remoteArchive = "$CampaignRoot/repo_sync_hosseini_vpop1000.tar"
        scp $tmpFile "${RemoteAlias}:$remoteArchive" | Out-Null
        if ($LASTEXITCODE -ne 0) { throw "scp upload failed with exit code $LASTEXITCODE" }

        ssh $RemoteAlias "mkdir -p '$RemoteRepo' && cd '$RemoteRepo' && tar -xf '$remoteArchive' && rm -f '$remoteArchive'"
        if ($LASTEXITCODE -ne 0) { throw "remote extract failed with exit code $LASTEXITCODE" }
    }
    finally {
        if (Test-Path $tmpFile) {
            Remove-Item -LiteralPath $tmpFile -Force
        }
    }
}
finally {
    Pop-Location
}

foreach ($relPath in $LargeSyncPaths) {
    $remoteRel = $relPath -replace "\\", "/"
    $requiredRel = $LargeSyncRequiredFiles[$relPath] -replace "\\", "/"
    ssh $RemoteAlias "test -f '$RemoteRepo/$remoteRel/$requiredRel'"
    if ($LASTEXITCODE -eq 0) {
        Write-Host "Remote already has $remoteRel/$requiredRel; skipping large sync."
        continue
    }

    $parentRel = (Split-Path -Parent $relPath) -replace "\\", "/"
    $localPath = (Resolve-Path (Join-Path $RepoRoot $relPath)).Path
    ssh $RemoteAlias "mkdir -p '$RemoteRepo/$parentRel'"
    if ($LASTEXITCODE -ne 0) { throw "remote large-sync parent mkdir failed with exit code $LASTEXITCODE" }

    Write-Host "Syncing large source root $relPath..."
    scp -r $localPath "${RemoteAlias}:$RemoteRepo/$parentRel/" | Out-Null
    if ($LASTEXITCODE -ne 0) { throw "large source sync failed for $relPath with exit code $LASTEXITCODE" }
}

if ($Bootstrap) {
    ssh $RemoteAlias "CAMPAIGN_ROOT='$CampaignRoot' bash '$RemoteRepo/scripts/server/bootstrap_emad_cpu_campaign.sh'"
    if ($LASTEXITCODE -ne 0) { throw "remote bootstrap failed with exit code $LASTEXITCODE" }
}

$LocalScript = [System.IO.Path]::ChangeExtension([System.IO.Path]::GetTempFileName(), ".sh")
try {
    $ScriptText = @"
#!/usr/bin/env bash
set -euo pipefail
while pgrep -f "$CampaignRoot/.+Pkg\.precompile" >/dev/null; do
  echo "Waiting for campaign Julia precompile to finish..."
  sleep 60
done
cd '$RemoteRepo'
export CAMPAIGN_ROOT='$CampaignRoot'
export JULIA_NUM_THREADS='$JuliaThreads'
export NICE_LEVEL='$NiceLevel'
export SUSILO_EFAST_HOSSEINI_VPOP_OUT_DIR='$OutDir'
export SUSILO_EFAST_VPOP_N_RESIM='$NResim'
export SUSILO_EFAST_VPOP_SEED='20260516'
export SUSILO_EFAST_HOSSEINI_VPOP_SHARD_SIZE='$ShardSize'
export SUSILO_EFAST_HOSSEINI_VPOP_RESIM_LIMIT='$ResimLimit'
export SUSILO_EFAST_VPOP_SELECT_N='1000'
export SUSILO_EFAST_VPOP_OUTPUT_STEM='selected_vpop1000'
export SUSILO_EFAST_VPOP_META_BASENAME='vpop1000_meta.json'
export SUSILO_EFAST_VPOP_PRUNE_SEED='20260516'
export SUSILO_EFAST_VPOP_PRUNE_INITIAL_CANDIDATES='600'
export SUSILO_EFAST_VPOP_PRUNE_STARTS='24'
export SUSILO_EFAST_VPOP_PRUNE_SWAP_ITERATIONS='3000'
export SUSILO_EFAST_VPOP_RANDOM_BASELINES='100'
export SUSILO_EFAST_VPOP_PRUNE_WORKERS='$PruneWorkers'
export SUSILO_EFAST_VPOP_USE_SINKHORN='1'
export SUSILO_EFAST_VPOP_SINKHORN_TARGET_N='501'
export SUSILO_EFAST_VPOP_SINKHORN_EPSILON='4.0'
export SUSILO_EFAST_VPOP_SINKHORN_TAU='0.65'
export SUSILO_EFAST_VPOP_SINKHORN_ITERATIONS='150'
export SUSILO_EFAST_VPOP_SINKHORN_STRENGTH='0.80'
exec bash '$RemoteRepo/scripts/server/launch_emad_cpu_susilo_efast_hosseini2020_fig5_vpop1000.sh'
"@

    $ScriptText = $ScriptText -replace "`r`n", "`n"
    [System.IO.File]::WriteAllText($LocalScript, $ScriptText, [System.Text.Encoding]::ASCII)

    scp $LocalScript "${RemoteAlias}:$RemoteScript" | Out-Null
    if ($LASTEXITCODE -ne 0) { throw "remote launcher upload failed with exit code $LASTEXITCODE" }

    ssh $RemoteAlias "chmod +x '$RemoteScript' && tmux has-session -t '$SessionName' 2>/dev/null && echo 'session-exists' || tmux new-session -d -s '$SessionName' 'bash $RemoteScript'"
    if ($LASTEXITCODE -ne 0) { throw "tmux launch failed with exit code $LASTEXITCODE" }
}
finally {
    if (Test-Path $LocalScript) {
        Remove-Item -LiteralPath $LocalScript -Force
    }
}

Write-Host "Launched or found ${RemoteAlias}:$SessionName"
Write-Host "Attach with: ssh $RemoteAlias `"tmux attach -t $SessionName`""
Write-Host "Output root: $OutDir"
