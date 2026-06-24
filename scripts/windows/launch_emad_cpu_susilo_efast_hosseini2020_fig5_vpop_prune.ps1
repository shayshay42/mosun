param(
    [string]$RemoteAlias = "emad-cpu",
    [string]$CampaignRoot = "",
    [string]$SessionName = "susilo_efast_hosseini_fig5_vpop_prune_20260505",
    [int]$Workers = 24,
    [int]$InitialCandidates = 240,
    [int]$Starts = 6,
    [int]$SwapIterations = 900,
    [int]$RandomBaselines = 32,
    [int]$NiceLevel = 5
)

$ErrorActionPreference = "Stop"
if ([string]::IsNullOrWhiteSpace($CampaignRoot)) {
    $remoteUser = if ([string]::IsNullOrWhiteSpace($env:MOSUN_REMOTE_USER)) { $env:USERNAME } else { $env:MOSUN_REMOTE_USER }
    $CampaignRoot = "/projects/$remoteUser/mosun_susilo_fig5_efast"
}
$RemoteRepo = "$CampaignRoot/repo"
$RemoteScript = "$CampaignRoot/run_susilo_efast_hosseini_fig5_vpop_prune_tmux.sh"

$LocalScript = [System.IO.Path]::ChangeExtension([System.IO.Path]::GetTempFileName(), ".sh")
try {
    $ScriptText = @"
#!/usr/bin/env bash
set -euo pipefail
cd '$RemoteRepo'
export CAMPAIGN_ROOT='$CampaignRoot'
export NICE_LEVEL='$NiceLevel'
export SUSILO_EFAST_VPOP_PRUNE_WORKERS='$Workers'
export SUSILO_EFAST_VPOP_PRUNE_INITIAL_CANDIDATES='$InitialCandidates'
export SUSILO_EFAST_VPOP_PRUNE_STARTS='$Starts'
export SUSILO_EFAST_VPOP_PRUNE_SWAP_ITERATIONS='$SwapIterations'
export SUSILO_EFAST_VPOP_RANDOM_BASELINES='$RandomBaselines'
exec bash '$RemoteRepo/scripts/server/launch_emad_cpu_susilo_efast_hosseini2020_fig5_vpop_prune.sh'
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
