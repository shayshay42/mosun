param(
    [string]$RemoteAlias = "emad-cpu",
    [string]$CampaignRoot = "",
    [string]$SessionName = "susilo_efast_hosseini_fig5_vpop_resim_20260505",
    [int]$JuliaThreads = 100,
    [int]$NiceLevel = 5,
    [int]$ShardSize = 100,
    [int]$ResimLimit = 0
)

$ErrorActionPreference = "Stop"
if ([string]::IsNullOrWhiteSpace($CampaignRoot)) {
    $remoteUser = if ([string]::IsNullOrWhiteSpace($env:MOSUN_REMOTE_USER)) { $env:USERNAME } else { $env:MOSUN_REMOTE_USER }
    $CampaignRoot = "/projects/$remoteUser/mosun_susilo_fig5_efast"
}
$RemoteRepo = "$CampaignRoot/repo"
$RemoteScript = "$CampaignRoot/run_susilo_efast_hosseini_fig5_vpop_resim_tmux.sh"

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
export SUSILO_EFAST_HOSSEINI_VPOP_SHARD_SIZE='$ShardSize'
export SUSILO_EFAST_HOSSEINI_VPOP_RESIM_LIMIT='$ResimLimit'
exec bash '$RemoteRepo/scripts/server/launch_emad_cpu_susilo_efast_hosseini2020_fig5_vpop_resim.sh'
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
