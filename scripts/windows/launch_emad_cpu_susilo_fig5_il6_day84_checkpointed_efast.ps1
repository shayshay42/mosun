param(
    [string]$RemoteAlias = "emad-cpu",
    [string]$CampaignRoot = "",
    [string]$SessionName = "susilo_fig5_il6_day84_ckpt_20260504",
    [int]$JuliaThreads = 16,
    [int]$NiceLevel = 10,
    [int]$SmokeN = 64,
    [string]$SmokeSeeds = "20260504",
    [int]$FullN = 1024,
    [string]$FullSeeds = "20260504,20260505",
    [switch]$NoSmoke,
    [switch]$NoFull
)

$ErrorActionPreference = "Stop"
if ([string]::IsNullOrWhiteSpace($CampaignRoot)) {
    $remoteUser = if ([string]::IsNullOrWhiteSpace($env:MOSUN_REMOTE_USER)) { $env:USERNAME } else { $env:MOSUN_REMOTE_USER }
    $CampaignRoot = "/projects/$remoteUser/mosun_susilo_fig5_efast"
}
$RemoteRepo = "$CampaignRoot/repo"
$RunSmoke = if ($NoSmoke) { "0" } else { "1" }
$RunFull = if ($NoFull) { "0" } else { "1" }
$RemoteScript = "$CampaignRoot/run_susilo_fig5_il6_day84_ckpt_efast_tmux.sh"

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
export RUN_SMOKE='$RunSmoke'
export RUN_FULL='$RunFull'
export SUSILO_FIG5_CKPT_EFAST_SMOKE_N='$SmokeN'
export SUSILO_FIG5_CKPT_EFAST_SMOKE_SEEDS='$SmokeSeeds'
export SUSILO_FIG5_CKPT_EFAST_FULL_N='$FullN'
export SUSILO_FIG5_CKPT_EFAST_FULL_SEEDS='$FullSeeds'
exec bash '$RemoteRepo/scripts/server/launch_emad_cpu_susilo_fig5_il6_day84_checkpointed_efast.sh'
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
