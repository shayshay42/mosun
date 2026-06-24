param(
    [string]$RemoteAlias = "emad-cpu",
    [string]$CampaignRoot = "",
    [switch]$Bootstrap
)

$ErrorActionPreference = "Stop"
if ([string]::IsNullOrWhiteSpace($CampaignRoot)) {
    $remoteUser = if ([string]::IsNullOrWhiteSpace($env:MOSUN_REMOTE_USER)) { $env:USERNAME } else { $env:MOSUN_REMOTE_USER }
    $CampaignRoot = "/projects/$remoteUser/mosun_susilo_fig5_efast"
}
$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$RemoteRepo = "$CampaignRoot/repo"

$CandidatePaths = @(
    "julia",
    "scripts",
    "README.md",
    "generated/model_tables",
    "generated/model_reconstruction_bundle",
    "generated/figures/reference/dlbcl_stack_parameter_ranges_all",
    "generated/figures/vpop_pruning/susilo_efast_hosseini2020_fig5_vpop250_20260505/high_importance_parameters.csv",
    "generated/figures/vpop_pruning/susilo_efast_hosseini2020_fig5_vpop250_20260505/trajectory_resim_candidate_ids.csv",
    "generated/figures/vpop_pruning/susilo_efast_hosseini2020_fig5_vpop250_20260505/trajectory_resim_candidate_parameters.csv",
    "generated/figures/vpop_pruning/susilo_efast_hosseini2020_fig5_vpop250_20260505/vpop250_prefilter_meta.json"
)

$Paths = @($CandidatePaths | Where-Object { Test-Path (Join-Path $RepoRoot $_) })
if ($Paths.Count -eq 0) {
    throw "No sync paths resolved from $RepoRoot"
}

ssh $RemoteAlias "mkdir -p '$RemoteRepo'"
if ($LASTEXITCODE -ne 0) { throw "remote mkdir failed with exit code $LASTEXITCODE" }

Push-Location $RepoRoot
try {
    $tmpFile = [System.IO.Path]::ChangeExtension([System.IO.Path]::GetTempFileName(), ".tar")
    try {
        & tar -cf $tmpFile @Paths
        if ($LASTEXITCODE -ne 0) { throw "tar archive creation failed with exit code $LASTEXITCODE" }

        $remoteArchive = "$CampaignRoot/repo_sync_hosseini_vpop.tar"
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

if ($Bootstrap) {
    ssh $RemoteAlias "CAMPAIGN_ROOT='$CampaignRoot' bash '$RemoteRepo/scripts/server/bootstrap_emad_cpu_campaign.sh'"
    if ($LASTEXITCODE -ne 0) { throw "remote bootstrap failed with exit code $LASTEXITCODE" }
}

Write-Host "Synced eFAST-derived Hosseini Fig.5 VPop workflow to ${RemoteAlias}:$RemoteRepo"
