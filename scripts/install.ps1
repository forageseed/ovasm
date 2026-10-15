[CmdletBinding()]
param(
    [string]$PackageRoot = $PSScriptRoot,
    [string]$InstallDir = (Join-Path $env:LOCALAPPDATA 'Programs\ovasm'),
    [switch]$NoPathUpdate
)
$ErrorActionPreference = 'Stop'
function Get-Sha256([string]$LiteralPath) {
    $stream = [IO.File]::OpenRead($LiteralPath)
    $sha = [Security.Cryptography.SHA256]::Create()
    try { return [BitConverter]::ToString($sha.ComputeHash($stream)).Replace('-', '') }
    finally { $stream.Dispose(); $sha.Dispose() }
}
if ($env:PROCESSOR_ARCHITECTURE -ne 'AMD64' -and $env:PROCESSOR_ARCHITEW6432 -ne 'AMD64') {
    throw 'This package requires Windows x86_64.'
}
$PackageRoot = (Resolve-Path -LiteralPath $PackageRoot).Path
$source = Join-Path $PackageRoot 'bin\ovasm.exe'
$manifest = Join-Path $PackageRoot 'FILES-SHA256.txt'
if (!(Test-Path -LiteralPath $source -PathType Leaf) -or !(Test-Path -LiteralPath $manifest -PathType Leaf)) {
    throw 'Extract the complete official ZIP before running install.cmd.'
}
$entries = @(Get-Content -LiteralPath $manifest | Where-Object { $_ -match '^([0-9a-fA-F]{64})\s+bin/ovasm\.exe$' })
if ($entries.Count -ne 1) { throw 'Missing or duplicate executable SHA-256 entry.' }
$expected = ($entries[0] -split '\s+')[0]
$actual = Get-Sha256 $source
if ($actual -ine $expected) { throw 'SHA-256 mismatch; nothing installed.' }
$version = & $source --version
if ($LASTEXITCODE -ne 0 -or $version -ne 'ovasm 0.2.0') { throw 'Unexpected executable version.' }
New-Item -ItemType Directory -Force -Path $InstallDir | Out-Null
$target = Join-Path $InstallDir 'ovasm.exe'
if ([IO.Path]::GetFullPath($source) -ine [IO.Path]::GetFullPath($target)) {
    Copy-Item -LiteralPath $source -Destination $target -Force
}
if ((Get-Sha256 $target) -ine $expected) {
    throw 'Installed executable hash mismatch.'
}
if (!$NoPathUpdate) {
    $userPath = [Environment]::GetEnvironmentVariable('PATH', 'User')
    $entries = @($userPath -split ';' | Where-Object { $_ })
    if ($entries -inotcontains $InstallDir) {
        $newPath = (@($entries) + @($InstallDir)) -join ';'
        [Environment]::SetEnvironmentVariable('PATH', $newPath, 'User')
    }
    if (@($env:PATH -split ';') -inotcontains $InstallDir) { $env:PATH += ";$InstallDir" }
}
Write-Output "Installed $version to $target"
if (!$NoPathUpdate) { Write-Output 'Open a new terminal, then run: ovasm --version' }
