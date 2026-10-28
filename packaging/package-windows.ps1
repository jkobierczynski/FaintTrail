# Usage: packaging/package-windows.ps1 -Version X.Y.Z -BuildDir build -OutDir dist
# Makes OutDir/fainttrail-X.Y.Z-windows-x64.zip from a Release build with the GUI:
# fainttrail.exe, fainttrail-gui.exe, the Qt and Visual C++ runtime (copied by the
# deploy_gui target) and the OpenMP runtime.
param(
  [Parameter(Mandatory = $true)][string]$Version,
  [Parameter(Mandatory = $true)][string]$BuildDir,
  [Parameter(Mandatory = $true)][string]$OutDir
)
$ErrorActionPreference = 'Stop'
$src = Split-Path -Parent $PSScriptRoot
$rel = Join-Path $BuildDir 'Release'
New-Item -ItemType Directory -Force -Path $OutDir | Out-Null
$name = "fainttrail-$Version-windows-x64"
$stage = Join-Path ([IO.Path]::GetTempPath()) $name
if (Test-Path $stage) { Remove-Item -Recurse -Force $stage }
New-Item -ItemType Directory -Path $stage | Out-Null

# programs, Qt DLLs and plugin folders; not the tests, import libraries or symbols
Get-ChildItem $rel | Where-Object {
  $_.Name -notmatch '^fainttrail_tests' -and $_.Extension -notin @('.lib', '.exp', '.pdb', '.ilk')
} | Copy-Item -Destination $stage -Recurse

# OpenMP runtime (vcomp140.dll) is not copied by windeployqt
$vs = & "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe" -latest -products * -property installationPath
$omp = Get-ChildItem -Path (Join-Path $vs 'VC\Redist\MSVC') -Recurse -Filter vcomp140.dll |
  Where-Object { $_.FullName -match '\\x64\\Microsoft\.VC\d+\.OpenMP\\' } | Select-Object -First 1
if (-not $omp) { throw 'vcomp140.dll not found in the Visual Studio redistributables' }
Copy-Item $omp.FullName $stage

foreach ($f in 'LICENSE', 'README.md', 'BUILDING.md', 'DEVELOPMENT.md') { Copy-Item (Join-Path $src $f) $stage }
Copy-Item -Recurse (Join-Path $src 'docs') (Join-Path $stage 'docs')
@"
FaintTrail $Version for Windows x64

fainttrail-gui.exe   the GUI
fainttrail.exe       command-line tool (run "fainttrail --help" in a Command Prompt)

Unpack the whole folder; the DLLs must stay next to the programs. The CUDA backend
is built in (NVIDIA GPUs from the GTX 16xx / RTX 20xx generation on, driver 570 or
newer) and is used automatically; without an NVIDIA GPU the CPU is used.
Optional: ASTAP (https://www.hnsky.org/astap.htm) for plate solving. curl, used for
the Gaia photometry, is part of Windows 10 and 11.
"@ | Set-Content -Encoding UTF8 (Join-Path $stage 'README-BINARY.txt')

$zip = Join-Path $OutDir "$name.zip"
if (Test-Path $zip) { Remove-Item $zip }
Compress-Archive -Path $stage -DestinationPath $zip
Get-Item $zip | Format-List Name, Length
