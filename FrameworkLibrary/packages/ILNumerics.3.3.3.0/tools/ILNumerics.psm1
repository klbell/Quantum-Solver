
$ScriptDir = Split-Path $MyInvocation.MyCommand.Path

Add-Type -Path $ScriptDir\ILNumerics32.dll
               
. $ScriptDir\ILAccelerators.ps1
. $ScriptDir\ILUsing.ps1