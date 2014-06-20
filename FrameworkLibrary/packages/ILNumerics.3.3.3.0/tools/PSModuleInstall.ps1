
$ScriptDir = Split-Path $MyInvocation.MyCommand.Path

# copy assemblies to this folder
Copy-Item "$ScriptDir\..\lib\*.*" -Destination $ScriptDir
Copy-Item "$ScriptDir\..\content\*.*" -Destination $ScriptDir
