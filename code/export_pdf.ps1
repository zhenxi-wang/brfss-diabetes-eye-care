param([string[]]$DocumentNames = @('Manuscript', 'Supplement'))
$ErrorActionPreference = 'Stop'
$projectRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$wordApp = $null
try {
    $wordApp = New-Object -ComObject Word.Application
    $wordApp.Visible = $false
    $wordApp.DisplayAlerts = 0
    $wordApp.AutomationSecurity = 3
    foreach ($documentName in $DocumentNames) {
        $inputFile = Join-Path $projectRoot ('deliverables\' + $documentName + '.docx')
        $outputFile = Join-Path $projectRoot ('deliverables\' + $documentName + '.pdf')
        $document = $wordApp.Documents.Open($inputFile, $false, $true, $false)
        try {
            $document.Repaginate()
            $document.ExportAsFixedFormat($outputFile, 17)
            Write-Output ($documentName + ' exported as PDF')
        }
        finally { $document.Close(0) }
    }
}
finally {
    if ($null -ne $wordApp) { $wordApp.Quit(0) }
}
