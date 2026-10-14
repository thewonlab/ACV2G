@echo off
Rscript -e "status <- ACV2G::acv2g_cli(commandArgs(trailingOnly = TRUE)); quit(save = 'no', status = status, runLast = FALSE)" %*
exit /b %ERRORLEVEL%
