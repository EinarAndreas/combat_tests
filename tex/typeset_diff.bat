@echo off

if "%1"=="" goto help

set typeset=D:\Projects\_Appl_\_TeX_\typeset.bat
set manus=manuscript
set diff=diff

latexdiff %1 %manus%.tex > %diff%.tex
%typeset% %diff%
pause

goto :EOF


:help
echo :
echo :
pause
goto :EOF
