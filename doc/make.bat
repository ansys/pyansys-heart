@ECHO OFF

pushd %~dp0

REM Command file for Sphinx documentation

if "%SPHINXBUILD%" == "" (
	set SPHINXOPTS=-j auto
	set SPHINXBUILD=sphinx-build
)
set SOURCEDIR=source
set BUILDDIR=_build
set APIDIR=api

if "%1" == "" goto help
if "%1" == "clean" goto clean
if "%1" == "pdf" goto pdf  REM Add this line to handle PDF builds
if "%1" == "html" goto html

%SPHINXBUILD% >NUL 2>NUL
if errorlevel 9009 (
	echo.
	echo.The 'sphinx-build' command was not found. Make sure you have Sphinx
	echo.installed, then set the SPHINXBUILD environment variable to point
	echo.to the full path of the 'sphinx-build' executable. Alternatively you
	echo.may add the Sphinx directory to PATH.
	echo.
	echo.If you don't have Sphinx installed, grab it from
	echo.http://sphinx-doc.org/
	exit /b 1
)

%SPHINXBUILD% -M %1 %SOURCEDIR% %BUILDDIR% %SPHINXOPTS% %O%
goto end

:html
%SPHINXBUILD% -M html %SOURCEDIR% %BUILDDIR% %SPHINXOPTS% %O%
goto end


:pdf
echo Building LaTeX files...
%SPHINXBUILD% -M latex %SOURCEDIR% %BUILDDIR% %SPHINXOPTS% %O%

cd "%BUILDDIR%\latex"
echo Compiling LaTeX into PDF...
for %%f in (*.tex) do (
pdflatex "%%f" --interaction=nonstopmode)

REM Check if any PDF was generated
set "pdfFound=false"
for %%F in (*.pdf) do (
    set "pdfFound=true"
    echo PDF successfully built: %%F
)

if "%pdfFound%" == "false" (
    echo No PDF generated!
    exit /b 1
)

goto end

:clean
rmdir /s /q %BUILDDIR% > /NUL 2>&1
for /d /r %SOURCEDIR% %%d in (%APIDIR) do @if exist "%%d" rmdir /s /q "%%d"
goto end

:help
%SPHINXBUILD% -M help %SOURCEDIR% %BUILDDIR% %SPHINXOPTS% %O%

:end
popd
