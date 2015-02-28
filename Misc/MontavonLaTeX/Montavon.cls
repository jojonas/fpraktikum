%LaTeX2e Klasse
%Benjamin Montavon
%13. Febaruar 2013
%Version 0.1
%-----------------------------------------------------------------------------------------------------------------------
% Klasseninformationen
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{Montavon}[2013/02/13]
%-----------------------------------------------------------------------------------------------------------------------
%Basiert auf Article
\LoadClass[a4paper,12pt]{article}
% PDF-Kompression
\pdfminorversion=5 \pdfcompresslevel=9 \pdfobjcompresslevel=9
%Seitenraender anpassen
\RequirePackage[a4paper, bottom=3.5cm]{geometry}
% Benoetigte Pakete
\RequirePackage[T1]{fontenc}
\RequirePackage[utf8]{inputenc}
\RequirePackage{lmodern}
\RequirePackage[ngerman]{babel}
\RequirePackage{fancyhdr}
\RequirePackage{graphicx}
\RequirePackage{amsmath}
\RequirePackage{amssymb}
\RequirePackage{upgreek}
\RequirePackage[pdfborder={0 0 0},breaklinks=true]{hyperref}
%-----------------------------------------------------------------------------------------------------------------------
% Seitenformatierungsoptionen
\pagestyle{fancy}
\numberwithin{equation}{section}
\renewcommand{\theequation}{\arabic{section}.\arabic{equation}}
%-----------------------------------------------------------------------------------------------------------------------
% Neu definierte Kurzbefehle
% Formatierung
\newcommand{\heads}[3]{
	\fancyhead[L]{#1}
	\fancyhead[C]{#2}
	\fancyhead[R]{#3}
}
\newcommand{\legendde}{\InputIfFileExists{legendde.tex}{}{\ClassError{Montavon}{Inputfile legendde.tex not found. Make sure it is placed in the directory of the class.Press Enter to continue without any input.}{}}}
% Mathe
\newcommand{\li}{\left(}
\newcommand{\re}{\right)}
\newcommand{\pdiff}[2]{\frac{\partial #1}{\partial #2}}
\newcommand{\diff}[2]{\frac{\mathrm d #1}{\mathrm d #2}}
\newcommand{\di}{\mathrm d} 
\newcommand{\poi}[1]{\left\{ #1 \right\}}
\newcommand{\ord}[1]{\mathcal O \left( #1 \right)}
\newcommand{\reell}{\mathbb R}
\newcommand{\svec}[3]{\begin{pmatrix} #1 \\ #2 \\ #3 \end{pmatrix}}
\newcommand{\tdvec}[2]{\begin{pmatrix} #1 \\ #2 \end{pmatrix} }
\newcommand{\Kern}{\mathrm{Kern} \, }
\newcommand{\Norm}[1]{\left| #1 \right| }
\newcommand{\Normm}[1]{\left| \left| #1 \right| \right| }
\newcommand{\vecb}{\boldsymbol}
\newcommand{\nullvec}{\boldsymbol{\it 0}} 
\newcommand{\matb}[1]{\mathbf #1} 
\newcommand{\vecbe}[1]{\vecb{\hat e}_{#1}}
\newcommand{\uuline}[1]{\underline{\underline{ #1}}}
\newcommand{\im}{\mathrm i}
\newcommand{\expe}{\mathrm e} 
\newcommand{\unit}[1]{\, \mathrm{#1} }
\newcommand{\unitb}[1]{\, \mathbf{#1} }
\newcommand{\pow}[1]{\cdot 10^{#1} } 
\renewcommand{\\}{\newline \noindent} 
\newcommand{\degree}{^{\circ}}