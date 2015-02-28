%LaTeX2e Klasse
%Benjamin Montavon // Jonas Lieb
%10. M�rz 2013
%Version 0.1
%-----------------------------------------------------------------------------------------------------------------------
% Klasseninformationen
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{MontavonBeamer}[2013/03/10]
%-----------------------------------------------------------------------------------------------------------------------
%Basiert auf Beamer (Vorlage: http://www.informatik.uni-freiburg.de/~frank/latex-kurs/latex-kurs-3/beamer-class-bsp-2011-6.html)
\LoadClass[hyperref={pdfpagelabels=false},xcolor=dvipsnames]{beamer}
% Benoetigte Pakete
\RequirePackage[T1]{fontenc}
\RequirePackage[utf8]{inputenc}
\RequirePackage{lmodern}
\RequirePackage[english,ngerman]{babel}
\RequirePackage{fancyhdr}
\RequirePackage{graphicx}
\RequirePackage{amsmath}
\RequirePackage{amssymb}
\RequirePackage{upgreek}
%-----------------------------------------------------------------------------------------------------------------------
% Seitenformatierungsoptionen
%\pagestyle{fancy}
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

\newtheorem*{bem}{Bemerkung}
\newtheorem*{bemen}{Note}