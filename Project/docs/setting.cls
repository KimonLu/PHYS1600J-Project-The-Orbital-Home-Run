% PHYS1600J Project Template based on RevTeX 4.2.
\documentclass[reprint,aps,pra,floatfix]{revtex4-2}

\usepackage{graphicx}
\usepackage{amsfonts,amsmath,amssymb}
\usepackage{bm}
\usepackage{booktabs}
\usepackage{array}
\usepackage{tabularx}
\usepackage{siunitx}
\usepackage{xcolor}
\usepackage{microtype}
\usepackage{placeins}
\usepackage{xurl}
\usepackage[colorlinks,citecolor=black,urlcolor=black,linkcolor=black]{hyperref}

\sisetup{
  separate-uncertainty=true,
  per-mode=symbol,
  range-phrase=--,
  range-units=single,
  exponent-product=\times
}
\graphicspath{{figures/}{../figures/}}
\renewcommand{\arraystretch}{1.08}
\newcommand{\dd}{\mathrm{d}}
\newcommand{\e}{\mathrm{e}}
\newcommand{\vect}[1]{\bm{#1}}

\hypersetup{
  pdftitle={Problem B: The Orbital Home Run - From a Measure-Zero Ideal Orbit to a Terrain-Aware Lunar Return},
  pdfauthor={Kemeng Lu; Jiahao Yin; Yuxin Wu},
  pdfkeywords={lunar orbit, two-body dynamics, sensitivity analysis, lunar gravity, terrain clearance}
}
