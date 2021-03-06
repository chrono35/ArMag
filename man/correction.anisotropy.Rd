% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ArMag.R
\name{correction.anisotropy}
\alias{correction.anisotropy}
\title{Corrige les mesures avec la matrice d'anisotropie donnée
Il faut une matrice 3 x3}
\usage{
correction.anisotropy(mesures.frame, ani.matrix)
}
\arguments{
\item{mesures.frame}{data.frame avec les mesures et variables X, Y, Z, i, D, F}

\item{ani.matric}{matrice 3 x 3 correspondant à la matrice (symétrique et normalisé) à utiliser pour la correction}
}
\value{
un data.frame du même type que mesures.frame
}
\description{
Corrige les mesures avec la matrice d'anisotropie donnée
Il faut une matrice 3 x3
}
