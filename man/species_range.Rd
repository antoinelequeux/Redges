% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/species_range.R
\name{species_range}
\alias{species_range}
\title{Trace et extrait les données de distribution d'une espèce}
\usage{
species_range(species, raster_data, species_points_list, return_plot = FALSE)
}
\arguments{
\item{species}{Caractère. Le nom de l'espèce.}

\item{raster_data}{Objet raster (par exemple, un RasterLayer) à convertir en data.frame.}

\item{species_points_list}{Liste. Une liste contenant pour chaque espèce ses points.}

\item{return_plot}{Logique. Si \code{TRUE}, la fonction renvoie une liste contenant le tableau
filtré et l'objet ggplot. Par défaut, \code{FALSE}.}
}
\value{
\itemize{
    \item Si \code{return_plot = FALSE}, un \code{data.frame} contenant le tableau filtré.
    \item Si \code{return_plot = TRUE}, une liste avec :
      \itemize{
        \item \code{data} : le \code{data.frame} filtré,
        \item \code{plot} : l'objet \code{ggplot} associé.
      }
  }
}
\description{
Cette fonction convertit un objet raster en data.frame, fusionne les données avec
les points d'une espèce donnée, filtre les cellules valides et trace la distribution
avec ggplot2.
}
