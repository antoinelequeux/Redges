% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/df_boundary.R
\name{df_boundary}
\alias{df_boundary}
\title{Calcule la distribution d'une espèce et identifie ses bordures sur un raster}
\usage{
df_boundary(
  input_raster_path,
  species_points_list,
  dbf_path,
  species_name,
  percentage_threshold = 20,
  difference_threshold = 10
)
}
\arguments{
\item{input_raster_path}{Caractère. Chemin vers le fichier raster d'entrée.}

\item{species_points_list}{Liste. Liste contenant pour chaque espèce ses points spatiaux.
Chaque élément doit être un objet spatialisé compatible avec \code{terra} (par exemple, un \code{SpatVector}).}

\item{dbf_path}{Caractère. Chemin vers le fichier .dbf fournissant des informations complémentaires sur l'espèce.}

\item{species_name}{Caractère. Nom de l'espèce à traiter.}

\item{percentage_threshold}{Numérique. Seuil de pourcentage pour identifier les anciennes bordures rouges.
Par défaut \code{20}.}

\item{difference_threshold}{Numérique. Seuil de différence entre le pourcentage d'eau et celui de l'espèce.
Par défaut \code{10}.}
}
\value{
Un \code{data.frame} contenant les informations calculées sur les bordures, incluant les coordonnées,
  les pourcentages (eau et espèce), la différence calculée, ainsi que la couleur associée à chaque cellule.
}
\description{
Cette fonction charge un raster, fusionne ses données avec une liste de points
correspondant à une espèce, et identifie les bordures en fonction de seuils de pourcentage.
Elle calcule également la différence entre le pourcentage d'eau et le pourcentage de présence
de l'espèce, afin de modifier l'affectation des couleurs pour la visualisation.
}
\details{
La fonction réalise les opérations suivantes :
\enumerate{
  \item Charge le raster et le convertit en \code{data.frame} tout en arrondissant les coordonnées.
  \item Récupère les points de l'espèce depuis la \code{species_points_list} et fusionne ces données
        avec le \code{data.frame} du raster à partir de la colonne \code{PageName}.
  \item Crée un objet \code{SpatRaster} à partir des valeurs de pourcentage.
  \item Identifie les anciennes bordures rouges et leurs cellules adjacentes (marquées en vert si elles sont vides).
  \item Fusionne ces informations avec celles du fichier \code{.dbf} pour enrichir les données.
  \item Calcule la différence entre le pourcentage d'eau et le pourcentage d'espèce,
        et ajuste la couleur en fonction de ce résultat.
}
}
\examples{
\dontrun{
  result <- df_boundary(
    input_raster_path = "path/to/raster.tif",
    species_points_list = species_points_list,
    dbf_path = "path/to/species.dbf",
    species_name = "Lynx lynx",
    percentage_threshold = 20,
    difference_threshold = 10
  )
  head(result)
}

}
