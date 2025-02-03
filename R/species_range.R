#' Trace et extrait les données de distribution d'une espèce
#'
#' Cette fonction convertit un objet raster en data.frame, fusionne les données avec
#' les points d'une espèce donnée, filtre les cellules valides et trace la distribution
#' avec ggplot2.
#'
#' @param species Caractère. Le nom de l'espèce.
#' @param raster_data Objet raster (par exemple, un RasterLayer) à convertir en data.frame.
#' @param species_points_list Liste. Une liste contenant pour chaque espèce ses points.
#' @param return_plot Logique. Si \code{TRUE}, la fonction renvoie une liste contenant le tableau
#'   filtré et l'objet ggplot. Par défaut, \code{FALSE}.
#'
#' @return
#'   \itemize{
#'     \item Si \code{return_plot = FALSE}, un \code{data.frame} contenant le tableau filtré.
#'     \item Si \code{return_plot = TRUE}, une liste avec :
#'       \itemize{
#'         \item \code{data} : le \code{data.frame} filtré,
#'         \item \code{plot} : l'objet \code{ggplot} associé.
#'       }
#'   }
#'
#' @import ggplot2
#'
#'
#'
#' @export
species_range <- function(species, raster_data, species_points_list, return_plot = FALSE) {
  # Convertir le raster en data.frame pour ggplot
  raster_df <- as.data.frame(raster_data, xy = TRUE, na.rm = FALSE)

  # Supprimer les NA uniquement pour les cellules à afficher (basé sur la colonne PageName)
  raster_df <- raster_df[!is.na(raster_df$PageName), ]

  # Obtenir les données des points pour l'espèce
  species_points <- species_points_list[[species]]

  # Fusionner les points de l'espèce avec les données raster via PageName
  merged_grid <- merge(raster_df, as.data.frame(species_points), by = "PageName", all.x = TRUE)

  # Filtrer les cellules valides (où PERCENTAGE est non-NA)
  filtered_grid <- merged_grid[!is.na(merged_grid$PERCENTAGE), ]

  # Définir les limites de la grille pour le graphique
  xmin <- min(filtered_grid$x)
  xmax <- max(filtered_grid$x)
  ymin <- min(filtered_grid$y)
  ymax <- max(filtered_grid$y)

  # Créer le graphique avec ggplot2
  p <- ggplot() +
    geom_tile(
      data = raster_df,
      aes(x = x, y = y),
      fill = "white",
      color = "black"  # Grilles visibles en fond
    ) +
    geom_tile(
      data = filtered_grid,
      aes(x = x, y = y, fill = PERCENTAGE),
      color = NA
    ) +
    scale_fill_gradient(
      low = "white", high = "red", name = "% de remplissage"
    ) +
    ggtitle(paste("Aire de répartition de l'espèce :", species)) +
    theme_minimal() +
    coord_cartesian(
      xlim = c(xmin - 1000000, xmax + 1000000),
      ylim = c(ymin - 1000000, ymax + 1000000)
    )

  # Retourner soit le tableau seul, soit le tableau et le graphique
  if (return_plot) {
    return(list(data = filtered_grid, plot = p))
  } else {
    return(filtered_grid)
  }
}
