#' Calcule la distribution d'une espèce et identifie ses bordures sur un raster
#'
#' Cette fonction charge un raster, fusionne ses données avec une liste de points
#' correspondant à une espèce, et identifie les bordures en fonction de seuils de pourcentage.
#' Elle calcule également la différence entre le pourcentage d'eau et le pourcentage de présence
#' de l'espèce, afin de modifier l'affectation des couleurs pour la visualisation.
#'
#' @param input_raster_path Caractère. Chemin vers le fichier raster d'entrée.
#' @param species_points_list Liste. Liste contenant pour chaque espèce ses points spatiaux.
#'   Chaque élément doit être un objet spatialisé compatible avec \code{terra} (par exemple, un \code{SpatVector}).
#' @param dbf_path Caractère. Chemin vers le fichier .dbf fournissant des informations complémentaires sur l'espèce.
#' @param species_name Caractère. Nom de l'espèce à traiter.
#' @param percentage_threshold Numérique. Seuil de pourcentage pour identifier les anciennes bordures rouges.
#'   Par défaut \code{20}.
#' @param difference_threshold Numérique. Seuil de différence entre le pourcentage d'eau et celui de l'espèce.
#'   Par défaut \code{10}.
#'
#' @return Un \code{data.frame} contenant les informations calculées sur les bordures, incluant les coordonnées,
#'   les pourcentages (eau et espèce), la différence calculée, ainsi que la couleur associée à chaque cellule.
#'
#' @details La fonction réalise les opérations suivantes :
#' \enumerate{
#'   \item Charge le raster et le convertit en \code{data.frame} tout en arrondissant les coordonnées.
#'   \item Récupère les points de l'espèce depuis la \code{species_points_list} et fusionne ces données
#'         avec le \code{data.frame} du raster à partir de la colonne \code{PageName}.
#'   \item Crée un objet \code{SpatRaster} à partir des valeurs de pourcentage.
#'   \item Identifie les anciennes bordures rouges et leurs cellules adjacentes (marquées en vert si elles sont vides).
#'   \item Fusionne ces informations avec celles du fichier \code{.dbf} pour enrichir les données.
#'   \item Calcule la différence entre le pourcentage d'eau et le pourcentage d'espèce,
#'         et ajuste la couleur en fonction de ce résultat.
#' }
#'
#' @import terra
#' @importFrom foreign read.dbf
#'
#' @examples
#' \dontrun{
#'   result <- df_boundary(
#'     input_raster_path = "path/to/raster.tif",
#'     species_points_list = species_points_list,
#'     dbf_path = "path/to/species.dbf",
#'     species_name = "Lynx lynx",
#'     percentage_threshold = 20,
#'     difference_threshold = 10
#'   )
#'   head(result)
#' }
#'
#' @export
df_boundary <- function(input_raster_path, species_points_list, dbf_path,
                                    species_name,
                                    percentage_threshold = 20, difference_threshold = 10) {
  # Charger le raster
  raster_data <- rast(input_raster_path)
  raster_df <- as.data.frame(raster_data, xy = TRUE, na.rm = FALSE)

  # Arrondir les coordonnées du raster
  raster_df$x <- round(raster_df$x)
  raster_df$y <- round(raster_df$y)

  # Obtenir les points de l'espèce
  species_points <- species_points_list[[species_name]]

  # Fusionner les données raster et espèce
  merged_grid <- merge(raster_df, as.data.frame(species_points), by = "PageName", all.x = TRUE)
  merged_area <- merge(raster_df, as.data.frame(species_points), by = "PageName", all.x = FALSE)
  merged_grid$PERCENTAGE[is.na(merged_grid$PERCENTAGE)] <- 0

  # Créer un SpatRaster à partir de merged_grid
  merged_raster <- rast(merged_grid[, c("x", "y", "PERCENTAGE")], type = "xyz", crs = crs(raster_data))

  # Identifier les anciennes bordures rouges (PERCENTAGE > threshold et <= 100)
  old_border_cells <- which(values(merged_raster) > percentage_threshold & values(merged_raster) <= 100)
  adj_cells <- adjacent(merged_raster, old_border_cells, directions = 8, pairs = FALSE)
  green_cells <- adj_cells[values(merged_raster)[adj_cells] == 0]
  green_coords <- xyFromCell(merged_raster, green_cells)

  # Identifier les pixels rouges initiaux
  potential_red_cells <- which(values(merged_raster) > 0 & values(merged_raster) <= percentage_threshold)
  adj_cells_red <- adjacent(merged_raster, potential_red_cells, directions = 8, pairs = TRUE)
  adj_empty_cells <- adj_cells_red[values(merged_raster)[adj_cells_red[, 2]] == 0, ]
  border_red_cells <- unique(adj_empty_cells[, 1])
  red_coords <- xyFromCell(merged_raster, border_red_cells)

  # Arrondir les coordonnées des bordures rouges et vertes
  red_coords <- data.frame(x = round(red_coords[, 1]), y = round(red_coords[, 2]))
  green_coords <- data.frame(x = round(green_coords[, 1]), y = round(green_coords[, 2]))

  # Préparer les données pour ggplot
  red_df <- data.frame(x = red_coords$x, y = red_coords$y, color = "brown")
  green_df <- data.frame(x = green_coords$x, y = green_coords$y, color = "brown")
  combined_df <- rbind(red_df, green_df)

  # Ajouter les informations de raster_df à combined_df
  bordures_infos <- merge(combined_df, raster_df, by = c("x", "y"), all.x = TRUE)
  dbf_data <- read.dbf(dbf_path)

  bordures_infos <- merge(bordures_infos, dbf_data, by = "PageName", all.x = TRUE)

  # Nettoyer les colonnes et mettre à jour les couleurs
  colnames(bordures_infos)[colnames(bordures_infos) == "PERCENTAGE"] <- "PERCENTAGEWATER"
  bordures_infos <- unique(bordures_infos)
  #bordures_infos$PERCENTAGEWATER[is.na(bordures_infos$PERCENTAGEWATER)] <- 0

  bordures_infos$color[is.na(bordures_infos$PageName) & is.na(bordures_infos$coast)] <- "darkblue"

  # Ajouter les colonnes PERCENTAGE et calculer les différences
  bordures_infos <- merge(bordures_infos, merged_grid[, c("x", "y", "PERCENTAGE")], by = c("x", "y"), all.x = TRUE)
  colnames(bordures_infos)[colnames(bordures_infos) == "PERCENTAGE"] <- "SPECIES_PERCENTAGE"

  bordures_infos$difference <- ifelse(
    !is.na(bordures_infos$PERCENTAGEWATER) & !is.na(bordures_infos$SPECIES_PERCENTAGE),
    abs(bordures_infos$PERCENTAGEWATER - bordures_infos$SPECIES_PERCENTAGE),
    NA
  )

  bordures_infos$color[!is.na(bordures_infos$difference) & bordures_infos$difference <= difference_threshold] <- "lightblue"

  # Étape 3 : Définir les brown (tout le reste)
  bordures_infos$color[is.na(bordures_infos$color)] <- "brown"

  return(bordures_infos)
}
