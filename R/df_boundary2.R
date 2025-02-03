#' Calcule la distribution d'une espèce et identifie ses bordures sur un raster (version 2)
#'
#' Cette fonction charge un raster principal (en conservant les coordonnées en nombres réels), fusionne ses données avec une liste de points
#' correspondant à une espèce, et identifie les bordures vertes et rouges selon des seuils définis. Les informations issues d'un fichier DBF sont ensuite
#' ajoutées pour enrichir les données. La fonction calcule enfin la différence entre le pourcentage d'eau et le pourcentage de présence de l'espèce,
#' et ajuste la couleur associée à chaque cellule en fonction d'un seuil de différence.
#'
#' @param input_raster_path Caractère. Chemin vers le fichier raster principal.
#' @param species_points_list Liste. Liste contenant pour chaque espèce ses points spatiaux (objets compatibles avec \code{terra}, par exemple des \code{SpatVector}).
#' @param dbf_path Caractère. Chemin vers le fichier .dbf fournissant des informations complémentaires sur l'espèce.
#' @param species_name Caractère. Nom de l'espèce à traiter.
#' @param percentage_threshold Numérique. Seuil de pourcentage pour identifier les bordures (cellules ayant une valeur supérieure à ce seuil et inférieure ou égale à 100). Par défaut, \code{20}.
#' @param difference_threshold Numérique. Seuil de différence entre le pourcentage d'eau et le pourcentage d'espèce pour ajuster la couleur. Par défaut, \code{10}.
#'
#' @return Un \code{data.frame} contenant les informations sur les bordures, incluant :
#' \itemize{
#'   \item les coordonnées (\code{x} et \code{y}),
#'   \item le pourcentage de terre (\code{PERCENTAGEWATER}),
#'   \item le pourcentage de l'espèce (\code{SPECIES_PERCENTAGE}),
#'   \item la différence entre ces deux pourcentages (\code{difference}),
#'   \item la couleur associée à la cellule (\code{color}).
#' }
#'
#' @details La fonction procède en dix étapes :
#' \enumerate{
#'   \item Charger le raster principal et le convertir en \code{data.frame} en conservant les coordonnées réelles.
#'   \item Fusionner les points de l'espèce avec le \code{data.frame} du raster à l'aide de la colonne \code{PageName}.
#'         \code{merged_grid} contient toutes les lignes du raster (même celles sans correspondance), tandis que \code{merged_area} ne conserve
#'         que les correspondances exactes.
#'   \item Remplacer les valeurs manquantes de la colonne \code{PERCENTAGE} par 0.
#'   \item Construire un objet \code{SpatRaster} à partir des colonnes \code{x}, \code{y} et \code{PERCENTAGE} pour faciliter la recherche des cellules adjacentes.
#'   \item Identifier les bordures vertes : les cellules dont la valeur est supérieure au seuil (mais inférieure ou égale à 100) et dont les cellules adjacentes ont une valeur de 0.
#'   \item Identifier les bordures rouges : les cellules dont la valeur est comprise entre 0 et le seuil, puis repérer, parmi leurs voisines, celles qui valent 0.
#'   \item Construire deux \code{data.frame} (pour les bordures rouges et vertes) et les combiner.
#'   \item Joindre le \code{data.frame} combiné au \code{data.frame} du raster afin de récupérer des informations complémentaires (par exemple, la colonne \code{PageName}).
#'   \item Enrichir ces données en les fusionnant avec celles issues du fichier DBF (via la colonne \code{PageName}), puis renommer la colonne \code{PERCENTAGE} en \code{PERCENTAGEWATER}.
#'   \item Fusionner avec \code{merged_grid} pour ajouter la colonne \code{SPECIES_PERCENTAGE}, calculer la différence absolue entre \code{PERCENTAGEWATER} et \code{SPECIES_PERCENTAGE}
#'         et ajuster la couleur en fonction du seuil de différence.
#' }
#'
#' @import terra
#' @importFrom foreign read.dbf
#'
#' @examples
#' \dontrun{
#'   result <- df_boundary2(
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
df_boundary2 <- function(
    input_raster_path,
    species_points_list,
    dbf_path,
    species_name,
    percentage_threshold = 20,
    difference_threshold = 10
) {

  # 1. Charger le raster principal
  raster_data <- rast(input_raster_path)

  # 2. Passer le SpatRaster en data.frame (sans arrondir x,y)
  #    On conserve x,y tels qu’ils sont, en flottants.
  raster_df <- as.data.frame(raster_data, xy = TRUE, na.rm = FALSE)

  # 3. Charger et fusionner les points de l'espèce par PageName
  species_points <- species_points_list[[species_name]]

  # merged_grid = TOUTES les lignes de raster_df (all.x=TRUE),
  # on y ajoute les colonnes de species_points (PERCENTAGE, etc.)
  merged_grid <- merge(
    raster_df,
    as.data.frame(species_points),
    by = "PageName",
    all.x = TRUE
  )
  # merged_area = uniquement les lignes qui matchent (all.x=FALSE)
  merged_area <- merge(
    raster_df,
    as.data.frame(species_points),
    by = "PageName",
    all.x = FALSE
  )

  # Remplacer NA dans PERCENTAGE par 0
  merged_grid$PERCENTAGE[is.na(merged_grid$PERCENTAGE)] <- 0

  # 4. Construire un SpatRaster "merged_raster" pour faire l’adjacence
  #    contenant la colonne "PERCENTAGE"
  merged_raster <- rast(
    merged_grid[, c("x", "y", "PERCENTAGE")],
    type = "xyz",
    crs  = crs(raster_data)
  )

  # 5. Identifier les bordures vertes (PERCENTAGE > threshold & <= 100)
  old_border_cells <- which(
    values(merged_raster) >  percentage_threshold &
      values(merged_raster) <= 100
  )
  adj_cells  <- adjacent(merged_raster, old_border_cells, directions = 8, pairs = FALSE)
  # On repère ceux où la valeur est 0 (vide) => "green" border
  green_cells <- adj_cells[ values(merged_raster)[adj_cells] == 0 ]
  green_coords <- xyFromCell(merged_raster, green_cells)

  # 6. Identifier les bordures rouges
  #    (ici : cells dont PERCENTAGE est >0 et <= threshold)
  potential_red_cells <- which(
    values(merged_raster) > 0 &
      values(merged_raster) <= percentage_threshold
  )
  adj_cells_red   <- adjacent(merged_raster, potential_red_cells, directions = 8, pairs = TRUE)
  # On cherche parmi leurs voisins ceux qui valent 0
  adj_empty_cells <- adj_cells_red[ values(merged_raster)[adj_cells_red[, 2]] == 0, ]
  border_red_cells <- unique(adj_empty_cells[, 1])
  red_coords       <- xyFromCell(merged_raster, border_red_cells)

  # 7. Construire data.frame (sans arrondir) pour rouge & vert
  red_df <- data.frame(
    x     = red_coords[, 1],
    y     = red_coords[, 2],
    color = "brown"
  )
  green_df <- data.frame(
    x     = green_coords[, 1],
    y     = green_coords[, 2],
    color = "brown"
  )
  combined_df <- rbind(red_df, green_df)

  # 8. Joindre à raster_df pour récupérer PageName, etc.
  #    Vu que raster_df$x,y sont flottants et
  #    red_coords/green_coords aussi, on espère un match exact.
  bordures_infos <- merge(
    combined_df,
    raster_df,
    by = c("x", "y"),
    all.x = TRUE
  )

  # 9. Ajouter DBF (par PageName) pour enrichir les infos
  dbf_data <- read.dbf(dbf_path)
  bordures_infos <- merge(
    bordures_infos,
    dbf_data,
    by = "PageName",
    all.x = TRUE
  )

  # Nettoyage de colonnes
  colnames(bordures_infos)[
    colnames(bordures_infos) == "PERCENTAGE"
  ] <- "PERCENTAGEWATER"

  bordures_infos <- unique(bordures_infos)

  # Colorer en "darkblue" si PageName ET coast sont NA (exemple)
  bordures_infos$color[
    is.na(bordures_infos$PageName) & is.na(bordures_infos$coast)
  ] <- "darkblue"

  # 10. Ajouter la PERCENTAGE de l’espèce (SPECIES_PERCENTAGE)
  bordures_infos <- merge(
    bordures_infos,
    merged_grid[, c("x", "y", "PERCENTAGE")],
    by = c("x", "y"),
    all.x = TRUE
  )
  colnames(bordures_infos)[
    colnames(bordures_infos) == "PERCENTAGE"
  ] <- "SPECIES_PERCENTAGE"

  # On calcule la différence absolute (PERCENTAGEWATER vs SPECIES_PERCENTAGE)
  bordures_infos$difference <- ifelse(
    !is.na(bordures_infos$PERCENTAGEWATER) &
      !is.na(bordures_infos$SPECIES_PERCENTAGE),
    abs(bordures_infos$PERCENTAGEWATER - bordures_infos$SPECIES_PERCENTAGE),
    NA
  )

  # On colorie en lightblue si la différence est <= difference_threshold
  bordures_infos$color[
    !is.na(bordures_infos$difference) &
      bordures_infos$difference <= difference_threshold
  ] <- "lightblue"

  # Tout ce qui reste en NA => brown
  bordures_infos$color[is.na(bordures_infos$color)] <- "brown"

  return(bordures_infos)
}
