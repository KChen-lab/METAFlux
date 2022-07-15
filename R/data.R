#' Excel data of human gem model
#'
#' Contains the reaction, gene, gene protein relationship,pathway information.
#'
#'
#' @format A data frame with 13082 rows and 16 variables:
#' \describe{
#'   \item{ID}{Reaction ID}
#'   \item{NAME}{Name of reaction}
#'   \item{EQUATION}{Equation of reaction}
#'   ...
#' }
#' @source \url{https://github.com/SysBioChalmers/Human-GEM/blob/main/model/Human-GEM.xlsx}
#' @examples
#' data(human_gem)
"human_gem"






#' This file contains the hams medium  nutrients information
#'
#' This default file indicates the nutrients available in the medium.
#'  There are total 44 metabolites available. This file can be changed depending on users.
#' We will allow users to define the medium profile based on their knowledge.
#'
#'
#' @format A data frame with 44 rows and 2 variables:
#' \describe{
#'   \item{metabolite}{stands for name of nutrients available in the medium. }
#'   \item{reaction_name}{V2 stands for reaction ID of corresponding nutreints}
#'
#'
#'
#' }
#'
#' @examples
#' data(cell_medium)
"cell_medium"




#' This file contains the hams medium  nutrients information in human blood
#'
#' This default file indicates the nutrients available in the medium.
#'  There are total 44 metabolites available. This file can be changed depending on users.
#' We will allow users to define the medium profile based on their knowledge.
#'
#'
#' @format A data frame with 44 rows and 2 variables:
#' \describe{
#'   \item{metabolite}{ stands for name of nutrients available in the human blood. }
#'   \item{reaction_name}{tands for reaction ID of corresponding nutreints}
#'
#'
#'
#' }
#'
#' @examples
#' data(human_blood)
"human_blood"



#' Bulk test example
#' @format A data frame with 58581 rows and 5 variables
#' @examples
#' data(bulk_test_example)
"bulk_test_example"


#' single cell test seurat object example
#' @format A toy example containing 350 cells(Tumor and T cells)
#' @examples
#' data(sc_test_example)
"sc_test_example"


#'1648 exchange reactions:mathematical representation of uptake/secrete
#' metabolites into the extracellular space

#' @format 1648 rows
#' @examples
#' data(nutrient_lookup_files")
"nutrient_lookup_files"

