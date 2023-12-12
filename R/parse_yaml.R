#' @title Parse YAML
#'
#' @description This function parses YAML files
#'
#' @param file The path to the YAML file
#'
#' @return A list containing the parameters specified in the YAML file
#'
#' @export
#'
#' @examples NULL
parse_yaml <- function(file) {
  return(yaml::read_yaml(file = file))
}
