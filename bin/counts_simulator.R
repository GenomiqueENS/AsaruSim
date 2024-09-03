#!/usr/bin/env Rscript

library(dplyr)
library(SPARSim)

'%!in%' <- function(x,y)!('%in%'(x,y))

calc.intensity <- function(mtx, ids, type) {
  cells <- ids %>%
           filter(cell_type == type) %>%
           rownames()
  mtx %>% `[`(, colnames(.) %in% cells) %>% rowMeans()
}

calc.library.size <- function(mtx, ids, type) {
  cells <- ids %>%
           filter(cell_type == type) %>%
           rownames()
  mtx %>% `[`(, colnames(.) %in% cells) %>% colSums()
}

calc_count_variability <- function(mtx, df, cell_type) {
  cells_of_type <- df$cell_type == cell_type
  mtx_filtered <- mtx[, cells_of_type]

  mean_normalized_counts <- rowMeans(mtx_filtered)
  count_variability <- rowMeans((mtx_filtered - mean_normalized_counts)^2)
  
  return(count_variability)
}


run_simulation <- function(matrix_file, cell_types_file, output_file) {

  preset <- c('Bacher', 'Camp', 'Chu', 'Engel', 'Horning', 'Tung', 'Zheng', 'Macosko', 'Saunders')
  preset <- paste0(preset, '_param_preset')

  if (matrix_file %in% preset){
    data(list=matrix_file)
    simulation_params <- get(matrix_file)
  }else{
    cat("Loading data ...\n")
    mtx <- read.csv(matrix_file, sep=",", header = TRUE, row.names = 1)
    cell_types <- read.csv(cell_types_file, sep=",", header = TRUE, row.names = 1)

    cat("Normalization ...\n")
    mtx.norm <- scran_normalization(mtx, positive = TRUE)

    cat("Parameters estimation ...\n")
    params_list <- list()
    for (type in unique(cell_types$cell_type)) {
      params_list[[type]] <- list(
        intensity = calc.intensity(mtx.norm, cell_types, type),
        variability = calc_count_variability(mtx.norm, cell_types, type),
        lib_size = calc.library.size(mtx, cell_types, type)
      )
    }
    
    simulation_params <- list()
    for (type in unique(cell_types$cell_type)) {
      param <- SPARSim_create_simulation_parameter(
        intensity = params_list[[type]]$intensity,
        variability = params_list[[type]]$variability,
        library_size = params_list[[type]]$lib_size,
        condition_name = type
      )
      simulation_params[[type]] <- param
    }
  }

  cat("Simulation ...\n")
  SPARSim_result <- SPARSim_simulation(simulation_params)
  sim.mtx <- SPARSim_result$count_matrix

  if (matrix_file %!in% preset){
      rownames(sim.mtx) <- names(params_list[[1]]$intensity)
      cols <- c()
      for (el in params_list) {
        cols <- c(cols, names(el$lib_size))
      }
      colnames(sim.mtx) <- cols
  }

  cat("Saving ...\n")
  write.csv(sim.mtx, file=output_file)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  cat("Usage: Rscript this_script.R matrix_file.csv cell_types_file.csv output_file.csv\n")
  stop("Incorrect number of arguments supplied.")
}

tryCatch({
  run_simulation(args[1], args[2], args[3])
}, error = function(e) {
  cat(paste0("An error occurred: ", e$message, "\n"))
})

