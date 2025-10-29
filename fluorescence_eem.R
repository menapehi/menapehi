# ============================================================================
# Módulo para procesamiento de matrices de espectros de fluorescencia
# de excitación-emisión (EEM)
#
# Este módulo proporciona funciones para aplanar y reconstruir matrices EEM,
# que son comúnmente utilizadas en análisis de fluorescencia de muestras
# químicas y ambientales.
#
# Author: @menapehi
# ============================================================================

# Clase S3 para matrices EEM
#' @title Crear objeto FluorescenceEEM
#' @description Crea un objeto de clase FluorescenceEEM para manejar matrices
#'              de espectros de fluorescencia de excitación-emisión
#'
#' @param data Matriz con los datos de fluorescencia (filas=excitación, columnas=emisión)
#' @param excitation_wavelengths Vector con longitudes de onda de excitación (opcional)
#' @param emission_wavelengths Vector con longitudes de onda de emisión (opcional)
#'
#' @return Objeto de clase FluorescenceEEM
#' @export
FluorescenceEEM <- function(data,
                           excitation_wavelengths = NULL,
                           emission_wavelengths = NULL) {

  # Validar que data es una matriz
  if (!is.matrix(data)) {
    stop("Los datos deben ser una matriz")
  }

  dims <- dim(data)
  n_excitation <- dims[1]
  n_emission <- dims[2]

  # Longitudes de onda por defecto si no se proporcionan
  if (is.null(excitation_wavelengths)) {
    excitation_wavelengths <- seq_len(n_excitation)
  }
  if (is.null(emission_wavelengths)) {
    emission_wavelengths <- seq_len(n_emission)
  }

  # Validar dimensiones
  if (length(excitation_wavelengths) != n_excitation) {
    stop("La longitud de excitation_wavelengths debe coincidir con las filas de data")
  }
  if (length(emission_wavelengths) != n_emission) {
    stop("La longitud de emission_wavelengths debe coincidir con las columnas de data")
  }

  # Crear objeto
  eem <- list(
    data = data,
    excitation_wavelengths = excitation_wavelengths,
    emission_wavelengths = emission_wavelengths,
    n_excitation = n_excitation,
    n_emission = n_emission
  )

  class(eem) <- "FluorescenceEEM"
  return(eem)
}


#' @title Aplanar matriz EEM
#' @description Convierte una matriz EEM 2D en un vector 1D
#'
#' @param eem Objeto FluorescenceEEM o matriz
#' @param byrow Lógico. Si TRUE, aplana por filas (default). Si FALSE, por columnas.
#'
#' @return Vector numérico con los datos aplanados
#' @export
flatten_eem <- function(eem, byrow = TRUE) {
  if (inherits(eem, "FluorescenceEEM")) {
    data <- eem$data
  } else if (is.matrix(eem)) {
    data <- eem
  } else {
    stop("El argumento debe ser un objeto FluorescenceEEM o una matriz")
  }

  # Aplanar la matriz
  if (byrow) {
    # Por filas (equivalente a order='C' en Python)
    flat <- as.vector(t(data))
  } else {
    # Por columnas (equivalente a order='F' en Python)
    flat <- as.vector(data)
  }

  return(flat)
}


#' @title Aplanar matriz EEM con etiquetas
#' @description Aplana una matriz EEM y retorna las etiquetas de cada punto
#'
#' @param eem Objeto FluorescenceEEM
#' @param byrow Lógico. Si TRUE, aplana por filas. Si FALSE, por columnas.
#'
#' @return Lista con 'data' (vector aplanado) y 'labels' (data.frame con excitación y emisión)
#' @export
flatten_eem_with_labels <- function(eem, byrow = TRUE) {
  if (!inherits(eem, "FluorescenceEEM")) {
    stop("El argumento debe ser un objeto FluorescenceEEM")
  }

  flat_data <- flatten_eem(eem, byrow = byrow)

  # Crear etiquetas
  if (byrow) {
    # Recorrer por filas
    labels <- expand.grid(
      emission = eem$emission_wavelengths,
      excitation = eem$excitation_wavelengths
    )[, c("excitation", "emission")]
  } else {
    # Recorrer por columnas
    labels <- expand.grid(
      excitation = eem$excitation_wavelengths,
      emission = eem$emission_wavelengths
    )
  }

  return(list(
    data = flat_data,
    labels = labels
  ))
}


#' @title Convertir EEM a data.frame
#' @description Convierte una matriz EEM aplanada a un data.frame
#'
#' @param eem Objeto FluorescenceEEM
#'
#' @return data.frame con columnas: excitation, emission, intensity
#' @export
eem_to_dataframe <- function(eem) {
  if (!inherits(eem, "FluorescenceEEM")) {
    stop("El argumento debe ser un objeto FluorescenceEEM")
  }

  result <- flatten_eem_with_labels(eem, byrow = TRUE)

  df <- data.frame(
    excitation_nm = result$labels$excitation,
    emission_nm = result$labels$emission,
    intensity = result$data
  )

  return(df)
}


#' @title Reconstruir matriz EEM desde datos aplanados
#' @description Reconstruye una matriz EEM 2D desde un vector 1D
#'
#' @param flat_data Vector con datos aplanados
#' @param n_excitation Número de longitudes de onda de excitación
#' @param n_emission Número de longitudes de onda de emisión
#' @param byrow Lógico. Orden usado en el aplanamiento
#' @param excitation_wavelengths Vector con longitudes de onda de excitación (opcional)
#' @param emission_wavelengths Vector con longitudes de onda de emisión (opcional)
#'
#' @return Objeto FluorescenceEEM con la matriz reconstruida
#' @export
reconstruct_eem <- function(flat_data,
                           n_excitation,
                           n_emission,
                           byrow = TRUE,
                           excitation_wavelengths = NULL,
                           emission_wavelengths = NULL) {

  expected_length <- n_excitation * n_emission
  if (length(flat_data) != expected_length) {
    stop(sprintf("El tamaño del vector aplanado (%d) no coincide con n_excitation * n_emission (%d)",
                length(flat_data), expected_length))
  }

  # Reconstruir matriz
  if (byrow) {
    # Por filas
    data <- matrix(flat_data, nrow = n_excitation, ncol = n_emission, byrow = TRUE)
  } else {
    # Por columnas
    data <- matrix(flat_data, nrow = n_excitation, ncol = n_emission, byrow = FALSE)
  }

  # Crear objeto EEM
  eem <- FluorescenceEEM(data, excitation_wavelengths, emission_wavelengths)
  return(eem)
}


#' @title Obtener estadísticas de matriz EEM
#' @description Calcula estadísticas básicas de la matriz EEM
#'
#' @param eem Objeto FluorescenceEEM
#'
#' @return Lista con estadísticas: min, max, mean, sd, median, shape, total_points
#' @export
get_eem_stats <- function(eem) {
  if (!inherits(eem, "FluorescenceEEM")) {
    stop("El argumento debe ser un objeto FluorescenceEEM")
  }

  data <- eem$data

  stats <- list(
    min = min(data, na.rm = TRUE),
    max = max(data, na.rm = TRUE),
    mean = mean(data, na.rm = TRUE),
    sd = sd(data, na.rm = TRUE),
    median = median(data, na.rm = TRUE),
    shape = dim(data),
    total_points = length(data)
  )

  return(stats)
}


#' @title Normalizar matriz EEM
#' @description Normaliza la matriz EEM usando diferentes métodos
#'
#' @param eem Objeto FluorescenceEEM
#' @param method Método de normalización: "max", "area", o "zscore"
#'
#' @return Nuevo objeto FluorescenceEEM normalizado
#' @export
normalize_eem <- function(eem, method = "max") {
  if (!inherits(eem, "FluorescenceEEM")) {
    stop("El argumento debe ser un objeto FluorescenceEEM")
  }

  data_norm <- eem$data

  if (method == "max") {
    max_val <- max(data_norm, na.rm = TRUE)
    if (max_val != 0) {
      data_norm <- data_norm / max_val
    }
  } else if (method == "area") {
    total_area <- sum(data_norm, na.rm = TRUE)
    if (total_area != 0) {
      data_norm <- data_norm / total_area
    }
  } else if (method == "zscore") {
    mean_val <- mean(data_norm, na.rm = TRUE)
    sd_val <- sd(data_norm, na.rm = TRUE)
    if (sd_val != 0) {
      data_norm <- (data_norm - mean_val) / sd_val
    }
  } else {
    stop(sprintf("Método de normalización '%s' no reconocido", method))
  }

  eem_norm <- FluorescenceEEM(data_norm,
                             eem$excitation_wavelengths,
                             eem$emission_wavelengths)
  return(eem_norm)
}


#' @title Eliminar dispersión de Rayleigh
#' @description Elimina la dispersión de Rayleigh (excitación ≈ emisión)
#'
#' @param eem Objeto FluorescenceEEM
#' @param bandwidth Ancho de banda en unidades de longitud de onda
#'
#' @return Nuevo objeto FluorescenceEEM con dispersión eliminada (valores = NA)
#' @export
remove_rayleigh_scatter <- function(eem, bandwidth = 10) {
  if (!inherits(eem, "FluorescenceEEM")) {
    stop("El argumento debe ser un objeto FluorescenceEEM")
  }

  data_cleaned <- eem$data

  for (i in seq_len(eem$n_excitation)) {
    for (j in seq_len(eem$n_emission)) {
      ex <- eem$excitation_wavelengths[i]
      em <- eem$emission_wavelengths[j]

      if (abs(em - ex) <= bandwidth) {
        data_cleaned[i, j] <- NA
      }
    }
  }

  eem_cleaned <- FluorescenceEEM(data_cleaned,
                                eem$excitation_wavelengths,
                                eem$emission_wavelengths)
  return(eem_cleaned)
}


#' @title Aplanar múltiples matrices EEM
#' @description Aplana múltiples matrices EEM en una matriz 2D
#'
#' @param eem_list Lista de objetos FluorescenceEEM
#' @param byrow Lógico. Orden de aplanamiento
#'
#' @return Matriz donde cada fila es una EEM aplanada
#' @export
flatten_eem_batch <- function(eem_list, byrow = TRUE) {
  if (length(eem_list) == 0) {
    stop("La lista de EEM está vacía")
  }

  # Validar que todas tengan la misma forma
  first_shape <- dim(eem_list[[1]]$data)

  for (i in seq_along(eem_list)) {
    if (!inherits(eem_list[[i]], "FluorescenceEEM")) {
      stop(sprintf("El elemento %d no es un objeto FluorescenceEEM", i))
    }

    current_shape <- dim(eem_list[[i]]$data)
    if (!all(current_shape == first_shape)) {
      stop(sprintf("EEM %d tiene forma (%d, %d), esperada (%d, %d)",
                  i, current_shape[1], current_shape[2],
                  first_shape[1], first_shape[2]))
    }
  }

  # Aplanar todas las EEM
  flat_list <- lapply(eem_list, flatten_eem, byrow = byrow)

  # Combinar en matriz
  batch_matrix <- do.call(rbind, flat_list)

  return(batch_matrix)
}


# Método print para objetos FluorescenceEEM
#' @export
print.FluorescenceEEM <- function(x, ...) {
  cat("FluorescenceEEM\n")
  cat(sprintf("  Dimensiones: %d x %d\n", x$n_excitation, x$n_emission))
  cat(sprintf("  Rango excitación: [%.1f, %.1f]\n",
             min(x$excitation_wavelengths), max(x$excitation_wavelengths)))
  cat(sprintf("  Rango emisión: [%.1f, %.1f]\n",
             min(x$emission_wavelengths), max(x$emission_wavelengths)))
  cat(sprintf("  Rango intensidad: [%.2f, %.2f]\n",
             min(x$data, na.rm = TRUE), max(x$data, na.rm = TRUE)))
}


# Método summary para objetos FluorescenceEEM
#' @export
summary.FluorescenceEEM <- function(object, ...) {
  cat("=== Resumen FluorescenceEEM ===\n\n")
  print(object)
  cat("\nEstadísticas:\n")
  stats <- get_eem_stats(object)
  cat(sprintf("  Media: %.2f\n", stats$mean))
  cat(sprintf("  Desviación estándar: %.2f\n", stats$sd))
  cat(sprintf("  Mediana: %.2f\n", stats$median))
  cat(sprintf("  Total de puntos: %d\n", stats$total_points))
}


# ============================================================================
# FUNCIONES AUXILIARES SIMPLES
# ============================================================================

#' @title Aplanar matriz simple
#' @description Función auxiliar simple para aplanar una matriz
#'
#' @param data Matriz 2D
#' @param byrow Lógico. Si TRUE, aplana por filas
#'
#' @return Vector aplanado
#' @export
flatten_matrix <- function(data, byrow = TRUE) {
  if (!is.matrix(data)) {
    stop("El argumento debe ser una matriz")
  }

  if (byrow) {
    return(as.vector(t(data)))
  } else {
    return(as.vector(data))
  }
}


#' @title Reconstruir matriz simple
#' @description Función auxiliar simple para reconstruir una matriz
#'
#' @param flat_data Vector aplanado
#' @param n_rows Número de filas
#' @param n_cols Número de columnas
#' @param byrow Lógico. Orden usado en el aplanamiento
#'
#' @return Matriz reconstruida
#' @export
reconstruct_matrix <- function(flat_data, n_rows, n_cols, byrow = TRUE) {
  if (byrow) {
    return(matrix(flat_data, nrow = n_rows, ncol = n_cols, byrow = TRUE))
  } else {
    return(matrix(flat_data, nrow = n_rows, ncol = n_cols, byrow = FALSE))
  }
}


# ============================================================================
# Mensaje de carga
# ============================================================================

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Módulo de análisis de matrices EEM de fluorescencia cargado")
  packageStartupMessage("Author: @menapehi")
}
