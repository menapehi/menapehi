# ============================================================================
# Ejemplos de uso del módulo fluorescence_eem.R
#
# Este archivo demuestra cómo usar las funciones para aplanar y procesar
# matrices de espectros de fluorescencia de excitación-emisión en R.
#
# Author: @menapehi
# ============================================================================

# Cargar el módulo
source("fluorescence_eem.R")

cat("\n")
cat("================================================================\n")
cat("EJEMPLOS DE USO: fluorescence_eem (R)\n")
cat("Procesamiento de matrices de espectros de fluorescencia\n")
cat("================================================================\n\n")


# ============================================================================
# EJEMPLO 1: Aplanamiento básico de matriz EEM
# ============================================================================

ejemplo_basico <- function() {
  cat("================================================================\n")
  cat("EJEMPLO 1: Aplanamiento básico de matriz EEM\n")
  cat("================================================================\n\n")

  # Crear una matriz EEM de ejemplo (10 excitaciones x 15 emisiones)
  set.seed(42)
  n_ex <- 10
  n_em <- 15
  data <- matrix(runif(n_ex * n_em) * 100, nrow = n_ex, ncol = n_em)

  cat(sprintf("Matriz original: %d x %d\n", n_ex, n_em))
  cat("Primeros valores:\n")
  print(data[1:3, 1:5])
  cat("\n")

  # Aplanar la matriz
  flat_data <- flatten_matrix(data, byrow = TRUE)
  cat(sprintf("Vector aplanado: %d elementos\n", length(flat_data)))
  cat("Primeros 10 valores:\n")
  print(flat_data[1:10])
  cat("\n")

  # Reconstruir la matriz
  reconstructed <- reconstruct_matrix(flat_data, n_ex, n_em, byrow = TRUE)
  cat(sprintf("Matriz reconstruida: %d x %d\n", nrow(reconstructed), ncol(reconstructed)))
  cat("Primeros valores:\n")
  print(reconstructed[1:3, 1:5])
  cat("\n")

  # Verificar que la reconstrucción es correcta
  is_equal <- all.equal(data, reconstructed)
  cat(sprintf("¿Reconstrucción exitosa? %s\n\n", isTRUE(is_equal)))
}


# ============================================================================
# EJEMPLO 2: EEM con longitudes de onda reales
# ============================================================================

ejemplo_con_longitudes_de_onda <- function() {
  cat("================================================================\n")
  cat("EJEMPLO 2: EEM con longitudes de onda reales\n")
  cat("================================================================\n\n")

  # Longitudes de onda típicas en espectroscopía de fluorescencia
  excitation_wl <- seq(250, 445, by = 5)  # 250-445 nm, cada 5 nm
  emission_wl <- seq(300, 595, by = 5)    # 300-595 nm, cada 5 nm

  # Simular datos de fluorescencia con un pico gaussiano
  n_ex <- length(excitation_wl)
  n_em <- length(emission_wl)
  data <- matrix(0, nrow = n_ex, ncol = n_em)

  # Crear un pico gaussiano en ex=350nm, em=450nm
  for (i in seq_along(excitation_wl)) {
    for (j in seq_along(emission_wl)) {
      ex <- excitation_wl[i]
      em <- emission_wl[j]
      # Pico gaussiano
      intensity <- 100 * exp(-((ex - 350)^2 / (2 * 20^2) + (em - 450)^2 / (2 * 30^2)))
      data[i, j] <- intensity
    }
  }

  # Crear objeto EEM
  eem <- FluorescenceEEM(data, excitation_wl, emission_wl)
  print(eem)
  cat("\n")

  # Obtener estadísticas
  stats <- get_eem_stats(eem)
  cat("Estadísticas:\n")
  cat(sprintf("  Min: %.4f\n", stats$min))
  cat(sprintf("  Max: %.4f\n", stats$max))
  cat(sprintf("  Media: %.4f\n", stats$mean))
  cat(sprintf("  Desv. estándar: %.4f\n", stats$sd))
  cat(sprintf("  Mediana: %.4f\n", stats$median))
  cat("\n")

  # Aplanar con etiquetas
  result <- flatten_eem_with_labels(eem)
  cat(sprintf("Datos aplanados: %d puntos\n", length(result$data)))
  cat("Primeras 5 etiquetas (Ex, Em):\n")
  print(head(result$labels, 5))
  cat("\nPrimeros 5 valores de intensidad:\n")
  print(result$data[1:5])
  cat("\n")
}


# ============================================================================
# EJEMPLO 3: Normalización de matrices EEM
# ============================================================================

ejemplo_normalizacion <- function() {
  cat("================================================================\n")
  cat("EJEMPLO 3: Normalización de matrices EEM\n")
  cat("================================================================\n\n")

  # Crear datos de ejemplo
  set.seed(123)
  excitation_wl <- seq(250, 440, by = 10)
  emission_wl <- seq(300, 590, by = 10)
  n_ex <- length(excitation_wl)
  n_em <- length(emission_wl)

  # Datos con valores entre 50 y 150
  data <- matrix(runif(n_ex * n_em) * 100 + 50, nrow = n_ex, ncol = n_em)

  eem <- FluorescenceEEM(data, excitation_wl, emission_wl)
  stats_orig <- get_eem_stats(eem)
  cat(sprintf("Datos originales - min: %.2f, max: %.2f\n",
             stats_orig$min, stats_orig$max))

  # Normalización por máximo
  eem_norm_max <- normalize_eem(eem, method = "max")
  stats_max <- get_eem_stats(eem_norm_max)
  cat(sprintf("Normalizado por max - min: %.2f, max: %.2f\n",
             stats_max$min, stats_max$max))

  # Normalización por área
  eem_norm_area <- normalize_eem(eem, method = "area")
  cat(sprintf("Normalizado por área - suma total: %.6f\n",
             sum(eem_norm_area$data)))

  # Normalización z-score
  eem_norm_z <- normalize_eem(eem, method = "zscore")
  stats_z <- get_eem_stats(eem_norm_z)
  cat(sprintf("Normalizado z-score - media: %.6f, sd: %.6f\n\n",
             stats_z$mean, stats_z$sd))
}


# ============================================================================
# EJEMPLO 4: Eliminación de dispersión de Rayleigh
# ============================================================================

ejemplo_eliminar_rayleigh <- function() {
  cat("================================================================\n")
  cat("EJEMPLO 4: Eliminación de dispersión de Rayleigh\n")
  cat("================================================================\n\n")

  # Crear datos de ejemplo
  set.seed(456)
  excitation_wl <- seq(250, 495, by = 5)
  emission_wl <- seq(250, 595, by = 5)
  n_ex <- length(excitation_wl)
  n_em <- length(emission_wl)

  data <- matrix(runif(n_ex * n_em) * 50, nrow = n_ex, ncol = n_em)

  # Añadir dispersión de Rayleigh (alta intensidad cuando ex ≈ em)
  for (i in seq_along(excitation_wl)) {
    for (j in seq_along(emission_wl)) {
      ex <- excitation_wl[i]
      em <- emission_wl[j]
      if (abs(em - ex) < 15) {
        data[i, j] <- 1000  # Alta intensidad de dispersión
      }
    }
  }

  eem <- FluorescenceEEM(data, excitation_wl, emission_wl)
  cat("Datos originales:\n")
  cat(sprintf("  Max: %.2f\n", max(eem$data, na.rm = TRUE)))
  cat(sprintf("  Puntos con datos válidos: %d\n", sum(!is.na(eem$data))))

  # Eliminar dispersión de Rayleigh
  eem_cleaned <- remove_rayleigh_scatter(eem, bandwidth = 15)
  cat("\nDatos después de eliminar Rayleigh:\n")
  cat(sprintf("  Max: %.2f\n", max(eem_cleaned$data, na.rm = TRUE)))
  cat(sprintf("  Puntos con datos válidos: %d\n", sum(!is.na(eem_cleaned$data))))
  cat(sprintf("  Puntos eliminados (NA): %d\n\n", sum(is.na(eem_cleaned$data))))
}


# ============================================================================
# EJEMPLO 5: Procesamiento por lotes de múltiples EEM
# ============================================================================

ejemplo_batch_processing <- function() {
  cat("================================================================\n")
  cat("EJEMPLO 5: Procesamiento por lotes de múltiples EEM\n")
  cat("================================================================\n\n")

  # Simular 5 muestras diferentes
  set.seed(789)
  n_samples <- 5
  n_ex <- 20
  n_em <- 30

  eem_list <- list()
  for (i in seq_len(n_samples)) {
    data <- matrix(runif(n_ex * n_em) * i * 20, nrow = n_ex, ncol = n_em)
    eem_list[[i]] <- FluorescenceEEM(data)
  }

  cat(sprintf("Número de muestras: %d\n", n_samples))
  cat(sprintf("Dimensiones de cada EEM: %d x %d\n", n_ex, n_em))

  # Aplanar todas las EEM en una matriz
  batch_matrix <- flatten_eem_batch(eem_list, byrow = TRUE)
  cat(sprintf("\nMatriz de lote: %d x %d\n", nrow(batch_matrix), ncol(batch_matrix)))
  cat("  (cada fila es una EEM aplanada)\n")

  # Mostrar estadísticas de cada muestra
  cat("\nEstadísticas por muestra:\n")
  for (i in seq_along(eem_list)) {
    stats <- get_eem_stats(eem_list[[i]])
    cat(sprintf("  Muestra %d: media=%.2f, max=%.2f\n", i, stats$mean, stats$max))
  }
  cat("\n")
}


# ============================================================================
# EJEMPLO 6: Convertir a data.frame
# ============================================================================

ejemplo_dataframe_export <- function() {
  cat("================================================================\n")
  cat("EJEMPLO 6: Exportación a data.frame\n")
  cat("================================================================\n\n")

  # Crear datos de ejemplo pequeños
  set.seed(101)
  excitation_wl <- c(250, 300, 350)
  emission_wl <- c(300, 350, 400, 450)
  data <- matrix(runif(3 * 4) * 100, nrow = 3, ncol = 4)

  eem <- FluorescenceEEM(data, excitation_wl, emission_wl)

  # Convertir a data.frame
  df <- eem_to_dataframe(eem)
  cat(sprintf("Data.frame generado (%d filas):\n", nrow(df)))
  print(head(df, 10))
  cat(sprintf("\nColumnas: %s\n", paste(names(df), collapse = ", ")))

  # Ejemplo de filtrado
  high_intensity <- df[df$intensity > 50, ]
  cat(sprintf("\nPuntos con intensidad > 50: %d puntos\n\n", nrow(high_intensity)))
}


# ============================================================================
# EJEMPLO 7: Comparación de órdenes de aplanamiento
# ============================================================================

ejemplo_comparacion_ordenes <- function() {
  cat("================================================================\n")
  cat("EJEMPLO 7: Comparación de órdenes de aplanamiento\n")
  cat("================================================================\n\n")

  # Matriz pequeña para visualizar la diferencia
  data <- matrix(c(
    1, 2, 3,
    4, 5, 6,
    7, 8, 9,
    10, 11, 12
  ), nrow = 4, ncol = 3, byrow = TRUE)

  cat("Matriz original (4 x 3):\n")
  print(data)
  cat("\n")

  eem <- FluorescenceEEM(data)

  # Orden byrow = TRUE (por filas)
  flat_byrow <- flatten_eem(eem, byrow = TRUE)
  cat("Aplanado byrow = TRUE (por filas):\n")
  print(flat_byrow)
  cat("Recorre: [1,2,3,4,5,6,7,8,9,10,11,12]\n\n")

  # Orden byrow = FALSE (por columnas)
  flat_bycol <- flatten_eem(eem, byrow = FALSE)
  cat("Aplanado byrow = FALSE (por columnas):\n")
  print(flat_bycol)
  cat("Recorre: [1,4,7,10,2,5,8,11,3,6,9,12]\n\n")
}


# ============================================================================
# EJEMPLO 8: Uso de métodos print y summary
# ============================================================================

ejemplo_metodos_s3 <- function() {
  cat("================================================================\n")
  cat("EJEMPLO 8: Uso de métodos print y summary\n")
  cat("================================================================\n\n")

  # Crear objeto EEM
  set.seed(202)
  excitation_wl <- seq(250, 400, by = 10)
  emission_wl <- seq(300, 500, by = 10)
  n_ex <- length(excitation_wl)
  n_em <- length(emission_wl)
  data <- matrix(runif(n_ex * n_em) * 200, nrow = n_ex, ncol = n_em)

  eem <- FluorescenceEEM(data, excitation_wl, emission_wl)

  # Usar método print
  cat("Método print():\n")
  print(eem)
  cat("\n")

  # Usar método summary
  cat("Método summary():\n")
  summary(eem)
  cat("\n")
}


# ============================================================================
# EJECUTAR TODOS LOS EJEMPLOS
# ============================================================================

if (interactive() || !interactive()) {
  ejemplo_basico()
  ejemplo_con_longitudes_de_onda()
  ejemplo_normalizacion()
  ejemplo_eliminar_rayleigh()
  ejemplo_batch_processing()
  ejemplo_dataframe_export()
  ejemplo_comparacion_ordenes()
  ejemplo_metodos_s3()

  cat("================================================================\n")
  cat("EJEMPLOS COMPLETADOS\n")
  cat("================================================================\n")
}
