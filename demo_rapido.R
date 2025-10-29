#!/usr/bin/env Rscript
# ============================================================================
# Demo rápido del módulo de matrices EEM de fluorescencia
# Ejecutar con: Rscript demo_rapido.R
# ============================================================================

source("fluorescence_eem.R")

cat("\n")
cat("================================================================\n")
cat("DEMO RÁPIDO: Análisis de matrices EEM de fluorescencia\n")
cat("================================================================\n\n")

# 1. Crear datos de ejemplo con un pico de fluorescencia
cat("1. Creando matriz EEM simulada con pico de fluorescencia...\n")
excitation_wl <- seq(250, 400, by = 5)  # 31 longitudes de onda
emission_wl <- seq(300, 500, by = 5)    # 41 longitudes de onda

n_ex <- length(excitation_wl)
n_em <- length(emission_wl)
data <- matrix(0, nrow = n_ex, ncol = n_em)

# Simular pico gaussiano en ex=320nm, em=410nm
for (i in seq_along(excitation_wl)) {
  for (j in seq_along(emission_wl)) {
    ex <- excitation_wl[i]
    em <- emission_wl[j]
    intensity <- 100 * exp(-((ex - 320)^2 / (2 * 15^2) + (em - 410)^2 / (2 * 25^2)))
    data[i, j] <- intensity + runif(1) * 5  # Añadir ruido
  }
}

# 2. Crear objeto EEM
eem <- FluorescenceEEM(data, excitation_wl, emission_wl)
cat("\n")
print(eem)
cat("\n")

# 3. Estadísticas
cat("2. Estadísticas de la matriz EEM:\n")
stats <- get_eem_stats(eem)
cat(sprintf("   - Intensidad mínima: %.2f\n", stats$min))
cat(sprintf("   - Intensidad máxima: %.2f\n", stats$max))
cat(sprintf("   - Intensidad media: %.2f\n", stats$mean))
cat(sprintf("   - Desviación estándar: %.2f\n", stats$sd))
cat(sprintf("   - Total de puntos: %d\n\n", stats$total_points))

# 4. Aplanar la matriz
cat("3. Aplanando la matriz EEM...\n")
flat_data <- flatten_eem(eem, byrow = TRUE)
cat(sprintf("   - Datos originales: %d x %d = %d puntos\n",
            n_ex, n_em, n_ex * n_em))
cat(sprintf("   - Datos aplanados: vector de %d elementos\n\n",
            length(flat_data)))

# 5. Normalización
cat("4. Normalizando por máximo...\n")
eem_norm <- normalize_eem(eem, method = "max")
stats_norm <- get_eem_stats(eem_norm)
cat(sprintf("   - Nuevo rango: [%.2f, %.2f]\n\n",
            stats_norm$min, stats_norm$max))

# 6. Eliminar dispersión de Rayleigh
cat("5. Eliminando dispersión de Rayleigh...\n")
eem_clean <- remove_rayleigh_scatter(eem, bandwidth = 10)
puntos_eliminados <- sum(is.na(eem_clean$data)) - sum(is.na(eem$data))
cat(sprintf("   - Puntos eliminados: %d\n", puntos_eliminados))
cat(sprintf("   - Puntos válidos restantes: %d\n\n",
            sum(!is.na(eem_clean$data))))

# 7. Convertir a data.frame
cat("6. Convirtiendo a data.frame...\n")
df <- eem_to_dataframe(eem)
cat(sprintf("   - Data.frame creado con %d filas y %d columnas\n",
            nrow(df), ncol(df)))
cat("   - Primeras 5 filas:\n")
print(head(df, 5))
cat("\n")

# 8. Reconstruir desde datos aplanados
cat("7. Reconstruyendo matriz desde datos aplanados...\n")
eem_reconstructed <- reconstruct_eem(flat_data, n_ex, n_em, byrow = TRUE,
                                    excitation_wl, emission_wl)
diferencia_max <- max(abs(eem$data - eem_reconstructed$data))
cat(sprintf("   - Diferencia máxima: %.10f\n", diferencia_max))
if (diferencia_max < 1e-10) {
  cat("   - ✓ Reconstrucción perfecta!\n\n")
} else {
  cat("   - ⚠ Hay diferencias en la reconstrucción\n\n")
}

# 9. Procesar múltiples muestras
cat("8. Procesamiento por lotes de 3 muestras...\n")
eem_list <- list()
for (i in 1:3) {
  data_sample <- matrix(runif(n_ex * n_em) * i * 30, nrow = n_ex, ncol = n_em)
  eem_list[[i]] <- FluorescenceEEM(data_sample, excitation_wl, emission_wl)
}

batch_matrix <- flatten_eem_batch(eem_list, byrow = TRUE)
cat(sprintf("   - Matriz de lote: %d muestras x %d características\n",
            nrow(batch_matrix), ncol(batch_matrix)))
cat("   - Esta matriz puede usarse para machine learning!\n\n")

cat("================================================================\n")
cat("DEMO COMPLETADO\n")
cat("================================================================\n\n")

cat("Para más ejemplos, ejecuta: Rscript examples_eem.R\n")
cat("Para tests unitarios, ejecuta: Rscript test_fluorescence_eem.R\n\n")
