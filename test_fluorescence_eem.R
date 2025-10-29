# ============================================================================
# Tests unitarios para el módulo fluorescence_eem.R
#
# Ejecutar con: Rscript test_fluorescence_eem.R
# O en R interactivo: source("test_fluorescence_eem.R")
#
# Si tienes testthat instalado, puedes usar:
# library(testthat)
# test_file("test_fluorescence_eem.R")
#
# Author: @menapehi
# ============================================================================

# Cargar el módulo
source("fluorescence_eem.R")

# Sistema simple de testing sin dependencias externas
test_count <- 0
test_passed <- 0
test_failed <- 0

run_test <- function(test_name, test_func) {
  test_count <<- test_count + 1
  cat(sprintf("Test %d: %s ... ", test_count, test_name))

  tryCatch({
    test_func()
    test_passed <<- test_passed + 1
    cat("OK\n")
  }, error = function(e) {
    test_failed <<- test_failed + 1
    cat(sprintf("FAILED\n  Error: %s\n", e$message))
  })
}

assert_true <- function(condition, message = "Assertion failed") {
  if (!isTRUE(condition)) {
    stop(message)
  }
}

assert_equal <- function(a, b, tolerance = 1e-7, message = "Values not equal") {
  if (!isTRUE(all.equal(a, b, tolerance = tolerance))) {
    stop(sprintf("%s: %s != %s", message, toString(a), toString(b)))
  }
}

assert_error <- function(expr, message = "Expected error did not occur") {
  result <- tryCatch({
    expr
    FALSE
  }, error = function(e) {
    TRUE
  })

  if (!result) {
    stop(message)
  }
}


cat("\n")
cat("================================================================\n")
cat("EJECUTANDO TESTS UNITARIOS: fluorescence_eem (R)\n")
cat("================================================================\n\n")


# ============================================================================
# Tests de inicialización
# ============================================================================

run_test("Inicialización básica", function() {
  data <- matrix(runif(150), nrow = 10, ncol = 15)
  eem <- FluorescenceEEM(data)

  assert_equal(dim(eem$data), c(10, 15))
  assert_equal(eem$n_excitation, 10)
  assert_equal(eem$n_emission, 15)
  assert_equal(length(eem$excitation_wavelengths), 10)
  assert_equal(length(eem$emission_wavelengths), 15)
})


run_test("Inicialización con longitudes de onda", function() {
  data <- matrix(runif(150), nrow = 10, ncol = 15)
  ex_wl <- seq(250, 295, by = 5)
  em_wl <- seq(300, 370, by = 5)

  eem <- FluorescenceEEM(data, ex_wl, em_wl)

  assert_equal(eem$excitation_wavelengths, ex_wl)
  assert_equal(eem$emission_wavelengths, em_wl)
})


run_test("Inicialización con dimensiones incorrectas (no matriz)", function() {
  assert_error({
    FluorescenceEEM(c(1, 2, 3))
  })
})


run_test("Inicialización con longitudes de onda incorrectas", function() {
  data <- matrix(runif(150), nrow = 10, ncol = 15)
  ex_wl_wrong <- seq(1, 5)  # Solo 5 elementos, debería ser 10

  assert_error({
    FluorescenceEEM(data, ex_wl_wrong, NULL)
  })
})


# ============================================================================
# Tests de aplanamiento
# ============================================================================

run_test("Aplanar matriz con byrow=TRUE", function() {
  data <- matrix(1:12, nrow = 3, ncol = 4, byrow = TRUE)
  eem <- FluorescenceEEM(data)
  flat <- flatten_eem(eem, byrow = TRUE)

  assert_equal(length(flat), 12)
  assert_equal(flat[1], 1)
  assert_equal(flat[4], 4)
  assert_equal(flat, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
})


run_test("Aplanar matriz con byrow=FALSE", function() {
  data <- matrix(1:12, nrow = 3, ncol = 4, byrow = TRUE)
  eem <- FluorescenceEEM(data)
  flat <- flatten_eem(eem, byrow = FALSE)

  assert_equal(length(flat), 12)
  assert_equal(flat[1], 1)
  assert_equal(flat[2], 5)
  assert_equal(flat, c(1, 5, 9, 2, 6, 10, 3, 7, 11, 4, 8, 12))
})


run_test("Aplanar con etiquetas", function() {
  data <- matrix(runif(20), nrow = 4, ncol = 5)
  ex_wl <- c(250, 260, 270, 280)
  em_wl <- c(300, 310, 320, 330, 340)

  eem <- FluorescenceEEM(data, ex_wl, em_wl)
  result <- flatten_eem_with_labels(eem, byrow = TRUE)

  assert_equal(length(result$data), 20)
  assert_equal(nrow(result$labels), 20)
  assert_equal(ncol(result$labels), 2)
  assert_equal(result$labels$excitation[1], 250)
  assert_equal(result$labels$emission[1], 300)
})


# ============================================================================
# Tests de reconstrucción
# ============================================================================

run_test("Reconstruir matriz desde datos aplanados", function() {
  data <- matrix(runif(150), nrow = 10, ncol = 15)
  eem <- FluorescenceEEM(data)
  flat <- flatten_eem(eem, byrow = TRUE)

  eem_reconstructed <- reconstruct_eem(flat, 10, 15, byrow = TRUE)

  assert_equal(eem$data, eem_reconstructed$data, tolerance = 1e-10)
})


run_test("Reconstruir con tamaño incorrecto", function() {
  flat <- c(1, 2, 3, 4, 5)

  assert_error({
    reconstruct_eem(flat, 3, 3, byrow = TRUE)
  })
})


run_test("Roundtrip: aplanar y reconstruir", function() {
  data <- matrix(runif(200), nrow = 20, ncol = 10)
  eem <- FluorescenceEEM(data)

  # Roundtrip con byrow=TRUE
  flat <- flatten_eem(eem, byrow = TRUE)
  reconstructed <- reconstruct_eem(flat, 20, 10, byrow = TRUE)
  assert_equal(data, reconstructed$data, tolerance = 1e-10)

  # Roundtrip con byrow=FALSE
  flat2 <- flatten_eem(eem, byrow = FALSE)
  reconstructed2 <- reconstruct_eem(flat2, 20, 10, byrow = FALSE)
  assert_equal(data, reconstructed2$data, tolerance = 1e-10)
})


# ============================================================================
# Tests de estadísticas
# ============================================================================

run_test("Calcular estadísticas básicas", function() {
  data <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3, ncol = 3)
  eem <- FluorescenceEEM(data)
  stats <- get_eem_stats(eem)

  assert_equal(stats$min, 1)
  assert_equal(stats$max, 9)
  assert_equal(stats$mean, 5)
  assert_equal(stats$median, 5)
  assert_equal(stats$total_points, 9)
})


# ============================================================================
# Tests de normalización
# ============================================================================

run_test("Normalización por máximo", function() {
  data <- matrix(runif(100) * 50 + 50, nrow = 10, ncol = 10)
  eem <- FluorescenceEEM(data)
  eem_norm <- normalize_eem(eem, method = "max")

  assert_true(max(eem_norm$data) <= 1.0 + 1e-10)
  assert_equal(max(eem_norm$data), 1.0, tolerance = 1e-10)
})


run_test("Normalización por área", function() {
  data <- matrix(runif(100) * 100, nrow = 10, ncol = 10)
  eem <- FluorescenceEEM(data)
  eem_norm <- normalize_eem(eem, method = "area")

  assert_equal(sum(eem_norm$data), 1.0, tolerance = 1e-10)
})


run_test("Normalización z-score", function() {
  data <- matrix(runif(100) * 100, nrow = 10, ncol = 10)
  eem <- FluorescenceEEM(data)
  eem_norm <- normalize_eem(eem, method = "zscore")

  assert_equal(mean(eem_norm$data), 0.0, tolerance = 1e-10)
  assert_equal(sd(eem_norm$data), 1.0, tolerance = 1e-10)
})


run_test("Normalización con método inválido", function() {
  data <- matrix(runif(100), nrow = 10, ncol = 10)
  eem <- FluorescenceEEM(data)

  assert_error({
    normalize_eem(eem, method = "invalid_method")
  })
})


# ============================================================================
# Tests de eliminación de dispersión de Rayleigh
# ============================================================================

run_test("Eliminar dispersión de Rayleigh", function() {
  n <- 20
  excitation_wl <- seq(300, 300 + (n - 1) * 10, by = 10)
  emission_wl <- seq(300, 300 + (n - 1) * 10, by = 10)
  data <- matrix(50, nrow = n, ncol = n)

  eem <- FluorescenceEEM(data, excitation_wl, emission_wl)
  eem_cleaned <- remove_rayleigh_scatter(eem, bandwidth = 5)

  # La diagonal debe ser NA
  for (i in seq_len(n)) {
    assert_true(is.na(eem_cleaned$data[i, i]))
  }

  # Valores alejados de la diagonal deben mantenerse
  assert_true(!is.na(eem_cleaned$data[1, n]))
})


# ============================================================================
# Tests de procesamiento por lotes
# ============================================================================

run_test("Aplanar múltiples EEM (batch)", function() {
  n_samples <- 5
  eem_list <- list()

  for (i in seq_len(n_samples)) {
    data <- matrix(runif(150), nrow = 10, ncol = 15)
    eem_list[[i]] <- FluorescenceEEM(data)
  }

  batch_matrix <- flatten_eem_batch(eem_list, byrow = TRUE)

  assert_equal(nrow(batch_matrix), 5)
  assert_equal(ncol(batch_matrix), 150)
})


run_test("Batch con lista vacía", function() {
  assert_error({
    flatten_eem_batch(list(), byrow = TRUE)
  })
})


run_test("Batch con formas diferentes", function() {
  eem1 <- FluorescenceEEM(matrix(runif(150), nrow = 10, ncol = 15))
  eem2 <- FluorescenceEEM(matrix(runif(200), nrow = 10, ncol = 20))

  assert_error({
    flatten_eem_batch(list(eem1, eem2), byrow = TRUE)
  })
})


# ============================================================================
# Tests de conversión a data.frame
# ============================================================================

run_test("Convertir a data.frame", function() {
  data <- matrix(runif(12), nrow = 3, ncol = 4)
  ex_wl <- c(250, 300, 350)
  em_wl <- c(300, 350, 400, 450)

  eem <- FluorescenceEEM(data, ex_wl, em_wl)
  df <- eem_to_dataframe(eem)

  assert_equal(nrow(df), 12)
  assert_equal(ncol(df), 3)
  assert_true("excitation_nm" %in% names(df))
  assert_true("emission_nm" %in% names(df))
  assert_true("intensity" %in% names(df))
})


# ============================================================================
# Tests de funciones auxiliares
# ============================================================================

run_test("Función auxiliar flatten_matrix", function() {
  data <- matrix(1:12, nrow = 3, ncol = 4, byrow = TRUE)
  flat <- flatten_matrix(data, byrow = TRUE)

  assert_equal(length(flat), 12)
  assert_equal(flat, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
})


run_test("Función auxiliar reconstruct_matrix", function() {
  flat <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
  reconstructed <- reconstruct_matrix(flat, 3, 4, byrow = TRUE)

  expected <- matrix(1:12, nrow = 3, ncol = 4, byrow = TRUE)
  assert_equal(reconstructed, expected)
})


# ============================================================================
# Tests de casos límite
# ============================================================================

run_test("Matriz de un solo elemento", function() {
  data <- matrix(42, nrow = 1, ncol = 1)
  eem <- FluorescenceEEM(data)
  flat <- flatten_eem(eem)

  assert_equal(length(flat), 1)
  assert_equal(flat[1], 42)
})


run_test("Matriz con ceros", function() {
  data <- matrix(0, nrow = 5, ncol = 7)
  eem <- FluorescenceEEM(data)
  stats <- get_eem_stats(eem)

  assert_equal(stats$min, 0)
  assert_equal(stats$max, 0)
  assert_equal(stats$mean, 0)
})


run_test("Matriz con valores NA", function() {
  data <- matrix(runif(100), nrow = 10, ncol = 10)
  data[5, 5] <- NA
  data[3, 8] <- NA

  eem <- FluorescenceEEM(data)
  flat <- flatten_eem(eem)

  # Debe preservar NA
  assert_true(sum(is.na(flat)) >= 2)
})


# ============================================================================
# Resumen de tests
# ============================================================================

cat("\n")
cat("================================================================\n")
cat("RESUMEN DE TESTS\n")
cat("================================================================\n")
cat(sprintf("Total de tests: %d\n", test_count))
cat(sprintf("Tests exitosos: %d\n", test_passed))
cat(sprintf("Tests fallidos: %d\n", test_failed))

if (test_failed == 0) {
  cat("\n¡Todos los tests pasaron exitosamente!\n")
} else {
  cat(sprintf("\nALERTA: %d tests fallaron\n", test_failed))
}

cat("================================================================\n\n")

# Retornar código de salida apropiado
if (test_failed > 0) {
  quit(status = 1)
}
