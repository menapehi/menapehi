# 👋 Hi, I'm @menapehi

- 👀 I'm interested in science, chemistry, environment and data science
- 🌱 I'm currently learning data science
- 💞️ I'm looking to collaborate on science projects
- 📫 How to reach me: melissa.perez.hincapie@gmail.com

## 🔬 Projects

### Fluorescence EEM Analysis

Herramientas para el procesamiento de matrices de espectros de fluorescencia de excitación-emisión (EEM - Excitation-Emission Matrix).

**Disponible en Python y R**

#### Características principales:

- **Aplanamiento de matrices EEM**: Convierte matrices 2D (excitación × emisión) en vectores 1D
- **Reconstrucción**: Recupera matrices 2D desde datos aplanados
- **Normalización**: Múltiples métodos (máximo, área, z-score)
- **Eliminación de dispersión de Rayleigh**: Limpia artefactos de dispersión
- **Procesamiento por lotes**: Maneja múltiples muestras simultáneamente
- **Exportación a data.frame/pandas**: Convierte a DataFrame para análisis avanzado
- **Tests completos**: 27+ tests unitarios en cada lenguaje

---

### 🐍 Versión Python

#### Instalación:

```bash
pip install numpy
pip install pandas  # opcional
```

#### Uso rápido:

```python
from fluorescence_eem import FluorescenceEEM
import numpy as np

# Crear datos de ejemplo
data = np.random.rand(20, 30) * 100  # 20 excitaciones × 30 emisiones

# Crear objeto EEM
eem = FluorescenceEEM(data)

# Aplanar la matriz
flat_data = eem.flatten()
print(f"Datos aplanados: {flat_data.shape}")

# Normalizar
eem_normalized = eem.normalize(method='max')

# Obtener estadísticas
stats = eem.get_stats()
print(stats)
```

#### Archivos Python:

- `fluorescence_eem.py`: Módulo principal
- `examples_eem.py`: Ejemplos completos de uso
- `test_fluorescence_eem.py`: Tests unitarios

#### Ejecutar ejemplos Python:

```bash
python examples_eem.py
```

#### Ejecutar tests Python:

```bash
python test_fluorescence_eem.py
# o con pytest:
pytest test_fluorescence_eem.py
```

---

### 📊 Versión R

#### Instalación:

No se requieren paquetes adicionales. El módulo usa solo funciones base de R.

```r
# Cargar el módulo
source("fluorescence_eem.R")
```

#### Uso rápido:

```r
source("fluorescence_eem.R")

# Crear datos de ejemplo
data <- matrix(runif(600) * 100, nrow = 20, ncol = 30)

# Crear objeto EEM
excitation_wl <- seq(250, 345, by = 5)
emission_wl <- seq(300, 445, by = 5)
eem <- FluorescenceEEM(data, excitation_wl, emission_wl)

# Mostrar información
print(eem)

# Aplanar la matriz
flat_data <- flatten_eem(eem, byrow = TRUE)
cat("Datos aplanados:", length(flat_data), "elementos\n")

# Normalizar
eem_normalized <- normalize_eem(eem, method = "max")

# Obtener estadísticas
stats <- get_eem_stats(eem)
print(stats)

# Convertir a data.frame
df <- eem_to_dataframe(eem)
head(df)
```

#### Archivos R:

- `fluorescence_eem.R`: Módulo principal
- `examples_eem.R`: Ejemplos completos de uso (8 ejemplos)
- `test_fluorescence_eem.R`: Tests unitarios (27 tests)

#### Ejecutar ejemplos R:

```bash
Rscript examples_eem.R
# o en R interactivo:
# source("examples_eem.R")
```

#### Ejecutar tests R:

```bash
Rscript test_fluorescence_eem.R
# o en R interactivo:
# source("test_fluorescence_eem.R")
```

---

### 📚 Funciones principales disponibles:

| Función | Python | R |
|---------|--------|---|
| Crear objeto EEM | `FluorescenceEEM(data)` | `FluorescenceEEM(data)` |
| Aplanar matriz | `eem.flatten()` | `flatten_eem(eem)` |
| Reconstruir matriz | `FluorescenceEEM.from_flattened()` | `reconstruct_eem()` |
| Normalizar | `eem.normalize(method)` | `normalize_eem(eem, method)` |
| Eliminar Rayleigh | `eem.remove_rayleigh_scatter()` | `remove_rayleigh_scatter(eem)` |
| Estadísticas | `eem.get_stats()` | `get_eem_stats(eem)` |
| A DataFrame | `eem.to_dataframe()` | `eem_to_dataframe(eem)` |
| Batch processing | `flatten_eem_batch(list)` | `flatten_eem_batch(list)` |

---

<!---
menapehi/menapehi is a ✨ special ✨ repository because its `README.md` (this file) appears on your GitHub profile.
You can click the Preview link to take a look at your changes.
--->
