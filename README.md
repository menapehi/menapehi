# 👋 Hi, I'm @menapehi

- 👀 I'm interested in science, chemistry, environment and data science
- 🌱 I'm currently learning data science
- 💞️ I'm looking to collaborate on science projects
- 📫 How to reach me: melissa.perez.hincapie@gmail.com

## 🔬 Projects

### Fluorescence EEM Analysis

Herramientas para el procesamiento de matrices de espectros de fluorescencia de excitación-emisión (EEM - Excitation-Emission Matrix).

#### Características principales:

- **Aplanamiento de matrices EEM**: Convierte matrices 2D (excitación × emisión) en vectores 1D
- **Reconstrucción**: Recupera matrices 2D desde datos aplanados
- **Normalización**: Múltiples métodos (máximo, área, z-score)
- **Eliminación de dispersión de Rayleigh**: Limpia artefactos de dispersión
- **Procesamiento por lotes**: Maneja múltiples muestras simultáneamente
- **Exportación a pandas**: Convierte a DataFrame para análisis avanzado

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

#### Archivos:

- `fluorescence_eem.py`: Módulo principal con todas las funciones
- `examples_eem.py`: Ejemplos completos de uso
- `test_fluorescence_eem.py`: Tests unitarios

#### Ejecutar ejemplos:

```bash
python examples_eem.py
```

#### Ejecutar tests:

```bash
python test_fluorescence_eem.py
# o con pytest:
pytest test_fluorescence_eem.py
```

---

<!---
menapehi/menapehi is a ✨ special ✨ repository because its `README.md` (this file) appears on your GitHub profile.
You can click the Preview link to take a look at your changes.
--->
