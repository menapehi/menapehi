"""
Ejemplos de uso del módulo fluorescence_eem

Este archivo demuestra cómo usar las funciones para aplanar y procesar
matrices de espectros de fluorescencia de excitación-emisión.
"""

import numpy as np
from fluorescence_eem import (
    FluorescenceEEM,
    flatten_eem_matrix,
    reconstruct_eem_matrix,
    flatten_eem_batch
)


def ejemplo_basico():
    """Ejemplo básico de aplanamiento de matriz EEM"""
    print("=" * 60)
    print("EJEMPLO 1: Aplanamiento básico de matriz EEM")
    print("=" * 60)

    # Crear una matriz EEM de ejemplo (10 excitaciones x 15 emisiones)
    n_ex = 10
    n_em = 15
    data = np.random.rand(n_ex, n_em) * 100

    print(f"\nMatriz original: {data.shape}")
    print(f"Primeros valores:\n{data[:3, :5]}\n")

    # Aplanar la matriz
    flat_data = flatten_eem_matrix(data, order='C')
    print(f"Vector aplanado: {flat_data.shape}")
    print(f"Primeros 10 valores: {flat_data[:10]}\n")

    # Reconstruir la matriz
    reconstructed = reconstruct_eem_matrix(flat_data, n_ex, n_em, order='C')
    print(f"Matriz reconstruida: {reconstructed.shape}")
    print(f"Primeros valores:\n{reconstructed[:3, :5]}\n")

    # Verificar que la reconstrucción es correcta
    is_equal = np.allclose(data, reconstructed)
    print(f"¿Reconstrucción exitosa? {is_equal}\n")


def ejemplo_con_longitudes_de_onda():
    """Ejemplo usando longitudes de onda reales"""
    print("=" * 60)
    print("EJEMPLO 2: EEM con longitudes de onda reales")
    print("=" * 60)

    # Longitudes de onda típicas en espectroscopía de fluorescencia
    excitation_wl = np.arange(250, 450, 5)  # 250-445 nm, cada 5 nm
    emission_wl = np.arange(300, 600, 5)    # 300-595 nm, cada 5 nm

    # Simular datos de fluorescencia con un pico
    n_ex = len(excitation_wl)
    n_em = len(emission_wl)
    data = np.zeros((n_ex, n_em))

    # Crear un pico gaussiano en ex=350nm, em=450nm
    for i, ex in enumerate(excitation_wl):
        for j, em in enumerate(emission_wl):
            # Pico gaussiano
            intensity = 100 * np.exp(-((ex - 350)**2 / (2 * 20**2) + (em - 450)**2 / (2 * 30**2)))
            data[i, j] = intensity

    # Crear objeto EEM
    eem = FluorescenceEEM(data, excitation_wl, emission_wl)
    print(f"\n{eem}")

    # Obtener estadísticas
    stats = eem.get_stats()
    print(f"\nEstadísticas:")
    for key, value in stats.items():
        print(f"  {key}: {value}")

    # Aplanar con etiquetas
    flat_data, labels = eem.flatten_with_labels()
    print(f"\nDatos aplanados: {len(flat_data)} puntos")
    print(f"Primeras 5 etiquetas (Ex, Em): {labels[:5]}")
    print(f"Primeros 5 valores de intensidad: {flat_data[:5]}\n")


def ejemplo_normalizacion():
    """Ejemplo de normalización de EEM"""
    print("=" * 60)
    print("EJEMPLO 3: Normalización de matrices EEM")
    print("=" * 60)

    # Crear datos de ejemplo
    excitation_wl = np.arange(250, 450, 10)
    emission_wl = np.arange(300, 600, 10)
    n_ex, n_em = len(excitation_wl), len(emission_wl)

    # Datos con valores entre 50 y 150
    data = np.random.rand(n_ex, n_em) * 100 + 50

    eem = FluorescenceEEM(data, excitation_wl, emission_wl)
    print(f"\nDatos originales - min: {eem.get_stats()['min']:.2f}, max: {eem.get_stats()['max']:.2f}")

    # Normalización por máximo
    eem_norm_max = eem.normalize(method='max')
    print(f"Normalizado por max - min: {eem_norm_max.get_stats()['min']:.2f}, max: {eem_norm_max.get_stats()['max']:.2f}")

    # Normalización por área
    eem_norm_area = eem.normalize(method='area')
    print(f"Normalizado por área - suma total: {np.sum(eem_norm_area.data):.4f}")

    # Normalización z-score
    eem_norm_z = eem.normalize(method='zscore')
    print(f"Normalizado z-score - media: {eem_norm_z.get_stats()['mean']:.4f}, std: {eem_norm_z.get_stats()['std']:.4f}\n")


def ejemplo_eliminar_rayleigh():
    """Ejemplo de eliminación de dispersión de Rayleigh"""
    print("=" * 60)
    print("EJEMPLO 4: Eliminación de dispersión de Rayleigh")
    print("=" * 60)

    # Crear datos de ejemplo
    excitation_wl = np.arange(250, 500, 5)
    emission_wl = np.arange(250, 600, 5)
    n_ex, n_em = len(excitation_wl), len(emission_wl)

    data = np.random.rand(n_ex, n_em) * 50

    # Añadir dispersión de Rayleigh (alta intensidad cuando ex ≈ em)
    for i, ex in enumerate(excitation_wl):
        for j, em in enumerate(emission_wl):
            if abs(em - ex) < 15:
                data[i, j] = 1000  # Alta intensidad de dispersión

    eem = FluorescenceEEM(data, excitation_wl, emission_wl)
    print(f"\nDatos originales:")
    print(f"  Max: {eem.get_stats()['max']:.2f}")
    print(f"  Puntos con datos válidos: {np.sum(~np.isnan(eem.data))}")

    # Eliminar dispersión de Rayleigh
    eem_cleaned = eem.remove_rayleigh_scatter(bandwidth=15)
    print(f"\nDatos después de eliminar Rayleigh:")
    print(f"  Max: {np.nanmax(eem_cleaned.data):.2f}")
    print(f"  Puntos con datos válidos: {np.sum(~np.isnan(eem_cleaned.data))}")
    print(f"  Puntos eliminados (NaN): {np.sum(np.isnan(eem_cleaned.data))}\n")


def ejemplo_batch_processing():
    """Ejemplo de procesamiento por lotes de múltiples EEM"""
    print("=" * 60)
    print("EJEMPLO 5: Procesamiento por lotes de múltiples EEM")
    print("=" * 60)

    # Simular 5 muestras diferentes
    n_samples = 5
    n_ex, n_em = 20, 30

    eem_list = []
    for i in range(n_samples):
        data = np.random.rand(n_ex, n_em) * (i + 1) * 20  # Diferentes intensidades
        eem = FluorescenceEEM(data)
        eem_list.append(eem)

    print(f"\nNúmero de muestras: {n_samples}")
    print(f"Dimensiones de cada EEM: {n_ex} x {n_em}")

    # Aplanar todas las EEM en una matriz
    batch_matrix = flatten_eem_batch(eem_list, order='C')
    print(f"\nMatriz de lote: {batch_matrix.shape}")
    print(f"  (cada fila es una EEM aplanada)")

    # Mostrar estadísticas de cada muestra
    print(f"\nEstadísticas por muestra:")
    for i, eem in enumerate(eem_list):
        stats = eem.get_stats()
        print(f"  Muestra {i+1}: mean={stats['mean']:.2f}, max={stats['max']:.2f}")

    print()


def ejemplo_pandas_export():
    """Ejemplo de exportación a pandas DataFrame"""
    print("=" * 60)
    print("EJEMPLO 6: Exportación a pandas DataFrame")
    print("=" * 60)

    try:
        import pandas as pd

        # Crear datos de ejemplo pequeños
        excitation_wl = np.array([250, 300, 350])
        emission_wl = np.array([300, 350, 400, 450])
        data = np.random.rand(3, 4) * 100

        eem = FluorescenceEEM(data, excitation_wl, emission_wl)

        # Convertir a DataFrame
        df = eem.to_dataframe()
        print(f"\nDataFrame generado ({len(df)} filas):")
        print(df.head(10))
        print(f"\nColumnas: {list(df.columns)}")

        # Ejemplo de filtrado
        high_intensity = df[df['intensity'] > 50]
        print(f"\nPuntos con intensidad > 50: {len(high_intensity)} puntos")

    except ImportError:
        print("\n[NOTA] pandas no está instalado. Instalar con: pip install pandas")

    print()


def ejemplo_comparacion_ordenes():
    """Ejemplo comparando órdenes de aplanamiento 'C' vs 'F'"""
    print("=" * 60)
    print("EJEMPLO 7: Comparación de órdenes de aplanamiento")
    print("=" * 60)

    # Matriz pequeña para visualizar la diferencia
    data = np.array([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
        [10, 11, 12]
    ])

    print(f"\nMatriz original (4 x 3):")
    print(data)

    eem = FluorescenceEEM(data)

    # Orden 'C' (por filas)
    flat_c = eem.flatten(order='C')
    print(f"\nAplanado orden 'C' (row-major):")
    print(flat_c)
    print("Recorre: [1,2,3,4,5,6,7,8,9,10,11,12]")

    # Orden 'F' (por columnas)
    flat_f = eem.flatten(order='F')
    print(f"\nAplanado orden 'F' (column-major):")
    print(flat_f)
    print("Recorre: [1,4,7,10,2,5,8,11,3,6,9,12]")

    print()


if __name__ == "__main__":
    print("\n" + "="*60)
    print("EJEMPLOS DE USO: fluorescence_eem")
    print("Procesamiento de matrices de espectros de fluorescencia")
    print("="*60 + "\n")

    ejemplo_basico()
    ejemplo_con_longitudes_de_onda()
    ejemplo_normalizacion()
    ejemplo_eliminar_rayleigh()
    ejemplo_batch_processing()
    ejemplo_pandas_export()
    ejemplo_comparacion_ordenes()

    print("="*60)
    print("EJEMPLOS COMPLETADOS")
    print("="*60)
