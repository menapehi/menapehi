"""
Módulo para procesamiento de matrices de espectros de fluorescencia de excitación-emisión (EEM).

Este módulo proporciona funciones para aplanar y reconstruir matrices EEM,
que son comúnmente utilizadas en análisis de fluorescencia de muestras químicas
y ambientales.

Author: @menapehi
"""

import numpy as np
from typing import Tuple, Optional, Dict, List


class FluorescenceEEM:
    """
    Clase para manejar matrices de espectros de fluorescencia de excitación-emisión.

    Una matriz EEM es una representación 2D donde:
    - Filas: longitudes de onda de excitación
    - Columnas: longitudes de onda de emisión
    - Valores: intensidad de fluorescencia
    """

    def __init__(self,
                 data: np.ndarray,
                 excitation_wavelengths: Optional[np.ndarray] = None,
                 emission_wavelengths: Optional[np.ndarray] = None):
        """
        Inicializa una matriz EEM.

        Args:
            data: Matriz 2D con los datos de fluorescencia
            excitation_wavelengths: Array 1D con las longitudes de onda de excitación
            emission_wavelengths: Array 1D con las longitudes de onda de emisión
        """
        if data.ndim != 2:
            raise ValueError("Los datos deben ser una matriz 2D")

        self.data = np.array(data)
        self.n_excitation, self.n_emission = data.shape

        # Longitudes de onda por defecto si no se proporcionan
        self.excitation_wavelengths = (
            excitation_wavelengths if excitation_wavelengths is not None
            else np.arange(self.n_excitation)
        )
        self.emission_wavelengths = (
            emission_wavelengths if emission_wavelengths is not None
            else np.arange(self.n_emission)
        )

        # Validar dimensiones
        if len(self.excitation_wavelengths) != self.n_excitation:
            raise ValueError("La longitud de excitation_wavelengths debe coincidir con las filas de data")
        if len(self.emission_wavelengths) != self.n_emission:
            raise ValueError("La longitud de emission_wavelengths debe coincidir con las columnas de data")

    def flatten(self, order: str = 'C') -> np.ndarray:
        """
        Aplana la matriz EEM en un vector 1D.

        Args:
            order: Orden de aplanamiento
                'C' = row-major (C-style, por filas) - por defecto
                'F' = column-major (Fortran-style, por columnas)

        Returns:
            Vector 1D con los datos aplanados
        """
        return self.data.flatten(order=order)

    def flatten_with_labels(self, order: str = 'C') -> Tuple[np.ndarray, List[Tuple[float, float]]]:
        """
        Aplana la matriz EEM y retorna las etiquetas de cada punto.

        Args:
            order: Orden de aplanamiento ('C' o 'F')

        Returns:
            Tupla con:
            - Vector 1D con los datos aplanados
            - Lista de tuplas (excitación, emisión) para cada punto
        """
        flat_data = self.flatten(order=order)

        # Crear etiquetas
        labels = []
        if order == 'C':
            # Recorrer por filas
            for i, ex in enumerate(self.excitation_wavelengths):
                for j, em in enumerate(self.emission_wavelengths):
                    labels.append((ex, em))
        else:
            # Recorrer por columnas
            for j, em in enumerate(self.emission_wavelengths):
                for i, ex in enumerate(self.excitation_wavelengths):
                    labels.append((ex, em))

        return flat_data, labels

    def to_dataframe(self):
        """
        Convierte la matriz EEM aplanada a un DataFrame de pandas.

        Returns:
            DataFrame con columnas: excitation, emission, intensity
        """
        try:
            import pandas as pd
        except ImportError:
            raise ImportError("pandas es requerido para esta función. Instalar con: pip install pandas")

        flat_data, labels = self.flatten_with_labels()

        return pd.DataFrame({
            'excitation_nm': [label[0] for label in labels],
            'emission_nm': [label[1] for label in labels],
            'intensity': flat_data
        })

    @classmethod
    def from_flattened(cls,
                      flat_data: np.ndarray,
                      n_excitation: int,
                      n_emission: int,
                      order: str = 'C',
                      excitation_wavelengths: Optional[np.ndarray] = None,
                      emission_wavelengths: Optional[np.ndarray] = None) -> 'FluorescenceEEM':
        """
        Reconstruye una matriz EEM desde datos aplanados.

        Args:
            flat_data: Vector 1D con los datos aplanados
            n_excitation: Número de longitudes de onda de excitación
            n_emission: Número de longitudes de onda de emisión
            order: Orden usado en el aplanamiento ('C' o 'F')
            excitation_wavelengths: Longitudes de onda de excitación
            emission_wavelengths: Longitudes de onda de emisión

        Returns:
            Objeto FluorescenceEEM con la matriz reconstruida
        """
        if len(flat_data) != n_excitation * n_emission:
            raise ValueError(
                f"El tamaño del vector aplanado ({len(flat_data)}) "
                f"no coincide con n_excitation * n_emission ({n_excitation * n_emission})"
            )

        data = flat_data.reshape((n_excitation, n_emission), order=order)
        return cls(data, excitation_wavelengths, emission_wavelengths)

    def get_stats(self) -> Dict[str, float]:
        """
        Calcula estadísticas básicas de la matriz EEM.

        Returns:
            Diccionario con estadísticas: min, max, mean, std, median
        """
        return {
            'min': float(np.min(self.data)),
            'max': float(np.max(self.data)),
            'mean': float(np.mean(self.data)),
            'std': float(np.std(self.data)),
            'median': float(np.median(self.data)),
            'shape': self.data.shape,
            'total_points': self.data.size
        }

    def remove_rayleigh_scatter(self, bandwidth: int = 10) -> 'FluorescenceEEM':
        """
        Elimina la dispersión de Rayleigh (excitación = emisión ± bandwidth).

        Args:
            bandwidth: Ancho de banda en unidades de longitud de onda

        Returns:
            Nueva instancia de FluorescenceEEM con dispersión eliminada
        """
        data_cleaned = self.data.copy()

        for i, ex in enumerate(self.excitation_wavelengths):
            for j, em in enumerate(self.emission_wavelengths):
                if abs(em - ex) <= bandwidth:
                    data_cleaned[i, j] = np.nan

        return FluorescenceEEM(data_cleaned,
                              self.excitation_wavelengths,
                              self.emission_wavelengths)

    def normalize(self, method: str = 'max') -> 'FluorescenceEEM':
        """
        Normaliza la matriz EEM.

        Args:
            method: Método de normalización
                'max': divide por el valor máximo
                'area': normalización por área (suma = 1)
                'zscore': normalización z-score

        Returns:
            Nueva instancia de FluorescenceEEM normalizada
        """
        data_norm = self.data.copy()

        if method == 'max':
            max_val = np.nanmax(data_norm)
            if max_val != 0:
                data_norm = data_norm / max_val
        elif method == 'area':
            total_area = np.nansum(data_norm)
            if total_area != 0:
                data_norm = data_norm / total_area
        elif method == 'zscore':
            mean = np.nanmean(data_norm)
            std = np.nanstd(data_norm)
            if std != 0:
                data_norm = (data_norm - mean) / std
        else:
            raise ValueError(f"Método de normalización '{method}' no reconocido")

        return FluorescenceEEM(data_norm,
                              self.excitation_wavelengths,
                              self.emission_wavelengths)

    def __repr__(self):
        return (f"FluorescenceEEM(shape={self.data.shape}, "
                f"ex_range=[{self.excitation_wavelengths[0]}, {self.excitation_wavelengths[-1]}], "
                f"em_range=[{self.emission_wavelengths[0]}, {self.emission_wavelengths[-1]}])")


def flatten_eem_batch(eem_list: List[FluorescenceEEM], order: str = 'C') -> np.ndarray:
    """
    Aplana múltiples matrices EEM en una matriz 2D.

    Args:
        eem_list: Lista de objetos FluorescenceEEM
        order: Orden de aplanamiento ('C' o 'F')

    Returns:
        Matriz 2D donde cada fila es una EEM aplanada
    """
    if not eem_list:
        raise ValueError("La lista de EEM está vacía")

    # Validar que todas tengan la misma forma
    first_shape = eem_list[0].data.shape
    for i, eem in enumerate(eem_list):
        if eem.data.shape != first_shape:
            raise ValueError(f"EEM {i} tiene forma {eem.data.shape}, esperada {first_shape}")

    # Aplanar todas
    flat_eems = [eem.flatten(order=order) for eem in eem_list]

    return np.array(flat_eems)


# Funciones auxiliares para compatibilidad y uso simple
def flatten_eem_matrix(data: np.ndarray, order: str = 'C') -> np.ndarray:
    """
    Función auxiliar simple para aplanar una matriz EEM.

    Args:
        data: Matriz 2D con datos de fluorescencia
        order: Orden de aplanamiento ('C' o 'F')

    Returns:
        Vector 1D aplanado
    """
    eem = FluorescenceEEM(data)
    return eem.flatten(order=order)


def reconstruct_eem_matrix(flat_data: np.ndarray,
                          n_excitation: int,
                          n_emission: int,
                          order: str = 'C') -> np.ndarray:
    """
    Función auxiliar simple para reconstruir una matriz EEM desde datos aplanados.

    Args:
        flat_data: Vector 1D con datos aplanados
        n_excitation: Número de longitudes de onda de excitación
        n_emission: Número de longitudes de onda de emisión
        order: Orden usado en el aplanamiento ('C' o 'F')

    Returns:
        Matriz 2D reconstruida
    """
    eem = FluorescenceEEM.from_flattened(flat_data, n_excitation, n_emission, order=order)
    return eem.data
