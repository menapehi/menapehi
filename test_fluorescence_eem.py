"""
Tests unitarios para el módulo fluorescence_eem

Ejecutar con: python -m pytest test_fluorescence_eem.py
o: python test_fluorescence_eem.py
"""

import numpy as np
import unittest
from fluorescence_eem import (
    FluorescenceEEM,
    flatten_eem_matrix,
    reconstruct_eem_matrix,
    flatten_eem_batch
)


class TestFluorescenceEEM(unittest.TestCase):
    """Tests para la clase FluorescenceEEM"""

    def setUp(self):
        """Configuración inicial para cada test"""
        self.n_ex = 10
        self.n_em = 15
        self.data = np.random.rand(self.n_ex, self.n_em) * 100
        self.excitation_wl = np.arange(250, 250 + self.n_ex * 5, 5)
        self.emission_wl = np.arange(300, 300 + self.n_em * 5, 5)

    def test_init_basic(self):
        """Test inicialización básica"""
        eem = FluorescenceEEM(self.data)
        self.assertEqual(eem.data.shape, (self.n_ex, self.n_em))
        self.assertEqual(len(eem.excitation_wavelengths), self.n_ex)
        self.assertEqual(len(eem.emission_wavelengths), self.n_em)

    def test_init_with_wavelengths(self):
        """Test inicialización con longitudes de onda"""
        eem = FluorescenceEEM(self.data, self.excitation_wl, self.emission_wl)
        np.testing.assert_array_equal(eem.excitation_wavelengths, self.excitation_wl)
        np.testing.assert_array_equal(eem.emission_wavelengths, self.emission_wl)

    def test_init_invalid_dimensions(self):
        """Test que falla con dimensiones incorrectas"""
        # Datos 1D
        with self.assertRaises(ValueError):
            FluorescenceEEM(np.array([1, 2, 3]))

        # Datos 3D
        with self.assertRaises(ValueError):
            FluorescenceEEM(np.random.rand(5, 5, 5))

        # Longitudes de onda incorrectas
        with self.assertRaises(ValueError):
            FluorescenceEEM(self.data, np.arange(5), self.emission_wl)

    def test_flatten_order_c(self):
        """Test aplanamiento con orden 'C' (por filas)"""
        eem = FluorescenceEEM(self.data)
        flat = eem.flatten(order='C')

        self.assertEqual(len(flat), self.n_ex * self.n_em)
        # Primer elemento debe ser data[0, 0]
        self.assertEqual(flat[0], self.data[0, 0])
        # Elemento n_em debe ser data[1, 0]
        self.assertEqual(flat[self.n_em], self.data[1, 0])

    def test_flatten_order_f(self):
        """Test aplanamiento con orden 'F' (por columnas)"""
        eem = FluorescenceEEM(self.data)
        flat = eem.flatten(order='F')

        self.assertEqual(len(flat), self.n_ex * self.n_em)
        # Primer elemento debe ser data[0, 0]
        self.assertEqual(flat[0], self.data[0, 0])
        # Elemento n_ex debe ser data[0, 1]
        self.assertEqual(flat[self.n_ex], self.data[0, 1])

    def test_flatten_with_labels(self):
        """Test aplanamiento con etiquetas"""
        eem = FluorescenceEEM(self.data, self.excitation_wl, self.emission_wl)
        flat_data, labels = eem.flatten_with_labels(order='C')

        self.assertEqual(len(flat_data), len(labels))
        self.assertEqual(len(labels), self.n_ex * self.n_em)

        # Verificar primera etiqueta
        self.assertEqual(labels[0], (self.excitation_wl[0], self.emission_wl[0]))

    def test_from_flattened(self):
        """Test reconstrucción desde datos aplanados"""
        eem = FluorescenceEEM(self.data)
        flat = eem.flatten(order='C')

        # Reconstruir
        eem_reconstructed = FluorescenceEEM.from_flattened(
            flat, self.n_ex, self.n_em, order='C'
        )

        np.testing.assert_array_almost_equal(eem.data, eem_reconstructed.data)

    def test_from_flattened_invalid_size(self):
        """Test reconstrucción con tamaño incorrecto"""
        flat = np.array([1, 2, 3, 4, 5])

        with self.assertRaises(ValueError):
            FluorescenceEEM.from_flattened(flat, 3, 3, order='C')

    def test_get_stats(self):
        """Test cálculo de estadísticas"""
        eem = FluorescenceEEM(self.data)
        stats = eem.get_stats()

        self.assertIn('min', stats)
        self.assertIn('max', stats)
        self.assertIn('mean', stats)
        self.assertIn('std', stats)
        self.assertIn('median', stats)

        # Verificar valores
        self.assertAlmostEqual(stats['min'], np.min(self.data))
        self.assertAlmostEqual(stats['max'], np.max(self.data))
        self.assertAlmostEqual(stats['mean'], np.mean(self.data))

    def test_normalize_max(self):
        """Test normalización por máximo"""
        eem = FluorescenceEEM(self.data)
        eem_norm = eem.normalize(method='max')

        self.assertAlmostEqual(np.max(eem_norm.data), 1.0)
        self.assertGreaterEqual(np.min(eem_norm.data), 0.0)

    def test_normalize_area(self):
        """Test normalización por área"""
        eem = FluorescenceEEM(self.data)
        eem_norm = eem.normalize(method='area')

        self.assertAlmostEqual(np.sum(eem_norm.data), 1.0)

    def test_normalize_zscore(self):
        """Test normalización z-score"""
        eem = FluorescenceEEM(self.data)
        eem_norm = eem.normalize(method='zscore')

        self.assertAlmostEqual(np.mean(eem_norm.data), 0.0, places=10)
        self.assertAlmostEqual(np.std(eem_norm.data), 1.0, places=10)

    def test_normalize_invalid_method(self):
        """Test normalización con método inválido"""
        eem = FluorescenceEEM(self.data)

        with self.assertRaises(ValueError):
            eem.normalize(method='invalid_method')

    def test_remove_rayleigh_scatter(self):
        """Test eliminación de dispersión de Rayleigh"""
        # Crear datos donde ex = em
        n = 20
        excitation_wl = np.arange(300, 300 + n * 10, 10)
        emission_wl = np.arange(300, 300 + n * 10, 10)
        data = np.ones((n, n)) * 50

        eem = FluorescenceEEM(data, excitation_wl, emission_wl)
        eem_cleaned = eem.remove_rayleigh_scatter(bandwidth=5)

        # La diagonal debe ser NaN
        for i in range(n):
            self.assertTrue(np.isnan(eem_cleaned.data[i, i]))

        # Valores alejados de la diagonal deben mantenerse
        self.assertFalse(np.isnan(eem_cleaned.data[0, n-1]))

    def test_repr(self):
        """Test representación en string"""
        eem = FluorescenceEEM(self.data, self.excitation_wl, self.emission_wl)
        repr_str = repr(eem)

        self.assertIn('FluorescenceEEM', repr_str)
        self.assertIn('shape', repr_str)


class TestAuxiliaryFunctions(unittest.TestCase):
    """Tests para funciones auxiliares"""

    def setUp(self):
        """Configuración inicial"""
        self.n_ex = 8
        self.n_em = 12
        self.data = np.random.rand(self.n_ex, self.n_em) * 100

    def test_flatten_eem_matrix(self):
        """Test función auxiliar flatten_eem_matrix"""
        flat = flatten_eem_matrix(self.data, order='C')

        self.assertEqual(len(flat), self.n_ex * self.n_em)
        self.assertEqual(flat[0], self.data[0, 0])

    def test_reconstruct_eem_matrix(self):
        """Test función auxiliar reconstruct_eem_matrix"""
        flat = flatten_eem_matrix(self.data, order='C')
        reconstructed = reconstruct_eem_matrix(flat, self.n_ex, self.n_em, order='C')

        np.testing.assert_array_almost_equal(self.data, reconstructed)

    def test_flatten_eem_batch(self):
        """Test procesamiento por lotes"""
        n_samples = 5
        eem_list = []

        for i in range(n_samples):
            data = np.random.rand(self.n_ex, self.n_em) * 100
            eem = FluorescenceEEM(data)
            eem_list.append(eem)

        batch_matrix = flatten_eem_batch(eem_list, order='C')

        self.assertEqual(batch_matrix.shape, (n_samples, self.n_ex * self.n_em))

    def test_flatten_eem_batch_empty(self):
        """Test procesamiento por lotes con lista vacía"""
        with self.assertRaises(ValueError):
            flatten_eem_batch([], order='C')

    def test_flatten_eem_batch_different_shapes(self):
        """Test procesamiento por lotes con formas diferentes"""
        eem1 = FluorescenceEEM(np.random.rand(10, 15))
        eem2 = FluorescenceEEM(np.random.rand(10, 20))  # Diferente forma

        with self.assertRaises(ValueError):
            flatten_eem_batch([eem1, eem2], order='C')


class TestRoundTrip(unittest.TestCase):
    """Tests de ida y vuelta (flatten -> reconstruct)"""

    def test_roundtrip_small_matrix(self):
        """Test ida y vuelta con matriz pequeña"""
        data = np.array([[1, 2, 3], [4, 5, 6]])
        eem = FluorescenceEEM(data)

        flat = eem.flatten(order='C')
        eem_reconstructed = FluorescenceEEM.from_flattened(flat, 2, 3, order='C')

        np.testing.assert_array_equal(data, eem_reconstructed.data)

    def test_roundtrip_large_matrix(self):
        """Test ida y vuelta con matriz grande"""
        data = np.random.rand(50, 80) * 1000
        eem = FluorescenceEEM(data)

        flat = eem.flatten(order='F')
        eem_reconstructed = FluorescenceEEM.from_flattened(flat, 50, 80, order='F')

        np.testing.assert_array_almost_equal(data, eem_reconstructed.data)

    def test_roundtrip_with_wavelengths(self):
        """Test ida y vuelta preservando longitudes de onda"""
        excitation = np.linspace(250, 450, 40)
        emission = np.linspace(300, 600, 60)
        data = np.random.rand(40, 60) * 100

        eem = FluorescenceEEM(data, excitation, emission)
        flat = eem.flatten(order='C')

        eem_reconstructed = FluorescenceEEM.from_flattened(
            flat, 40, 60, order='C',
            excitation_wavelengths=excitation,
            emission_wavelengths=emission
        )

        np.testing.assert_array_almost_equal(data, eem_reconstructed.data)
        np.testing.assert_array_equal(excitation, eem_reconstructed.excitation_wavelengths)
        np.testing.assert_array_equal(emission, eem_reconstructed.emission_wavelengths)


class TestEdgeCases(unittest.TestCase):
    """Tests de casos límite"""

    def test_single_element_matrix(self):
        """Test con matriz de un solo elemento"""
        data = np.array([[42.0]])
        eem = FluorescenceEEM(data)

        flat = eem.flatten()
        self.assertEqual(len(flat), 1)
        self.assertEqual(flat[0], 42.0)

    def test_matrix_with_zeros(self):
        """Test con matriz de ceros"""
        data = np.zeros((5, 7))
        eem = FluorescenceEEM(data)

        stats = eem.get_stats()
        self.assertEqual(stats['min'], 0.0)
        self.assertEqual(stats['max'], 0.0)
        self.assertEqual(stats['mean'], 0.0)

    def test_matrix_with_nan(self):
        """Test con valores NaN"""
        data = np.random.rand(10, 15)
        data[5, 5] = np.nan
        data[3, 8] = np.nan

        eem = FluorescenceEEM(data)
        flat = eem.flatten()

        # Debe preservar NaN
        self.assertTrue(np.isnan(flat).sum() >= 2)

    def test_matrix_with_negative_values(self):
        """Test con valores negativos"""
        data = np.random.randn(8, 10)  # Puede tener negativos
        eem = FluorescenceEEM(data)

        stats = eem.get_stats()
        # Debe funcionar incluso con negativos
        self.assertIsNotNone(stats['min'])
        self.assertIsNotNone(stats['max'])


def run_tests():
    """Ejecutar todos los tests"""
    unittest.main(verbosity=2)


if __name__ == '__main__':
    print("="*60)
    print("EJECUTANDO TESTS UNITARIOS: fluorescence_eem")
    print("="*60)
    print()
    run_tests()
